#' Fit zero-inflated random forests given a pre-specified set of
#' tuning parameters.
#'
#' @param x  Matrix or data.frame of covariates for the count part
#'             of the model.
#' @param z  Matrix or data.frame of covariates for the zero-inflation
#'            part of the model.
#' @param y  vector of count outcomes.
#' @param rounds  how many iterations the EM algorithm should be run.
#' @param mtry  which value of mtry to use in each forest.
#' @param ntree  number of trees
#' @param nodesize  controls depth of regression trees
#' @param perm_frac what percentage of observations should be used for fitting the model
#'                  versus estimating the vims
#' @param count_coef values for initial estimate of count model.
#' @param zero_coef values for initial estimate of zero-inflation model.
#' @param newx count model covariate matrix for new data.
#' @param newz zero inflation model covariate matrix for new data.
#' @param iter_keep_pi boolean, if TRUE values of pi are kept across EM iterations
#' @param iter_keep_preds boolean, if TRUE prediction values are kept across EM iterations
#' @param iter_keep_xsi boolean, if TRUE zero-inflation model coefficients
#'                      are kept across EM iterations
#' @param nlambda number of folds of cross-validation
#' @param dimension_reduction If TRUE, dimension reduction will be applied.  In this case,
#'                            permutation vims will be calculated.
#' @param MAST_pval matrix of pvalues sizes from MAST.
#' @param pval_threshold threshold for pvalues for screening covariates based on marginal
#'                       effect sizes.
#' @param MAST_fc matrix of effect sizes from MAST.
#' @param fc_threshold threshold for effect sizes for screening covariates based on marginal
#'                     effect sizes.
#' @param pca_genes genes used pca
#' @param inputTF input transcription factors
#' @return A list.
#' @export
hyper_zirf_fit <- function(x, z, y, rounds, mtry,
                           ntree=100, nodesize=5,
                           perm_frac=.75,
                           count_coef=NULL, zero_coef=NULL,
                           newx=NULL, newz=NULL,
                           iter_keep_pi=F, iter_keep_preds=F,
                           iter_keep_xsi=F, nlambda=10,
                           dimension_reduction=T,
                           MAST_pval=NULL,
                           pval_threshold=NULL,
                           MAST_fc=NULL,
                           fc_threshold=NULL,
                           pca_genes=NULL,
                           inputTF=NULL){
  cgene <- colnames(y)
  yname <- deparse(substitute(y))
  n <- dim(x)[1]; p <- dim(x)[2]
  perm_ind <- rep(FALSE, n)
  if(!is.null(inputTf)){
    if(cgene %in% inputTF){
      isTF <- TRUE
      cgene_ind <- which(inputTF %in% cgene)
    } else {
      isTF <- FALSE
    }
    inputTF <- setdiff(inputTF, cgene)
  }
  #split data into two parts:one for fitting the model, the other for calculating vims
  perm_ind[sample(1:n, floor(perm_frac*n))] <- TRUE
  xperm <- x[perm_ind, ]; zperm <- z[perm_ind, , drop=F]; yperm <- y[perm_ind]
  x <- x[!perm_ind, ]; z <- z[!perm_ind, , drop=F]; y <- y[!perm_ind]
  x_ind <- 1:dim(x)[2]
  z_ind <- (dim(x)[2] + 1):(dim(x)[2] + dim(z)[2])

  #convert matrices to data.frames
  if(class(x) == "matrix"){
      x <- as.data.frame(x, stringsAsFactors=FALSE)
  }
  if(class(z) == "matrix"){
    z <- as.data.frame(z, stringsAsFactors=FALSE)
  }
  if(!is.null(newx)){
    if(class(newx) == "matrix"){
        newx <- as.data.frame(newx, stringsAsFactors=FALSE)
       }
    names(newx) <- names(x)
  }
  if(class(y) == "matrix"){y <- as.data.frame(y, stringsAsFactors=FALSE)}
  #

  #set up test data set
  if(is.null(newx)){
    newx <- as.data.frame(x)
    newz <- as.matrix(z)
    newz <- as.matrix(cbind(rep(1, nrow(newz)), newz))
  }
  else {
    newx <- as.data.frame(newx)
    newz <- as.matrix(newz)
    newz <- as.matrix(cbind(rep(1, nrow(newz)), newz))
  }

  #create data.frame that will be used for fitting inital zip model
  zilm_dat <- data.frame(x, z, y)
  names(zilm_dat)[dim(zilm_dat)[2]] <- yname

  #create formula
  mod_formula <- XZ_to_form(yname, x, z)

  #calculate PCA scores
  pca_dat <- x[, pca_genes]
  pca_fit <- corpcor::fast.svd(scale(pca_dat))
  PCs <- as.matrix(pca_fit$v[, 1:2])
  colnames(PCs) <- c("PC1", "PC2")
  PC_scores <- as.matrix(pca_dat)%*%PCs

  #calculate PCA scores for the new data set
  new_pca_dat <- newx[, pca_genes]
  new_PC_scores <- as.matrix(new_pca_dat)%*%PCs

  #initial estimate of count model
  pois_mod <- mpath::cv.zipath(mod_formula, data = zilm_dat,  nlambda=nlambda,
                        penalty.factor.zero=0,
                        family="poisson",
                        plot.it = F)

  #now filter out TFs based on marginal effects
  if(!is.null(MAST_pval)){
    pval_TFs <- colnames(as.matrix(MAST_pval))[MAST_pval[cgene, ] < pval_threshold]
    keep_TFs <- pval_TFs
  }
  if(!is.null(MAST_fc)){
    fc_TFs <- colnames(as.matrix(MAST_fc))[MAST_fc[cgene, ] < fc_threshold]
    keep_TFs <- intersect(keep_TFs, fc_TFs)
  }

  ##apply initial E-step to get initial probabilities of inclusion
  pois_coef <- stats::coef(pois_mod)
  zero_coef <- pois_coef$zero
  count_coef <- pois_coef$count
  z_mat_int <- as.matrix(cbind(rep(1, nrow(z)), z))
  #x is a data.frame without the intercept.
  #This is appropriate for using with randomForest, but it won't work for
  #matrix multiplication.  We need a matrix instead of a data.frame and
  #we need the intercept.
  x_mat_int <- as.matrix(cbind(rep(1, nrow(x)), x))
  linear_cut <- x_mat_int%*%count_coef

  #filter out TFs screened by MAST
  filter_TFs <- function(x){
    #keep all nonTF columns
    ind_noninput_TF <- which(!colnames(x) %in% inputTF)
    #keep the TFs selected by MAST
    ind_keep_TF <- which(colnames(x) %in% keep_TFs)
    x <- x[, c(ind_noninput_TF,  ind_keep_TF)]
  }
  x <- data.frame(x, PC_scores, linear_cut=linear_cut)
  x <- filter_TFs(x)
  newx_linear_cut <- as.matrix(cbind(rep(1, dim(newx)[1]), newx))%*%count_coef
  newx <- data.frame(newx, new_PC_scores, linear_cut=newx_linear_cut)
  newx <- filter_TFs(newx)

  #create data.frame that will be used for ranger fit
  ranger_dat <- data.frame(y, x)

  #Next apply formula (9) of Wang et al (2014) from Stats in Medicine
  #We know z_{i}=0 for these observations
  mtry <- dim(x)[2]
  y1 <- y > 0
  logit_pi <- z_mat_int%*%zero_coef
  log_init_pois <- x_mat_int%*%count_coef
  expit <- function(x){1/(1 + exp(-x))}
  probi <- expit(logit_pi)
  mui <- exp(log_init_pois)
  probi <- probi/(probi + (1 - probi)*stats::dpois(0, mui))
  probi[y1] <- 0
  #edit this now
  n <- dim(x)[1]
  nnew <- dim(newx)[1]
  keep_preds <- matrix(NA, nrow=rounds, ncol=n)
  keep_pi <- matrix(NA, nrow=rounds, ncol=n)
  keep_xsi <- matrix(NA, nrow=rounds, ncol=length(zero_coef))
  for(i in 1:rounds){
      inbag_arg <- vector("list", length=ntree)
      for(k in 1:ntree){
        inbag_arg[[k]] <- stats::rbinom(n, 1, 1 - probi)
      }
      rf_fit <- ranger::ranger(y ~ ., data=ranger_dat, mtry=mtry,
                               write.forest = T, importance="none",
                               num.trees = ntree, inbag=inbag_arg)
      rf_preds_z0 <- stats::predict(rf_fit, data=x)$predictions
      if(i == rounds){
        new_rf_preds_z0 <- stats::predict(rf_fit, data=newx)$predictions
      }
      mui <- rf_preds_z0
      mui <- round(mui, 12)
      model_zero <- suppressWarnings(stats::glm.fit(z_mat_int, probi,
                              family = stats::binomial(link = "logit"),
                              start = zero_coef,
                              control=list(maxit=1000)))
      zero_coef <-  stats::coef(model_zero)
      logit_pi <- z_mat_int%*%zero_coef
      probi <- expit(logit_pi)
      probi <- probi/(probi + (1 - probi)*stats::dpois(0, mui))
      probi[y1] <- 0
      if(iter_keep_pi){
        keep_pi[i, ] <- probi
      }
      if(iter_keep_preds){
        keep_preds[i, ] <- mui
      }
      if(iter_keep_xsi){
        keep_xsi[i, ] <- zero_coef
      }
    }
  new_log_struc_zero <- newz%*%zero_coef
  new_probi <- 1/(1 + exp(-new_log_struc_zero))
  new_fit <- data.frame(lambda_preds=new_rf_preds_z0,
                        count_preds=new_rf_preds_z0*(1 - new_probi),
                        unconditional_pi=new_probi)

  unconditional_probi <- 1/(1 + exp(-logit_pi))
  data_fit <- data.frame(lambda_preds=rf_preds_z0,
                         count_preds=rf_preds_z0*(1 - unconditional_probi),
                         fitted_pi=probi,
                         unconditional_pi=unconditional_probi)
  if(!iter_keep_preds){
    keep_preds = NULL
  }
  if(!iter_keep_pi){
    keep_pi = NULL
  }

  #now compute permutation variable importance measures
  df_mult <- function(x, coef){
      out <- as.matrix(cbind(rep(1, dim(x)[1]), x))%*%coef
      }
  m <- 2000
  xperm_perturb <- apply(xperm, 2,
                         function(x){sample(x, m, replace=T)})
  xperm <- xperm[sample(1:dim(xperm)[1], m, replace = T), ]
  xperm_linear_cut <- df_mult(xperm, count_coef)
  xperm_PCs <- as.matrix(xperm[, pca_genes])%*%PCs
  colnames(xperm_PCs) <- c("PC1", "PC2")
  xperm <- data.frame(xperm, xperm_PCs, linear_cut=xperm_linear_cut)
  ind_noninput_TF <- which(!colnames(xperm) %in% inputTF)
  ind_keep_TF <- which(colnames(xperm) %in% keep_TFs)
  ranger_ind <- c(ind_noninput_TF,  ind_keep_TF)
  l <- 1
  vim_ests <- rep(NA, length(inputTF))
  for(l in 1:length(inputTF)){
    perturbed_covariate <- xperm
    perturbed_covariate[, inputTF[l]] <- xperm_perturb[, inputTF[l]]
    perturbed_covariate$linear_cut <- df_mult(perturbed_covariate[, names(count_coef[-1])],
                                              count_coef)
    #form PCs for perturbed covariate
    perturbed_PCs <- as.matrix(perturbed_covariate[, pca_genes])%*%PCs
    colnames(perturbed_PCs) <- c("PC1", "PC2")
    perturbed_covariate[, c("PC1", "PC2")] <- perturbed_PCs

    ranger_perm <- xperm[, ranger_ind]
    ranger_perturbed <- perturbed_covariate[, ranger_ind]

    perturbed_fit <- stats::predict(rf_fit, data=ranger_perturbed)$predictions
    vim_fit <- stats::predict(rf_fit, data=ranger_perm)$predictions
    vim_ests[l] <- mean((perturbed_fit - vim_fit)^2)
  }

  if(!isTF){
    names(vim_ests) <- inputTF
  }
  if(isTF){
    vim_ests_new <- rep(NA, length(inputTF) + 1)
    vim_ests_new[cgene_ind] <- 0
    names(vim_ests_new)[cgene_ind] <- cgene
    vim_ests_new[-cgene_ind] <- vim_ests
    names(vim_ests_new)[-cgene_ind] <- inputTF
    vim_ests <- vim_ests_new
  }
  out <- list(new_fit=new_fit, data_fit=data_fit,
              importance_measures=vim_ests,
              zero_coef=zero_coef,
              iter_keep_preds=keep_preds,
              iter_keep_pi=keep_pi,
              iter_keep_xsi=keep_xsi)
}











