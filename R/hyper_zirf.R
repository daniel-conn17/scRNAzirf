#' Fit zero-inflated random forests given a pre-specified set of
#' tuning parameters.
#'
#' @param x  Matrix or data.frame of covariates for the count part
#'             of the model.
#' @param z  Matrix or data.frame of covariates for the zero-inflation
#'            part of the model.
#' @param y  vector of count outcomes.
#' @param rounds  how many iterations the EM algorithm should be run.
#' @param mtry which value of mtry to use in each forest.
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
#' @param mtry_factor multiplier of mtry used when dimension_reduction=T.  mtry is set to
#'                    the number of screened variables times mtry_factor.
#' @param pvals vector of pvalues sizes
#' @param log_fcs vector of log-fold changes
#' @param pval_threshold threshold for pvalues for screening covariates based on marginal
#'                       effect sizes.
#' @param log_fc_threshold threshold for log fold-changes
#' @return A list.
#' @export
hyper_zirf_fit <- function(x, z, y, rounds, mtry,
                           ntree=100, nodesize=5,
                           perm_frac=.15,
                           count_coef=NULL, zero_coef=NULL,
                           newx=NULL, newz=NULL,
                           iter_keep_pi=F, iter_keep_preds=F,
                           iter_keep_xsi=F, nlambda=10,
                           dimension_reduction=T,
                           mtry_factor=1,
                           pvals=NULL,
                           log_fcs=NULL,
                           pval_threshold=NULL,
                           log_fc_threshold=NULL){
  cgene <- colnames(y)
  yname <- deparse(substitute(y))
  n <- dim(x)[1]; p <- dim(x)[2]
  perm_ind <- rep(FALSE, n)
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

  #we can add pca_genes as argument if desired later
  #calculate PCA scores
  #pca_dat <- x[, pca_genes]
  pca_dat <- x
  pca_fit <- corpcor::fast.svd(scale(pca_dat))
  PCs <- as.matrix(pca_fit$v[, 1:2])
  colnames(PCs) <- c("PC1", "PC2")
  PC_scores <- as.matrix(pca_dat)%*%PCs
  #calculate PCA scores for the new data set
  #new_pca_dat <- newx[, pca_genes]
  new_pca_dat <- newx
  new_PC_scores <- as.matrix(new_pca_dat)%*%PCs

  #now filter out TFs based on marginal effects
  if(!is.null(pvals)){
    pval_keep_ind <- which(pvals < pval_threshold)
  }
  if(!is.null(log_fcs)){
    log_fcs_keep_ind <- which(abs(log_fcs) > log_fcs)
  }
  if(!is.null(pvals) & !is.null(log_fcs)){
    cov_keep <- intersect(pval_keep_ind, log_fcs_keep_ind)
  }
  if(!is.null(pvals) & is.null(log_fcs)){
    cov_keep <- pval_keep_ind
  }
  if(!is.null(log_fcs) & is.null(pvals)){
    cov_keep <- log_fcs_keep_ind
  }

  #create formula
  mod_formula <- XZ_to_form(yname, x, z)
  zilm_dat <- data.frame(x, z, y)
  names(zilm_dat)[dim(zilm_dat)[2]] <- yname
  #initial estimate of count model
  pois_mod <- mpath::cv.zipath(mod_formula, data = zilm_dat,  nlambda=nlambda,
                               penalty.factor.zero=0,
                               family="poisson",
                               plot.it = F)
  ##apply initial E-step to get initial probabilities of inclusion
  pois_coef <- stats::coef(pois_mod)
  zero_coef <- pois_coef$zero
  count_coef <- pois_coef$count
  #this is just for simulations; must delete when done with testing.
  #what happens if we use the right values for linear cuts
  #count_coef <-  c(-1, 1, .5, rep(0, 20 - 3))
  #zero_coef <- c(rep(-.25, 0), rep(0, 20 - 1))
  count_coef_full <- count_coef
  z_mat_int <- as.matrix(cbind(rep(1, nrow(z)), z))
  #x is a data.frame without the intercept.
  #This is appropriate for using with randomForest, but it won't work for
  #matrix multiplication.  We need a matrix instead of a data.frame and
  #we need the intercept.
  x_mat_int <- as.matrix(cbind(rep(1, nrow(x)), x))
  linear_cut <- x_mat_int%*%count_coef
  log_init_pois <- linear_cut
  z_mat_int <- as.matrix(cbind(rep(1, nrow(z)), z))
  logit_pi <- z_mat_int%*%zero_coef

  x <- x[, cov_keep, drop=F]
  x <- data.frame(PC_scores, linear_cut=linear_cut)
  #x <- data.frame(x, PC_scores, linear_cut=linear_cut)
  newx_linear_cut <- as.matrix(cbind(rep(1, dim(newx)[1]), newx))%*%count_coef
  newx <- newx[, cov_keep, drop=F]
  newx <- data.frame(newx, new_PC_scores, linear_cut=newx_linear_cut)

  ##get new fit to poisson model to
  ##create formula
  #x <- x[, c("PC1", "PC2", "linear_cut")]
  #mod_formula <- XZ_to_form(yname, x, z)
  #zilm_dat <- data.frame(x, z, y)
  #names(zilm_dat)[dim(zilm_dat)[2]] <- yname
  ##initial estimate of count model
  #pois_mod <- mpath::cv.zipath(mod_formula, data = zilm_dat,  nlambda=nlambda,
  #                             penalty.factor.zero=0,
  #                             family="poisson",
  #                             plot.it = F)
  ##apply initial E-step to get initial probabilities of inclusion
  #pois_coef <- stats::coef(pois_mod)
  #zero_coef <- pois_coef$zero
  #count_coef <- pois_coef$count

  #Next apply formula (9) of Wang et al (2014) from Stats in Medicine
  #We know z_{i}=0 for these observations
  mtry <- dim(x)[2]
  mtry <- floor(mtry_factor*mtry)
  y1 <- y > 0
  expit <- function(x){1/(1 + exp(-x))}
  probi <- expit(logit_pi)
  mui <- exp(log_init_pois)
  probi <- probi/(probi + (1 - probi)*stats::dpois(0, mui))
  probi[y1] <- 0
  n <- dim(x)[1]
  nnew <- dim(newx)[1]
  keep_preds <- matrix(NA, nrow=rounds, ncol=n)
  keep_pi <- matrix(NA, nrow=rounds, ncol=n)
  keep_xsi <- matrix(NA, nrow=rounds, ncol=length(zero_coef))
  #create data.frame that will be used for ranger fit
  ranger_dat <- data.frame(y, x)
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
  df_mult <- function(x, coef, return_mat=T){
    out <- as.matrix(cbind(rep(1, dim(x)[1]), x))%*%coef
    if(!return_mat){
      out <- out[, 1]
    }
    out
  }
  m <- 2000
  xperm_perturb <- apply(xperm, 2,
                         function(x){sample(x, m, replace=T)})
  xperm <- xperm[sample(1:dim(xperm)[1], m, replace = T), ]
  xperm_linear_cut <- df_mult(xperm, count_coef_full)
  #xperm_PCs <- as.matrix(xperm[, pca_genes])%*%PCs
  xperm_PCs <- as.matrix(xperm)%*%PCs
  colnames(xperm_PCs) <- c("PC1", "PC2")
  vim_ests <- rep(NA, dim(xperm)[2])
  for(l in 1:dim(xperm)[2]){
    perturbed_covariate <- xperm
    perturbed_covariate[, l] <- xperm_perturb[, l]
    perturbed_covariate_linear_cut <- df_mult(perturbed_covariate[, names(count_coef_full[-1])],
                                              count_coef_full, return_mat = FALSE)
    #form PCs for perturbed covariate
    #perturbed_PCs <- as.matrix(perturbed_covariate[, pca_genes])%*%PCs
    perturbed_PCs <- as.matrix(perturbed_covariate)%*%PCs
    colnames(perturbed_PCs) <- c("PC1", "PC2")
    ranger_perturbed <- perturbed_covariate[, cov_keep, drop=F]
    ranger_perturbed <- data.frame(ranger_perturbed,
                                   perturbed_PCs,
                                   linear_cut=perturbed_covariate_linear_cut)
    #ranger_perturbed <- data.frame(perturbed_PCs,
    #                               linear_cut=perturbed_covariate_linear_cut)

    ranger_perm <- xperm[, cov_keep, drop=F]
    ranger_perm <- data.frame(ranger_perm,
                              xperm_PCs,
                              linear_cut=xperm_linear_cut)
    #ranger_perm <- data.frame(xperm_PCs,
    #                          linear_cut=xperm_linear_cut)

    perturbed_fit <- stats::predict(rf_fit, data=ranger_perturbed)$predictions
    vim_fit <- stats::predict(rf_fit, data=ranger_perm)$predictions
    vim_ests[l] <- mean((perturbed_fit - vim_fit)^2)
  }
  out <- list(new_fit=new_fit, data_fit=data_fit,
              importance_measures=vim_ests,
              zero_coef=zero_coef,
              iter_keep_preds=keep_preds,
              iter_keep_pi=keep_pi,
              iter_keep_xsi=keep_xsi)
}
