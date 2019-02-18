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
#' @return A list.
#' @export
hyper_zirf_fit <- function(x, z, y, rounds, mtry,
                           ntree=100, nodesize=5,
                           count_coef=NULL, zero_coef=NULL,
                           newx=NULL, newz=NULL,
                           iter_keep_pi=F, iter_keep_preds=F,
                           iter_keep_xsi=F, nlambda=10, dimension_reduction=T){
  yname <- deparse(substitute(y))
  zero_coef_og <- zero_coef
  x_ind <- 1:dim(x)[2]
  z_ind <- (dim(x)[2] + 1):(dim(x)[2] + dim(z)[2])
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
  if(class(y) == "matrix"){y <- as.data.fame(y, stringsAsFactors=FALSE)}
  zilm_dat <- data.frame(x, z, y)
  names(zilm_dat)[dim(zilm_dat)[2]] <- yname
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
  mod_formula <- XZ_to_form(yname, x, z)

  pois_mod <- mpath::cv.zipath(mod_formula, data = zilm_dat,  nlambda=nlambda,
                        penalty.factor.zero=0,
                        family="poisson",
                        plot.it = F)

  ##apply initial E-step to get initial probabilities of inclusion
  pois_coef <- stats::coef(pois_mod)
  zero_coef <- pois_coef$zero
  count_coef <- pois_coef$count
  z_mat_int <- zilm_dat[, z_ind, drop=F]
  z_mat_int <- as.matrix(cbind(rep(1, nrow(z_mat_int)), z_mat_int))
  #x is a data.frame without the intercept.
  #This is appropriate for using with randomForest, but it won't work for
  #matrix multiplication.  We need a matrix instead of a data.frame and
  #we need the intercept.
  x_mat_int <- zilm_dat[, x_ind, drop=F]
  x_mat_int <- as.matrix(cbind(rep(1, nrow(x_mat_int)), x_mat_int))
  linear_cut <- x_mat_int%*%count_coef
  x <- data.frame(x, linear_cut=linear_cut)
  newx_linear_cut <- as.matrix(cbind(rep(1, dim(newx)[1]), newx))%*%count_coef
  newx <- data.frame(newx, linear_cut=newx_linear_cut)
  #Next apply formula (9) of Wang et al (2014) from Stats in Medicine
  #We know z_{i}=0 for these observations
  sure_nonzero_ind <- which(y != 0)
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
  #ptm <- proc.time()
  tree_preds <- matrix(NA, nrow=n, ncol=ntree)
  newdat_tree_preds <- matrix(NA, nrow=nnew, ncol=ntree)
  keep_preds <- matrix(NA, nrow=rounds, ncol=n)
  keep_pi <- matrix(NA, nrow=rounds, ncol=n)
  keep_xsi <- matrix(NA, nrow=rounds, ncol=length(zero_coef))
  for(i in 1:rounds){
      #rf_list <- vector("list", length=ntree)
      if(i == rounds){
        importance_measures <- matrix(NA, dim(x)[2], ntree)
      }
      rf_out <- for(j in 1:ntree){
        #first generate values of Z for each individual
        czi <- stats::rbinom(n, 1, probi)
        cz0 <- which(czi == 0)
        cx <- x[cz0, , drop=F]
        cy <- y[cz0]
        cdat <- data.frame(cy=cy, cx)
        rf_fit <- randomForest::randomForest(cx, cy, mtry=mtry, ntree=1,
                                             nodesize=nodesize)
        if(i == rounds){
          newdat_tree_preds[, j] <- stats::predict(rf_fit, newdata=newx)
          importance_measures[, j] <- randomForest::importance(rf_fit, type=2)
        } else {
          tree_preds[, j] <- stats::predict(rf_fit, newdata=x)
        }
      }
      #print(proc.time() - ptm)
      #we have predictions e.g. estimated values of \lambda_{i}
      #we now update the coefficients xsi
      rf_preds_z0 <- rowSums(tree_preds)/ntree
      mui <- rf_preds_z0
      mui <- round(mui, 12)
      if(sum(is.na(probi)) > 0){browser()}
      if(sum(is.nan(probi)) > 0){browser()}
      if(any(probi < 0 | probi > 1)){browser()}
      model_zero <- suppressWarnings(stats::glm.fit(z_mat_int, probi,
                              family = stats::binomial(link = "logit"),
                              start = zero_coef,
                              control=list(maxit=1000)))
      zero_coef <-  stats::coef(model_zero)
      logit_pi <- z_mat_int%*%zero_coef
      probi <- expit(logit_pi)
      probi <- probi/(probi + (1 - probi)*stats::dpois(0, mui))
      probi[y1] <- 0
      if(sum(is.na(probi)) > 0){browser()}
      if(sum(is.nan(probi)) > 0){browser()}
      if(any(probi < 0 | probi > 1)){browser()}
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
  new_rf_preds_z0 <- rowSums(newdat_tree_preds)/ntree
  new_probi <- 1/(1 + exp(-new_log_struc_zero))
  new_fit <- data.frame(lambda_preds=new_rf_preds_z0,
                        count_preds=new_rf_preds_z0*(1 - new_probi),
                        unconditional_pi=new_probi)

  unconditional_probi <- 1/(1 + exp(-logit_pi))
  data_fit <- data.frame(lambda_preds=rf_preds_z0,
                         count_preds=rf_preds_z0*(1 - unconditional_probi),
                         fitted_pi=probi,
                         unconditional_pi=unconditional_probi)
  importance_measures <- rowMeans(importance_measures)
  #should prob return variable importance measures eventually
  if(!iter_keep_preds){
    keep_preds = NULL
  }
  if(!iter_keep_pi){
    keep_pi = NULL
  }
  out <- list(new_fit=new_fit, data_fit=data_fit,
              importance_measures=importance_measures,
              zero_coef=zero_coef,
              iter_keep_preds=keep_preds,
              iter_keep_pi=keep_pi,
              iter_keep_xsi=keep_xsi)
}
