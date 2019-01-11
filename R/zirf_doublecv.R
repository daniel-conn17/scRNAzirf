#' Fit zero-inflated random forests using double cross-validation.
#' In other words, nested cross-validation is used to estimate the
#' mse of a cross-validated zero-inflated random forest model.
#'
#' @param x  Matrix or data.frame of covariates for the count part
#'             of the model.
#' @param z  Matrix or data.frame of covariates for the zero-inflation
#'            part of the model.
#' @param y  vector of count outcomes.
#' @param rounds  how many iterations the EM algorithm should be run.
#' @param nlambda  number folds used in estimating initial zip model
#' @param ntree  number of trees
#' @param count_coef values for initial estimate of count model.
#' @param zero_coef values for initial estimate of zero-inflation model.
#' @param nfolds  number of folds used in outer cv
#' @param n_innerfolds number of folds used in inner cv rounds
#' @param tuning_mat matrix with tuning parameters to carry out cross-validation over
#' @param OuterFolds Index of which observation is in which outer cv fold.
#' This is important because it easily allows us to compare predictions from
#' multiple runs zero-inflated random forests or regular random forests.
#' @param iter_keep_pi boolean, if TRUE values of pi are kept across EM iterations
#' @param iter_keep_preds boolean, if TRUE prediction values are kept across EM iterations
#' @param iter_keep_xsi boolean, if TRUE zero-inflation model coefficients
#'                      are kept across EM iterations
#' @return A list
#' @export
zirf_doublecv <- function(x, z, y,
                          rounds=rounds,
                          nlambda=20,
                          ntree=100,
                          count_coef=NULL,
                          zero_coef=NULL,
                          nfolds=2,
                          n_innerfolds=4,
                          tuning_mat=NULL,
                          OuterFolds=NULL,
                          iter_keep_preds=NULL,
                          iter_keep_pi=NULL,
                          iter_keep_xsi=NULL){
  n <- dim(x)[1]
  p <- dim(x)[2]
  if(is.null(tuning_mat)){
    tuning_mat <- matrix(c(5, 50, 100, rep(ceiling(p/3), 3)), 3, 2)
    colnames(tuning_mat) <- c("nodesize", "mtry")
  }
  num_tuning_params <- dim(tuning_mat)[1]
  if(is.null(OuterFolds)){
    folds <- cvTools::cvFolds(n=n, K=nfolds)
    folds <- as.data.frame(cbind(folds$which, folds$subsets))
    names(folds) <- c("Fold", "Index")
  }
  if(!is.null(OuterFolds)){
    folds <- OuterFolds
  }
  fit_list <- vector("list", length=nfolds)
  rf_fit_list <- vector("list", length=nfolds)
  obs_preds <- vector("list", length=nfolds)
  inner_nfold_mses <- vector("list", length=nfolds)
  outer_nfold_mses <- matrix(NA, nfolds, dim(tuning_mat)[2] + 1)
  outer_zirf_preds <- vector("list", length=nfolds)
  zero_coef_mat <- matrix(NA, nfolds, dim(z)[2] + 1)
  colnames(outer_nfold_mses) <- c("nodesize", "mtry", "mse")
  for(i in 1:nfolds){
    ctrain_ind <- folds$Index[folds$Fold != i]
    ctest_ind <- folds$Index[folds$Fold == i]
    cx <- x[ctrain_ind, , drop=F]
    cz <- z[ctrain_ind, , drop=F]
    cy <- y[ctrain_ind]
    newx <- x[ctest_ind, , drop=F]
    newz <- z[ctest_ind, , drop=F]
    newy <- y[ctest_ind]
    #carry out inner round of cross-validation for nodesize
    #matrix with rows equal to nodesize,mtry,mse
    innerfolds <- cvTools::cvFolds(n=dim(cx)[1], K=n_innerfolds)
    innerfolds <- as.data.frame(cbind(innerfolds$which,
                                      innerfolds$subsets))
    names(innerfolds) <- c("Fold", "Index")
    for(j in 1:n_innerfolds){
      #set up data set in inner-fold
      inner_ctrain_ind <- innerfolds$Index[innerfolds$Fold != j]
      inner_ctest_ind <- innerfolds$Index[innerfolds$Fold == j]
      inner_cx <- cx[inner_ctrain_ind, , drop=F]
      inner_cz <- cz[inner_ctrain_ind, , drop=F]
      inner_cy <- cy[inner_ctrain_ind]
      inner_newx <- cx[inner_ctest_ind, , drop=F]
      inner_newz <- cz[inner_ctest_ind, , drop=F]
      inner_newy <- cy[inner_ctest_ind]
      #inner_nodesize <- tuning_mat[j, 1]
      #inner_mtry <- tuning_mat[j, 2]
      tuning_mat_mse <- matrix(NA, dim(tuning_mat)[1], dim(tuning_mat)[2]+1)
      colnames(tuning_mat_mse) <- c("nodesize", "mtry", "mse")
      for(k in 1:dim(tuning_mat)[1]){
          cnodesize <- tuning_mat[k, 1]
          cmtry <- tuning_mat[k, 2]
          inner_fit <- zirf_fit(inner_cx, inner_cz, inner_cy, rounds=rounds,
                                mtry=cmtry,
                                nlambda=nlambda,
                                ntree=ntree, nodesize=cnodesize,
                                count_coef=count_coef, zero_coef=zero_coef,
                                newx=inner_newx, newz=inner_newz,
                                iter_keep_preds=iter_keep_preds,
                                iter_keep_pi=iter_keep_pi,
                                iter_keep_xsi=iter_keep_xsi)
          inner_fitnew <- inner_fit$new_fit
          inner_preds <- inner_fitnew$count_preds
          tuning_mat_mse[k, ]<- c(tuning_mat[k, ],
                               mean((inner_preds - inner_newy)^2))
      }
      inner_nfold_mses[[i]][[j]] <- tuning_mat_mse
    }
    #fit best zirf model to the data set
    cbest_setting <- tuning_mat_mse[which.min(tuning_mat_mse[, 3]), ]
    cnodesize_star <- cbest_setting[1]
    cmtry_star <- cbest_setting[2]
    outer_zirf_fit <- zirf_fit(cx, cz, cy, rounds,
                               mtry=cmtry_star, nlambda=20, ntree=ntree,
                               nodesize=cnodesize_star,
                               newx=newx, newz=newz,
                               count_coef=count_coef,
                               zero_coef=zero_coef,
                               iter_keep_preds=iter_keep_preds,
                               iter_keep_pi=iter_keep_pi,
                               iter_keep_xsi=iter_keep_xsi)
    zero_coef_mat[i, ] <- outer_zirf_fit$zero_coef
    outer_zirf_preds[[i]] <- outer_zirf_fit
    outer_nfold_mses[i, 3] <- mean((outer_zirf_fit$new_fit$count_preds - newy)^2)
    outer_nfold_mses[i, 1:2] <- c(cnodesize_star, cmtry_star)
  }
  tuning_params <- list(tuning_mat=tuning_mat, ntree=ntree, rounds=rounds,
                        count_coef=count_coef, zero_coef=zero_coef,
                        nfolds=nfolds, n_innerfolds=n_innerfolds)

  fold_names <- paste("OuterFold", 1:nfolds, sep="")
  names(outer_zirf_preds) <- fold_names
  out <- list(inner=inner_nfold_mses, outer=outer_nfold_mses,
              tuning_params=tuning_params, preds_vims=outer_zirf_preds,
              folds=folds, zero_coef_mat=zero_coef_mat)
}
