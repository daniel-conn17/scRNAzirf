#' Fit zero-inflated random forests using double cross-validation.
#' In other words, nested cross-validation is used to estimate the
#' mse of a cross-validated zero-inflated random forest model.
#'
#' @param x  Matrix or data.frame of covariates for the count part
#'             of the model.
#' @param z  Matrix or data.frame of covariates for the zero-inflation
#'            part of the model.
#' @param y  vector of count outcomes.
#' @param ntree  number of trees
#' @param nfolds  number of folds used in outer cv
#' @param n_innerfolds number of folds used in inner cv rounds
#' @param tuning_mat matrix with tuning parameters to carry out cross-validation over
#' @param OuterFolds Index of which observation is in which outer cv fold.
#' This is important because it easily allows us to compare predictions from
#' multiple runs zero-inflated random forests or regular random forests.
#' @return A list
#' @export
rf_doublecv <- function(x, z, y,
                        ntree=100,
                        nfolds=2, n_innerfolds=4,
                        tuning_mat=NULL, OuterFolds=NULL){
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
  rf_preds_vims  <- vector("list", length=nfolds)
  inner_nfold_mses <- vector("list", length=nfolds)
  outer_nfold_mses <- matrix(NA, nfolds, dim(tuning_mat)[2] + 1)
  outer_preds_vims <- vector("list", length=nfolds)
  colnames(outer_nfold_mses) <- c("nodesize", "mtry", "mse")
  for(i in 1:nfolds){
    ctrain_ind <- folds$Index[folds$Fold != i]
    ctest_ind <- folds$Index[folds$Fold == i]
    cx <- x[ctrain_ind, , drop=F]
    cz <- z[ctrain_ind, , drop=F]
    cxz <- cbind(cx, cz)
    cy <- y[ctrain_ind]
    newx <- x[ctest_ind, , drop=F]
    newz <- z[ctest_ind, , drop=F]
    newxz <- cbind(newx, newz)
    newy <- y[ctest_ind]
    #carry out inner round of cross-validation for nodesize
    #matrix with rows equal to nodesize,mtry,mse
    innerfolds <- cvTools::cvFolds(n=dim(cxz)[1], K=n_innerfolds)
    innerfolds <- as.data.frame(cbind(innerfolds$which,
                                      innerfolds$subsets))
    names(innerfolds) <- c("Fold", "Index")
    for(j in 1:n_innerfolds){
      #set up data set in inner-fold
      inner_ctrain_ind <- innerfolds$Index[innerfolds$Fold != j]
      inner_ctest_ind <- innerfolds$Index[innerfolds$Fold == j]
      inner_cx <- cx[inner_ctrain_ind, , drop=F]
      inner_cz <- cz[inner_ctrain_ind, , drop=F]
      inner_cxz <- cbind(inner_cx, inner_cz)
      inner_cy <- cy[inner_ctrain_ind]
      inner_newx <- cx[inner_ctest_ind, , drop=F]
      inner_newz <- cz[inner_ctest_ind, , drop=F]
      inner_newxz <- cbind(inner_newx, inner_newz)
      inner_newy <- cy[inner_ctest_ind]
      #inner_nodesize <- tuning_mat[j, 1]
      #inner_mtry <- tuning_mat[j, 2]
      tuning_mat_mse <- matrix(NA, dim(tuning_mat)[1], dim(tuning_mat)[2]+1)
      colnames(tuning_mat_mse) <- c("nodesize", "mtry", "mse")
      for(k in 1:dim(tuning_mat)[1]){
          cnodesize <- tuning_mat[k, 1]
          cmtry <- tuning_mat[k, 2]
          inner_fit <- randomForest::randomForest(inner_cxz, inner_cy,
                                mtry=cmtry,
                                ntree=ntree, nodesize=cnodesize)
          inner_preds <- stats::predict(inner_fit, newdata=inner_newxz)
          tuning_mat_mse[k, ]<- c(tuning_mat[k, ],
                               mean((inner_preds - inner_newy)^2))
      }
      inner_nfold_mses[[i]][[j]] <- tuning_mat_mse
    }
    #fit best zirf model to the data set
    cbest_setting <- tuning_mat_mse[which.min(tuning_mat_mse[, 3]), ]
    cnodesize_star <- cbest_setting[1]
    cmtry_star <- cbest_setting[2]
    outer_rf_fit <- randomForest::randomForest(cxz, cy,
                                 mtry=cmtry_star, ntree=100,
                                 nodesize=cnodesize_star)
    outer_vims <- randomForest::importance(outer_rf_fit, type=2)
    outer_preds <- stats::predict(outer_rf_fit, newdata=newxz)
    outer_preds_vims <- list(preds=outer_preds, vims=outer_vims)
    rf_preds_vims[[i]] <- outer_preds_vims
    outer_nfold_mses[i, 3] <- mean((outer_preds - newy)^2)
    outer_nfold_mses[i, 1:2] <- c(cnodesize_star, cmtry_star)
  }
  tuning_params <- list(tuning_mat=tuning_mat, ntree=ntree,
                        nfolds=nfolds, n_innerfolds=n_innerfolds)
  fold_names <- paste("OuterFold", 1:nfolds, sep="")
  names(rf_preds_vims) <- fold_names
  out <- list(inner=inner_nfold_mses, outer=outer_nfold_mses,
              tuning_params=tuning_params,
              preds_vims=rf_preds_vims, folds=folds)
}
