#' Fit random forests in multi-outcome given a pre-specified set of
#' tuning parameters.
#'
#' @param x  matrix or data.frame of covariates for the count part
#'             of the model.
#' @param y  matrix of outcome values.  each row should correspond to a sample in the data set.
#' @param mtry  which value of mtry to use in each forest.
#' @param ntree  number of trees
#' @param min_leaf  controls depth of regression trees
#' @param perm_frac  controls percentage of observations used for fitting model and estimatings vims
#' @export
mvrf_vim <- function(x, y, mtry=dim(x)[2],
                     ntree=100, min_leaf=5,
                     perm_frac=.75){
  rf_fit <- MultivariateRandomForest::build_forest_predict(x, y, ntree, mtry,
                                                           min_leaf, )
  browser()
}
