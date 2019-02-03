#' Produce formula from X and Z
#'
#' The function takes data.frames X and Z, and it produces
#' a formula that could be used by other functions such
#' as pscl::zeroinfl.
#'
#' @param out_name outcome name
#' @param X design matrix for count model, should be a data.frame
#' @param Z  design matrix for zero-inflation model, should be a data.frame
#' @return A formula
#' @export
XZ_to_form <- function(out_name, X, Z){
  X <- as.data.frame(X)
  Z <- as.data.frame(Z)
  x_names <- names(X)
  z_names <- names(Z)
  x_term <- paste(x_names, collapse = "+")
  z_term <- paste(z_names, collapse = "+")
  xz_term <- paste(x_term, z_term, sep=" | ")
  out_form <- paste(out_name, "~ ", xz_term)
  out_form <- stats::as.formula(out_form)
}

#' ZIP model pdf
#'
#' Evaluates zip model at particular value k.
#'
#' @param x Value the pdf is to be evaluated at
#' @param lambda The expected value for the Poisson model
#' @param pi The probability of being a structural zero
#' @return The probability
#' @export
dzip <- function(x, lambda, pi){
   out <- rep(NA, length(x))
   zero_ind <- x == 0
   out[zero_ind] <-  pi + (1-pi)*stats::dpois(x[zero_ind], lambda)
   out[!zero_ind] <- (1 - pi)*stats::dpois(x[!zero_ind], lambda)
   out
}

#' Produce predictions for Poisson model
#'
#' Takes two design matrices and coefficients and produces predicted
#' values.
#'
#' @param beta Count model coefficients
#' @param xi Zero-inflation model coefficients
#' @param X Count model design matrix
#' @param Z Zero-inflation model design matrix
#' @param type Type of prediction to return
#' @param add_intercept Should intercept be added to X and Z
#' @return A vector of predictions
#' @export
get_preds <- function(beta, xi, X, Z, type='count', add_intercept=T){
  n <- dim(X)[1]
  X <- as.matrix(X)
  Z <- as.matrix(Z)
  if(add_intercept){
    X <- cbind(rep(1, n), X)
    Z <- cbind(rep(1, n), Z)
  }
  lambda_preds <- exp(X%*%beta)
  expit <- function(x){exp(x)/(1 + exp(x))}
  pi <- expit(Z%*%xi)
  count_preds <- lambda_preds*(1 - pi)
  if(type == "count")
    out <- count_preds
  if(type == "lambda")
    out <- lambda_preds
  if(type == "pi")
    out <- pi
  out <- as.numeric(out)
}




















