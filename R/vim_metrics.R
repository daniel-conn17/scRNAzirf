#' Fit zero-inflated random forests given a pre-specified set of
#' tuning parameters.
#'
#' @param coefs Matrix or data.frame of covariates for the count part
#'              of the model.
#' @param vims  Matrix or data.frame of covariates for the zero-inflation
#'              part of the model.
#' @param type  Precision-recall curve or ROC curve
#' @return A list.
#' @export
calc_auc <- function(coefs, vims, type){
  sig_ind <- which(abs(coefs) > 0)
  zero_ind <- which(abs(coefs) == 0)
  if(type == "pr"){
    out <- PRROC::pr.curve(scores.class0=vims[sig_ind],
                    scores.class1=vims[zero_ind],
                    curve=T)
  }
  if(type == "roc"){
    out <- PRROC::roc.curve(scores.class0=vims[sig_ind],
                     scores.class1=vims[zero_ind],
                     curve=T)
  }
  out
}
