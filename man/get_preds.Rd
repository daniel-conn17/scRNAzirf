% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/count_model_utilities.R
\name{get_preds}
\alias{get_preds}
\title{Produce predictions for Poisson model}
\usage{
get_preds(beta, xi, X, Z, type = "count", add_intercept = T)
}
\arguments{
\item{beta}{Count model coefficients}

\item{xi}{Zero-inflation model coefficients}

\item{X}{Count model design matrix}

\item{Z}{Zero-inflation model design matrix}

\item{type}{Type of prediction to return}

\item{add_intercept}{Should intercept be added to X and Z}
}
\value{
A vector of predictions
}
\description{
Takes two design matrices and coefficients and produces predicted
values.
}
