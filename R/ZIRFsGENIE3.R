#' Fits ZIRFs for all regulators using regulators and covariates.
#'
#' @param exprMatrix  expression matrix (genes x samples). It must be a matrix.
#' @param countMatrix  count matrix (genes x samples). It must be a matrix.
#' @param regulators  a vector of rownames representing putative regulators.
#'                    ZIRFs will be run with all of these regulators as covariates.
#' @param targets     a vector of rownames representing subset of genes to be used as targets in GENIE3.
#'                    These genes serve as outcomes each regression model.
#' @param z           matrix holding the covariates for the zero-inflation model.
#' @param rounds      how many iterations the EM algorithm should be run.
#' @param mtry        value of mtry used in random forest.
#' @param ntree       number of trees.
#' @param nodesize    controls depth of regression trees.
#' @param nlambda     number of lambda values to use when fitting initial ZIP-LASSO model for model
#'                    initialization.
#' @param nCores      number of cores to use for parallelization.
#' @return a list of vims to be used in SCENIC
#' @export
zirf_genie3 <- function(exprMatrix, countMatrix, regulators = NULL, targets = NULL, z,
                        rounds = 10, mtry = 5, ntree=100, nodesize=5, nlambda=10,
                        nCores = 1){
  if(is.null(regulators))
    regulators <- rownames(exprMatrix)
  if(is.null(targets))
    targets <- rownames(exprMatrix)
  importance_list <- vector("list", length(targets))
  zirf_gene_fit <- function(gene) {
    y <- as.numeric(countMatrix[gene, ,drop = F])
    if( (gene %in% regulators)  )
      regulators <- setdiff(regulators, gene)
    x <- t(exprMatrix[regulators, ,drop = F])
    gene_fit <- zirf_fit(x, z, y, rounds = rounds, mtry = mtry, ntree = ntree,
                         nodesize = nodesize, nlambda = nlambda)
    vims <- gene_fit$importance
  }

  cl <- parallel::makeCluster(as.integer(nCores))
  parallel::clusterExport(cl,
                          varlist = c("z", "rounds", "mtry", "exprMatrix", "countMatrix",
                                      "ntree", "nodesize", "nlambda"),
                          envir = environment())
  vims <- parallel::parLapply(cl, targets, zirf_gene_fit)
  parallel::stopCluster(cl)
  names(vims) <- targets

  reg_and_targ <- intersect(regulators, targets)
  for(i in reg_and_targ) {
    vims[[i]] <- c(0, vims[[i]])
    names(vims[[i]])[1] <- i
  }
  for(i in 1:length(vims)) {
    vims[[i]] <- as.matrix(vims[[i]])
    colnames(vims[[i]]) <- targets[i]
  }
  vims
}
