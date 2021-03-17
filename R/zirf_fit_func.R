#' Fit zero-inflated random forests given a pre-specified set of
#' tuning parameters.
#'
#' @param x  Matrix or data.frame of covariates for the count part
#'             of the model.
#' @param z  Matrix or data.frame of covariates for the zero-inflation
#'            part of the model.
#' @param y  vector of count outcomes.
#' @param rounds  how many iterations the EM algorithm should be run at maximum.
#' @param mtry  which value of mtry to use in each forest.
#' @param ntree  number of trees.
#' @param nodesize  controls depth of regression trees.
#' @param count_coef values for initial estimate of count model.
#' @param zero_coef values for initial estimate of zero-inflation model.
#' @param newx A list of count model covariate matrices for new data.
#' @param newz A list of zero inflation model covariate matrices for new data.
#' @param iter_keep_pi boolean, if TRUE values of pi are kept across EM iterations.
#' @param iter_keep_preds boolean, if TRUE prediction values are kept across EM iterations.
#' @param iter_keep_xsi boolean, if TRUE zero-inflation model coefficients
#'                      are kept across EM iterations.
#' @param nlambda number of folds of cross-validation.
#' @param threshold_prob boolean, if TRUE.
#' @param epsilon tolerance used for stopping iterations.
#'                 The maximum difference in probability of inclusion.
#'                 between iterations is the criterion.

#' @return A list.
#' @export
zirf_fit <- function(x, z, y, rounds = 25, mtry,
                     ntree=100, nodesize=5,
                     count_coef=NULL, zero_coef=NULL,
                     newx=NULL, newz=NULL,
                     iter_keep_pi=T, iter_keep_preds=T,
                     iter_keep_xsi=T, nlambda=10,
                     threshold_prob = T, epsilon = .01){
  #newx and newz are processed as lists of test sets. this just makes them a list
  if(!identical(class(x),class(z))) {
    stop("class(x) must equal class(z)")
  }
  if(is.data.frame(newx) | is.matrix(newx)){
     newx = list(newx)
  }
  if(is.data.frame(newz) | is.matrix(newz)){
     newz = list(newz)
  }
  #add a character to start of every variable so that we don't run into other issues
  # to get the original names run the following command:
  # substring(names(zilm_dat), 5)
  zero_coef_og <- zero_coef
  x_ind <- 1:dim(x)[2]
  z_ind <- (dim(x)[2] + 1):(dim(x)[2] + dim(z)[2])
  nb_part01 <- paste("y", "~ ", collapse="")
  if(is.data.frame(x)){
    names(x) <- paste0("X", names(x))
    #this is to deal withh the fact that some gene names start with numbers
    #we are concatenating X to all of the gene names so that it is easier to
    #undo this mild transformation of the gene names.
    x_names <- names(x)
    nb_part02 <- paste(x_names, collapse="+")
    z_names <- names(z)
    zero_part <- paste0("|", paste0(z_names, collapse="+"))
    if(!is.null(newx)){
      for(i in 1:length(newx)){
        names(newx[[i]]) <- names(x)
      }
    }
  }
  if(is.matrix(x)){
    colnames(x) <- paste0("X", colnames(x))
    x_names <- colnames(x)
    nb_part02 <- paste(x_names, collapse=" + ")
    z_names <- colnames(z)
    zero_part <- paste0("|", paste0(z_names, collapse=" + "))
    if(!is.null(newx)){
      for(i in 1:length(newx)){
        colnames(newx[[i]]) <- colnames(x)
      }
    }
  }
  mod_string <- paste(nb_part01, nb_part02, zero_part, sep="")
  mod_formula <- stats::as.formula(mod_string)
  zilm_dat <- data.frame(x, z, y=y)
  if(is.null(newx)){
    newx <- as.data.frame(x)
    newz <- as.matrix(z)
    newz <- as.matrix(cbind(rep(1, nrow(newz)), newz))
    newx <- list(newx)
    newz <- list(newz)
  }
  else {
    for(i in 1:length(newx)){
      newx[[i]] <- as.data.frame(newx[[i]])
      newz[[i]] <- as.matrix(newz[[i]])
      newz[[i]] <- as.matrix(cbind(rep(1, nrow(newz[[i]])), newz[[i]]))
    }
  }
  #fit initial zip model, comment out the next line, if desired
  #normally this should be set to 100-3000
  #i've commented this out because the additional source of randomness makes
  #predictions harder to compare
  #quick_zilm_dat <- zilm_dat[sample(1:dim(zilm_dat)[1], dim(zilm_dat)[1]), ]
  quick_zilm_dat <- zilm_dat
  #ptm <- proc.time()
  names(quick_zilm_dat)[dim(quick_zilm_dat)[2]] <- "y"
  if(is.null(count_coef)){
    pois_mod <- tryCatch({mpath::zipath(mod_formula, data = quick_zilm_dat,  nlambda=nlambda,
                          family="poisson",
                          plot.it = F)},
                         error = function(err) {
                           print(err)
                           out <- "error"
                         })
    if(class(pois_mod) == "zipath"){
      best_ind <- which.min(pois_mod$bic)
      pois_coef <- stats::coef(pois_mod)
      zero_coef <- pois_coef$zero[, best_ind]
      count_coef <- pois_coef$count[, best_ind]
    } else {
      mean_y <- mean(y)
      zero_prop <- sum(y == 0)/length(y)
      zero_coef <- c(log(zero_prop/(1 - zero_prop)), rep(0, dim(z)[2]))
      count_coef <- c(exp(mean_y), rep(0, dim(x)[2]))
    }
  }
  zdat <- zilm_dat[, z_ind, drop=F]
  zdat <- as.matrix(cbind(rep(1, nrow(zdat)), zdat))
  xdat <- zilm_dat[, x_ind, drop=F]
  xdat <- as.matrix(cbind(rep(1, nrow(xdat)), xdat))
  xdat_df <- as.data.frame(xdat[, -1, drop=F])

  #Next apply formula (9) of Wang et al (2014) from Stats in Medicine
  #We know z_{i}=0 for these observations
  sure_nonzero_ind <- which(y != 0)
  y1 <- y > 0
  logit_pi <- zdat%*%zero_coef
  log_init_pois <- xdat%*%count_coef
  expit <- function(x){1/(1 + exp(-x))}
  probi <- expit(logit_pi)
  mui <- exp(log_init_pois)
  probi <- probi/(probi + (1 - probi)*stats::dpois(0, mui))
  probi[y1] <- 0
  #edit this now
  n <- dim(x)[1]
  #ptm <- proc.time()
  tree_preds <- matrix(NA, nrow=dim(xdat_df)[1], ncol=ntree)
  newdat_tree_preds <- vector("list", length(newx))
  for(i in 1:length(newx)){
    newdat_tree_preds[[i]] <- matrix(NA, nrow=dim(newx[[i]])[1], ncol=ntree)
  }
  keep_preds <- matrix(NA, nrow=rounds, ncol=dim(xdat_df)[1])
  keep_pi <- matrix(NA, nrow=rounds, ncol=dim(xdat_df)[1])
  keep_xsi <- matrix(NA, nrow=rounds, ncol=length(zero_coef))
  for(i in 1:rounds){
      #rf_list <- vector("list", length=ntree)
      if(i == rounds){
        importance_measures <- matrix(NA, dim(x)[2], ntree)
      }
      for(j in 1:ntree){
        if(!threshold_prob){
          #generate values of Z for each individual
          czi <- stats::rbinom(n, 1, probi)
          cz0 <- czi == 0
        }
        if(threshold_prob){
          #generate values of Z for each individual
          cz0 <- probi < .5
        }
        cx <- x[cz0, , drop=F]
        cy <- y[cz0]
        cdat <- data.frame(cy=cy, cx)
        rf_fit <- suppressWarnings(randomForest::randomForest(cx, cy, mtry=mtry, ntree=1,
                                                              nodesize=nodesize))
        if(i == rounds){
          importance_measures[, j] <- randomForest::importance(rf_fit, type=2)
          for(rr in 1:length(newx)){
            newdat_tree_preds[[rr]][, j] <- stats::predict(rf_fit, newdata=newx[[rr]])
          }
        } else {
          tree_preds[, j] <- stats::predict(rf_fit, newdata=xdat_df)
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
      model_zero <- suppressWarnings(stats::glm.fit(zdat, probi,
                              family = stats::binomial(link = "logit"),
                              start = zero_coef,
                              control=list(maxit=1000)))
      zero_coef <-  stats::coef(model_zero)
      logit_pi <- zdat%*%zero_coef
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
      if(i > 1){

        current_prob_diff <- max(abs(keep_pi[i, ] - keep_pi[i - 1, ]))
        if(!is.logical(current_prob_diff < epsilon)){browser()}
        if(current_prob_diff < epsilon){
          converged = TRUE
        }
      }
  }
  total_iters <- i - 1

  new_fit <- vector("list", length(newx))
  for(i in 1:length(newx)){
    new_log_struc_zero <- newz[[i]]%*%zero_coef
    new_rf_preds_z0 <- rowSums(newdat_tree_preds[[i]])/ntree
    new_probi <- 1/(1 + exp(-new_log_struc_zero))
    new_prob_zero <- 1 - (new_probi + (1 - new_probi)*exp(-new_rf_preds_z0))
    new_fit[[i]] <- data.frame(lambda_preds=new_rf_preds_z0,
                          count_preds=new_rf_preds_z0*(1 - new_probi),
                          unconditional_pi=new_probi,
                          prob_zero = new_prob_zero)
  }
  if(length(new_fit) == 1){
    new_fit <- new_fit[[1]]
  }
  #conditional probi is the probability a structural zero given that y=0
  #unconditiona probi is the probability that we have a zero count given
  #the covariates
  unconditional_probi <- 1/(1 + exp(-logit_pi))
  prob_zero <- 1 - (unconditional_probi + (1 - unconditional_probi)*exp(-rf_preds_z0))
  data_fit <- data.frame(lambda_preds=rf_preds_z0,
                         count_preds=rf_preds_z0*(1 - unconditional_probi),
                         fitted_pi=probi,
                         unconditional_pi=unconditional_probi,
                         prob_zero = prob_zero)
  importance_measures <- rowMeans(importance_measures)
  if(is.data.frame(x)) {
    imp_names <- names(x)
  } else {
    imp_names <- colnames(x)
  }
  names(importance_measures) <- imp_names
  names(importance_measures) <- gsub('^.', '', names(importance_measures))
  if(!iter_keep_preds){
    keep_preds = NULL
  }
  if(!iter_keep_pi){
    keep_pi = NULL
  }
  iter_keep_preds <- keep_preds[1:total_iters, ]
  iter_keep_pi <- keep_pi[1:total_iters, ]
  iter_keep_xsi <- keep_xsi[1:total_iters, ]
  out <- list(new_fit=new_fit, data_fit=data_fit,
              importance_measures=importance_measures,
              zero_coef=zero_coef,
              iter_keep_preds=keep_preds,
              iter_keep_pi=keep_pi,
              iter_keep_xsi=keep_xsi,
              iterations = total_iters)
  class(out) <- "zirf_fit"
  out
}
