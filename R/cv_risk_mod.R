#' Cross-Validation to set lambda0
#'
#' Runs k-fold cross-validation and records class accuracy and deviance
#' @param X input matrix with dimension n x p, every row is an observation
#' @param y numeric vector for the response variable (binomial)
#' @param weights numeric vector of length n with weights for each
#' observation (defult NULL - will give equal weights)
#' @param a integer lower bound for betas (default -10)
#' @param b integer upper bound for betas (default 10)
#' @param max_iters maximum number of iterations (default 100)
#' @param tol tolerance for convergence
#' @param nlambda number of lambda values to try (default 10)
#' @param lambda_min_ratio smallest value for lambda, as a fraction of
#' lambda_max, the (data derived) entry value (i.e. the smallest value
#' for which all coefficients are zero). The default depends on the sample size
#' (n) relative to the number of variables (p). If n > p, the default is 0.0001,
#' close to zero.  If n < p, the default is 0.01.
#' @param lambda0 optional sequence of lambda values (default NULL)
#' @param nfolds number of folds, implied if foldids provided (default 10)
#' @param foldids optional vector of values between 1 and nfolds (default NULL)
#' @param parallel If \code{TRUE}, use parallel \code{foreach} to fit each fold.
#' Must register parallel before hand, such as \code{doParallel} or others.
#' @param seed An integer that is used as argument by the set.seed() for
#' offsetting the random number generator. Default is to leave the random number
#' generator alone.
#' @return class of cv_risk_mod with a list containing a data.frame of results
#' along with the lambda_min and lambda_1se
#' @importFrom foreach %dopar%
#' @export
cv_risk_mod <- function(X, y, weights = NULL, a = -10, b = 10, max_iters = 100,
                        tol= 1e-5, nlambda = 25,
                        lambda_min_ratio = ifelse(nrow(X) < ncol(X), 0.01, 1e-04),
                        lambda0 = NULL, nfolds = 10, foldids = NULL, parallel=F,
                        seed = NULL) {


  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Get folds
  if (is.null(foldids) & is.null(nfolds)) stop("Must provide foldids or nfolds")
  if (is.null(foldids)){
    foldids <- sample(rep(seq(nfolds), length = nrow(X)))
  } else {
    nfolds <- max(foldids)
  }

  # Check at least 3 folds
  if (nfolds <= 3) stop("Must have more than 3 folds")

  # Get lambda sequence
  if (is.null(lambda0)){
    sd_n <- function(y) sqrt(sum((y-mean(y))^2)/length(y))

    X_scaled <- scale(X[,-1], scale=apply(X[,-1], 2, sd_n))
    X_scaled <- as.matrix(X_scaled, ncol = ncol(X[,-1]), nrow = nrow(X[,-1]))
    y_weighted <- ifelse(y==0, -mean(y == 1), mean(y == 0))

    lambda_max <- max(abs(colSums(X_scaled*y_weighted)))/length(y_weighted)
    lambda0 <- exp(seq(log(lambda_max), log(lambda_max * lambda_min_ratio),
                       length.out=nlambda))
  }

  num_lambda0 <- length(lambda0)
  if (num_lambda0 < 2) stop("Need at least two values for lambda0")

  # Results data frame
  res_df <- data.frame(lambda0 = rep(lambda0, nfolds),
                       fold = sort(rep(1:nfolds, num_lambda0)),
                       dev = rep(0, nfolds*num_lambda0),
                       acc = rep(0, nfolds*num_lambda0),
                       non_zeros = rep(0, nfolds*num_lambda0))


  # Function to run for single fold and lambda0
  fold_fcn <- function(l0, foldid){
    X_train <- X[foldids != foldid, ]
    y_train <- y[foldids != foldid]
    weight_train <- weights[foldids != foldid]
    mod <- risk_mod(X_train, y_train, gamma = NULL, beta = NULL,
                    weights = weight_train, lambda0 = l0, a = a, b = b,
                    max_iters = max_iters, tol= 1e-5)
    res <- get_metrics(mod, X[foldids == foldid,], y[foldids == foldid])
    non_zeros <- sum(mod$beta != 0)
    return(c(res$dev, res$acc, non_zeros))
  }

  # Run through all folds
  # Parallel Programming. ! Must register parallel beforehand
  i = NULL # set global variable
  if (parallel) {

    outlist = foreach::foreach(i = 1:nrow(res_df)) %dopar%
      {
        fold_fcn(res_df[i,1],res_df[i,2])
      }
    res_df[,3:5] <- t(sapply(1:nrow(res_df), function(i) res_df[i,3:5] <- outlist[[i]]))
  } else {

    res_df[,3:5] <- t(sapply(1:nrow(res_df),
                             function(i) fold_fcn(res_df$lambda0[i],
                                                  res_df$fold[i])))
  }


  # Summarize
  dev = acc = NULL # set global variables
  res_df <- res_df %>%
    dplyr::group_by(lambda0) %>%
    dplyr::summarize(mean_dev = mean(dev), sd_dev = stats::sd(dev),
              mean_acc = mean(acc), sd_acc = stats::sd(acc))

  # Find number of nonzero coefficients when fit on full data
  full_fcn <- function(l0) {
    mod <- risk_mod(X, y,  gamma = NULL, beta = NULL,
                    weights = weights, lambda0 = l0, a = a, b = b,
                    max_iters = max_iters, tol= 1e-5)
    non_zeros <- sum(mod$beta[-1] != 0)
    return(non_zeros)
  }

  res_df$nonzero <- sapply(1:nrow(res_df),
                           function(i) full_fcn(res_df$lambda0[i]))

  # Find lambda_min and lambda1_se for deviance
  lambda_min_ind <- which.min(res_df$mean_dev)
  lambda_min <- res_df$lambda0[lambda_min_ind]
  min_dev_1se <- res_df$mean_dev[lambda_min_ind] +
    res_df$sd_dev[lambda_min_ind]
  lambda_1se <- res_df$lambda0[max(which(res_df$mean_dev <= min_dev_1se))]

  cv_obj <- list(results = res_df, lambda_min = lambda_min,
                 lambda_1se =lambda_1se)
  class(cv_obj) <- "cv_risk_mod"
  return(cv_obj)
}
