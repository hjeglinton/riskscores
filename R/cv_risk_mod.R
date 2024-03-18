#' Run Cross-Validation to Tune Lambda0
#'
#' Runs k-fold cross-validation on a grid of \eqn{\lambda_0} values. Records
#'  class accuracy and deviance for each \eqn{\lambda_0}. Returns an object of
#'  class "cv_risk_mod".
#' @inheritParams risk_mod
#' @param nlambda Number of lambda values to try (default: 25).
#' @param lambda_min_ratio Smallest value for lambda, as a fraction of
#'  lambda_max (the smallest value for which all coefficients are zero).
#'  The default depends on the sample size (\eqn{n}) relative to the number of
#'  variables (\eqn{p}). If \eqn{n > p}, the default is 0.0001, close to zero.
#'  If \eqn{n < p}, the default is 0.01.
#' @param lambda0 Optional sequence of lambda values. By default, the function
#'  will derive the lambda0 sequence based on the data (see `lambda_min_ratio`).
#' @param nfolds Number of folds, implied if `foldids` provided (default: 10).
#' @param foldids Optional vector of values between 1 and `nfolds`.
#' @param parallel If `TRUE`, parallel processing (using [foreach]) is implemented
#'    during cross-validation to increase efficiency (default: `FALSE`).
#'    User must first register parallel backend with a function such as
#'    [doParallel::registerDoParallel].
#' @return An object of class "cv_risk_mod" with the following attributes:
#'  \item{results}{Dataframe containing a summary of deviance and accuracy for each
#'    value of `lambda0` (mean and SD). Also includes the number of nonzero
#'    coefficients that are produced by each `lambda0` when fit on the full data.}
#'  \item{lambda_min}{Numeric value indicating the `lambda0` that resulted in the
#'    lowest mean deviance.}
#'  \item{lambda_1se}{Numeric value indicating the largest `lamdba0` that
#'    had a mean deviance within one standard error of `lambda_min`.}
#' @importFrom foreach %dopar%
#' @export
cv_risk_mod <- function(X, y, weights = NULL, beta = NULL, a = -10, b = 10,
                        max_iters = 100, tol= 1e-5, nlambda = 25,
                        lambda_min_ratio = ifelse(nrow(X) < ncol(X), 0.01, 1e-04),
                        lambda0 = NULL, nfolds = 10, foldids = NULL, parallel = FALSE,
                        seed = NULL) {

  # Set seed
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Check that X is a matrix
  if (!is.matrix(X)) stop ("X must be a matrix")

  # Add intercept column
  if (!all(X[,1] == rep(1, nrow(X)))) {
    X <- cbind(rep(1, nrow(X)), X)
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

  # Check no numeric issues
  if (length(y) != nrow(X)) stop("y and X non-compatible")

  # Check that outcome is 0/1
  if (!all(sort(unique(y)) == c(0,1))) stop ("y must contain values of 0 and 1")

  # Get lambda sequence
  if (is.null(lambda0)){
    sd_n <- function(y) sqrt(sum((y-mean(y))^2)/length(y))

    X_scaled <- scale(X[,-1], scale=apply(X[,-1], 2, sd_n))
    X_scaled <- as.matrix(X_scaled, ncol = ncol(X[,-1]), nrow = nrow(X[,-1]))
    y_weighted <- ifelse(y==0, -mean(y == 1), mean(y == 0))

    lambda_max <- max(abs(colSums(X_scaled*y_weighted)), na.rm = TRUE)/length(y_weighted)
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

    mod <- risk_mod(X_train, y_train, gamma = NULL, beta = beta,
                    weights = weight_train, lambda0 = l0, a = a, b = b,
                    max_iters = max_iters, tol= 1e-5)
    res <- get_metrics_internal(mod, X[foldids == foldid,], y[foldids == foldid])
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
    res_df[,3:5] <- base::t(sapply(1:nrow(res_df), function(i) res_df[i,3:5] <- outlist[[i]]))
  } else {

    res_df[,3:5] <- base::t(sapply(1:nrow(res_df),
                             function(i) fold_fcn(res_df$lambda0[i],
                                                  res_df$fold[i])))
  }


  # Summarize
  dev = acc = NULL # set global variables
  res_df_summary <- res_df %>%
    dplyr::group_by(lambda0) %>%
    dplyr::summarize(mean_dev = mean(dev), sd_dev = stats::sd(dev),
              mean_acc = mean(acc), sd_acc = stats::sd(acc))

  # Find number of nonzero coefficients when fit on full data
  full_fcn <- function(l0) {
    mod <- risk_mod(X, y,  gamma = NULL, beta = NULL,
                    weights = weights, lambda0 = l0, a = a, b = b,
                    max_iters = max_iters, tol= 1e-5)
    non_zeros <- sum(mod$beta[-1] != 0)
    return(c(non_zeros))
  }

  res_df_summary$nonzero <- sapply(1:nrow(res_df_summary),
                           function(i) full_fcn(res_df_summary$lambda0[i]))

  # Find lambda_min and lambda1_se for deviance
  lambda_min_ind <- which.min(res_df_summary$mean_dev)
  lambda_min <- res_df_summary$lambda0[lambda_min_ind]
  min_dev_1se <- res_df_summary$mean_dev[lambda_min_ind] +
    res_df_summary$sd_dev[lambda_min_ind]
  lambda_1se <- res_df_summary$lambda0[max(which(res_df_summary$mean_dev <= min_dev_1se))]

  cv_obj <- list(results = res_df_summary, lambda_min = lambda_min,
                 lambda_1se =lambda_1se)
  class(cv_obj) <- "cv_risk_mod"
  return(cv_obj)
}
