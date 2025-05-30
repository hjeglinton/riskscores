#' Get Model Metrics
#'
#' Calculates a risk model's accuracy, sensitivity, and specificity
#' given a set of data.
#' @param mod An object of class `risk_mod`, usually a result of a call to
#'  [risk_mod()].
#' @param threshold Numeric vector of classification threshold values used to
#'  calculate the accuracy, sensitivity, and specificity of the model. Defaults
#'  to a range of risk probability thresholds from 0.1 to 0.9 by 0.1.
#' @param threshold_type Defines whether the `threshold` vector contains
#'  risk probability values ("response") or threshold values expressed as scores
#'  from the risk score model ("score"). Default: "response".
#' @inheritParams risk_mod
#' @return Data frame with accuracy, sensitivity, and specificity for each threshold.
#' @examples
#' y <- breastcancer[[1]]
#' X <- as.matrix(breastcancer[,2:ncol(breastcancer)])
#'
#' mod <- risk_mod(X, y)
#' get_metrics(mod, X, y)
#'
#' get_metrics(mod, X, y, threshold = c(150, 175, 200), threshold_type = "score")
#' @export
get_metrics <- function(mod, X = NULL, y = NULL, weights = NULL,
                        threshold = NULL, threshold_type = c("response", "score")) {

  threshold_type <- match.arg(threshold_type)

  if (is.null(threshold)) {
    threshold <- seq(0.1, 0.9, 0.1)
    threshold_type <- "response"
  }

  metrics <- data.frame(threshold_risk = threshold,
                        threshold_score = threshold,
                        accuracy = rep(NA, length(threshold)),
                        sensitivity = rep(NA, length(threshold)),
                        specificity = rep(NA, length(threshold)))

  for (i in 1:length(threshold)) {

    res <- get_metrics_internal(mod, X, y, weights, threshold[i], threshold_type)
    metrics[i,] <- c(threshold[i], threshold[i], res$acc, res$sens, res$spec)

  }

  if (threshold_type == "score") {
    metrics$threshold_risk <- round(get_risk(mod, metrics$threshold_score),3)
  } else if (threshold_type == "response") {
    metrics$threshold_score <- round(get_score(mod, metrics$threshold_risk),1)
  }

  return(metrics)

}

#' Calculate Risk Probability from Score
#'
#' Returns the risk probabilities for the provided score value(s).
#' @param object An object of class "risk_mod", usually a result of a call to
#'  [risk_mod()].
#' @param score Numeric vector with score value(s).
#' @return Numeric vector with the same length as `score`.
#' @examples
#' y <- breastcancer[[1]]
#' X <- as.matrix(breastcancer[,2:ncol(breastcancer)])
#'
#' mod <- risk_mod(X, y)
#' get_risk(mod, score = c(1, 10, 20))
#'
#' @export
get_risk <- function(object, score) {

  # Check that object is "risk_mod"
  if (!inherits(object, "risk_mod"))
    stop("'object' must be of class 'risk_mod'")

  risk <- exp(object$gamma*(object$beta[[1]] + score))/
    (1+exp(object$gamma*(object$beta[[1]] + score)))

  return(risk)

}

#' Calculate Score from Risk Probability
#'
#' Returns the score(s) for the provided risk probabilities.
#' @param object An object of class "risk_mod", usually a result of a call to
#'  [risk_mod()].
#' @param risk Numeric vector with probability value(s).
#' @return Numeric vector with the same length as `risk`.
#' @examples
#' y <- breastcancer[[1]]
#' X <- as.matrix(breastcancer[,2:ncol(breastcancer)])
#'
#' mod <- risk_mod(X, y)
#' get_score(mod, risk = c(0.25, 0.50, 0.75))
#'
#' @export
get_score <- function(object, risk) {

  # Check that object is "risk_mod"
  if (!inherits(object, "risk_mod"))
    stop("'object' must be of class 'risk_mod'")

  # Check that risk is between 0 and 1
  if (any(risk <=0) | any(risk >= 1))
    stop("'risk' must contain values between 0 and 1")

  logit_p <- log(risk/(1-risk))

  score <- (1/object$gamma) * logit_p - object$beta[[1]]

  return(score)

}

#' Generate Stratified Fold IDs
#'
#' Returns a vector of fold IDs that preserves class proportions.
#' @inheritParams cv_risk_mod
#' @param nfolds Number of folds (default: 10).
#' @return Numeric vector with the same length as `y`.
#' @examples
#' y <- rbinom(100, 1, 0.3)
#' foldids <- stratify_folds(y, nfolds = 5)
#' table(y, foldids)
#' @export
stratify_folds <- function(y, nfolds = 10, seed = NULL) {

  # Set seed
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Check at least 3 folds
  if (nfolds <= 3) stop("Must have more than 3 folds")

  # Check that y is binomial
  if (length(unique(y)) != 2) stop("'y' must be binomial")

  # Create folds
  y_vals <- unique(y)

  index_y0 <- which(y == y_vals[1])
  index_y1 <- which(y == y_vals[2])
  folds_y0 <- sample(rep(seq(nfolds), length = length(index_y0)))
  folds_y1 <- sample(rep(seq(nfolds), length = length(index_y1)))

  foldids <- rep(NA, length(y))
  foldids[index_y0] <- folds_y0
  foldids[index_y1] <- folds_y1

  return(foldids)

}

#' Run risk model with random start
#'
#' Runs `nstart` iterations of `risk_mod()`, each with a different
#' warm start, and selects the best model. Each coefficient start is
#' randomly selected as -1, 0, or 1.
#' @inheritParams risk_mod
#' @param nstart Number of different random starts to try
#' (default: 5).
#' @export
risk_mod_random_start <- function(X, y, weights = NULL,
                                  lambda0 = 0, a = -10, b = 10,
                                  max_iters = 100, tol= 1e-5,
                                  seed = NULL, nstart = 5) {
  # Set seed
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Initialize best model and objective function value
  best_mod <- NULL
  best_objective <- 99999999

  for (i in 1:nstart) {

    # Randomly choose coefficient values
    beta <- sample(c(-1,0,1), ncol(X), replace = TRUE)

    # Run risk_mod
    mod <- risk_mod(X, y, gamma = NULL, beta = beta, weights = weights,
                    lambda0 = lambda0, a = a, b = b, max_iters = max_iters,
                    tol = tol)

    # Calculate objective function
    objective <- obj_fcn(X, y, gamma = mod$gamma, beta = mod$beta,
                         weights = mod$weights, lambda0 = mod$lambda0)

    # Save model if lowest objective value
    if (objective < best_objective) {
      best_mod <- mod
      best_dev <- objective
    }
  }
  return(best_mod)
}

