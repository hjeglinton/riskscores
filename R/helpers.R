#' Get Model Metrics
#'
#' Calculates a risk model's deviance, accuracy, sensitivity, and specificity
#' given a set of data.
#' @param mod An object of class `risk_mod`, usually a result of a call to
#'  [risk_mod()].
#' @inheritParams risk_mod
#' @return List with deviance (dev), accuracy (acc), sensitivity (sens), and
#'  specificity (spec).
#' @examples
#' y <- breastcancer[[1]]
#' X <- as.matrix(breastcancer[,2:ncol(breastcancer)])
#'
#' mod <- risk_mod(X, y)
#' get_metrics(mod, X, y)
#' @export
get_metrics <- function(mod, X = NULL, y = NULL, weights = NULL){

  # Check if new data
  if (is.null(X)+is.null(y) == 1) stop("Must provide both X and y")
  if (is.null(X) & is.null(y)){
    X = mod$X
    y = mod$y
  }

  # Add intercept column
  if (!all(X[,1] == rep(1, nrow(X)))) {
    X <- cbind(rep(1, nrow(X)), X)
  }

  # Check compatibility
  if (nrow(X) != length(y)) stop("X and y must match in number of observations")
  if (ncol(X) != length(mod$beta)) stop("X is incompatible with the model")
  if (sum(! (y %in% c(0,1)))) stop("y must be 0/1 valued")

  # Get predicted probs and classes
  v <- mod$gamma * X %*% mod$beta
  v <- clip_exp_vals(v)
  p <- exp(v)/(1+exp(v))
  pred <- ifelse(p>=0.5, 1, 0)

  # Deviance
  p[p == 1] <- 0.99999
  p[p == 0] <- 0.00001
  dev <- -2*sum(y*log(p)+(1-y)*log(1-p))

  # Confusion matrix
  tp <- sum(pred == 1 & y == 1)
  tn <- sum(pred == 0 & y == 0)
  fp <- sum(pred == 1 & y == 0)
  fn <- sum(pred == 0 & y == 1)

  # Accuracy values
  acc <- (tp+tn)/(tp+tn+fp+fn)
  sens <- tp/(tp+fn)
  spec <- tn/(tn+fp)

  return(list(dev = dev, acc=acc, sens=sens, spec=spec))
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


