#' Get metrics from beta and gamma
#'
#' Calculates deviance, accuracy, sensitivity, and specificity
#' @param mod riskMod object
#' @param X input matrix with dimensions n x p, must match dimensions of beta
#' in mod (default NULL)
#' @param y numeric vector for the response variable (binomial) of length n,
#' (default NULL)
#' @param weights numeric vector of length n with weights for each
#' observation (defult NULL - will give equal weights)
#' @return list with deviance (dev), accuracy (acc), sensitivity (sens), and
#' specificity (spec)
get_metrics <- function(mod, X = NULL, y = NULL, weights = NULL){

  # Check if new data
  if (is.null(X)+is.null(y) == 1) stop("Must provide both X and y")
  if (is.null(X) & is.null(y)){
    X = mod$X
    y = mod$y
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

#' Assign stratified fold ids
#'
#' Returns a vector of fold ids with an equal percentage of samples for each class
#' @param y numeric vector for the response variable (binomial)
#' @param nfolds number of folds (default 10)
#' @param seed An integer that is used as argument by the set.seed() for
#' offsetting the random number generator. Default is to leave the random number
#' generator alone.
#' @return vector with the same length as y
#' @export
stratify_folds <- function(y, nfolds = 10, seed = NULL) {

  if (!is.null(seed)) {
    set.seed(seed)
  }

  index_y0 <- which(y == 0)
  index_y1 <- which(y == 1)
  folds_y0 <- sample(rep(seq(nfolds), length = length(index_y0)))
  folds_y1 <- sample(rep(seq(nfolds), length = length(index_y1)))

  foldids <- rep(NA, length(y))
  foldids[index_y0] <- folds_y0
  foldids[index_y1] <- folds_y1

  return(foldids)
}

#' Get risk from score
#'
#' Returns the risk probabilities for the provided score values
#' @param object an object of class "risk_mod", usually a result of a call to
#' risk_mod()
#' @param score numeric vector with score value(s)
#' @return numeric vector with the same length as `score`
#' @export
get_risk <- function(object, score) {

  risk <- exp(object$gamma*(object$beta[[1]] + score))/
    (1+exp(object$gamma*(object$beta[[1]] + score)))

  return(risk)

}


