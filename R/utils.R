#' Clip Values
#'
#' Clip values prior to exponentiation to avoid numeric errors.
#' @param x Numeric vector.
#' @return Input vector `x` with all values between -709.78 and 709.78.
#' @examples
#' clip_exp_vals(710)
#' @export
clip_exp_vals <- function(x){

  return(pmax(pmin(x, 709.78), -709.78))

}


#' Define Method for Converting to "glm" Object
#'
#' Creates the "glm_fit_risk" method that will be used to convert the risk score
#'  model to a "glm" object. Modified from [glm.fit()] code.
#' @param x Input covariate matrix with dimension \eqn{n \times p};
#'  every row is an observation.
#' @param y Numeric vector for the (binomial) response variable.
#' @param weights An optional vector of ‘prior weights’ to be used in the fitting
#' process. Should be NULL or a numeric vector.
#' @param start Starting values for the parameters in the linear predictor.
#' @param etastart Starting values for the linear predictor.
#' @param mustart Starting values for the vector of means.
#' @param offset This can be used to specify an a priori known component to be
#' included in the linear predictor during fitting. This should be NULL or a
#' numeric vector of length equal to the number of cases. One or more \code{offset}
#' terms can be included in the formula instead or as well, and if more than one
#' is specified their sum is used.
#' @param family A description of the error distribution and link function to be
#' used in the model. For glm this can be a character string naming a family
#' function, a family function or the result of a call to a family function. For
#' glm.fit only the third option is supported.
#' @param control A list of parameters for controlling the fitting process.
#' @param intercept Logical. Should an intercept be included in the null model?
#' @param singular.ok Logical; if FALSE a singular fit is an error.
#' @return Object of class "glm".
#' @noRd
glm_fit_risk <- function (x, y, weights = rep(1, nobs), start = NULL,
                          etastart = NULL, mustart = NULL, offset = rep(0, nobs),
                          family = stats::gaussian(), control = list(),
                          intercept = TRUE, singular.ok=TRUE){

  # dimensions and names
  x <- as.matrix(x)
  xnames <- dimnames(x)[[2L]]
  ynames <- if(is.matrix(y)) rownames(y) else names(y)
  nobs <- NROW(y)
  nvars <- ncol(x)
  EMPTY <- nvars == 0

  ## define weights and offset if needed
  if (is.null(weights))
    weights <- rep.int(1, nobs)
  if (is.null(offset))
    offset <- rep.int(0, nobs)

  ## get family functions:
  variance <- family$variance
  linkinv  <- family$linkinv
  dev.resids <- family$dev.resids
  aic <- family$aic
  mu.eta <- family$mu.eta
  unless.null <- function(x, if.null) if(is.null(x)) if.null else x
  valideta <- unless.null(family$valideta, function(eta) TRUE)
  validmu  <- unless.null(family$validmu,  function(mu) TRUE)
  eval(family$initialize)

  # starting coefficients, eta, and mu
  coef <- start
  eta <- offset + as.vector(if (NCOL(x) == 1L) x * start else x %*% start)
  mu <- linkinv(eta)
  dev <- sum(dev.resids(y, mu, weights))
  boundary <- FALSE
  conv <- TRUE

  varmu <- variance(mu)
  mu.eta.val <- mu.eta(eta)
  good <- (weights > 0) & (mu.eta.val != 0)
  z <- (eta - offset)[good] + (y - mu)[good]/mu.eta.val[good]
  w <- sqrt((weights[good] * mu.eta.val[good]^2)/variance(mu)[good])

  residuals <-  (y - mu)/mu.eta(eta)
  names(coef) <- xnames

  names(residuals) <- ynames
  names(mu) <- ynames
  names(eta) <- ynames

  # for compatibility with lm, which has a full-length weights vector
  wt <- rep.int(0, nobs)
  wt[good] <- w^2
  names(wt) <- ynames
  names(weights) <- ynames
  names(y) <- ynames

  ## calculate null deviance -- corrected in glm() if offset and intercept
  wtdmu <-
    if (intercept) sum(weights * y)/sum(weights) else linkinv(offset)
  nulldev <- sum(dev.resids(y, wtdmu, weights))

  ## calculate df
  n.ok <- nobs - sum(weights==0)
  nulldf <- n.ok - as.integer(intercept)
  rank <- dim(x)[2]
  resdf  <- n.ok - rank

  ## calculate AIC
  aic.model <- aic(y, n.ok, mu, weights, dev) + 2*rank

  # return list
  list(coefficients = coef, residuals = residuals, fitted.values = mu,
       effects = c("(Intercept)", xnames),
       rank = rank,
       family = family,
       linear.predictors = eta, deviance = dev,
       aic = aic.model,
       null.deviance = nulldev, iter = 0, weights = wt,
       prior.weights = weights,
       df.residual = resdf,
       df.null = nulldf,
       y = y, converged = conv, boundary = boundary)
}


#' Get Model Metrics for a Single Threshold
#'
#' Calculates a risk model's deviance, accuracy, sensitivity, and specificity
#' given a set of data and a threshold value.
#' @param mod An object of class `risk_mod`, usually a result of a call to
#'  [risk_mod()].
#' @inheritParams get_metrics
#' @return List with deviance (dev), accuracy (acc), sensitivity (sens),
#'  specificity (spec), and auc.
get_metrics_internal <- function(mod, X = NULL, y = NULL, weights = NULL,
                                 threshold = 0.50, threshold_type = c("response", "score")){
  
  threshold_type <- match.arg(threshold_type)
  
  # Check threshold value against type
  if (!is.null(threshold) & threshold_type == "response") {
    if (threshold < 0 | threshold > 1) stop("threshold must be between 0 and 1 when threshold_type = 'response'")
  } else if (!is.null(threshold) & threshold_type == "score") {
    if(threshold > 0 & threshold < 1) warning("Threshold input is being interpreted as a score but it may be a probability. Use `threshold_type = 'response'` to interpret input as a probability.")
  }
  
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
  
  # Define threshold
  if (threshold_type == "response") {
    prob_cutoff <- threshold
  } else if (threshold_type == "score") {
    prob_cutoff <- get_risk(mod, threshold)
  }
  
  # Get predicted probs and classes
  v <- mod$gamma * X %*% mod$beta
  v <- clip_exp_vals(v)
  p <- exp(v)/(1+exp(v))
  pred <- ifelse(p >= prob_cutoff, 1, 0)
  
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
  
  # AUC values
  roc <-roc(y, predict(mod, X, type = "response")[,1], quiet = TRUE)
  return(list(dev = dev, acc=acc, sens=sens, spec=spec, auc=roc$auc))
}

#' Randomly round the initialized coefficients before coordinate descent
#'
#' Round each LR coefficient based on its decimal value. The decimal is the probability
#' of rounding the coefficient up to the next integer
#' @param beta Numeric vector or  logistic regression coefficients initialized
#' before cyclical coordinate descent in `risk_mod()`. The first element is the
#' intercept and is not modified.
#' @return A numeric vector with randomized rounding (apart from the first element).
randomized_rounding <- function(beta) {

  # Extract the decimals of the coefficients, excluding the intercept
  beta_dec <- beta[-1] %% 1
  # Binary rounding outcome
  decision <- rbinom(n = length(beta_dec), size = 1, prob = beta_dec)

  # Apply the randomized rounding for columns 2 to n
  beta[-1] <- floor(beta[-1]) + decision

  return(beta)
}

