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
#' @useDynLib riskscores, .registration = TRUE
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

  # call Fortran code via C wrapper
  fit <- .Call("Cdqrls", x[good, , drop = FALSE] * w, z * w,
               min(1e-7, control$epsilon/1000), check=FALSE)
  xxnames <- xnames[fit$pivot]

  # r matrix
  residuals <-  (y - mu)/mu.eta(eta)
  fit$qr <- as.matrix(fit$qr)
  nr <- min(sum(good), nvars)
  if (nr < nvars) {
    Rmat <- diag(nvars)
    Rmat[1L:nr, 1L:nvars] <- fit$qr[1L:nr, 1L:nvars]
  }
  else Rmat <- fit$qr[1L:nvars, 1L:nvars]

  Rmat <- as.matrix(Rmat)
  Rmat[row(Rmat) > col(Rmat)] <- 0
  names(coef) <- xnames
  colnames(fit$qr) <- xxnames
  dimnames(Rmat) <- list(xxnames, xxnames)

  names(residuals) <- ynames
  names(mu) <- ynames
  names(eta) <- ynames

  # for compatibility with lm, which has a full-length weights vector
  wt <- rep.int(0, nobs)
  wt[good] <- w^2
  names(wt) <- ynames
  names(weights) <- ynames
  names(y) <- ynames
  if(!EMPTY)
    names(fit$effects) <-
    c(xxnames[seq_len(fit$rank)], rep.int("", sum(good) - fit$rank))

  ## calculate null deviance -- corrected in glm() if offset and intercept
  wtdmu <-
    if (intercept) sum(weights * y)/sum(weights) else linkinv(offset)
  nulldev <- sum(dev.resids(y, wtdmu, weights))

  ## calculate df
  n.ok <- nobs - sum(weights==0)
  nulldf <- n.ok - as.integer(intercept)
  rank <- if(EMPTY) 0 else fit$rank
  resdf  <- n.ok - rank

  ## calculate AIC
  aic.model <- aic(y, n.ok, mu, weights, dev) + 2*rank

  # return list
  list(coefficients = coef, residuals = residuals, fitted.values = mu,
       effects = if(!EMPTY) fit$effects, R = if(!EMPTY) Rmat, rank = rank,
       qr = if(!EMPTY) structure(fit[c("qr", "rank", "qraux", "pivot", "tol")], class = "qr"),
       family = family,
       linear.predictors = eta, deviance = dev, aic = aic.model,
       null.deviance = nulldev, iter = 0, weights = wt,
       prior.weights = weights, df.residual = resdf, df.null = nulldf,
       y = y, converged = conv, boundary = boundary)
}


