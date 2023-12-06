#' Partial derivative of the negative log-likelihood
#'
#' Calculates the partial derivative for beta_j of the objective function
#' @param X input matrix with dimension n x p, every row is an observation
#' @param y numeric vector for the response variable (binomial)
#' @param gamma scalar to rescale betas for prediction
#' @param beta numeric vector with p coefficients
#' @param weights numeric vector of length n with weights for each
#' observation
#' @param j index of beta
#' @return the numeric partial derivative value
par_deriv <- function(X, y, gamma, beta, weights, j) {

  # Calculate partial derivative for NLL
  pd_1 <- sum(gamma * weights * (y * X[,j]))
  exp_pred <- exp(clip_exp_vals(gamma * (X %*% beta)))
  pd_2 <- sum(gamma*weights * X[,j] * (exp_pred / (1.0 + exp_pred)))
  nll_pd <- (-1/nrow(X))*(pd_1-pd_2)

  return(nll_pd)
}

#' Objective function for NLL+penalty
#'
#' Calculates the objective function for gamma, beta (NLL+penalty)
#' @param X input matrix with dimension n x p, every row is an observation
#' @param y numeric vector for the response variable (binomial)
#' @param gamma scalar to rescale betas for prediction
#' @param beta numeric vector with p coefficients
#' @param weights numeric vector of length n with weights for each
#' observation
#' @param lambda0 penalty coefficient for L0 term (default 0)
#' @return numeric objective function value
obj_fcn <- function(X, y, gamma, beta, weights, lambda0=0) {

  # Calculate partial derivative for NLL
  v <- gamma * (X %*% beta)
  v <- clip_exp_vals(v) # avoids numeric errors
  nll_fcn <- (-1/nrow(X))*sum(weights * (y * v - log(1+exp(v))))

  # Penalty term for lambda0*||beta||_0
  pen_fcn <- lambda0*sum(beta[-1] != 0)
  return (nll_fcn + pen_fcn)
}

#' Bisection search for coordinate descent
#'
#' Returns optimal value on beta_j using bisection search
#' @param X input matrix with dimension n x p, every row is an observation
#' @param y numeric vector for the response variable (binomial)
#' @param gamma scalar to rescale betas for prediction
#' @param beta numeric vector with p coefficients
#' @param weights numeric vector of length n with weights for each
#' observation
#' @param j index of beta
#' @param lambda0 penalty coefficient for L0 term (default 0)
#' @param a integer lower bound for beta_j (default -10)
#' @param b integer upper bound for beta_j (default 10)
#' @return numeric vector beta with optimal value for beta\[j\] updated
bisec_search <- function(X, y, gamma, beta, weights, j, lambda0 = 0,
                         a = -10, b = 10) {

  # Initial betas to compare
  beta_a <- beta
  beta_a[j] <- a
  beta_b <- beta
  beta_b[j] <- b
  beta_0 <- beta
  beta_0[j] <- 0

  # If no zero derivative in range, skip while loop
  der_a <- par_deriv(X, y, gamma, beta_a, weights, j)
  der_b <- par_deriv(X, y, gamma, beta_b, weights, j)
  search <- TRUE
  if (sign(der_a) == sign(der_b)) search <- FALSE

  while (((b - a) > 1) & search){
    # Find partial derivative at midpoint
    c <- floor((a+b)/2)
    beta_c <- beta
    beta_c[j] <- c
    der_c <- par_deriv(X, y, gamma, beta_c, weights, j)

    # Update interval
    if (der_c == 0)
    {
      # If partial derivative is zero then break loop
      beta_a <- beta_c
      beta_b <- beta_c
      break
    } else if (sign(der_c) == sign(der_a)){
      # Move to right
      a <- c
      beta_a <- beta_c
      der_a <- der_c
    } else {
      # Move to left
      b <- c
      beta_b <- beta_c
      der_b <- der_c
    }
  }

  # Find best of a, b, and 0 in objective function
  obj_a <- obj_fcn(X, y, gamma, beta_a, weights, lambda0)
  obj_b <- obj_fcn(X, y, gamma, beta_b, weights, lambda0)
  obj_0 <- obj_fcn(X, y, gamma, beta_0, weights, lambda0)

  # Return optimal solution
  if ((obj_0 <= obj_a) & (obj_0 <= obj_b)){
    return (beta_0)
  } else if (obj_a <= obj_b) {
    return (beta_a)
  }
  return (beta_b)
}

#' Update gamma and beta\[1\]
#'
#' Finds optimal gamma value and intercept
#' @param X input matrix with dimension n x p, every row is an observation
#' @param y numeric vector for the response variable (binomial)
#' @param beta numeric vector with p coefficients
#' @param weights numeric vector of length n with weights for each
#' observation
#' @return gamma and beta (with updated beta\[1\] value) in a list
update_gamma_intercept <- function(X, y, beta, weights) {

  # Calculate current integer scores and run logistic regression
  z <- X %*% beta - beta[1]*X[,1]
  lr_mod <- stats::glm(y ~ z, weights = weights, family="binomial")

  # Find gamma and beta[1]
  coef_vec <- unname(stats::coef(lr_mod))
  gamma <- coef_vec[2]
  if (is.na(gamma)){
    gamma <- 1
  }
  beta[1] <- coef_vec[1] / gamma

  return (list(gamma=gamma, beta=beta))
}

#' Coordinate descent to find the optimal risk model
#'
#' Returns the estimated gamma and beta for the risk score model
#' @param X input matrix with dimension n x p, every row is an observation
#' @param y numeric vector for the response variable (binomial)
#' @param gamma scalar to rescale betas for prediction
#' @param beta numeric vector with p coefficients
#' @param weights numeric vector of length n with weights for each
#' observation
#' @param lambda0 penalty coefficient for L0 term (default 0)
#' @param a integer lower bound for betas (default -10)
#' @param b integer upper bound for betas (default 10)
#' @param max_iters maximum number of iterations (default 100)
#' @param tol tolerance for convergence
#' @return optimal gamma (numeric) and beta (numeric vector) as a list
risk_coord_desc <- function(X, y, gamma, beta, weights, lambda0 = 0,
                            a = -10, b = 10, max_iters = 100, tol= 1e-5) {

  # Run for maximum number of iterations
  iters <- 1
  while (iters < max_iters)
  {
    # Keep track of old value to check convergence
    old_beta <- beta

    # Iterate through all variables and update intercept/gamma after each
    for (j in (2:ncol(X))){
      beta <- bisec_search(X, y, gamma, beta, weights, j, lambda0, a, b)
      upd <- update_gamma_intercept(X, y, beta, weights)
      gamma <- upd$gamma
      beta <- upd$beta

      # Check for NaN
      if (is.nan(gamma) | sum(is.nan(beta)) > 0){
        stop("Algorithm did not converge - encountered NaN")
      }
    }

    # Check if change in beta is within tolerance to converge
    if (max(abs(old_beta - beta)) < tol){
      break
    }
    iters <- iters+1
  }

  # Check if max iterations
  if(iters >= max_iters) warning("Algorithm reached maximum number of
                                 iterations")
  return(list(gamma=gamma, beta=beta))
}

#' Risk model estimation
#'
#' Returns the estimated gamma and beta for the risk score model along
#' with a glm object with the corresponding coefficients
#' @param X input matrix with dimension n x p, every row is an observation
#' @param y numeric vector for the response variable (binomial)
#' @param gamma starting value to rescale betas for prediction (default NULL)
#' @param beta starting numeric vector with p coefficients (default NULL)
#' @param weights numeric vector of length n with weights for each
#' observation (defult NULL - will give equal weights)
#' @param lambda0 penalty coefficient for L0 term (default 0)
#' @param a integer lower bound for betas (default -10)
#' @param b integer upper bound for betas (default 10)
#' @param max_iters maximum number of iterations (default 100)
#' @param tol tolerance for convergence
#' @return optimal gamma (numeric) and beta (numeric vector) and corresponding
#' glm object (mod) as a list
#' @export
risk_mod <- function(X, y, gamma = NULL, beta = NULL, weights = NULL,
                     lambda0 = 0, a = -10, b = 10, max_iters = 100, tol= 1e-5) {


  # Add intercept column
  if (!all(X[,1] == rep(1, nrow(X)))) {
    X <- cbind(rep(1, nrow(X)), X)
  }

  # Weights
  if (is.null(weights))
    weights <- rep(1, nrow(X))

  # If initial gamma is null but have betas then use update function
  if (is.null(gamma) & (!is.null(beta))){
    upd <- update_gamma_intercept(X, y, beta, weights)
    gamma <- upd$gamma
    beta <- upd$beta
  }

  # Initial beta is null then round LR coefficients using median
  if (is.null(beta)){
    # Initial model
    df <- data.frame(X, y)
    init_mod <- stats::glm(y~.-1, family = "binomial", weights = weights, data = df)

    # Replace NA's with 0's
    coef_vals <- unname(stats::coef(init_mod))
    coef_vals[is.na(coef_vals)] <- 0

    # Round so betas within range
    gamma <- min(abs(a), abs(b))/max(abs(coef_vals[-1]))
    beta <- coef_vals*gamma
    beta[-1] <- round(beta[-1])
  }

  # Check no numeric issues
  if (is.nan(gamma) | sum(is.nan(beta)) > 0){
    stop("Initial gamma or beta is NaN - check starting value for beta")
  }
  if (is.na(gamma) | sum(is.na(beta)) > 0){
    stop("Initial gamma or beta is NA - check starting value for beta")
  }
  if (length(beta) != ncol(X)) stop("beta and X non-compatible")
  if (length(y) != nrow(X)) stop("y and X non-compatible")
  #print(beta)
  # Run coordinate descent from initial solution
  res <- risk_coord_desc(X, y, gamma, beta, weights, lambda0, a, b, max_iters,
                         tol)
  gamma <- res$gamma
  beta <- res$beta



  # Convert to GLM object
  glm_mod <- stats::glm(y~.-1, family = "binomial", weights = weights,
                 start = gamma*beta, method=glm_fit_risk, data = df)
  names(beta) <- names(stats::coef(glm_mod))
  names(beta)[1] <- "(Intercept)"

  # save model score card
  nonzero_beta <- beta[beta != 0][-1]
  if (length(nonzero_beta) <= 1) {
    model_card <- NULL
    score_map <- NULL
  } else {
    model_card <- data.frame(Points = nonzero_beta)

    # get range of possible scores
    X_nonzero <- X[,which(beta != 0)]
    X_nonzero <- X_nonzero[,-1]
    min_pts <- rep(NA, length(nonzero_beta))
    max_pts <- rep(NA, length(nonzero_beta))
    for (i in 1:ncol(X_nonzero)) {
      temp <- nonzero_beta[i] * c(min(X_nonzero[,i]), max(X_nonzero[,i]))
      min_pts[i] <- min(temp)
      max_pts[i] <- max(temp)
    }

    score_range <- seq(sum(min_pts), sum(max_pts))

    # map scores to risk
    v <- gamma*(beta[1] + score_range)
    p <- exp(v)/(1+exp(v))

    # save score map
    score_map <- data.frame(Score = score_range,
                            Risk = round(100*p,1))
  }


  # Return risk_mod object
  mod <- list(gamma=gamma, beta=beta, glm_mod=glm_mod, X=X, y=y, weights=weights,
              lambda0 = lambda0, model_card = model_card, score_map = score_map)
  class(mod) <- "risk_mod"
  return(mod)
}
