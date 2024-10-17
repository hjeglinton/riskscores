#' Partial Derivative of the Negative Log-Likelihood
#'
#' Calculates the partial derivative of the objective function for `beta[j]`.
#' @inheritParams bisec_search
#' @return Numeric partial derivative value.
#' @noRd
par_deriv <- function(X, y, gamma, beta, weights, j) {

  # Calculate partial derivative for NLL
  pd_1 <- sum(gamma * weights * (y * X[,j]))
  exp_pred <- exp(clip_exp_vals(gamma * (X %*% beta)))
  pd_2 <- sum(gamma*weights * X[,j] * (exp_pred / (1.0 + exp_pred)))
  nll_pd <- (-1/nrow(X))*(pd_1-pd_2)

  return(nll_pd)
}

#' Objective Function for NLL+Penalty
#'
#' Calculates the objective function for gamma, beta (NLL+penalty).
#' @inheritParams risk_coord_desc
#' @return Numeric objective function value.
#' @noRd
obj_fcn <- function(X, y, gamma, beta, weights, lambda0=0) {

  # Calculate partial derivative for NLL
  v <- gamma * (X %*% beta)
  v <- clip_exp_vals(v) # avoids numeric errors
  nll_fcn <- (-1/nrow(X))*sum(weights * (y * v - log(1+exp(v))))

  # Penalty term for lambda0*||beta||_0
  pen_fcn <- lambda0*sum(beta[-1] != 0)
  return (nll_fcn + pen_fcn)
}

#' Run Bisection Search
#'
#' Returns optimal value on `beta[j]` using bisection search. For use in each
#'  iteration of the coordinate descent algorithm.
#' @inheritParams risk_coord_desc
#' @param j Index of `beta`.
#' @return Numeric vector `beta` with optimal value for `beta[j]` updated.
#' @noRd
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

#' Update Gamma and Intercept
#'
#' Finds optimal gamma value and intercept (`beta[1]`).
#' @inheritParams risk_coord_desc
#' @return  A list containing the optimal gamma (numeric) and
#'  beta (numeric vector), with gamma and `beta[1]` updated.
#' @noRd
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

#' Run Coordinate Descent
#'
#' Find the optimal risk score model through a coordinate descent algorithm.
#' At each iteration, the algorithm updates the model intercept and scalar value
#' (gamma) based on results from bisection search. Coordinate descent algorithm
#' runs until it converges or reaches the maximum number of iterations.
#' @inheritParams risk_mod
#' @param gamma Scalar to rescale coefficients for prediction.
#' @param beta Numeric vector with \eqn{p} coefficients.
#' @return A list containing the optimal gamma (numeric) and
#'  beta (numeric vector).
#' @noRd
risk_coord_desc <- function(X, y, gamma, beta, weights, lambda0 = 0,
                            a = -10, b = 10, max_iters = 100, tol= 1e-5,
                            shuffle = TRUE) {

  # Run for maximum number of iterations
  iters <- 1
  while (iters < max_iters)
  {
    # Keep track of old value to check convergence
    old_beta <- beta

    # Shuffle order of variables
    if (shuffle == TRUE) {
      variable_index <- sample(2:ncol(X), length(2:ncol(X)), replace = FALSE)
    } else {
      variable_index <- 2:ncol(X)
    }

    # Iterate through all variables and update intercept/gamma after each
    for (j in variable_index){
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

#' Fit an Integer Risk Score Model
#'
#' Fits an optimized integer risk score model using a cyclical coordinate descent
#'  algorithm. Returns an object of class "risk_mod".
#'
#' @details
#'
#' This function uses a cyclical coordinate descent algorithm to solve the
#' following optimization problem.
#'
#'  \deqn{\min_{\alpha,\beta} \quad \frac{1}{n} \sum_{i=1}^{n} (\gamma y_i x_i^T \beta - log(1 + exp(\gamma x_i^T \beta))) + \lambda_0 \sum_{j=1}^{p} 1(\beta_{j} \neq 0)}
#'
#'  \deqn{l \le \beta_j \le u \; \; \; \forall j = 1,2,...,p}
#'  \deqn{\beta_j \in \mathbb{Z} \; \; \; \forall j = 1,2,...,p }
#'  \deqn{\beta_0, \gamma \in \mathbb{R}}
#'
#' These constraints ensure that the model will be sparse and include
#' only integer coefficients.
#'
#' @param X Input covariate matrix with dimension \eqn{n \times p};
#'  every row is an observation.
#' @param y Numeric vector for the (binomial) response variable.
#' @param gamma Starting value to rescale coefficients for prediction (optional).
#' @param beta Starting numeric vector with \eqn{p} coefficients.
#'  Default starting coefficients are rounded coefficients from a
#'  logistic regression model.
#' @param weights Numeric vector of length \eqn{n} with weights for each
#'  observation. Unless otherwise specified, default will give equal weight to
#'  each observation.
#'  @param n_train_runs A positive integer representing the number of times to
#'  train the model, returning the run with the highest accuracy for the
#'  training data.
#' @param lambda0 Penalty coefficient for L0 term (default: 0).
#'  See [cv_risk_mod()] for `lambda0` tuning.
#' @param a Integer lower bound for coefficients (default: -10).
#' @param b Integer upper bound for coefficients (default: 10).
#' @param max_iters Maximum number of iterations (default: 100).
#' @param tol Tolerance for convergence (default: 1e-5).
#' @param shuffle Whether order of coefficients is shuffled during coordinate
#' descent (default: TRUE).
#' @param seed An integer that is used as argument by `set.seed()` for
#'    offsetting the random number generator. Default is to not set a
#'    particular randomization seed.
#' @return An object of class "risk_mod" with the following attributes:
#'  \item{gamma}{Final scalar value.}
#'  \item{beta}{Vector of integer coefficients.}
#'  \item{glm_mod}{Logistic regression object of class "glm" (see [stats::glm]).}
#'  \item{X}{Input covariate matrix.}
#'  \item{y}{Input response vector.}
#'  \item{weights}{Input weights.}
#'  \item{lambda0}{Imput `lambda0` value.}
#'  \item{model_card}{Dataframe displaying the nonzero integer coefficients
#'    (i.e. "points") of the risk score model.}
#'  \item{score_map}{Dataframe containing a column of possible scores and a column
#'    with each score's associated risk probability.}
#' @examples
#' y <- breastcancer[[1]]
#' X <- as.matrix(breastcancer[,2:ncol(breastcancer)])
#'
#' mod1 <- risk_mod(X, y)
#' mod1$model_card
#'
#' mod2 <- risk_mod(X, y, lambda0 = 0.01)
#' mod2$model_card
#'
#' mod3 <- risk_mod(X, y, lambda0 = 0.01, a = -5, b = 5)
#' mod3$model_card
#' @export
risk_mod <- function(X, y, gamma = NULL, beta = NULL, weights = NULL,
                     n_train_runs = 5, lambda0 = 0, a = -10, b = 10,
                     max_iters = 100,  tol = 1e-5, shuffle = TRUE, seed = NULL) {

  # Set seed
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Check that X is a matrix
  if (!is.matrix(X)) stop ("X must be a matrix")

  # Add intercept column
  if (!all(X[,1] == rep(1, nrow(X)))) {
    X <- cbind(Intercept = rep(1, nrow(X)), X)
  }

  # Convert beta to integers within range
  if (any(!(beta%%1==0)) | any(beta < a) | any(beta > b)) {
    if (max(abs(beta[-1])) == 0) {
      scalar <- 1
    } else {
      scalar <- max(abs(beta[-1]))/min(abs(a + 0.5), abs(b + 0.5))
    }
    beta <- beta/scalar
    beta[-1] <- round(beta[-1])
  }

  # Weights
  if (is.null(weights)) {
    weights <- rep(1, nrow(X))}

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
    gamma <- max(abs(coef_vals[-1]))/min(abs(a + 0.5), abs(b + 0.5))
    beta <- coef_vals/gamma
    beta <- randomized_rounding(beta)
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

  if (!is.integer(n_train_runs) | n_train_runs < 0) {
    stop("n_train_runs must be a positive integer")
  }

  # Function to run coordinate descent
  run_risk_mod <- function() {
    # Run coordinate descent from initial solution
    res <- risk_coord_desc(X, y, gamma, beta, weights, lambda0, a, b, max_iters,
                           tol, shuffle)

    gamma <- res$gamma
    beta <- res$beta

    # Convert to GLM object
    glm_mod <- stats::glm(y~.-1, family = "binomial", weights = weights,
                          start = gamma*beta, method=glm_fit_risk, data = data.frame(X, y))
    names(beta) <- names(stats::coef(glm_mod))

    # Save model score card
    nonzero_beta <- beta[beta != 0][-1]
    if (length(nonzero_beta) <= 1) {
      model_card <- NULL
      score_map <- NULL
    } else {
      model_card <- data.frame(Points = nonzero_beta)

      # Get range of possible scores
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

      # Map scores to risk
      v <- gamma*(beta[1] + score_range)
      p <- exp(v)/(1+exp(v))

      # Save score map
      score_map <- data.frame(Score = score_range,
                              Risk = round(p,4))
    }

    # Return risk_mod object
    mod <- list(gamma=gamma, beta=beta, glm_mod=glm_mod, X=X, y=y, weights=weights,
                lambda0 = lambda0, model_card = model_card, score_map = score_map)
    class(mod) <- "risk_mod"
    return(mod)
  }

  # Return the model with the highest accuracy after n_train_runs trains
  highest_acc <- -Inf
  highest_acc_mod <- NULL

  for (i in 1:n_train_runs) {
    mod <- run_risk_mod()
    mod_pred_der <- predict(mod, type = "response")[,1]

    if (mod_pred_der > highest_acc) {
      highest_acc <- mod_pred_der
      highest_acc_mod <- mod
    }
  }

  return(highest_acc_mod)
}
