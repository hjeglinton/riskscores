#' Initial Temperature for Simulated Annealing
#'
#' Finds initial temperature so that will accept 90% of current obj with 95% probability
getInitTemp <- function(X, y, gamma, beta, weights, lambda0) {
  obj <- obj_fcn(X, y, gamma, beta, weights, lambda0)
  max_obj = 1.1*obj
  T <- (obj - max_obj) / log(0.95)
  return(T)
}

#' Alpha for Simulated Annealing
#'
#' Calculates alpha value so that will accept initial obj with 15% probability
getAlpha <- function(p) {
  # Calculate alpha using the corrected formula
  alpha <- (log(0.95) / log(0.15))^(1 / p)
  return(max(alpha, 0.95))
}

#' @param e1 current objective
#' @param e2 candidate objective
getAcceptanceProb <- function(e1, e2, T) {
  if (e2 < e1) {
    return(1)
  }
  return(exp(-(e2 - e1) / T))
}

#' Run Simulated Annealing
#'
#' Find the optimal risk score model through a simulated annealing algorithm.
#' At each iteration, the algorithm updates the model intercept and scalar value
#' (gamma) based on a randomly generated neighbor. Algorithm runs until it
#' converges or reaches the maximum number of iterations.
#' @inheritParams risk_mod
#' @param gamma Scalar to rescale coefficients for prediction.
#' @param beta Numeric vector with \eqn{p} coefficients.
#' @return A list containing the optimal gamma (numeric) and
#'  beta (numeric vector).
#' @noRd
simulated_annealing <- function(X, y, a, b, gamma, beta, weights, tol=1e-6,
                                lambda0 = 0, max_iters = 1000) {
  # Getting initial objective function and temperature
  obj <- obj_fcn(X, y, gamma, beta, weights, lambda0)
  T <- getInitTemp(X, y, gamma, beta, weights, lambda0)
  p <- ncol(X)-1
  p_try <- max(floor(p/2),1)
  
  alpha <- getAlpha(p)
  best_beta <- beta
  best_obj <- obj
  best_gamma <- gamma

  for (iter in 1:max_iters) {
    # Cool down
    T <- alpha * T

    # Selecting subset of indices
    indices <- sample(1:p, p_try)
    
    X_sub <- X[,c(indices+1)]
    score = as.vector(X %*% beta)
    X_sub <- cbind(X_sub, score = score)
    
    # Building an LR model for residuals
    mod_scores <- glm(y ~ score, data = data.frame(X, y = y, score = score))
    resid <- y-predict(mod_scores, type = "response")
    mod_lm <- lm(resid ~ ., data = as.data.frame(X_sub), weights = weights)
    coef_lm <- coef(mod_lm)
    
    # Replace NA's with 0's
    coef_lm[is.na(coef_lm)] <- 0
    
    selected_coefs <- coef_lm[-c(1, length(coef_lm))]
    
    # If all coefficients are zero, then skip this iteration
    if (max(abs(selected_coefs)) == 0) {
      next
    }
    
    scaled_coefs <- abs(selected_coefs) / max(abs(selected_coefs))
    try_beta <- beta 

    # Randomly round
    for (i in 1:p_try) {
      var_index <- indices[i]
      
      # beta[i] +1 or -1 if p-value is significant
      if (runif(1) < scaled_coefs[i]) {
        try_beta[var_index + 1] <- try_beta[var_index + 1] + sign(selected_coefs[i])
      }
    }
    try_beta <- pmax(try_beta, a)
    try_beta <- pmin(try_beta, b)
    
    upd <- update_gamma_intercept(X, y, try_beta, weights)
    try_beta <- upd$beta
    try_gamma <- upd$gamma

    try_obj <- obj_fcn(X, y, try_gamma, try_beta, weights, lambda0)
    
    # Accept the new solution based on probability
    acc_prob <- getAcceptanceProb(best_obj, try_obj, T)
    
    if (runif(1) < acc_prob ) {
      obj <- try_obj
      beta <- try_beta
      gamma <- try_gamma
    } 
    
    # Update best seen objective function value
    if (obj < best_obj) {
      best_obj <- obj
      best_beta <- beta
      best_gamma <- gamma
    }
    
    if (T < tol) {
      break
    }
  }

  # Run one iteration of coordinate descent
  result <- risk_coord_desc(X, y, best_gamma, best_beta, weights, lambda0 = lambda0, a = a, b = b, max_iters = 2)
  return (list(gamma=result$gamma, beta=result$beta))
}

#' Fit an Integer Risk Score Model
#'
#' Fits an optimized integer risk score model using a simulated annealing
#'  algorithm. Returns an object of class "risk_mod".
#'
#' @details
#'
#' This function uses a simulated annealing algorithm to solve the
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
#'  initialize and train the model, returning the run with the lowest objective
#'  function for the training data.
#' @param lambda0 Penalty coefficient for L0 term (default: 0).
#'  See [cv_annealscore()] for `lambda0` tuning.
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
#' mod1 <- annealscore(X, y)
#' mod1$model_card
#'
#' mod2 <- annealscore(X, y, lambda0 = 0.01)
#' mod2$model_card
#'
#' mod3 <- annealscore(X, y, lambda0 = 0.01, a = -5, b = 5)
#' mod3$model_card
#' @export
annealscore <- function(X, y, gamma = NULL, beta = NULL, weights = NULL,
                        a = -10, b = 10,
                        lambda0 = 0, max_iters = 10000, tol= 1e-5, seed = NULL) {
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
  
  # Weights
  if (is.null(weights)) {
    weights <- rep(1, nrow(X))
  }
  
  # If initial gamma is null but have betas then use update function
  if (is.null(gamma) & (!is.null(beta))){
    upd <- update_gamma_intercept(X, y, beta, weights)
    gamma <- upd$gamma
    beta <- upd$beta
  }
  
  # Initial beta is null then round LR coefficients using median and randomized rounding
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
  
  # Run simulated annealing algorithm from initial solution
  res <- simulated_annealing(X, y, a, b, gamma, beta, weights, tol=tol, lambda0=lambda0, max_iters=max_iters)
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
