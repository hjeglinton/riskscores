#' Summarize Risk Model Fit
#'
#' Prints text that summarizes "risk_mod" objects.
#' @param object An object of class "risk_mod", usually a result of a call to
#'  [risk_mod()].
#' @param ... Additional arguments affecting the summary produced.
#' @return Printed text with intercept, nonzero coefficients, gamma, lambda,
#'  and deviance
#' @examples
#' y <- breastcancer[[1]]
#' X <- as.matrix(breastcancer[,2:ncol(breastcancer)])
#'
#' mod <- risk_mod(X, y, lambda0 = 0.01)
#' summary(mod)
#'
#' @export
summary.risk_mod <- function(object, ...) {

  coef <- object$beta
  res_metrics <- get_metrics(object)

  cat("\nIntercept: ", coef[1], "\n", sep = "")
  cat("\n")

  nonzero_beta <- coef[coef != 0][-1] %>%
    as.data.frame()
  cat("Non-zero coefficients:")
  stats::printCoefmat(nonzero_beta)
  cat("\n")

  cat("Gamma (multiplier): ", object$gamma, "\n")
  cat("Lambda (regularizer): ", object$lambda0, "\n\n")
  cat("Deviance: ", object$glm_mod$deviance, "\n")
  cat("AIC: ", object$glm_mod$aic, "\n\n")

}


#' Extract Model Coefficients
#'
#' Extracts a vector of model coefficients (both nonzero and zero) from a
#'  "risk_mod" object. Equivalent to accessing the `beta` attribute of a
#'  "risk_mod" object.
#' @param object An object of class "risk_mod", usually a result of a call to
#' [risk_mod()].
#' @return Numeric vector with coefficients.
#' @examples
#' y <- breastcancer[[1]]
#' X <- as.matrix(breastcancer[,2:ncol(breastcancer)])
#'
#' mod <- risk_mod(X, y, lambda0 = 0.01)
#' coef(mod)
#'
#' @export
coef.risk_mod <- function(object) {

  return(object$beta)

}


#' Predict Method for Risk Model Fits
#'
#' Obtains predictions from risk score models.
#' @param object An object of class "risk_mod", usually a result of a call to
#'  [risk_mod()].
#' @param newx Optional matrix of new values for `X` for which predictions are
#'  to be made. If ommited, the fitted values are used.
#' @param type The type of prediction required. The default ("link") is on the
#'  scale of the predictors (i.e. log-odds); the "response" type is on the scale
#'  of the response variable (i.e. risk probabilities); the "score" type returns
#'  the risk score calculated from the integer model.
#' @return Numeric vector of predicted values.
#' @examples
#' y <- breastcancer[[1]]
#' X <- as.matrix(breastcancer[,2:ncol(breastcancer)])
#' mod <- risk_mod(X, y, lambda0 = 0.01)

#' predict(mod, type = "link")[1]
#' predict(mod, type = "response")[1]
#' predict(mod, type = "score")[1]
#' @export
predict.risk_mod <- function(object, newx = NULL,
                             type = c("link", "response", "score")) {

  if (is.null(newx)) {

    X <- object$X

  } else {

    X <- newx

    # Add intercept column
    if (!all(X[,1] == rep(1, nrow(X)))) {
      X <- cbind(rep(1, nrow(X)), X)
    }
  }

  # Remove non-zero coefficients
  beta_new <- object$beta[object$beta != 0]
  X_new <- X[,c(1, which(dimnames(X)[[2]] %in% names(beta_new)))]

  if (!all(names(beta_new)[-1] %in% dimnames(X_new)[[2]][-1])) {
    stop(paste0("newx must contain all non-zero covariates (", paste(names(beta_new)[-1], collapse = ", "), ")"))
  }

  v <- object$gamma * X_new %*% beta_new
  v <- clip_exp_vals(v)
  p <- exp(v)/(1+exp(v))

  type = match.arg(type)

  if (type == "link") {
    return(v)
  } else if (type == "response") {
    return(p)
  } else if (type == "score") {
    return(X[,-1] %*% object$beta[-1])
  }

}

#' Plot Risk Score Cross-Validation Results
#'
#' Plots the mean deviance for each \eqn{lambda_0} tested during cross-validation.
#' @param x An object of class "cv_risk_mod", usually a result of a call to
#' [cv_risk_mod()].
#' @param ... Additional arguments affecting the plot produced
#' @return Object of class "ggplot".
#' @export
plot.cv_risk_mod <- function(x, ...) {

  # get mean/sd deviance of lambda_min
  min_mean <- x$results$mean_dev[x$results$lambda0 == x$lambda_min]
  min_sd <- x$results$sd_dev[x$results$lambda0 == x$lambda_min]

  # define x axis breaks
  lambda_grid <- log(x$results$lambda0)
  nlambda <- length(lambda_grid)
  nonzero_seq <- x$results$nonzero
  if (nlambda > 25) {
    new_n <- ceiling(nlambda/25)
    lambda_grid <- lambda_grid[seq(1, nlambda, new_n)]
    nonzero_seq[-seq(1, nlambda, new_n)] <- ""
  }

  # create plot
  lambda0 = mean_dev = sd_dev = NULL # set global variables
  cv_plot <- ggplot2::ggplot(x$results, ggplot2::aes(x = log(lambda0), y = mean_dev)) +
    ggplot2::geom_point() +
    ggplot2::geom_linerange(ggplot2::aes(ymin = mean_dev - sd_dev, ymax= mean_dev + sd_dev)) +
    ggplot2::geom_point(ggplot2::aes(x = log(x$lambda_min), y = min_mean), color = "red") +
    ggplot2::geom_linerange(ggplot2::aes(x = log(x$lambda_min), ymin = min_mean - min_sd,
                       ymax= min_mean + min_sd), color = "red", inherit.aes = FALSE) +
    ggplot2::geom_hline(yintercept = min_mean + min_sd, linetype = "dashed", color = "red") +

    ggplot2::geom_text(ggplot2::aes(x = log(lambda0), label = nonzero_seq,
                  y = (max(mean_dev) + max(sd_dev))*1.01),
              size = 3, col = 'grey30') +

    ggplot2::scale_x_continuous(breaks = lambda_grid, labels = round(lambda_grid, 1)) +
    ggplot2::labs(x = "Log Lambda", y = "Deviance") +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
          axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust=1))


  return(cv_plot)
}

#' Plot Risk Score Model Curve
#'
#' Plots the linear regression equation associated with the integer risk score
#'  model. Plots the scores on the x-axis and risk on the y-axis.
#' @param x An object of class "risk_mod", usually a result of a call to
#' [risk_mod()].
#' @param score_min The minimum score displayed on the x-axis. The default is the
#' minimum score predicted from model's training data.
#' @param score_max The maximum score displayed on the x-axis. The default is the
#' maximum score predicted from model's training data.
#' @param ... Additional arguments affecting the plot produced
#' @return Object of class "ggplot".
#' @examples
#' y <- breastcancer[[1]]
#' X <- as.matrix(breastcancer[,2:ncol(breastcancer)])
#' mod <- risk_mod(X, y, lambda0 = 0.01)
#'
#'plot(mod)
#' @export
plot.risk_mod <- function(x, score_min = NULL, score_max = NULL, ...) {

  if (is.null(score_min)) {
    score_min <- min(predict.risk_mod(x, type = "score"))
  }

  if (is.null(score_max)) {
    score_max <- max(predict.risk_mod(x, type = "score"))
  }


  ggplot2::ggplot() +
    ggplot2::geom_function(data = data.frame(x = seq(score_min, score_max)),
                           ggplot2::aes(x),
                  fun = function(i) get_risk(x, i)) +
    ggplot2::labs(x = "Score", y = "Risk") +
    ggplot2::theme_bw()

}

