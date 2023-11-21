#' Summarize Risk Model Fits
#'
#' Prints text that summarizes risk_mod objects
#' @param object an object of class "risk_mod", usually a result of a call to risk_mod()
#' @param ... additional arguments affecting the summary produced
#' @return text with intercept, nonzero coefficients, gamma, lambda, and deviance
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
#' Extracts model coefficients (both nonzero and zero) from risk_mod object
#' @param object an object of class "risk_mod", usually a result of a call to risk_mod()
#' @return numeric vector with coefficients
#' @export
coef.risk_mod <- function(object) {

  names(object$beta)[1] <- "(Intercept)"
  object$beta

}

#' Model Predictions
#'
#' Obtains predictions from a risk score model object.
#' @param object an object of class "risk_mod", usually a result of a call to
#' risk_mod()
#' @param newdata optionally, a data frame in which to look for variables with
#' which to predict. If omitted, the fitted linear predictors are used.
#' @param type the type of prediction required. The default is on the scale of
#' the linear predictors; the alternative "response is on the scale of the
#' response variable. Thus for a risk score model, the default predictions are of
#' log-odds (probabilities on logit scale) and type = "response" gives the
#' predicted risk probabilities. The "score" option returns integer risk scores.
#' @return array with predictions
#' @export
predict.risk_mod <- function(object, newdata = NULL,
                             type = c("link", "response", "score")) {

  if (is.null(newdata)) {
    X <- object$X
  } else {
    X <- newdata
  }

  v <- object$gamma * X %*% object$beta
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

#' Plot cross-validation results
#'
#' Plots the mean deviance for each lambda tested during cross-validation
#' @param x an object of class "cv_risk_mod", usually a result of a call to
#' cv_risk_mod()
#' @param ... additional arguments affecting the plot produced
#' @return ggplot object
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


