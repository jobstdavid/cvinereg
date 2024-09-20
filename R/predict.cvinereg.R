#' Predict conditional mean and quantiles from a C-vine regression model
#'
#' @param object an object of class \code{cvinereg}.
#' @param newdata matrix of covariate values for which to predict the quantile.
#' @param alpha vector of quantile levels; `NA` predicts the mean based on an
#'   average of the `1:10 / 11`-quantiles.
#' @param cores integer; the number of cores to use for computations.
#' @param ... unused.
#'
#' @return A data.frame of quantiles where each column corresponds to one
#' value of `alpha`.
#'
#' @examples
#' # simulate data
#' x <- matrix(rnorm(200), 100, 2)
#' y <- x %*% c(1, -2)
#' data <- data.frame(y = y, x = x)
#'
#' # fit vine regression model
#' fit <- cvinereg(y ~ ., data)
#'
#' # model predictions
#' mu_hat <- predict(fit, newdata = data, alpha = NA) # mean
#' med_hat <- predict(fit, newdata = data, alpha = 0.5) # median
#'
#' @seealso \code{\link{cvinereg}}
#'
#' @export
#'
#' @importFrom kde1d pkde1d qkde1d
#' @importFrom stats predict
predict.cvinereg <- function(object, newdata, alpha = 0.5, cores = 1, ...) {
  if (missing(newdata)) {
    return(fitted.cvinereg(object, alpha = alpha))
  }

  stopifnot(length(alpha) > 0)
  if (any(is.na(alpha)) & inherits(object$model_frame[[1]], "ordered")) {
    stop("cannot predict mean for ordered response.")
  }

  # predict the conditional mean if alpha contains NA
  if (any(is.na(alpha))) {
    alpha <- alpha[!is.na(alpha)] # remove NA for quantile estimation
    preds_mean <- predict_mean(object, newdata)
  } else {
    preds_mean <- NULL
  }

  ## computation of conditional quantiles
  if (length(alpha) > 0) {
    stopifnot(is.numeric(alpha), all(alpha > 0), all(alpha < 1))

    ## preprocessing
    newdata <- prepare_newdata(newdata, object)
    newdata <- to_uscale(newdata, object$margins[-1], add_response = TRUE)
    preds <- qcvine(newdata, alpha, vine = object$vine, cores)

    ## actual predictions on original scale
    preds <- to_yscale(preds, object)
    if (!is.null(preds_mean))
      preds <- cbind(preds_mean, preds)
  } else {
    preds <- preds_mean
  }

  preds
}

#' @rdname predict.cvinereg
#' @importFrom stats fitted
#' @export
fitted.cvinereg <- function(object, alpha = 0.5, ...) {
  predict.cvinereg(object, newdata = object$model_frame, alpha = alpha, ...)
}


#' predicts the conditional mean as the average of quantiles.
#' @noRd
predict_mean <- function(object, newdata) {
  preds <- predict.cvinereg(object, newdata, alpha = 1:10 / 11)
  data.frame(mean = rowMeans(preds))
}

qcvine <- function(u, alpha, vine, cores) {
  vine$var_types[1] <- "c"
  q_hat <- as.data.frame(cond_cvine_quantile_cpp(alpha, as.matrix(u), vine, cores))
  names(q_hat) <- alpha
  q_hat
}
