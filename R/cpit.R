#' Conditional probability integral transform
#'
#' Calculates the conditional distribution of the response given the covariates.
#'
#' @param object an object of class \code{cvinereg}.
#' @param newdata matrix of response and covariate values for which to compute
#'   the conditional distribution.
#' @param cores integer; the number of cores to use for computations.
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
#' # PIT
#' u <- cpit(fit, data)
#'
#' @export
#'

cpit <- function(object, newdata, cores = 1) {
  newdata <- prepare_newdata(newdata, object, use_response = TRUE)
  newdata <- to_uscale(newdata, object$margins)
  cond_cvine_dist_cpp(newdata, object$vine, cores)
}

#' Conditional log-likelihood
#'
#' Calculates the conditional log-likelihood of the response given the covariates.
#'
#' @param object an object of class \code{cvinereg}.
#' @param newdata matrix of response and covariate values for which to compute
#'   the conditional distribution.
#' @param cores integer; the number of cores to use for computations.
#'
#' @examples
#' # simulate data
#' x <- matrix(rnorm(200), 100, 2)
#' y <- x %*% c(1, -2)
#' data <- data.frame(y = y, x = x)
#'
#' # fit vine regression model
#' fit <- cvinereg(y ~ ., data)
#' # cll
#' cll(fit, data)
#'
#' @export
#'

cll <- function(object, newdata, cores = 1) {
  newdata <- prepare_newdata(newdata, object, use_response = TRUE)
  ll_marg <- if (inherits(object$margins[[1]], "kde1d")) {
    sum(log(kde1d::dkde1d(newdata[, 1], object$margins[[1]])))
  } else {
    0
  }
  newdata <- to_uscale(newdata, object$margins)
  ll_cop <- cond_cvine_loglik_cpp(newdata, object$vine, cores)
  ll_cop + ll_marg
}

#' Conditional PDF
#'
#' Calculates the conditional density of the response given the covariates.
#'
#' @param object an object of class \code{cvinereg}.
#' @param newdata matrix of response and covariate values for which to compute
#'   the conditional density
#' @param cores integer; the number of cores to use for computations.
#'
#' @examples
#' # simulate data
#' x <- matrix(rnorm(200), 100, 2)
#' y <- x %*% c(1, -2)
#' data <- data.frame(y = y, x = x)
#'
#' # fit vine regression model
#' fit <- cvinereg(y ~ ., data)
#' # pdf
#' d <- cpdf(fit, data)
#' # cll
#' sum(log(d))
#'
#' @export
#'

cpdf <- function(object, newdata, cores = 1) {
  newdata <- prepare_newdata(newdata, object, use_response = TRUE)
  dens_marg <- if (inherits(object$margins[[1]], "kde1d")) {
    kde1d::dkde1d(newdata[, 1], object$margins[[1]])
  } else {
    1
  }
  newdata <- to_uscale(newdata, object$margins)
  cond_cvine_dens_cpp(newdata, object$vine, cores) * dens_marg
}
