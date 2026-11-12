#' Transform endpoint values for model fitting
#'
#' Applies variance-stabilizing transformations to endpoint values and their
#' standard errors using the delta method.
#'
#' @param y Numeric vector of endpoint values
#' @param se Numeric vector of standard errors
#' @param type Transformation type: "log", "logit", "sqrt", or "identity"
#'
#' @return List with components:
#'   \item{y}{Transformed endpoint values}
#'   \item{se}{Transformed standard errors (via delta method)}
#'
#' @details
#' The delta method is used to transform standard errors:
#' \itemize{
#'   \item log: SE_new = SE / y
#'   \item logit: SE_new = SE / (y * (1 - y))
#'   \item sqrt: SE_new = SE / (2 * sqrt(y))
#'   \item identity: no transformation
#' }
#'
#' @export
#' @examples
#' # Log-transform rupture rates
#' y <- c(0.5, 1.2, 0.8)
#' se <- c(0.1, 0.2, 0.15)
#' slim_transform(y, se, type = "log")
slim_transform <- function(y, se, type = "identity") {

  if (type == "log") {
    y_out  <- log(y)
    se_out <- se / y
  } else if (type == "logit") {
    y_out  <- log(y / (1 - y))
    se_out <- se / (y * (1 - y))
  } else if (type == "sqrt") {
    y_out  <- sqrt(y)
    se_out <- se / (2 * sqrt(y))
  } else {
    y_out  <- y
    se_out <- se
  }

  list(y = y_out, se = se_out)
}


#' Back-transform estimates to original scale
#'
#' Reverses transformations applied by \code{\link{slim_transform}}.
#'
#' @param est Numeric vector of estimates
#' @param ci.lb Numeric vector of CI lower bounds
#' @param ci.ub Numeric vector of CI upper bounds
#' @param type Transformation type that was applied: "log", "logit", "sqrt", or "identity"
#'
#' @return Data frame with columns:
#'   \item{est}{Back-transformed estimates}
#'   \item{ci.lb}{Back-transformed CI lower bounds}
#'   \item{ci.ub}{Back-transformed CI upper bounds}
#'
#' @export
#' @examples
#' # Back-transform log-scale estimates
#' slim_backtransform(est = -0.5, ci.lb = -1.0, ci.ub = 0.0, type = "log")
slim_backtransform <- function(est, ci.lb, ci.ub, type = "identity") {

  if (type == "log") {
    data.frame(
      est   = exp(est),
      ci.lb = exp(ci.lb),
      ci.ub = exp(ci.ub)
    )
  } else if (type == "logit") {
    data.frame(
      est   = plogis(est),
      ci.lb = plogis(ci.lb),
      ci.ub = plogis(ci.ub)
    )
  } else if (type == "sqrt") {
    data.frame(
      est   = est^2,
      ci.lb = ci.lb^2,
      ci.ub = ci.ub^2
    )
  } else {
    data.frame(est = est, ci.lb = ci.lb, ci.ub = ci.ub)
  }
}


#' Apply display scaling to output
#'
#' Applies a multiplier or transformation for display purposes only.
#' Useful for converting units (e.g., per 1000 patient-years).
#'
#' @param est Numeric vector of estimates
#' @param ci.lb Numeric vector of CI lower bounds
#' @param ci.ub Numeric vector of CI upper bounds
#' @param spec Scaling specification: numeric multiplier (e.g., 1000) or "log"
#'
#' @return Data frame with scaled values (same structure as input)
#'
#' @export
#' @examples
#' # Convert to per 1000 patient-years
#' apply_display_scale(est = 0.005, ci.lb = 0.003, ci.ub = 0.008, spec = 1000)
apply_display_scale <- function(est, ci.lb, ci.ub, spec) {

  if (is.numeric(spec)) {
    data.frame(
      est   = est * spec,
      ci.lb = ci.lb * spec,
      ci.ub = ci.ub * spec
    )
  } else if (is.character(spec) && spec == "log") {
    data.frame(
      est   = log(est),
      ci.lb = log(ci.lb),
      ci.ub = log(ci.ub)
    )
  } else {
    data.frame(est = est, ci.lb = ci.lb, ci.ub = ci.ub)
  }
}
