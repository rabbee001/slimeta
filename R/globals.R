# Global variable declarations to avoid R CMD check NOTEs
# These are column names used in ggplot2 aes() calls and data.frame operations

utils::globalVariables(c(
  "Model_Est", "CI_lb", "CI_ub", "Observed",
  "X_obs", "Y_obs",
  "x", "y", "y_lo", "y_hi",
  "yi", "Endpoint",
  "Original", "Adjusted", "SE_orig", "SE_adj",
  "Study", "ci_lo", "ci_hi", "se"
))

#' @importFrom stats aggregate as.formula coef confint pchisq plogis pt qnorm vcov cor.test
NULL
