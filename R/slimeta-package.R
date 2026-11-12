#' slimeta: Stacked Linear Mixed-Effects Framework for Multivariate Meta-Analysis
#'
#' The slimeta package implements the Stacked Linear Mixed-Effects (SLIM)
#' framework for multivariate meta-analysis with correlated endpoints.
#'
#' @section Main functions:
#' \itemize{
#'   \item \code{\link{slim_fit}}: Fit a SLIM model via metafor::rma.mv()
#'   \item \code{\link{slim2EP_plot_by_covariate}}: Panel plot of model vs observed by covariate
#'   \item \code{\link{slim2EP_plot_by_study}}: Scatter plot with E(Y2|Y1) line
#'   \item \code{\link{slim_eta_by_stratum}}: Test eta heterogeneity (Cochran's Q)
#'   \item \code{\link{slim_eta_dunnett}}: Test eta heterogeneity (Dunnett test)
#'   \item \code{\link{slim_sample_size}}: Sample size calculation
#' }
#'
#' @section Helper functions:
#' \itemize{
#'   \item \code{\link{slim_transform}}: Transform y and SE (log, logit, sqrt)
#'   \item \code{\link{slim_backtransform}}: Back-transform to original scale
#'   \item \code{\link{apply_display_scale}}: Apply multiplier for display
#' }
#'
#' @section Mathematical framework:
#' The SLIM model is:
#' \deqn{Y_{ij} = \mu_j + u_{ij} + \varepsilon_{ij}}
#'
#' where \eqn{(u_{i1}, u_{i2})' \sim N(0, G)} and \eqn{G} is the between-study
#' covariance matrix with elements \eqn{\tau_{11}, \tau_{12}, \tau_{22}}.
#'
#' Key quantities:
#' \itemize{
#'   \item \eqn{R^2_{between} = \tau_{12}^2 / (\tau_{11} \times \tau_{22})}
#'   \item \eqn{\eta = \tau_{12} / \tau_{11}} (translation coefficient)
#'   \item \eqn{E(Y_2 | Y_1) = \mu_2 + \eta \times (Y_1 - \mu_1)}
#' }
#'
#' @references
#' Rabbee N. (2026). Stacked Linear Mixed-Effects Models for Evidence
#' Synthesis with Multiple Correlated Endpoints.
#'
#' @docType package
#' @name slimeta-package
#' @aliases slimeta
"_PACKAGE"
