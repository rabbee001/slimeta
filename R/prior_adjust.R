#' Apply Bayesian prior adjustment to treatment effects
#'
#' Pre-processes treatment effects (e.g., rate ratios near zero) by applying
#' Bayesian shrinkage with a prior distribution. This helps stabilize extreme
#' estimates before stacking in the SLIM model.
#'
#' @importFrom stats median
#'
#' @param data Data frame containing treatment effects
#' @param y_col Column name for treatment effect estimates
#' @param se_col Column name for standard errors
#' @param prior_type Type of prior: "normal" (default), "flat", or "halfnormal"
#' @param prior_mean Prior mean (default: 0, i.e., no effect)
#' @param prior_sd Prior standard deviation. If NULL, estimated empirically
#'   from the data using the median SE.
#' @param scale Scale of the treatment effects: "log" (default) or "original".
#'   If "original" and values are rate/risk ratios, they will be log-transformed

#'   before applying the prior, then back-transformed.
#' @param min_value Minimum allowed value before transformation (default: 1e-6).
#'   Values below this are set to min_value to avoid log(0).
#' @param suffix Suffix for new columns (default: "_adj")
#'
#' @return Data frame with additional columns for adjusted estimates:
#'   \code{{y_col}{suffix}} (posterior mean) and \code{{se_col}{suffix}}
#'   (posterior SD).
#'
#' @details
#' For a normal prior N(prior_mean, prior_sd^2) and observed effect y with
#' standard error se, the posterior is:
#'
#' \deqn{posterior\_mean = \frac{y/se^2 + prior\_mean/prior\_sd^2}{1/se^2 + 1/prior\_sd^2}}
#' \deqn{posterior\_var = \frac{1}{1/se^2 + 1/prior\_sd^2}}
#'
#' This shrinks extreme estimates toward the prior mean, with shrinkage
#' proportional to uncertainty (larger SE = more shrinkage).
#'
#' A flat prior (prior_type = "flat") returns the original estimates unchanged,
#' which is useful as a baseline comparison.
#'
#' The half-normal prior (prior_type = "halfnormal") uses |N(0, prior_sd^2)|,
#' appropriate when effects are expected to be non-negative (e.g., hazard ratios > 1).
#'
#' @examples
#' # Shrink extreme or imprecise log rate ratios toward the prior mean (0)
#' d <- data.frame(
#'   study     = paste0("S", 1:8),
#'   log_RR    = c(-1.2, -0.4, -0.1, 0.05, 0.3, -0.7, -2.0, 0.6),
#'   SE_log_RR = c(0.50, 0.30, 0.20, 0.90, 0.25, 0.40, 0.80, 0.35)
#' )
#'
#' d_adj <- slim_prior_adjust(
#'   data     = d,
#'   y_col    = "log_RR",
#'   se_col   = "SE_log_RR",
#'   prior_sd = 1.0
#' )
#'
#' # Adjusted columns are added with the suffix "_adj"
#' d_adj[, c("study", "log_RR", "log_RR_adj", "SE_log_RR", "SE_log_RR_adj")]
#'
#' @export
slim_prior_adjust <- function(data,
                              y_col,
                              se_col,
                              prior_type = c("normal", "flat", "halfnormal"),
                              prior_mean = 0,
                              prior_sd = NULL,
                              scale = c("log", "original"),
                              min_value = 1e-6,
                              suffix = "_adj") {


  prior_type <- match.arg(prior_type)
  scale <- match.arg(scale)


  # Extract values

  y <- data[[y_col]]
  se <- data[[se_col]]


  if (is.null(y) || is.null(se)) {
    stop("Columns '", y_col, "' and '", se_col, "' must exist in data")
  }


  # Handle original scale (rate/risk ratios)
  if (scale == "original") {
    # Ensure positive values
    y[y < min_value] <- min_value
    # Log transform
    y_log <- log(y)
    # SE on log scale via delta method: SE_log = SE / y
    se_log <- se / y
  } else {
    y_log <- y
    se_log <- se
  }

  # Estimate prior_sd empirically if not provided
  if (is.null(prior_sd)) {
    prior_sd <- median(se_log, na.rm = TRUE) * 2
    message("Using empirically estimated prior_sd = ", round(prior_sd, 3))
  }

  # Apply prior adjustment
  if (prior_type == "flat") {
    # Flat prior: posterior = likelihood
    y_adj <- y_log
    se_adj <- se_log

  } else if (prior_type == "normal") {
    # Normal prior: conjugate update
    prior_prec <- 1 / prior_sd^2
    data_prec <- 1 / se_log^2

    # Posterior precision and variance
    post_prec <- prior_prec + data_prec
    post_var <- 1 / post_prec

    # Posterior mean (precision-weighted average)
    y_adj <- (y_log * data_prec + prior_mean * prior_prec) / post_prec
    se_adj <- sqrt(post_var)

  } else if (prior_type == "halfnormal") {
    # Half-normal: approximate via truncated normal update
    # For simplicity, use normal prior but constrain posterior mean >= 0
    prior_prec <- 1 / prior_sd^2
    data_prec <- 1 / se_log^2

    post_prec <- prior_prec + data_prec
    post_var <- 1 / post_prec

    y_adj <- (y_log * data_prec + 0 * prior_prec) / post_prec
    se_adj <- sqrt(post_var)

    # Truncate at 0 for half-normal (effects should be >= 0 on log scale means RR >= 1)
    y_adj <- pmax(y_adj, 0)
  }

  # Back-transform if original scale was provided
  if (scale == "original") {
    # Back to original scale
    y_final <- exp(y_adj)
    # SE via delta method: SE_original = y_original * SE_log
    se_final <- y_final * se_adj
  } else {
    y_final <- y_adj
    se_final <- se_adj
  }

  # Add adjusted columns to data
  data[[paste0(y_col, suffix)]] <- y_final
  data[[paste0(se_col, suffix)]] <- se_final

  # Add shrinkage diagnostics
  shrinkage <- 1 - (se_final^2 / se^2)
  shrinkage[shrinkage < 0] <- 0  # numerical cleanup

  attr(data, "prior_adjustment") <- list(
    prior_type = prior_type,
    prior_mean = prior_mean,
    prior_sd = prior_sd,
    scale = scale,
    mean_shrinkage = mean(shrinkage, na.rm = TRUE),
    n_adjusted = sum(!is.na(y_final))
  )

  data
}


#' Visualize prior adjustment effect
#'
#' Creates a diagnostic plot comparing original vs prior-adjusted estimates.
#'
#' @param data Data frame output from \code{\link{slim_prior_adjust}}
#' @param y_col Original y column name
#' @param se_col Original SE column name
#' @param suffix Suffix used in slim_prior_adjust (default: "_adj")
#' @param title Plot title
#'
#' @return A ggplot object
#'
#' @examples
#' d <- data.frame(
#'   study     = paste0("S", 1:8),
#'   log_RR    = c(-1.2, -0.4, -0.1, 0.05, 0.3, -0.7, -2.0, 0.6),
#'   SE_log_RR = c(0.50, 0.30, 0.20, 0.90, 0.25, 0.40, 0.80, 0.35)
#' )
#'
#' d_adj <- slim_prior_adjust(d, y_col = "log_RR", se_col = "SE_log_RR", prior_sd = 1.0)
#' slim_prior_plot(d_adj, y_col = "log_RR", se_col = "SE_log_RR")
#'
#' @export
slim_prior_plot <- function(data,
                            y_col,
                            se_col,
                            suffix = "_adj",
                            title = "Prior Adjustment: Original vs Adjusted") {

  y_orig <- data[[y_col]]
  y_adj <- data[[paste0(y_col, suffix)]]
  se_orig <- data[[se_col]]
  se_adj <- data[[paste0(se_col, suffix)]]

  if (is.null(y_adj)) {
    stop("Adjusted column '", paste0(y_col, suffix), "' not found. ",
         "Run slim_prior_adjust() first.")
  }

  plot_df <- data.frame(
    Original = y_orig,
    Adjusted = y_adj,
    SE_orig = se_orig,
    SE_adj = se_adj
  )

  # Get prior info if available
  prior_info <- attr(data, "prior_adjustment")
  subtitle <- if (!is.null(prior_info)) {
    sprintf("Prior: %s(%.2f, %.2f), Mean shrinkage: %.1f%%",
            prior_info$prior_type,
            prior_info$prior_mean,
            prior_info$prior_sd,
            prior_info$mean_shrinkage * 100)
  } else {
    NULL
  }

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = Original, y = Adjusted)) +
    ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed",
                         color = "gray50") +
    ggplot2::geom_point(ggplot2::aes(size = 1/SE_orig), alpha = 0.6) +
    ggplot2::geom_segment(ggplot2::aes(xend = Adjusted, yend = Adjusted),
                          alpha = 0.3, color = "steelblue") +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = paste("Original", y_col),
      y = paste("Adjusted", y_col),
      size = "Precision\n(1/SE)"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "right")

  p
}
