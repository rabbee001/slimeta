#' Panel plot of two endpoints by covariate levels
#'
#' Creates a side-by-side panel plot comparing model estimates to observed
#' means across levels of a covariate (e.g., diameter strata).
#'
#' @param result Output from \code{\link{slim_fit}}
#' @param data Original data frame
#' @param endpoints Character vector of length 2: c("EP1", "EP2"). Order determines panel order.
#' @param covariate Column name for x-axis covariate (e.g., "Diameter")
#' @param title Plot title
#' @param output What to return: "plot", "table", or "both" (default)
#' @param save_to File path to save plot (NULL = display only)
#' @param width Plot width in inches
#' @param height Plot height in inches
#'
#' @return A list with components `plot` (ggplot object) and `table`
#'   (data frame with observed and model estimates per covariate level).
#'
#' @export
#' @examples
#' \donttest{
#' data(aaa)
#' # RuptureRisk is stored on the log scale in aaa, but slim_fit() applies the
#' # log transform itself. Supply the rate and its SE on the original scale.
#' aaa_fit <- aaa
#' rup <- aaa_fit$Endpoint == "RuptureRisk"
#' aaa_fit$y[rup]  <- aaa_fit$Mean[rup]
#' aaa_fit$se[rup] <- aaa_fit$Mean[rup] * aaa_fit$se[rup]
#'
#' result <- slim_fit(
#'   data          = aaa_fit,
#'   formula       = ~ 0 + Endpoint:factor(Diameter),
#'   trial_col     = "Study",
#'   endpoint_col  = "Endpoint",
#'   y_col         = "y",
#'   se_col        = "se",
#'   transform     = list(RuptureRisk = "log"),
#'   display_scale = list(RuptureRisk = 1000)
#' )
#'
#' out <- slim2EP_plot_by_covariate(
#'   result,
#'   data      = aaa_fit,
#'   endpoints = c("GrowthRate", "RuptureRisk"),
#'   covariate = "Diameter",
#'   title     = "AAA: Model vs Observed by Diameter"
#' )
#' out$plot
#' }
slim2EP_plot_by_covariate <- function(
    result,
    data,
    endpoints,
    covariate,
    title   = NULL,
    output  = "both",
    save_to = NULL,
    width   = 10,
    height  = 5
) {

  # --- Require covariate ---
  if (missing(covariate) || is.null(covariate)) {
    stop("covariate is required for this plot. Use slim_forest() or slim2EP_plot_by_study() for unstratified data.")
  }

  # --- Extract metadata from result ---
  endpoint_col  <- result$endpoint_col
  trial_col     <- result$trial_col
  y_col         <- result$y_col
  display_scale <- result$display_scale
  pooled        <- result$pooled

  # --- Validate endpoints ---
  if (length(endpoints) != 2) {
    stop("Must specify exactly 2 endpoints")
  }

  # --- Compute observed means ---
  obs_list <- list()
  for (ep in endpoints) {
    ep_data <- data[data[[endpoint_col]] == ep, ]

    obs_agg <- aggregate(
      ep_data[[y_col]],
      by = list(ep_data[[covariate]]),
      FUN = mean,
      na.rm = TRUE
    )
    names(obs_agg) <- c(covariate, "Observed")
    obs_agg$Endpoint <- ep

    if (ep %in% names(display_scale) && is.numeric(display_scale[[ep]])) {
      obs_agg$Observed <- obs_agg$Observed * display_scale[[ep]]
    }

    obs_list[[ep]] <- obs_agg
  }
  observed <- do.call(rbind, obs_list)

  # --- Parse pooled estimates to extract endpoint and covariate ---
  pooled$Endpoint <- NA_character_
  pooled[[covariate]] <- NA

  for (i in seq_len(nrow(pooled))) {
    term <- pooled$term[i]
    for (ep in endpoints) {
      if (grepl(ep, term, fixed = TRUE)) {
        pooled$Endpoint[i] <- ep
        pattern <- paste0(".*factor\\(", covariate, "\\)")
        cov_val <- gsub(pattern, "", term)
        pooled[[covariate]][i] <- as.numeric(cov_val)
        break
      }
    }
  }

  pooled <- pooled[pooled$Endpoint %in% endpoints, ]

  # --- Merge observed and model estimates ---
  tbl <- merge(
    observed,
    pooled[, c("Endpoint", covariate, "est", "ci.lb", "ci.ub")],
    by = c("Endpoint", covariate),
    all = TRUE
  )
  names(tbl)[names(tbl) == "est"]   <- "Model_Est"
  names(tbl)[names(tbl) == "ci.lb"] <- "CI_lb"
  names(tbl)[names(tbl) == "ci.ub"] <- "CI_ub"

  tbl$Endpoint <- factor(tbl$Endpoint, levels = endpoints)

  # --- Build plot ---
  p <- ggplot2::ggplot(tbl, ggplot2::aes_string(x = covariate)) +
    ggplot2::geom_point(ggplot2::aes(y = Model_Est), color = "red", size = 3) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = CI_lb, ymax = CI_ub), color = "red", width = 0.3) +
    ggplot2::geom_point(ggplot2::aes(y = Observed), shape = 8, size = 4, color = "black") +
    ggplot2::facet_wrap(~ Endpoint, scales = "free_y") +
    ggplot2::labs(
      title    = title,
      subtitle = "Red = Model (95% CI), Black star = Observed mean",
      x        = covariate,
      y        = "Estimate"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      strip.background = ggplot2::element_rect(fill = "gray90", color = "gray50"),
      strip.text       = ggplot2::element_text(face = "bold")
    )

  if (!is.null(save_to)) {
    ggplot2::ggsave(save_to, p, width = width, height = height)
    message("Plot saved to: ", save_to)
  }

  if (output == "plot") {
    return(list(plot = p))
  } else if (output == "table") {
    return(list(table = tbl))
  } else {
    return(list(plot = p, table = tbl))
  }
}


#' Scatter plot of two endpoints by study
#'
#' Creates a scatter plot of one endpoint vs another, with each point
#' representing a study. Includes the conditional expectation line
#' E(Y2|Y1) = mu2 + eta*(Y1-mu1) with confidence band.
#'
#' @param result Output from \code{\link{slim_fit}}
#' @param data Original data frame
#' @param endpoints Character vector of length 2: c("X", "Y"). First is X-axis, second is Y-axis.
#' @param title Plot title
#' @param show_line Show E(Y|X) regression line (default: TRUE)
#' @param show_band Show conditional CI band (default: TRUE)
#' @param color_by Column name to color points by (optional)
#' @param facet_by Column name to facet panels by (optional, e.g., "Diameter")
#' @param facet_scales "free" or "fixed" axis scales across facets
#' @param eta_by_stratum Output from \code{\link{slim_eta_by_stratum}} for stratum-specific eta lines
#' @param scale "original" (back-transformed, default) or "transformed"
#' @param y_scale Y-axis scale: "linear" (default) or "log"
#' @param pct_changes Percent changes for translation table (default: -10 percent, -20 percent, -30 percent)
#' @param output What to return: "plot", "table", or "both" (default)
#' @param save_to File path to save plot
#' @param width Plot width in inches
#' @param height Plot height in inches
#'
#' @return A list with components `plot` (ggplot object), `table`
#'   (data frame with observed X, Y and model E(Y|X) per study), `translation`
#'   (data frame showing X percent change to Y percent change), and `eta`
#'   (translation coefficient tau_XY/tau_XX).
#'
#' @export
#' @examples
#' \donttest{
#' data(aaa)
#' # RuptureRisk is stored on the log scale in aaa, but slim_fit() applies the
#' # log transform itself. Supply the rate and its SE on the original scale.
#' aaa_fit <- aaa
#' rup <- aaa_fit$Endpoint == "RuptureRisk"
#' aaa_fit$y[rup]  <- aaa_fit$Mean[rup]
#' aaa_fit$se[rup] <- aaa_fit$Mean[rup] * aaa_fit$se[rup]
#'
#' result <- slim_fit(
#'   data          = aaa_fit,
#'   formula       = ~ 0 + Endpoint:factor(Diameter),
#'   trial_col     = "Study",
#'   endpoint_col  = "Endpoint",
#'   y_col         = "y",
#'   se_col        = "se",
#'   transform     = list(RuptureRisk = "log"),
#'   display_scale = list(RuptureRisk = 1000)
#' )
#'
#' out <- slim2EP_plot_by_study(
#'   result,
#'   data         = aaa_fit,
#'   endpoints    = c("GrowthRate", "RuptureRisk"),
#'   title        = "AAA: Growth vs Rupture",
#'   facet_by     = "Diameter",
#'   facet_scales = "fixed",
#'   y_scale      = "log"
#' )
#' out$plot
#' out$translation
#' }
slim2EP_plot_by_study <- function(
    result,
    data,
    endpoints,
    title           = NULL,
    show_line       = TRUE,
    show_band       = TRUE,
    color_by        = NULL,
    facet_by        = NULL,
    facet_scales    = "free",
    eta_by_stratum  = NULL,
    scale           = "original",
    y_scale         = "linear",
    pct_changes     = c(-0.10, -0.20, -0.30),
    output          = "both",
    save_to         = NULL,
    width           = 7,
    height          = 6
) {

  # --- Extract metadata ---
  endpoint_col  <- result$endpoint_col
  trial_col     <- result$trial_col
  y_col         <- result$y_col
  G             <- result$G
  fit           <- result$fit
 transform     <- result$transform
  display_scale <- result$display_scale

  if (length(endpoints) != 2) {
    stop("Must specify exactly 2 endpoints")
  }

  ep_x <- endpoints[1]
  ep_y <- endpoints[2]

  # --- Get G matrix indices ---
  G_names <- rownames(G)
  idx_x <- which(G_names == ep_x)
  idx_y <- which(G_names == ep_y)

  if (length(idx_x) == 0 || length(idx_y) == 0) {
    stop("Endpoints not found in G matrix. Available: ", paste(G_names, collapse = ", "))
  }

  tau_xx <- G[idx_x, idx_x]
  tau_yy <- G[idx_y, idx_y]
  tau_xy <- G[idx_x, idx_y]

  eta <- tau_xy / tau_xx
  cond_var <- tau_yy - (tau_xy^2) / tau_xx
  cond_sd  <- sqrt(max(0, cond_var))

  # --- Get overall means from fixed effects ---
  beta <- coef(fit)

  get_mu <- function(ep) {
    idx <- grepl(ep, names(beta), fixed = TRUE)
    if (sum(idx) > 0) mean(beta[idx]) else NA
  }

  mu_x <- get_mu(ep_x)
  mu_y <- get_mu(ep_y)

  # --- Pivot data to wide format ---
  keep_cols_x <- c(trial_col, y_col, color_by, facet_by)
  keep_cols_x <- unique(keep_cols_x[!is.null(keep_cols_x)])
  keep_cols_y <- c(trial_col, y_col, facet_by)
  keep_cols_y <- unique(keep_cols_y[!is.null(keep_cols_y)])

  data_x <- data[data[[endpoint_col]] == ep_x, keep_cols_x]
  data_y <- data[data[[endpoint_col]] == ep_y, keep_cols_y]
  names(data_x)[names(data_x) == y_col] <- "X_obs"
  names(data_y)[names(data_y) == y_col] <- "Y_obs"

  merge_by <- trial_col
  if (!is.null(facet_by)) merge_by <- c(merge_by, facet_by)
  study_data <- merge(data_x, data_y, by = merge_by, all = TRUE)

  # --- Apply display scale to observed ---
  if (ep_x %in% names(display_scale) && is.numeric(display_scale[[ep_x]])) {
    study_data$X_obs <- study_data$X_obs * display_scale[[ep_x]]
  }
  if (ep_y %in% names(display_scale) && is.numeric(display_scale[[ep_y]])) {
    study_data$Y_obs <- study_data$Y_obs * display_scale[[ep_y]]
  }

  # --- Compute E[Y|X] for each study ---
  if (scale == "original") {
    x_trans <- if (ep_x %in% names(transform)) {
      slim_transform(study_data$X_obs / ifelse(ep_x %in% names(display_scale), display_scale[[ep_x]], 1),
                     rep(1, nrow(study_data)), transform[[ep_x]])$y
    } else {
      study_data$X_obs / ifelse(ep_x %in% names(display_scale), display_scale[[ep_x]], 1)
    }

    y_cond_trans <- mu_y + eta * (x_trans - mu_x)

    if (ep_y %in% names(transform)) {
      bt <- slim_backtransform(y_cond_trans,
                               y_cond_trans - 1.96 * cond_sd,
                               y_cond_trans + 1.96 * cond_sd,
                               transform[[ep_y]])
      study_data$Y_model <- bt$est
      study_data$Y_CI_lb <- bt$ci.lb
      study_data$Y_CI_ub <- bt$ci.ub
    } else {
      study_data$Y_model <- y_cond_trans
      study_data$Y_CI_lb <- y_cond_trans - 1.96 * cond_sd
      study_data$Y_CI_ub <- y_cond_trans + 1.96 * cond_sd
    }

    if (ep_y %in% names(display_scale) && is.numeric(display_scale[[ep_y]])) {
      study_data$Y_model <- study_data$Y_model * display_scale[[ep_y]]
      study_data$Y_CI_lb <- study_data$Y_CI_lb * display_scale[[ep_y]]
      study_data$Y_CI_ub <- study_data$Y_CI_ub * display_scale[[ep_y]]
    }
  } else {
    study_data$Y_model <- mu_y + eta * (study_data$X_obs - mu_x)
    study_data$Y_CI_lb <- study_data$Y_model - 1.96 * cond_sd
    study_data$Y_CI_ub <- study_data$Y_model + 1.96 * cond_sd
  }

  # --- Build translation table ---
  translation <- data.frame(
    X_pct_change = pct_changes * 100,
    stringsAsFactors = FALSE
  )

  for (i in seq_along(pct_changes)) {
    pct <- pct_changes[i]
    x_new <- mu_x * (1 + pct)
    y_cond <- mu_y + eta * (x_new - mu_x)

    if (ep_y %in% names(transform) && transform[[ep_y]] == "log") {
      y_abs <- exp(y_cond)
      y_base <- exp(mu_y)
    } else {
      y_abs <- y_cond
      y_base <- mu_y
    }

    if (ep_y %in% names(display_scale) && is.numeric(display_scale[[ep_y]])) {
      y_abs  <- y_abs * display_scale[[ep_y]]
      y_base <- y_base * display_scale[[ep_y]]
    }

    translation$Y_abs[i] <- y_abs
    translation$Y_pct_change[i] <- ((y_abs - y_base) / y_base) * 100
  }

  # --- Build regression line(s) ---
  if (!is.null(facet_by)) {
    strata <- unique(study_data[[facet_by]])
    strata <- strata[!is.na(strata)]
    line_df_list <- list()
    global_x_range <- range(study_data$X_obs, na.rm = TRUE)

    for (s in strata) {
      mu_x_s <- beta[grep(paste0(ep_x, ".*", facet_by, "\\)", s, "$"), names(beta))]
      mu_y_s <- beta[grep(paste0(ep_y, ".*", facet_by, "\\)", s, "$"), names(beta))]

      if (length(mu_x_s) == 0 || length(mu_y_s) == 0) next
      mu_x_s <- as.numeric(mu_x_s)
      mu_y_s <- as.numeric(mu_y_s)

      if (!is.null(eta_by_stratum) && !is.null(eta_by_stratum$eta_table)) {
        eta_row <- eta_by_stratum$eta_table[eta_by_stratum$eta_table$stratum == s, ]
        if (nrow(eta_row) > 0 && !is.na(eta_row$eta[1])) {
          eta_s <- eta_row$eta[1]
        } else {
          eta_s <- eta
        }
      } else {
        eta_s <- eta
      }

      if (facet_scales == "fixed") {
        x_range_s <- global_x_range
      } else {
        stratum_data <- study_data[study_data[[facet_by]] == s, ]
        x_range_s <- range(stratum_data$X_obs, na.rm = TRUE)
      }
      if (any(is.infinite(x_range_s))) next
      x_seq <- seq(x_range_s[1], x_range_s[2], length.out = 50)

      if (scale == "original" && ep_x %in% names(transform)) {
        x_seq_trans <- slim_transform(
          x_seq / ifelse(ep_x %in% names(display_scale), display_scale[[ep_x]], 1),
          rep(1, length(x_seq)),
          transform[[ep_x]]
        )$y
      } else {
        x_seq_trans <- x_seq / ifelse(ep_x %in% names(display_scale), display_scale[[ep_x]], 1)
      }

      y_line_trans <- mu_y_s + eta_s * (x_seq_trans - mu_x_s)

      if (scale == "original" && ep_y %in% names(transform)) {
        line_data <- slim_backtransform(y_line_trans,
                                        y_line_trans - 1.96 * cond_sd,
                                        y_line_trans + 1.96 * cond_sd,
                                        transform[[ep_y]])
        y_line <- line_data$est
        y_band_lo <- line_data$ci.lb
        y_band_hi <- line_data$ci.ub
      } else {
        y_line <- y_line_trans
        y_band_lo <- y_line_trans - 1.96 * cond_sd
        y_band_hi <- y_line_trans + 1.96 * cond_sd
      }

      if (ep_y %in% names(display_scale) && is.numeric(display_scale[[ep_y]])) {
        y_line    <- y_line * display_scale[[ep_y]]
        y_band_lo <- y_band_lo * display_scale[[ep_y]]
        y_band_hi <- y_band_hi * display_scale[[ep_y]]
      }

      line_df_list[[as.character(s)]] <- data.frame(
        x = x_seq, y = y_line, y_lo = y_band_lo, y_hi = y_band_hi,
        stratum = s
      )
    }
    line_df <- do.call(rbind, line_df_list)
    names(line_df)[names(line_df) == "stratum"] <- facet_by

  } else {
    x_range <- range(study_data$X_obs, na.rm = TRUE)
    x_seq   <- seq(x_range[1], x_range[2], length.out = 100)

    if (scale == "original" && ep_x %in% names(transform)) {
      x_seq_trans <- slim_transform(
        x_seq / ifelse(ep_x %in% names(display_scale), display_scale[[ep_x]], 1),
        rep(1, length(x_seq)),
        transform[[ep_x]]
      )$y
    } else {
      x_seq_trans <- x_seq / ifelse(ep_x %in% names(display_scale), display_scale[[ep_x]], 1)
    }

    y_line_trans <- mu_y + eta * (x_seq_trans - mu_x)

    if (scale == "original" && ep_y %in% names(transform)) {
      line_data <- slim_backtransform(y_line_trans,
                                      y_line_trans - 1.96 * cond_sd,
                                      y_line_trans + 1.96 * cond_sd,
                                      transform[[ep_y]])
      y_line <- line_data$est
      y_band_lo <- line_data$ci.lb
      y_band_hi <- line_data$ci.ub
    } else {
      y_line <- y_line_trans
      y_band_lo <- y_line_trans - 1.96 * cond_sd
      y_band_hi <- y_line_trans + 1.96 * cond_sd
    }

    if (ep_y %in% names(display_scale) && is.numeric(display_scale[[ep_y]])) {
      y_line    <- y_line * display_scale[[ep_y]]
      y_band_lo <- y_band_lo * display_scale[[ep_y]]
      y_band_hi <- y_band_hi * display_scale[[ep_y]]
    }

    line_df <- data.frame(x = x_seq, y = y_line, y_lo = y_band_lo, y_hi = y_band_hi)
  }

  # --- Build plot ---
  if (!is.null(color_by)) {
    p <- ggplot2::ggplot(study_data, ggplot2::aes_string(x = "X_obs", y = "Y_obs",
                                                          color = paste0("factor(", color_by, ")")))
  } else {
    p <- ggplot2::ggplot(study_data, ggplot2::aes(x = X_obs, y = Y_obs))
  }

  if (show_band) {
    p <- p + ggplot2::geom_ribbon(data = line_df, ggplot2::aes(x = x, ymin = y_lo, ymax = y_hi),
                                   inherit.aes = FALSE, fill = "gray80", alpha = 0.5)
  }

  if (show_line) {
    p <- p + ggplot2::geom_line(data = line_df, ggplot2::aes(x = x, y = y),
                                 inherit.aes = FALSE, color = "blue", linewidth = 1)
  }

  p <- p + ggplot2::geom_point(size = 3, alpha = 0.8)

  y_label <- ep_y
  if (y_scale == "log") {
    y_label <- paste0(ep_y, " (log scale)")
  }

  if (!is.null(result$ci_eta) && !any(is.na(result$ci_eta))) {
    subtitle_text <- sprintf("R2 = %.2f, eta = %.3f (95%% CI: %.3f, %.3f)",
                             result$R2_between, eta, result$ci_eta[1], result$ci_eta[2])
  } else {
    subtitle_text <- sprintf("R2 = %.2f, eta = %.3f", result$R2_between, eta)
  }

  p <- p +
    ggplot2::labs(
      title    = title,
      subtitle = subtitle_text,
      x        = ep_x,
      y        = y_label,
      color    = color_by
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "right")

  if (y_scale == "log") {
    p <- p + ggplot2::scale_y_log10()
  }

  if (!is.null(facet_by)) {
    p <- p + ggplot2::facet_wrap(as.formula(paste0("~ ", facet_by)), scales = facet_scales)

    if (facet_scales == "fixed") {
      x_lim <- range(study_data$X_obs, na.rm = TRUE)
      y_lim <- range(c(study_data$Y_obs, line_df$y, line_df$y_lo, line_df$y_hi), na.rm = TRUE)
      p <- p + ggplot2::coord_cartesian(xlim = x_lim, ylim = y_lim)
    }
  }

  if (!is.null(save_to)) {
    ggplot2::ggsave(save_to, p, width = width, height = height)
    message("Plot saved to: ", save_to)
  }

  tbl <- study_data[, c(trial_col, "X_obs", "Y_obs", "Y_model", "Y_CI_lb", "Y_CI_ub", color_by)]

  if (output == "plot") {
    return(list(plot = p, eta = eta))
  } else if (output == "table") {
    return(list(table = tbl, translation = translation, eta = eta))
  } else {
    return(list(plot = p, table = tbl, translation = translation, eta = eta))
  }
}


#' Forest plot of endpoints by study
#'
#' Creates a forest plot showing effect estimates with confidence intervals
#' for each study, faceted by endpoint. Studies are ordered by the overall
#' mean of the first endpoint for consistent comparison across panels.
#'
#' @param data Data frame in stacked format
#' @param trial_col Column name for trial/study identifier
#' @param endpoint_col Column name for endpoint identifier
#' @param y_col Column name for effect estimate
#' @param se_col Column name for standard error
#' @param endpoints Character vector of endpoints to include (default: all)
#' @param strata_col Column name for stratification (optional). If NULL,
#'   aggregates across strata using inverse-variance weighting.
#' @param transform Named list of transformations to apply for display.
#'   E.g., \code{list(RuptureRisk = "log")} shows RuptureRisk on log scale.
#' @param title Plot title
#' @param save_to File path to save plot
#' @param width Plot width in inches
#' @param height Plot height in inches
#'
#' @return A list of ggplot objects (one per endpoint, or grid if stratified)
#'
#' @export
#' @examples
#' data(aaa)
#' # RuptureRisk is stored on the log scale in aaa, but slim_fit() applies the
#' # log transform itself. Supply the rate and its SE on the original scale.
#' aaa_fit <- aaa
#' rup <- aaa_fit$Endpoint == "RuptureRisk"
#' aaa_fit$y[rup]  <- aaa_fit$Mean[rup]
#' aaa_fit$se[rup] <- aaa_fit$Mean[rup] * aaa_fit$se[rup]
#'
#' # Overall (aggregated across strata)
#' slim_forest(
#'   data         = aaa_fit,
#'   trial_col    = "Study",
#'   endpoint_col = "Endpoint",
#'   y_col        = "y",
#'   se_col       = "se",
#'   endpoints    = c("GrowthRate", "RuptureRisk")
#' )
#'
#' # By stratum
#' slim_forest(
#'   data         = aaa_fit,
#'   trial_col    = "Study",
#'   endpoint_col = "Endpoint",
#'   y_col        = "y",
#'   se_col       = "se",
#'   strata_col   = "Diameter"
#' )
slim_forest <- function(data,
                        trial_col,
                        endpoint_col,
                        y_col,
                        se_col,
                        endpoints = NULL,
                        strata_col = NULL,
                        transform = list(),
                        title = NULL,
                        save_to = NULL,
                        width = 10,
                        height = 8) {

  if (is.null(endpoints)) {
    endpoints <- unique(data[[endpoint_col]])
  }

  df <- data[data[[endpoint_col]] %in% endpoints, ]

  # Apply transforms for display
  for (ep in names(transform)) {
    idx <- df[[endpoint_col]] == ep
    if (sum(idx) == 0) next

    if (transform[[ep]] == "log") {
      y_orig <- df[[y_col]][idx]
      se_orig <- df[[se_col]][idx]
      # Transform to log scale
      df[[y_col]][idx] <- log(pmax(y_orig, 1e-10))
      df[[se_col]][idx] <- se_orig / pmax(y_orig, 1e-10)
    }
  }

  # Compute overall mean for first endpoint to get center ordering
  ep1 <- endpoints[1]
  ep1_data <- df[df[[endpoint_col]] == ep1, ]
  center_means <- stats::aggregate(ep1_data[[y_col]],
                                   by = list(ep1_data[[trial_col]]),
                                   FUN = mean, na.rm = TRUE)
  names(center_means) <- c("Center", "mean_ep1")
  center_order <- center_means$Center[order(center_means$mean_ep1)]

  if (is.null(strata_col)) {
    # Aggregate across strata - one row per study per endpoint
    agg_list <- lapply(split(df, list(df[[trial_col]], df[[endpoint_col]]), drop = TRUE), function(sub) {
      if (nrow(sub) == 0) return(NULL)
      yi <- sub[[y_col]]
      se <- sub[[se_col]]
      vi <- se^2
      w <- 1/vi
      w[!is.finite(w)] <- 0
      if (sum(w) == 0) return(NULL)
      yi_pooled <- sum(w * yi) / sum(w)
      se_pooled <- sqrt(1 / sum(w))
      data.frame(
        Study = sub[[trial_col]][1],
        Endpoint = sub[[endpoint_col]][1],
        yi = yi_pooled,
        se = se_pooled,
        stringsAsFactors = FALSE
      )
    })
    plot_df <- do.call(rbind, agg_list)

    plot_df$ci_lo <- plot_df$yi - 1.96 * plot_df$se
    plot_df$ci_hi <- plot_df$yi + 1.96 * plot_df$se
    plot_df$Study <- factor(plot_df$Study, levels = center_order)

    plot_list <- lapply(endpoints, function(ep) {
      ep_df <- plot_df[plot_df$Endpoint == ep, ]

      x_label <- if (ep %in% names(transform) && transform[[ep]] == "log") {
        paste0(ep, " (log scale)")
      } else {
        ep
      }

      ggplot2::ggplot(ep_df, ggplot2::aes(x = yi, y = Study)) +
        ggplot2::geom_point(size = 2) +
        ggplot2::geom_errorbar(ggplot2::aes(xmin = ci_lo, xmax = ci_hi),
                               width = 0.2, orientation = "y") +
        ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
        ggplot2::labs(title = ep, x = "Effect (95% CI)", y = NULL) +
        ggplot2::theme_minimal() +
        ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
                       axis.text.y = ggplot2::element_text(size = 8))
    })

  } else {
    # By stratum
    plot_df <- data.frame(
      Study = df[[trial_col]],
      Endpoint = df[[endpoint_col]],
      Stratum = df[[strata_col]],
      yi = df[[y_col]],
      se = df[[se_col]],
      stringsAsFactors = FALSE
    )
    plot_df$ci_lo <- plot_df$yi - 1.96 * plot_df$se
    plot_df$ci_hi <- plot_df$yi + 1.96 * plot_df$se
    plot_df$Study <- factor(plot_df$Study, levels = center_order)

    strata_levels <- sort(unique(plot_df$Stratum))

    # Build grid: rows = strata, cols = endpoints
    # First row has endpoint names as titles, subsequent rows have stratum labels
    plot_list <- list()
    for (i in seq_along(strata_levels)) {
      stratum <- strata_levels[i]
      for (j in seq_along(endpoints)) {
        ep <- endpoints[j]
        ep_df <- plot_df[plot_df$Endpoint == ep & plot_df$Stratum == stratum, ]

        # First row: show endpoint as title; other rows: show stratum
        if (i == 1) {
          panel_title <- ep
        } else {
          panel_title <- NULL
        }

        # First column: show stratum as y-axis label
        if (j == 1) {
          y_label <- as.character(stratum)
        } else {
          y_label <- NULL
        }

        p <- ggplot2::ggplot(ep_df, ggplot2::aes(x = yi, y = Study)) +
          ggplot2::geom_point(size = 2) +
          ggplot2::geom_errorbar(ggplot2::aes(xmin = ci_lo, xmax = ci_hi),
                                 width = 0.2, orientation = "y") +
          ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
          ggplot2::labs(title = panel_title, x = NULL, y = y_label) +
          ggplot2::theme_minimal() +
          ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = 10),
                         axis.text.y = ggplot2::element_text(size = 6),
                         axis.title.y = ggplot2::element_text(face = "bold", size = 9))

        plot_list <- c(plot_list, list(p))
      }
    }
  }

  # Arrange plots into grid
  if (is.null(strata_col)) {
    g <- gridExtra::arrangeGrob(grobs = plot_list, ncol = length(endpoints), top = title)
  } else {
    strata_levels <- sort(unique(plot_df$Stratum))
    g <- gridExtra::arrangeGrob(grobs = plot_list, ncol = length(endpoints),
                                nrow = length(strata_levels), top = title)
  }

  if (!is.null(save_to)) {
    ggplot2::ggsave(save_to, g, width = width, height = height)
    message("Saved to: ", save_to)
  }

  # Display the plot
  grid::grid.draw(g)

  invisible(list(plot = g, panels = plot_list))
}
