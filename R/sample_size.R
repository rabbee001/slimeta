#' Calculate sample size for a trial based on SLIM estimates
#'
#' Uses pooled estimates from \code{\link{slim_fit}} to calculate required
#' sample size for a trial. Supports stratified designs with enrollment
#' weights and multi-arm trials with unequal allocation.
#'
#' @param result Output from \code{\link{slim_fit}}
#' @param endpoint Primary endpoint to power on
#' @param effect_size Treatment effect: numeric (absolute delta_mu) or string like "-20 percent" (relative)
#' @param sd Standard deviation of the endpoint. If NULL, estimated from G matrix.
#' @param strata_col Column name for stratification (NULL = unstratified)
#' @param weights Enrollment weights by stratum (NULL = equal weights). Must sum to 1.
#' @param adjustments Named list of adjustments to baseline, e.g., list(statin = -0.63)
#' @param power Target power (default: 0.90)
#' @param alpha Significance level (default: 0.05, two-sided)
#' @param n_arms Number of arms (default: 2)
#' @param arm_props Proportion of patients in each arm (default: equal). Must sum to 1.
#'
#' @return A list containing: `n_per_arm` (required N per arm), `n_total`
#'   (total sample size), `n_arms`, `arm_props`, `effect_size_input`,
#'   `effect_size_weighted` (weighted delta_mu), `sd_pooled`, `mu_baseline`
#'   (if unstratified), `power`, `alpha`, `strata_table` (if stratified),
#'   `adjustments`, and `translation` (implied effect on other endpoint via eta).
#'
#' @details
#' The sample size formula for a two-sample comparison is:
#' \deqn{N = var\_factor \times (z_\alpha + z_\beta)^2 \times SD^2 / \Delta\mu^2}
#'
#' where var_factor = 1/p1 + 1/p2 for arm proportions p1, p2.
#'
#' For stratified designs, the weighted treatment effect is:
#' \deqn{\Delta\mu_{weighted} = \sum_d w_d \times \Delta\mu_d}
#'
#' The translation via eta provides the implied effect on the secondary endpoint:
#' \deqn{\Delta\mu_2 = \eta \times \Delta\mu_1}
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
#' # 2-arm, 1:1 allocation
#' ss <- slim_sample_size(
#'   result       = result,
#'   endpoint     = "GrowthRate",
#'   effect_size  = -0.67,              # absolute
#'   sd           = 3.5,
#'   strata_col   = "Diameter",
#'   weights      = c(0.05, 0.45, 0.45, 0.03, 0.02),
#'   adjustments  = list(statin = -0.63),
#'   power        = 0.90,
#'   alpha        = 0.05
#' )
#'
#' ss$n_total        # total required
#' ss$translation    # implied effect on other endpoint
#'
#' # 2-arm, 2:1 allocation (treatment:control)
#' ss2 <- slim_sample_size(
#'   result       = result,
#'   endpoint     = "GrowthRate",
#'   effect_size  = -0.67,
#'   sd           = 3.5,
#'   strata_col   = "Diameter",
#'   weights      = c(0.05, 0.45, 0.45, 0.03, 0.02),
#'   arm_props    = c(1/3, 2/3)         # control, treatment
#' )
#' ss2$n_total
#' }
slim_sample_size <- function(
    result,
    endpoint,
    effect_size,
    sd              = NULL,
    strata_col      = NULL,
    weights         = NULL,
    adjustments     = list(),
    power           = 0.90,
    alpha           = 0.05,
    n_arms          = 2,
    arm_props       = NULL
) {

  pooled        <- result$pooled
  G             <- result$G
  eta           <- result$eta
  endpoint_col  <- result$endpoint_col
  transform     <- result$transform
  display_scale <- result$display_scale

  # --- Get endpoint estimates ---
  ep_rows <- grepl(endpoint, pooled$term, fixed = TRUE)
  if (sum(ep_rows) == 0) {
    stop("Endpoint '", endpoint, "' not found in pooled estimates")
  }

  ep_data <- pooled[ep_rows, ]

  # --- Parse effect size ---
  if (is.character(effect_size) && grepl("%$", effect_size)) {
    pct <- as.numeric(gsub("%", "", effect_size)) / 100
    relative_effect <- TRUE
  } else {
    delta_mu <- as.numeric(effect_size)
    relative_effect <- FALSE
  }

  # --- Stratified vs unstratified ---
  if (!is.null(strata_col)) {
    ep_data$stratum <- NA
    for (i in seq_len(nrow(ep_data))) {
      term <- ep_data$term[i]
      pattern <- paste0(".*factor\\(", strata_col, "\\)")
      stratum_val <- gsub(pattern, "", term)
      ep_data$stratum[i] <- stratum_val
    }

    strata <- ep_data$stratum
    k <- length(strata)

    if (is.null(weights)) {
      weights <- rep(1/k, k)
    } else {
      if (length(weights) != k) {
        stop("Length of weights (", length(weights), ") must equal number of strata (", k, ")")
      }
      weights <- weights / sum(weights)
    }

    mu_baseline <- ep_data$est

    for (adj_name in names(adjustments)) {
      mu_baseline <- mu_baseline + adjustments[[adj_name]]
    }

    if (relative_effect) {
      delta_mu_strata <- mu_baseline * pct
    } else {
      delta_mu_strata <- rep(delta_mu, k)
    }

    delta_mu_weighted <- sum(weights * delta_mu_strata)

    strata_table <- data.frame(
      stratum      = strata,
      weight       = weights,
      mu_baseline  = mu_baseline,
      delta_mu     = delta_mu_strata,
      mu_treatment = mu_baseline + delta_mu_strata,
      stringsAsFactors = FALSE
    )

  } else {
    mu_baseline <- mean(ep_data$est)

    for (adj_name in names(adjustments)) {
      mu_baseline <- mu_baseline + adjustments[[adj_name]]
    }

    if (relative_effect) {
      delta_mu_weighted <- mu_baseline * pct
    } else {
      delta_mu_weighted <- delta_mu
    }

    strata_table <- NULL
  }

  # --- Standard deviation ---
  if (is.null(sd)) {
    G_names <- rownames(G)
    ep_idx <- which(G_names == endpoint)
    if (length(ep_idx) > 0) {
      tau_ep <- G[ep_idx, ep_idx]
      sd_pooled <- sqrt(tau_ep)
    } else {
      stop("Cannot estimate SD: endpoint not found in G matrix")
    }
  } else {
    sd_pooled <- sd
  }

  # --- Arm proportions ---
  if (is.null(arm_props)) {
    arm_props <- rep(1/n_arms, n_arms)
  } else {
    if (length(arm_props) != n_arms) {
      stop("Length of arm_props (", length(arm_props), ") must equal n_arms (", n_arms, ")")
    }
    if (abs(sum(arm_props) - 1) > 1e-6) {
      stop("arm_props must sum to 1, got ", sum(arm_props))
    }
  }

  # --- Sample size formula ---
  z_alpha <- qnorm(1 - alpha/2)
  z_beta  <- qnorm(power)

  if (n_arms == 2) {
    var_factor <- 1/arm_props[1] + 1/arm_props[2]
  } else {
    var_factor <- 1/arm_props[1] + 1/arm_props[2]
  }

  n_total <- var_factor * (z_alpha + z_beta)^2 * sd_pooled^2 / delta_mu_weighted^2
  n_total <- ceiling(n_total)

  n_per_arm <- ceiling(n_total * arm_props)

  # --- Translation to other endpoint ---
  translation <- NULL
  G_names <- rownames(G)
  other_endpoints <- setdiff(G_names, endpoint)

  if (length(other_endpoints) == 1 && !is.na(eta)) {
    other_ep <- other_endpoints[1]
    delta_mu_other <- eta * delta_mu_weighted

    if (other_ep %in% names(transform) && transform[[other_ep]] == "log") {
      # delta_mu_other is a change in log(endpoint): report it on the log
      # scale and as a percent change. display_scale (e.g. x1000) applies to
      # levels, not to a log-scale change, so it is not applied here.
      pct_change_other       <- (exp(delta_mu_other) - 1) * 100
      delta_mu_other_display <- delta_mu_other
    } else {
      pct_change_other <- NA
      if (other_ep %in% names(display_scale) && is.numeric(display_scale[[other_ep]])) {
        delta_mu_other_display <- delta_mu_other * display_scale[[other_ep]]
      } else {
        delta_mu_other_display <- delta_mu_other
      }
    }

    translation <- list(
      endpoint        = other_ep,
      eta             = eta,
      implied_effect  = delta_mu_other_display,
      pct_change      = pct_change_other
    )
  }

  list(
    n_per_arm          = n_per_arm,
    n_total            = sum(n_per_arm),
    n_arms             = n_arms,
    arm_props          = arm_props,
    effect_size_input  = effect_size,
    effect_size_weighted = delta_mu_weighted,
    sd_pooled          = sd_pooled,
    mu_baseline        = if(is.null(strata_table)) mu_baseline else NULL,
    power              = power,
    alpha              = alpha,
    strata_table       = strata_table,
    adjustments        = adjustments,
    translation        = translation
  )
}
