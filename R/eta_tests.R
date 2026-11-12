#' Test whether eta varies by stratum (Cochran's Q)
#'
#' Fits separate SLIM models per stratum, extracts stratum-specific eta with CI,
#' and performs Cochran's Q test for heterogeneity.
#'
#' @param data Data frame in stacked format
#' @param strata_col Column name for stratification variable (e.g., "Diameter")
#' @param endpoints Character vector of length 2: c("EP1", "EP2")
#' @param trial_col Column name for trial/study/center identifier
#' @param endpoint_col Column name for endpoint identifier
#' @param y_col Column name for effect estimate
#' @param se_col Column name for standard error
#' @param transform Named list of transformations
#' @param display_scale Named list of display scaling
#' @param struct Variance structure passed to \code{\link{slim_fit}}:
#'   \code{"UN"} (default) or \code{"CS"}.
#'
#' @return A list with `eta_table` (data frame with stratum, eta, se_eta,
#'   ci_lb, ci_ub), `Q_test` (list with Q, df, p_value for Cochran's Q),
#'   `pooled_eta` (inverse-variance weighted pooled eta), and `se_pooled`
#'   (standard error of pooled eta).
#'
#' @details
#' Cochran's Q tests H0: eta is constant across all strata. The test statistic
#' is Q = sum(w_d * (eta_d - eta_pooled)^2), which follows chi-squared(k-1)
#' under H0, where k is the number of strata.
#'
#' Limitation: Assumes eta estimates are independent. When same centers appear
#' in multiple strata, this assumption is violated. Use \code{\link{slim_eta_dunnett}}
#' for a more powerful test that accounts for this correlation.
#'
#' @seealso \code{\link{slim_eta_dunnett}} for a more powerful test
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
#' eta_test <- slim_eta_by_stratum(
#'   data          = aaa_fit,
#'   strata_col    = "Diameter",
#'   endpoints     = c("GrowthRate", "RuptureRisk"),
#'   trial_col     = "Study",
#'   endpoint_col  = "Endpoint",
#'   y_col         = "y",
#'   se_col        = "se",
#'   transform     = list(RuptureRisk = "log")
#' )
#'
#' eta_test$eta_table   # stratum-specific eta
#' eta_test$Q_test      # Cochran's Q test of a common eta
#' }
slim_eta_by_stratum <- function(
    data,
    strata_col,
    endpoints,
    trial_col     = "trial",
    endpoint_col  = "endpoint",
    y_col         = "y",
    se_col        = "se",
    transform     = list(),
    display_scale = list(),
    struct        = c("UN", "CS")
) {

  struct <- match.arg(struct)
  strata <- sort(unique(data[[strata_col]]))
  ep_x <- endpoints[1]
  ep_y <- endpoints[2]

  eta_list <- list()

  for (s in strata) {
    data_s <- data[data[[strata_col]] == s, ]

    eps_present <- unique(data_s[[endpoint_col]])
    if (!all(endpoints %in% eps_present)) {
      warning("Stratum ", s, " missing one or both endpoints, skipping")
      next
    }

    fit_s <- try(slim_fit(
      data          = data_s,
      formula       = as.formula(paste0("~ 0 + ", endpoint_col)),
      random        = as.formula(paste0("~ ", endpoint_col, " | ", trial_col)),
      trial_col     = trial_col,
      endpoint_col  = endpoint_col,
      y_col         = y_col,
      se_col        = se_col,
      transform     = transform,
      display_scale = display_scale,
      struct        = struct,
      satterthwaite = FALSE
    ), silent = TRUE)

    if (inherits(fit_s, "try-error")) {
      warning("Stratum ", s, " failed to fit, skipping")
      next
    }

    G <- fit_s$G
    G_names <- rownames(G)
    idx_x <- which(G_names == ep_x)
    idx_y <- which(G_names == ep_y)

    if (length(idx_x) == 0 || length(idx_y) == 0) {
      warning("Stratum ", s, ": endpoints not found in G matrix, skipping")
      next
    }

    tau_xx <- G[idx_x, idx_x]
    tau_xy <- G[idx_x, idx_y]
    eta_s <- tau_xy / tau_xx

    if (struct == "UN") {
      # Delta-method CI for eta = tau_xy / tau_xx
      d        <- slim_eta_delta(fit_s$fit, x = idx_x)
      se_eta_s <- d$se
      ci_eta_s <- d$ci
    } else {
      # CS: eta = rho, no rescaling needed
      se_eta_s <- fit_s$se_eta
      ci_eta_s <- fit_s$ci_eta
    }

    eta_list[[as.character(s)]] <- data.frame(
      stratum = s,
      eta     = eta_s,
      se_eta  = se_eta_s,
      ci_lb   = ci_eta_s[1],
      ci_ub   = ci_eta_s[2],
      stringsAsFactors = FALSE
    )
  }

  if (length(eta_list) == 0) {
    stop("No strata could be fit")
  }

  eta_table <- do.call(rbind, eta_list)
  rownames(eta_table) <- NULL

  # Cochran's Q test
  valid_idx <- !is.na(eta_table$se_eta) & eta_table$se_eta > 0
  if (any(!valid_idx)) {
    message("Strata dropped from Cochran's Q (no usable SE for eta): ",
            paste(eta_table$stratum[!valid_idx], collapse = ", "))
  }
  if (sum(valid_idx) >= 2) {
    eta_valid <- eta_table$eta[valid_idx]
    se_valid <- eta_table$se_eta[valid_idx]
    w <- 1 / se_valid^2

    pooled_eta <- sum(w * eta_valid) / sum(w)
    se_pooled <- sqrt(1 / sum(w))

    Q <- sum(w * (eta_valid - pooled_eta)^2)
    df <- sum(valid_idx) - 1
    p_value <- 1 - pchisq(Q, df = df)

    Q_test <- list(
      Q       = Q,
      df      = df,
      p_value = p_value
    )
  } else {
    pooled_eta <- eta_table$eta[1]
    se_pooled <- NA
    Q_test <- list(Q = NA, df = NA, p_value = NA)
  }

  list(
    eta_table  = eta_table,
    Q_test     = Q_test,
    pooled_eta = pooled_eta,
    se_pooled  = se_pooled
  )
}


#' Dunnett test for eta heterogeneity across strata
#'
#' A more powerful alternative to Cochran's Q that accounts for correlation
#' between stratum-specific eta estimates due to shared centers.
#'
#' @param result Output from \code{\link{slim_fit}} (global model)
#' @param eta_by_stratum Output from \code{\link{slim_eta_by_stratum}}
#' @param data Original data frame
#' @param strata_col Column name for stratification variable
#' @param trial_col Column name for trial/study/center identifier
#' @param endpoint_col Column name for endpoint identifier
#'
#' @return A list with `test_stat` (Dunnett max-abs-t), `p_value` (joint p-value),
#'   `df` (degrees of freedom = n_centers - 3), `n_complete` (complete centers),
#'   `n_centers`, `n_strata`, `contrasts` (data frame with eta_d - eta_global,
#'   SE, t-stat, adjusted p-value), `corr_matrix`, `qc_pass` (logical), and
#'   `qc_message` (status message).
#'
#' @details
#' The Dunnett test computes contrasts delta_d = eta_d - eta_global and tests
#' whether any delta significantly differs from zero using the multivariate t
#' distribution with the correlation structure induced by shared centers.
#'
#' QC requirement: n_complete > k (number of complete centers must exceed
#' number of strata). A warning is issued if n_complete < 2k due to low power.
#'
#' The variance of each contrast is:
#' Var(delta_d) = Var(eta_d) + Var(eta_global) - 2*w_d*Var(eta_d)
#'
#' where w_d is the inverse-variance weight of stratum d in the pooled estimate.
#'
#' @references
#' Dunnett CW (1955). A Multiple Comparison Procedure for Comparing Several
#' Treatments with a Control. Journal of the American Statistical Association.
#'
#' @seealso \code{\link{slim_eta_by_stratum}} for Cochran's Q test
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
#' eta_test <- slim_eta_by_stratum(
#'   data          = aaa_fit,
#'   strata_col    = "Diameter",
#'   endpoints     = c("GrowthRate", "RuptureRisk"),
#'   trial_col     = "Study",
#'   endpoint_col  = "Endpoint",
#'   y_col         = "y",
#'   se_col        = "se",
#'   transform     = list(RuptureRisk = "log")
#' )
#'
#' dunnett <- slim_eta_dunnett(
#'   result         = result,
#'   eta_by_stratum = eta_test,
#'   data           = aaa_fit,
#'   strata_col     = "Diameter",
#'   trial_col      = "Study",
#'   endpoint_col   = "Endpoint"
#' )
#'
#' dunnett$p_value    # p > 0.05 -> common eta appropriate
#' dunnett$contrasts  # see which strata differ from global
#' }
slim_eta_dunnett <- function(
    result,
    eta_by_stratum,
    data,
    strata_col,
    trial_col    = "trial",
    endpoint_col = "endpoint"
) {

  if (!requireNamespace("mvtnorm", quietly = TRUE)) {
    stop("Package 'mvtnorm' required for Dunnett test. Install with: install.packages('mvtnorm')")
  }

  eta_table <- eta_by_stratum$eta_table
  eta_global <- result$eta
  se_global <- result$se_eta
  strata <- eta_table$stratum
  k <- length(strata)

  # --- QC: Count complete centers ---
  endpoints <- unique(data[[endpoint_col]])
  centers <- unique(data[[trial_col]])
  n_centers <- length(centers)

  complete_count <- 0
  for (ctr in centers) {
    ctr_data <- data[data[[trial_col]] == ctr, ]
    is_complete <- TRUE
    for (s in strata) {
      for (ep in endpoints) {
        has_obs <- any(ctr_data[[strata_col]] == s & ctr_data[[endpoint_col]] == ep &
                       !is.na(ctr_data[[result$y_col]]))
        if (!has_obs) {
          is_complete <- FALSE
          break
        }
      }
      if (!is_complete) break
    }
    if (is_complete) complete_count <- complete_count + 1
  }

  n_complete <- complete_count

  if (n_complete <= k) {
    return(list(
      test_stat  = NA,
      p_value    = NA,
      df         = NA,
      n_complete = n_complete,
      n_strata   = k,
      contrasts  = NULL,
      qc_pass    = FALSE,
      qc_message = sprintf("QC FAIL: Need > %d complete centers, have %d", k, n_complete)
    ))
  }

  qc_warn <- ""
  if (n_complete < 2 * k) {
    qc_warn <- sprintf(" (Warning: low power with only %d complete centers for %d strata)",
                       n_complete, k)
  }

  # --- Build contrasts ---
  eta_d <- eta_table$eta
  se_d <- eta_table$se_eta

  valid <- !is.na(se_d) & se_d > 0 & !is.na(eta_d)
  if (sum(valid) < 2) {
    return(list(
      test_stat  = NA,
      p_value    = NA,
      df         = NA,
      n_complete = n_complete,
      n_strata   = k,
      contrasts  = NULL,
      qc_pass    = FALSE,
      qc_message = "QC FAIL: Fewer than 2 strata with valid eta estimates"
    ))
  }

  eta_d <- eta_d[valid]
  se_d <- se_d[valid]
  strata_valid <- strata[valid]
  k_valid <- length(eta_d)

  w_d <- 1 / se_d^2
  w_d <- w_d / sum(w_d)

  delta <- eta_d - eta_global

  var_delta <- se_d^2 + se_global^2 - 2 * w_d * se_d^2
  se_delta <- sqrt(pmax(var_delta, 0))

  # --- Correlation matrix ---
  R <- matrix(0, k_valid, k_valid)
  for (i in 1:k_valid) {
    for (j in 1:k_valid) {
      if (i == j) {
        R[i, j] <- 1
      } else {
        cov_ij <- se_global^2 - w_d[i] * se_d[i]^2 - w_d[j] * se_d[j]^2
        R[i, j] <- cov_ij / (se_delta[i] * se_delta[j])
      }
    }
  }

  R <- (R + t(R)) / 2
  diag(R) <- 1
  R[R > 1] <- 1
  R[R < -1] <- -1

  df <- n_centers - 3
  df <- max(df, 1)

  t_stats <- delta / se_delta

  p_marginal <- 2 * pt(-abs(t_stats), df = df)

  max_t <- max(abs(t_stats))

  p_joint <- 1 - mvtnorm::pmvt(
    lower = rep(-max_t, k_valid),
    upper = rep(max_t, k_valid),
    delta = rep(0, k_valid),
    corr = R,
    df = df
  )[1]

  p_adjusted <- numeric(k_valid)
  for (i in 1:k_valid) {
    t_i <- abs(t_stats[i])
    p_adjusted[i] <- 1 - mvtnorm::pmvt(
      lower = rep(-t_i, k_valid),
      upper = rep(t_i, k_valid),
      delta = rep(0, k_valid),
      corr = R,
      df = df
    )[1]
  }

  contrasts_df <- data.frame(
    stratum     = strata_valid,
    eta_d       = eta_d,
    eta_global  = eta_global,
    difference  = delta,
    se_diff     = se_delta,
    t_stat      = t_stats,
    p_marginal  = p_marginal,
    p_adjusted  = p_adjusted,
    stringsAsFactors = FALSE
  )

  list(
    test_stat  = max_t,
    p_value    = p_joint,
    df         = df,
    n_complete = n_complete,
    n_centers  = n_centers,
    n_strata   = k_valid,
    contrasts  = contrasts_df,
    corr_matrix = R,
    qc_pass    = TRUE,
    qc_message = paste0("QC PASS: ", n_complete, " complete centers, ", k_valid, " strata", qc_warn)
  )
}
