#' Fit a SLIM (Stacked Linear Mixed-Effects) model
#'
#' Fits a multivariate meta-analysis model for correlated endpoints using
#' \code{metafor::rma.mv()}. The SLIM framework enables joint modeling of

#' multiple endpoints across studies with estimation of the between-study
#' covariance structure.
#'
#' @param data Data frame in stacked (long) format with one row per study-endpoint
#' @param formula Formula for fixed effects (passed to mods argument of rma.mv).
#'   Typically \code{~ 0 + Endpoint} or \code{~ 0 + Endpoint:factor(Stratum)}
#' @param random Formula for random effects. Default: \code{~ endpoint | trial}
#' @param trial_col Column name for trial/study/center identifier
#' @param endpoint_col Column name for endpoint identifier
#' @param y_col Column name for effect estimate (mean, proportion, etc.)
#' @param se_col Column name for standard error
#' @param transform Named list of transformations by endpoint.
#'   E.g., \code{list(RuptureRisk = "log")}. Options: "log", "logit", "sqrt"
#' @param display_scale Named list of display scaling by endpoint.
#'   E.g., \code{list(RuptureRisk = 1000)} for per-1000 display
#' @param method Estimation method (default: "REML")
#' @param satterthwaite If TRUE (default), the CIs for the pooled means use
#'   cluster-robust (CR2) variance estimation with Satterthwaite degrees of
#'   freedom from \code{clubSandwich}, clustered on \code{trial_col}. If FALSE,
#'   or if that calculation fails (with a warning), plain Wald z intervals are
#'   used. The CIs are back-transformed and display-scaled like the estimates.
#' @param struct Variance structure of G, passed to \code{metafor::rma.mv()}.
#'   \code{"UN"} (default): unstructured, with free tau11, tau12, tau22, so
#'   eta = tau12/tau11 as in the paper. \code{"CS"}: compound symmetry, which
#'   is metafor's own default; it forces tau11 = tau22 and therefore eta = rho.
#'
#' @return A list with class "slim_fit" containing: `fit` (rma.mv model object),
#'   `pooled` (data frame of pooled estimates with CIs), `G` (between-study
#'   covariance matrix), `R2_between` (between-study R-squared), `eta`
#'   (translation coefficient), `se_eta` (standard error of eta), `ci_eta`
#'   (95 percent CI for eta: delta method when \code{struct = "UN"}, profile CI
#'   for rho when \code{struct = "CS"}), `ci_rho` (95 percent profile CI for
#'   the between-study correlation rho), and metadata columns (`transform`, `display_scale`,
#'   `endpoint_col`, `trial_col`, `y_col`, `se_col`).
#'
#' @details
#' The SLIM model is:
#' \deqn{Y_{ij} = \mu_j + u_{ij} + \varepsilon_{ij}}
#'
#' where:
#' \itemize{
#'   \item \eqn{\mu_j} = population mean for endpoint j
#'   \item \eqn{u_{ij}} = study-level random effect
#'   \item \eqn{\varepsilon_{ij} \sim N(0, V_{ij})} = sampling error (known)
#' }
#'
#' The random effects follow a bivariate normal:
#' \deqn{(u_{i1}, u_{i2})' \sim N(0, G)}
#'
#' where G is the between-study covariance matrix with elements
#' tau11, tau12, tau22.
#'
#' The translation coefficient eta = tau12/tau11 enables prediction of one
#' endpoint from another via E(Y2|Y1) = mu2 + eta * (Y1 - mu1).
#'
#' @references
#' Rabbee N. (2026). Stacked Linear Mixed-Effects Models for Evidence
#' Synthesis with Multiple Correlated Endpoints.
#'
#' Viechtbauer W. (2010). Conducting meta-analyses in R with the metafor
#' package. Journal of Statistical Software, 36(3), 1-48.
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
#' result$eta         # translation coefficient
#' result$ci_eta      # 95% CI for eta
#' result$R2_between  # between-study R-squared
#' result$pooled      # pooled estimates
#' }
slim_fit <- function(
    data,
    formula,
    random        = NULL,
    trial_col     = "trial",
    endpoint_col  = "endpoint",
    y_col         = "y",
    se_col        = "se",
    transform     = list(),
    display_scale = list(),
    method        = "REML",
    struct        = c("UN", "CS"),
    satterthwaite = TRUE
) {

  struct <- match.arg(struct)

  # --- Validate inputs ---
  required_cols <- c(trial_col, endpoint_col, y_col, se_col)
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop("Missing columns: ", paste(missing_cols, collapse = ", "))
  }

  endpoints <- unique(data[[endpoint_col]])
  if (length(endpoints) < 2) {
    stop("Need at least 2 endpoints, found: ", length(endpoints))
  }

  # --- Build default random effects if not specified ---
  if (is.null(random)) {
    random <- as.formula(paste0("~ ", endpoint_col, " | ", trial_col))
  }

  # --- Prepare yi and vi columns ---
  data$yi <- data[[y_col]]
  data$vi <- data[[se_col]]^2

  # --- Apply transformations ---
  for (ep in names(transform)) {
    idx <- data[[endpoint_col]] == ep
    if (sum(idx) == 0) {
      warning("Endpoint '", ep, "' not found in data, skipping transform")
      next
    }

    trans <- slim_transform(
      y    = data[[y_col]][idx],
      se   = data[[se_col]][idx],
      type = transform[[ep]]
    )
    data$yi[idx] <- trans$y
    data$vi[idx] <- trans$se^2
  }

  # --- Check for sparse-data induced dependence ---
  for (ep in names(transform)) {
    if (!(transform[[ep]] %in% c("log", "logit"))) next

    idx <- data[[endpoint_col]] == ep
    yi_ep <- data$yi[idx]
    se_ep <- sqrt(data$vi[idx])

    # Need at least 5 studies to estimate correlation
    if (sum(!is.na(yi_ep) & !is.na(se_ep)) < 5) next

    # Test correlation between |yi| and SE
    cor_test <- try(cor.test(abs(yi_ep), se_ep, method = "pearson"), silent = TRUE)

    if (!inherits(cor_test, "try-error") &&
        !is.na(cor_test$p.value) &&
        cor_test$p.value < 0.05 &&
        cor_test$estimate > 0) {
      warning(
        sprintf(
          "Endpoint '%s' shows sparse-data induced dependence (cor(|yi|, SE) = %.2f, p = %.3f).\n",
          ep, cor_test$estimate, cor_test$p.value
        ),
        sprintf(
          "Consider: data_adj <- slim_prior_adjust(data, y_col='%s', se_col='%s', prior_sd=1.0)\n",
          y_col, se_col
        ),
        "Then re-run slim_fit() with the adjusted columns.",
        call. = FALSE
      )
    }
  }

  # --- Fit rma.mv model ---
  fit <- metafor::rma.mv(
    yi     = yi,
    V      = diag(data$vi),
    mods   = formula,
    random = random,
    struct = struct,
    data   = data,
    method = method,
    cvvc   = "varcor"
  )

  # --- Extract G matrix ---
  G <- fit$G

  # --- Compute R2_between, eta, and CI(eta) ---
  if (nrow(G) >= 2 && ncol(G) >= 2) {
    tau11 <- G[1, 1]
    tau22 <- G[2, 2]
    tau12 <- G[1, 2]
    R2_between <- (tau12^2) / (tau11 * tau22)

    # eta = tau12 / tau11
    eta <- tau12 / tau11

    ci_rho <- slim_rho_ci(fit)

    if (struct == "UN") {
      # Delta-method CI for eta = tau12 / tau11
      d      <- slim_eta_delta(fit, x = 1)
      se_eta <- d$se
      ci_eta <- d$ci
    } else {
      # CS: tau11 = tau22, so eta = rho and the CI for rho is the CI for eta
      ci_eta <- ci_rho
      se_eta <- if (any(is.na(ci_eta))) NA else (ci_eta[2] - ci_eta[1]) / (2 * qnorm(0.975))
    }
  } else {
    R2_between <- NA
    eta <- NA
    se_eta <- NA
    ci_eta <- c(NA, NA)
    ci_rho <- c(NA, NA)
  }

  # --- Extract pooled estimates ---
  beta    <- coef(fit)
  se_beta <- sqrt(diag(vcov(fit)))
  z       <- qnorm(0.975)

  ci_lb  <- as.numeric(beta - z * se_beta)
  ci_ub  <- as.numeric(beta + z * se_beta)
  df_mu  <- rep(NA_real_, length(beta))

  if (satterthwaite) {
    # CR2 cluster-robust variance + Satterthwaite d.f. (small-K correction)
    ok <- !is.na(data$yi) & !is.na(data$vi)
    cs <- try({
      V_cr <- clubSandwich::vcovCR(fit, cluster = data[[trial_col]][ok], type = "CR2")
      as.data.frame(clubSandwich::conf_int(fit, vcov = V_cr, test = "Satterthwaite"))
    }, silent = TRUE)
    if (inherits(cs, "try-error") || nrow(cs) != length(beta)) {
      warning("Satterthwaite CIs for the pooled means failed; using Wald z intervals.",
              call. = FALSE)
    } else {
      # conf_int columns: Coef, Estimate, SE, d.f., lower, upper
      ci_lb <- as.numeric(cs[[5]])
      ci_ub <- as.numeric(cs[[6]])
      df_mu <- as.numeric(cs[[4]])
      if (any(df_mu < 4, na.rm = TRUE)) {
        warning("Satterthwaite d.f. < 4 for some pooled means; those CIs are unreliable.",
                call. = FALSE)
      }
    }
  }

  pooled <- data.frame(
    term  = names(beta),
    est   = as.numeric(beta),
    ci.lb = ci_lb,
    ci.ub = ci_ub,
    df    = df_mu,
    stringsAsFactors = FALSE
  )

  # --- Back-transform and apply display scale ---
  for (ep in endpoints) {
    idx <- grepl(ep, pooled$term, fixed = TRUE)
    if (sum(idx) == 0) next

    if (ep %in% names(transform)) {
      bt <- slim_backtransform(
        pooled$est[idx],
        pooled$ci.lb[idx],
        pooled$ci.ub[idx],
        type = transform[[ep]]
      )
      pooled$est[idx]   <- bt$est
      pooled$ci.lb[idx] <- bt$ci.lb
      pooled$ci.ub[idx] <- bt$ci.ub
    }

    if (ep %in% names(display_scale)) {
      ds <- apply_display_scale(
        pooled$est[idx],
        pooled$ci.lb[idx],
        pooled$ci.ub[idx],
        spec = display_scale[[ep]]
      )
      pooled$est[idx]   <- ds$est
      pooled$ci.lb[idx] <- ds$ci.lb
      pooled$ci.ub[idx] <- ds$ci.ub
    }
  }

  # --- Return result ---
  result <- list(
    fit           = fit,
    pooled        = pooled,
    G             = G,
    R2_between    = R2_between,
    eta           = eta,
    se_eta        = se_eta,
    ci_eta        = ci_eta,
    ci_rho        = ci_rho,
    transform     = transform,
    display_scale = display_scale,
    endpoint_col  = endpoint_col,
    trial_col     = trial_col,
    y_col         = y_col,
    se_col        = se_col
  )

  class(result) <- c("slim_fit", "list")
  result
}


#' @export
print.slim_fit <- function(x, ...) {
  cat("SLIM Model Fit\n")
  cat("--------------\n")
  cat("R-squared (between):", round(x$R2_between, 3), "\n")
  cat("Eta (translation coefficient):", round(x$eta, 3))
  if (!any(is.na(x$ci_eta))) {
    cat(" (95% CI:", round(x$ci_eta[1], 3), ",", round(x$ci_eta[2], 3), ")\n")
  } else {
    cat("\n")
  }
  cat("\nG matrix (between-study covariance):\n")
  print(round(x$G, 4))
  cat("\nPooled estimates:\n")
  print(x$pooled, row.names = FALSE)
  invisible(x)
}


#' Profile-likelihood CI for the between-study correlation rho
#'
#' Internal helper. Looks through \code{confint(fit)} for the component whose
#' name starts with "rho" (its position depends on \code{struct}).
#'
#' @param fit An \code{rma.mv} object.
#' @return Numeric vector \code{c(lower, upper)}, or \code{c(NA, NA)}.
#' @noRd
slim_rho_ci <- function(fit) {
  ci_fit <- try(confint(fit), silent = TRUE)
  if (inherits(ci_fit, "try-error") || !is.list(ci_fit)) return(c(NA_real_, NA_real_))
  for (comp in ci_fit) {
    r <- try(comp$random, silent = TRUE)
    if (is.matrix(r) || is.data.frame(r)) {
      hit <- grep("^rho", rownames(r))
      if (length(hit) == 1) return(c(r[hit, "ci.lb"], r[hit, "ci.ub"]))
    }
  }
  c(NA_real_, NA_real_)
}


#' Delta-method CI for eta = tau_xy / tau_xx
#'
#' Internal helper for an unstructured fit. With theta = (tau11, tau22, rho)
#' and V = \code{fit$vvc} (requested with \code{cvvc = "varcor"}), and x the
#' index (1 or 2) of the endpoint that plays the role of Y1,
#' eta = rho * sqrt(tau_yy / tau_xx) and Var(eta) is approximately g' V g.
#' Returns NA with a warning when the fit is on the boundary, V is not finite
#' and positive definite, or the model does not have exactly two endpoints.
#'
#' @param fit An \code{rma.mv} object fitted with \code{struct = "UN"}.
#' @param x Index (1 or 2) of the predictor endpoint.
#' @param level Confidence level.
#' @return List with \code{se} and \code{ci}.
#' @noRd
slim_eta_delta <- function(fit, x = 1, level = 0.95) {

  na_out <- list(se = NA_real_, ci = c(NA_real_, NA_real_))
  warn   <- function(msg) {
    warning("Delta-method CI for eta unavailable: ", msg, call. = FALSE)
    na_out
  }

  tau2 <- fit$tau2
  rho  <- fit$rho
  V    <- fit$vvc

  if (length(tau2) != 2 || length(rho) != 1) {
    return(warn("needs exactly 2 endpoints with struct = 'UN'."))
  }
  if (is.null(V) || !is.matrix(V) || any(dim(V) != 3) || any(!is.finite(V))) {
    return(warn("covariance matrix of the variance components is missing or not finite."))
  }
  if (min(eigen(V, symmetric = TRUE, only.values = TRUE)$values) <= 0) {
    return(warn("covariance matrix of the variance components is not positive definite (fit likely on the boundary)."))
  }
  if (min(tau2) <= 1e-8 || abs(rho) > 0.999) {
    return(warn("estimate is on the boundary (a tau^2 near 0 or |rho| near 1)."))
  }

  tx  <- tau2[x]
  ty  <- tau2[3 - x]
  eta <- rho * sqrt(ty / tx)

  # V is assumed ordered (tau2[1], tau2[2], rho)
  g_tx <- -eta / (2 * tx)
  g_ty <-  eta / (2 * ty)
  g    <- if (x == 1) c(g_tx, g_ty, sqrt(ty / tx)) else c(g_ty, g_tx, sqrt(ty / tx))

  se <- sqrt(drop(t(g) %*% V %*% g))
  z  <- qnorm(1 - (1 - level) / 2)

  list(se = se, ci = eta + c(-1, 1) * z * se)
}
