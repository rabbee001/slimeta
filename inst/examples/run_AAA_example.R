# =============================================================================
# run_AAA_example.R - Complete SLIMETA Analysis for AAA Data
# =============================================================================

library(slimeta)

data(aaa)

# RuptureRisk is stored on the log scale in aaa, but slim_fit() applies the
# log transform itself. Supply the rate and its SE on the original scale.
aaa_data <- aaa
rup <- aaa_data$Endpoint == "RuptureRisk"
aaa_data$y[rup]  <- aaa_data$Mean[rup]
aaa_data$se[rup] <- aaa_data$Mean[rup] * aaa_data$se[rup]

# =============================================================================
# 1. Fit SLIM Model
# =============================================================================

result <- slim_fit(
  data          = aaa_data,
  formula       = ~ 0 + Endpoint:factor(Diameter),
  random        = ~ Endpoint | Study,
  trial_col     = "Study",
  endpoint_col  = "Endpoint",
  y_col         = "y",
  se_col        = "se",
  transform     = list(RuptureRisk = "log"),
  display_scale = list(RuptureRisk = 1000)
)

# Access results
result$pooled       # estimates with CIs
result$G            # covariance matrix
result$R2_between   # between-study R-squared
result$eta          # translation coefficient
result$se_eta       # SE of eta
result$ci_eta       # 95% CI for eta
result$fit          # rma.mv object for diagnostics

cat("\n=== Global Model Results ===\n")
cat(sprintf("R² = %.3f\n", result$R2_between))
cat(sprintf("η = %.3f (95%% CI: %.3f, %.3f)\n",
            result$eta, result$ci_eta[1], result$ci_eta[2]))


# =============================================================================
# 2. Forest Plot - Overall
# =============================================================================

slim_forest(
  data         = aaa_data,
  trial_col    = "Study",
  endpoint_col = "Endpoint",
  y_col        = "y",
  se_col       = "se",
  endpoints    = c("GrowthRate", "RuptureRisk"),
  title        = "AAA Forest Plot - Overall"
)


# =============================================================================
# 3. Forest Plot - By Diameter
# =============================================================================

slim_forest(
  data         = aaa_data,
  trial_col    = "Study",
  endpoint_col = "Endpoint",
  y_col        = "y",
  se_col       = "se",
  endpoints    = c("GrowthRate", "RuptureRisk"),
  strata_col   = "Diameter",
  title        = "AAA Forest Plot - By Diameter"
)


# =============================================================================
# 4. Panel Plot by Diameter
# =============================================================================

out1 <- slim2EP_plot_by_covariate(
  result,
  data      = aaa_data,
  endpoints = c("GrowthRate", "RuptureRisk"),
  covariate = "Diameter",
  title     = "AAA: Model vs Observed by Diameter"
)

out1$plot   # view
out1$table  # data


# =============================================================================
# 5. Scatter Plot by Study (faceted by Diameter)
# =============================================================================

out2 <- slim2EP_plot_by_study(
  result,
  data         = aaa_data,
  endpoints    = c("GrowthRate", "RuptureRisk"),
  title        = "AAA Center Data: Growth vs Rupture",
  facet_by     = "Diameter",
  facet_scales = "fixed",
  y_scale      = "log"
)

out2$plot         # scatter with E[Y|X] curve
out2$table        # observed vs model per study
out2$translation  # -10%, -20%, -30% table
out2$eta          # translation coefficient


# =============================================================================
# 6. Test for Stratum-Specific Eta (Cochran's Q)
# =============================================================================

eta_test <- slim_eta_by_stratum(
  data          = aaa_data,
  strata_col    = "Diameter",
  endpoints     = c("GrowthRate", "RuptureRisk"),
  trial_col     = "Study",
  endpoint_col  = "Endpoint",
  y_col         = "y",
  se_col        = "se",
  transform     = list(RuptureRisk = "log"),
  display_scale = list(RuptureRisk = 1000)
)

# View stratum-specific eta with CI
eta_test$eta_table

# Cochran's Q test for heterogeneity
# H0: eta is constant across strata
# p > 0.05 --> common eta sufficient
# p < 0.05 --> stratum-specific eta warranted
cat("\n=== Cochran's Q Test for Eta Heterogeneity ===\n")
cat(sprintf("Q = %.3f, df = %d, p = %.4f\n",
            eta_test$Q_test$Q, eta_test$Q_test$df, eta_test$Q_test$p_value))

if (eta_test$Q_test$p_value > 0.05) {
  cat("Conclusion: No significant heterogeneity. Common eta is appropriate.\n")
} else {
  cat("Conclusion: Significant heterogeneity. Consider stratum-specific eta.\n")
}

# Pooled eta (inverse-variance weighted)
cat(sprintf("\nPooled η = %.3f\n", eta_test$pooled_eta))


# =============================================================================
# 7. Dunnett Test for Eta Heterogeneity (more powerful)
# =============================================================================

dunnett <- slim_eta_dunnett(
  result         = result,
  eta_by_stratum = eta_test,
  data           = aaa_data,
  strata_col     = "Diameter",
  trial_col      = "Study",
  endpoint_col   = "Endpoint"
)

cat("\n=== Dunnett Test for Eta Heterogeneity ===\n")
cat(sprintf("QC: %s\n", dunnett$qc_message))

if (dunnett$qc_pass) {
  cat(sprintf("Test statistic (max|t|) = %.3f, df = %d, p = %.4f\n",
              dunnett$test_stat, dunnett$df, dunnett$p_value))

  if (dunnett$p_value > 0.05) {
    cat("Conclusion: No significant heterogeneity. Common eta is appropriate.\n")
  } else {
    cat("Conclusion: Significant heterogeneity. Consider stratum-specific eta.\n")
  }

  cat("\nContrasts (eta_d - eta_global):\n")
  print(dunnett$contrasts[, c("stratum", "eta_d", "difference", "t_stat", "p_adjusted")])
}


# =============================================================================
# 8. Stratum-Specific Eta Table
# =============================================================================

cat("\n=== Stratum-Specific Eta ===\n")
print(eta_test$eta_table)


# =============================================================================
# 9. Plot with Stratum-Specific Eta (only if heterogeneity detected)
# =============================================================================

# Only plot stratum-specific eta if Dunnett test rejects H0
if (dunnett$qc_pass && dunnett$p_value < 0.05) {
  cat("\n=== Stratum-specific η plot (heterogeneity detected) ===\n")

  out3 <- slim2EP_plot_by_study(
    result,
    data           = aaa_data,
    endpoints      = c("GrowthRate", "RuptureRisk"),
    title          = "AAA Center Data: Growth vs Rupture (stratum-specific η)",
    facet_by       = "Diameter",
    facet_scales   = "fixed",
    eta_by_stratum = eta_test,
    y_scale        = "log"
  )

  out3$plot  # each panel now has different slope
} else {
  cat("\n=== Stratum-specific η plot skipped (no heterogeneity detected) ===\n")
  cat("Global η is appropriate for all strata.\n")
}


# =============================================================================
# 10. Sample Size Calculation
# =============================================================================

# Example: Power a trial on GrowthRate
# Assume treatment effect of -0.67 mm/year (from HDL drug effect)
# Enrollment weights: 5% at 30mm, 45% at 35mm, 45% at 40mm, 3% at 45mm, 2% at 50mm

ss <- slim_sample_size(
  result       = result,
  endpoint     = "GrowthRate",
  effect_size  = -0.67,              # absolute effect: 0.67 mm/year reduction
  sd           = 3.5,                # from literature
  strata_col   = "Diameter",
  weights      = c(0.05, 0.45, 0.45, 0.03, 0.02),
  adjustments  = list(statin = -0.63),  # statin adjustment
  power        = 0.90,
  alpha        = 0.05,
  n_arms       = 2,                  # 2-arm trial
  arm_props    = c(0.5, 0.5)         # 1:1 allocation
)

cat("\n=== Sample Size Calculation ===\n")
cat(sprintf("Primary endpoint: GrowthRate\n"))
cat(sprintf("Treatment effect: %.3f mm/year\n", ss$effect_size_weighted))
cat(sprintf("Pooled SD: %.2f mm/year\n", ss$sd_pooled))
cat(sprintf("Power: %.0f%%, Alpha: %.2f (two-sided)\n", ss$power * 100, ss$alpha))
cat(sprintf("Number of arms: %d\n", ss$n_arms))
cat(sprintf("Arm proportions: %s\n", paste(ss$arm_props, collapse = ", ")))
cat(sprintf("\nRequired N per arm: %s\n", paste(ss$n_per_arm, collapse = ", ")))
cat(sprintf("Required N total:   %d\n", ss$n_total))

# Stratum details
cat("\nStratum-specific details:\n")
print(ss$strata_table)

# Translation to RuptureRisk
if (!is.null(ss$translation)) {
  cat(sprintf("\n=== Translation via η ===\n"))
  cat(sprintf("η = %.3f\n", ss$translation$eta))
  if (!is.na(ss$translation$pct_change)) {
    cat(sprintf("Implied change in %s: %.1f%% (change in log scale: %.3f)\n",
                ss$translation$endpoint, ss$translation$pct_change,
                ss$translation$implied_effect))
  } else {
    cat(sprintf("Implied effect on %s: %.3f\n",
                ss$translation$endpoint, ss$translation$implied_effect))
  }
}
