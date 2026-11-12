## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, message = FALSE, warning = FALSE, comment = "")
old_opts <- options(digits = 4)

## ----aaa-fit------------------------------------------------------------------
library(slimeta)
data(aaa)

# View structure
head(aaa)

# aaa stores RuptureRisk on the log scale (y, se), but slim_fit() applies the
# log transform itself. Give it the rate (per patient-year) and its SE on the
# original scale (delta method: SE(rate) = rate x SE(log rate)).
aaa_fit <- aaa
rup <- aaa_fit$Endpoint == "RuptureRisk"
aaa_fit$y[rup]  <- aaa_fit$Mean[rup]
aaa_fit$se[rup] <- aaa_fit$Mean[rup] * aaa_fit$se[rup]

# Fit SLIM model
result <- slim_fit(
  data          = aaa_fit,
  formula       = ~ 0 + Endpoint:factor(Diameter),
  trial_col     = "Study",
  endpoint_col  = "Endpoint",
  y_col         = "y",
  se_col        = "se",
  transform     = list(RuptureRisk = "log"),
  display_scale = list(RuptureRisk = 1000)
)

# Key results
result$R2_between  # Between-study R-squared
result$eta         # Translation coefficient
result$ci_eta      # 95% CI for eta
result$pooled      # pooled means by endpoint and diameter, on the display scale

## ----aaa-plot, fig.width=10, fig.height=5-------------------------------------
library(ggplot2)
library(gridExtra)

# Get pooled estimates from the model
res <- result$pooled
res$Endpoint <- ifelse(grepl("GrowthRate", res$term), "GrowthRate", "RuptureRisk")
res$Diameter <- as.numeric(gsub(".*Diameter)", "", res$term))

# Calculate observed means from data
obs_growth <- aggregate(Mean ~ Diameter,
                        data = subset(aaa, Endpoint == "GrowthRate"),
                        FUN = mean)
obs_growth$Endpoint <- "GrowthRate"
names(obs_growth)[2] <- "Observed"

obs_rupture <- aggregate(Mean1000 ~ Diameter,
                         data = subset(aaa, Endpoint == "RuptureRisk"),
                         FUN = mean)
obs_rupture$Endpoint <- "RuptureRisk"
names(obs_rupture)[2] <- "Observed"

obs <- rbind(obs_growth, obs_rupture)
obs$Diameter <- as.numeric(as.character(obs$Diameter))

# Merge model estimates with observed
plotdata <- merge(res, obs, by = c("Diameter", "Endpoint"))
plotdata$Model <- "SLIM"

# Growth Rate plot
p1 <- ggplot(subset(plotdata, Endpoint == "GrowthRate"),
             aes(x = Diameter, y = est)) +
  geom_point(aes(color = Model), size = 3) +
  geom_errorbar(aes(ymin = ci.lb, ymax = ci.ub, color = Model), width = 1) +
  geom_point(aes(x = Diameter, y = Observed), shape = 8, size = 3, color = "black") +
  labs(x = "Diameter (mm)", y = "Estimate / Observed Mean",
       title = "GrowthRate (mm/yr)") +
  scale_color_manual(values = c("SLIM" = "red")) +
  theme_minimal() +
  theme(legend.position = "bottom")

# Rupture Risk plot
p2 <- ggplot(subset(plotdata, Endpoint == "RuptureRisk"),
             aes(x = Diameter, y = est)) +
  geom_point(aes(color = Model), size = 3) +
  geom_errorbar(aes(ymin = ci.lb, ymax = ci.ub, color = Model), width = 1) +
  geom_point(aes(x = Diameter, y = Observed), shape = 8, size = 3, color = "black") +
  labs(x = "Diameter (mm)", y = "Estimate / Observed Mean",
       title = "RuptureRisk (per 1000 p.y.)") +
  scale_color_manual(values = c("SLIM" = "red")) +
  theme_minimal() +
  theme(legend.position = "bottom")

# Combine panels
grid.arrange(p1, p2, ncol = 2,
  top = "SLIM Estimates (red, with 95% CI) vs Observed Means (black star)")

## ----eval=FALSE---------------------------------------------------------------
# # Scatter plot by study with conditional expectation line
# out <- slim2EP_plot_by_study(
#   result,
#   data         = aaa_fit,
#   endpoints    = c("GrowthRate", "RuptureRisk"),
#   title        = "AAA: Growth vs Rupture",
#   facet_by     = "Diameter",
#   facet_scales = "fixed",
#   y_scale      = "log"
# )
# out$plot

## ----aaa-heterogeneity--------------------------------------------------------
# Cochran's Q test for eta heterogeneity
eta_test <- slim_eta_by_stratum(
  data          = aaa_fit,
  strata_col    = "Diameter",
  endpoints     = c("GrowthRate", "RuptureRisk"),
  trial_col     = "Study",
  endpoint_col  = "Endpoint",
  y_col         = "y",
  se_col        = "se",
  transform     = list(RuptureRisk = "log")
)

eta_test$Q_test  # Q statistic and p-value

## ----aaa-samplesize-----------------------------------------------------------
ss <- slim_sample_size(
  result       = result,
  endpoint     = "GrowthRate",
  effect_size  = -0.67,              # 0.67 mm/year reduction
  sd           = 3.5,
  strata_col   = "Diameter",
  weights      = c(0.05, 0.45, 0.45, 0.03, 0.02),
  adjustments  = list(statin = -0.63),
  power        = 0.90,
  alpha        = 0.05
)

ss$n_total        # Total N required
ss$translation    # Implied effect on RuptureRisk via eta

## ----eval=FALSE---------------------------------------------------------------
# # Fit SLIM with drug class as stratum
# result <- slim_fit(
#   data         = obesity_data,
#   formula      = ~ 0 + Endpoint:DrugClass,
#   random       = ~ Endpoint | Study,
#   trial_col    = "Study",
#   endpoint_col = "Endpoint",
#   y_col        = "Effect",
#   se_col       = "SE"
# )
# 
# # Test if eta differs by drug class
# eta_test <- slim_eta_by_stratum(
#   data       = obesity_data,
#   strata_col = "DrugClass",
#   endpoints  = c("WeightLoss", "HbA1c"),
#   ...
# )
# 
# dunnett <- slim_eta_dunnett(
#   result         = result,
#   eta_by_stratum = eta_test,
#   data           = obesity_data,
#   strata_col     = "DrugClass",
#   ...
# )
# 
# # If p < 0.05: GLP-1 and GIP have different translation coefficients
# # This would mean the weight loss -> glycemic control relationship
# # differs between drug mechanisms

## ----eval=FALSE---------------------------------------------------------------
# # Original data has some rate ratios near zero
# data_adj <- slim_prior_adjust(
#   data       = mydata,
#   y_col      = "log_RR",
#   se_col     = "SE_log_RR",
#   prior_type = "normal",
#   prior_mean = 0,         # no effect (RR = 1 on original scale)
#   prior_sd   = 1.0        # weakly informative
# )
# 
# # Visualize the shrinkage
# slim_prior_plot(data_adj, "log_RR", "SE_log_RR")
# 
# # Use adjusted values in SLIM
# result <- slim_fit(
#   data     = data_adj,
#   y_col    = "log_RR_adj",
#   se_col   = "SE_log_RR_adj",
#   ...
# )

## ----reset-options, include=FALSE---------------------------------------------
# Restore the user's options
options(old_opts)

