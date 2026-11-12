# slimeta

Stacked Linear Mixed-Effects Framework for Multivariate Meta-Analysis

## Installation

```r
# Install from GitHub
install.packages("remotes")
remotes::install_github("rabbee001/slimeta")
```

## Quick Start

```r
library(slimeta)

# Fit SLIM model
result <- slim_fit(
  data          = your_data,
  formula       = ~ 0 + Endpoint:factor(Stratum),
  random        = ~ Endpoint | Study,
  trial_col     = "Study",
  endpoint_col  = "Endpoint",
  y_col         = "Mean",
  se_col        = "SE",
  transform     = list(Endpoint2 = "log"),
  display_scale = list(Endpoint2 = 1000)
)

# Key results
result$eta         # translation coefficient
result$R2_between  # between-study R-squared
result$pooled      # pooled estimates

# Visualization
slim2EP_plot_by_study(result, data, endpoints = c("EP1", "EP2"))

# Test heterogeneity
eta_test <- slim_eta_by_stratum(data, strata_col = "Stratum", ...)
dunnett  <- slim_eta_dunnett(result, eta_test, data, ...)

# Sample size
slim_sample_size(result, endpoint = "EP1", effect_size = -0.5, sd = 2.0)
```

## Main Functions

| Function | Purpose |
|----------|---------|
| `slim_fit()` | Fit SLIM model |
| `slim2EP_plot_by_covariate()` | Panel plot by covariate |
| `slim2EP_plot_by_study()` | Scatter plot with E(Y2\|Y1) |
| `slim_eta_by_stratum()` | Cochran's Q test for eta |
| `slim_eta_dunnett()` | Dunnett test for eta |
| `slim_sample_size()` | Sample size calculation |

## References

Rabbee N. (2026). Translating treatment effects between correlated endpoints.
*BMC Medical Research Methodology*, 26, 126. https://doi.org/10.1186/s12874-026-02852-x

## License

MIT
