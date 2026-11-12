# slimeta 0.1.1

   * Added the method reference to the Description field.
   * Examples now run on the packaged data; slower ones use \donttest.
   * The vignette restores the user's options at the end.

# slimeta 0.1.0

* First release.
* `slim_fit()` fits the Stacked Linear Mixed-Effects (SLIM) model for two
  correlated endpoints with `metafor::rma.mv()`, with an unstructured
  between-study covariance by default (`struct = "UN"`; `"CS"` is also
  available).
* The translation coefficient `eta` is returned with a delta-method confidence
  interval, and pooled means can use CR2 intervals with Satterthwaite degrees
  of freedom (`satterthwaite = TRUE`).
* `slim_eta_by_stratum()` and `slim_eta_dunnett()` test whether `eta` differs
  across strata (Cochran's Q and Dunnett).
* `slim_prior_adjust()` and `slim_prior_plot()` apply and display mild
  regularization of effects close to zero.
* `slim_sample_size()`, forest and scatter plots, and the `aaa` example data.
