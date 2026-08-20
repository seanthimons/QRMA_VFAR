# singlehit 0.1.0

- Converted the legacy CAMRA exponential and beta-Poisson workflow into a
  side-effect-free R package API.
- Added grouped-binomial bootstrapping, diagnostics, effective doses, tidy
  summaries, prediction curves, and ggplot2 output.
- Added optional mirai-backed parallel model refits with backend-independent
  reproducible bootstrap samples.
- Reported chi-squared, AIC, and BIC model preferences separately, with the
  legacy chi-squared deviance comparison driving the primary result.
- Added machine-readable fit assessments and human-readable conclusions that
  keep absolute model adequacy separate from relative model preference.
- Added the experimental `consensus_model_decision()` vote parser for the
  chi-squared, AIC, and BIC model selections.
