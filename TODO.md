# TODO

## Owner-requested dose-response workflow

This checklist records the numbered CAMRA/owner-requested workflow items without
presenting them as completed package documentation. Keep package behavior in
`README.md`; use this file to track remaining alignment work and compatibility
tradeoffs.

1. Do data meet criteria: quantal response, more than one nonzero response, and minimum `N = 3`?
2. Run the Cochrane-Armitage test of trend, evaluated by chi-squared.
3. Fit exponential and typically appropriate beta-Poisson models via MLE, with absolute fit evaluated by chi-squared.
4. Calculate bootstrapped 95th percentile confidence intervals from the MLE model.
5. Compare relative fit: which is better, by difference in deviances compared with the chi-squared distribution?
6. Pool if desired, typically by analysts manually running code for pooled versus unpooled models.
7. Test sufficient fit by difference in deviances compared with the base case.
8. Output statistics: deviances and chi-squared p-values.
9. Output plots: exponential and beta-Poisson.
10. Output the best-fitting model choice.

## Current implementation notes

- `dose_trend_test()` currently preserves the package's existing log-dose,
  one-sided normal/Z trend gate. The owner-requested checklist item names
  chi-squared evaluation; if needed for reporting, that can be added
  compatibility-preservingly as `Z^2` with 1 degree of freedom without changing
  the current `p_value` or pass/fail rule.
- The `comparison` result keeps decision rules separate: `preferred` follows the
  legacy chi-squared deviance test, while `preferred_by_AIC` and
  `preferred_by_BIC` report the corresponding information-criterion choices.
- The combined `assessment` result is intended for programmatic parsing. Its
  `appropriate` flag comes only from absolute goodness-of-fit, while
  `recommendation` is one of `recommended`, `acceptable_alternative`,
  `preferred_but_inadequate`, or `not_recommended`. A model is only
  `recommended` when it both fits the data adequately and is preferred by the
  chi-squared deviance comparison. The `conclusion` column provides the matching
  human-readable statement.
- `consensus_model_decision()` is experimental. It treats the chi-squared
  comparison, AIC, and BIC as three model-selection votes and reports each vote,
  the total per model, and an agreement level of `unanimous`, `majority`, or
  `no_consensus`. This is a relative selection summary only; use
  `analysis$assessment` to determine whether the selected model also fits the
  data adequately.
- Boosterpak remains optional repository-development tooling. It is not a
  package dependency and is not needed to install or use `qrmavfar`.
