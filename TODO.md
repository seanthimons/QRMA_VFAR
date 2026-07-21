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

## Verified implementation status (code audit, 2026-07-20)

Each item below was independently verified against the working source in `R/*.R`
and `tests/testthat/*.R` — not against this file's own self-description. Verdicts:
IMPLEMENTED (code does what the item asks), PARTIAL (present but incomplete),
MISMATCH (code does something methodologically different), MISSING (absent).

| # | Verdict | Summary |
|---|---------|---------|
| 1 | IMPLEMENTED | Quantal structure, ≥3 distinct dose groups, and >1 positive response total (reading B) enforced as hard errors with explanatory messages. |
| 2 | IMPLEMENTED | Directional one-sided pre-screen (`Zca` on log dose) that reproduces the Weir gate exactly; see clarification below. |
| 3 | IMPLEMENTED | Binomial MLE (`optim`/BFGS) + chi-squared deviance goodness-of-fit. |
| 4 | IMPLEMENTED* | Percentile bootstrap CIs, 95% default, MLE refit per replicate. Caveats below. |
| 5 | IMPLEMENTED | Δdeviance vs chi-squared; hard-wired to exactly two models. |
| 6 | MISSING (by design) | No pooling; item states this is done manually / out of scope. |
| 7 | IMPLEMENTED | Residual deviance vs saturated base case, `qchisq(df_residual)`. |
| 8 | IMPLEMENTED | Deviances + chi-squared p-values surfaced as tibble columns. |
| 9 | PARTIAL | Exponential + approximate BP fully plotted; exact model's CI bands and bootstrap-parameter plot break. |
| 10 | IMPLEMENTED* | `recommendation` label + chi-squared-driven `preferred`; two-model only. |

### Detail and evidence

1. **Data criteria — PARTIAL.** Quantal/binomial structure is enforced
   (`R/data.R:115-131`: finite, non-negative whole-number counts, `positive +
   negative > 0` per group). Two gaps: (a) the `nrow(data) < 3` guard
   (`R/data.R:119`) runs *before* the group-by-dose aggregation
   (`R/data.R:86-92`), so it counts raw input rows, not distinct dose groups —
   3 rows collapsing to 2 doses passes. (b) There is **no check** that more than
   one dose group has `positive > 0`; a single-responding-dose or all-zero
   dataset validates clean. **→ Resolved:** both gaps fixed in
   `validate_dose_response_groups()` (see the DONE note below); the
   nonzero-response rule uses reading B (>1 positive response in total).

2. **Cochran-Armitage trend — IMPLEMENTED (requirement clarified).** The item's
   wording "evaluated by chi-squared" was a conflation. The stakeholder clarified
   (2026-07-20): item 2 is a **pre-screening test for a monotonic increasing
   trend** — response rising with dose — run before any model is fit, exactly as
   the Weir CAMRA code does it (`if (Zca > 1.644)` in
   `dev/CAMRA_bootstrap_v10_5_kh.R:57-62`). `dose_trend_test()` (`R/fit.R`)
   reproduces that reference **line for line**: log-dose scores, the same `Zca`
   formula, and a one-sided gate (`pnorm(Zca, lower.tail = FALSE) < alpha`, i.e.
   `Zca > 1.645`). The test is directional by design — a decreasing trend must
   not pass — so a chi-squared (two-sided `Z^2`) gate would be *incorrect* here.
   The chi-squared distribution is used elsewhere in the workflow (goodness-of-fit
   and best-model comparison, items 3/5/7/10), which is where the reproduction
   target's `chisq.*` columns come from — there is no trend statistic in that
   target at all. Reporting-only chi-squared columns briefly added to the trend
   test were reverted once the intent was clarified.

3. **Exp + beta-Poisson MLE, absolute fit by chi-squared — IMPLEMENTED.**
   Negative binomial log-likelihood minimized by `stats::optim(method="BFGS")`
   (`R/fit.R:61-73`); absolute fit via residual deviance vs `stats::qchisq` /
   `stats::pchisq` on `df_residual` (`R/fit.R:231-238`). Note "beta-Poisson" in
   the automatic two-model set (`fit_dose_response_models`, `R/fit.R:136-139`)
   means the **approximate** model; the exact beta-Poisson is fittable but
   opt-in only via `fit_dose_response(data, "exact_beta_poisson")`.

4. **Bootstrapped 95% CIs — IMPLEMENTED, with caveats.** Grouped binomial
   bootstrap refitting the MLE each replicate (`R/bootstrap.R:59`,
   `bootstrap_refit` at `R/bootstrap.R:105-109`); percentile CIs with 95% default
   (`bootstrap_confint`, `R/bootstrap.R:187-208`). Works for all three models
   including the exact one (model-agnostic via `model_parameters()` +
   `effective_dose()`). Caveats: (a) default `resample = "observed"` draws from
   observed group proportions, not the fitted-MLE curve — strict parametric
   bootstrap requires `resample = "fitted"`; (b) intervals cover parameters +
   ED10/ED50 but not a predicted-curve band; (c) no dedicated unit test exercises
   `bootstrap_confint` (only `test-bootstrap-plot.R` exists).

5. **Relative fit by Δdeviance vs chi-squared — IMPLEMENTED.**
   `compare_dose_response_models()` (`R/fit.R`) computes Δdeviance =
   simpler − fuller, compared to `qchisq(1-alpha, df = Δparameters)`; `preferred`
   follows this test (AIC/BIC reported separately, not driving `preferred`). It is
   hard-wired to exactly two fits (`if (length(fits) != 2L) stop(...)`), so the
   exact beta-Poisson third model cannot enter — see the three-model decision
   note below.

6. **Pooling — MISSING by design (deferred).** No pooling / multi-dataset
   aggregation / pooled-vs-unpooled comparison anywhere in `R/`. The only
   "combining" is `as_dose_response()` deduplicating rows at the same dose within
   one dataset (`R/data.R:86-92`). The item itself states pooling is done manually
   by analysts, so the absence is consistent, not a defect.

   **Validation observation (2026-07-20, deferred):** the QMRA-wiki Excel target
   contains ~22 pre-pooled experiments (comma-separated ids, e.g. `253, 254`).
   Feeding their rows through `as_dose_response()` *aggregates* repeated doses
   (sums the counts), whereas the reference pooled fits keep each sub-study's row
   as a separate binomial observation. That changes the saturated model, the
   residual dof, and the deviance — e.g. `253, 254` has 7 rows collapsing to 4
   distinct doses, deviance 3.05 vs the target's 8.85. This is the expected
   consequence of not having a pooling feature; reproducing the pooled targets
   would require stacking rather than aggregating same-dose rows. Left for the
   future pooling work, not a fit-correctness bug.

7. **Sufficient fit vs base case — IMPLEMENTED.** The base case is the saturated
   model: `saturated_log_lik` (`R/fit.R:78-83`), `deviance = 2*(saturated_log_lik
   - log_lik)` (`R/fit.R:99`), tested against `qchisq(1-alpha, df_residual)` in
   `goodness_of_fit()` (`R/fit.R:231-238`). Genuinely distinct from item 5
   (which compares the two fitted models to each other on Δparameters df).

8. **Output deviances + chi-squared p-values — IMPLEMENTED.** Surfaced as real
   columns, not internal scalars: `goodness_of_fit()` → `deviance`, `p_value`
   (`R/fit.R:236-237`); `compare_dose_response_models()` → `deviance`,
   `deviance_difference`, `chi_square_p_value`; `glance.qdr_fit()` → `deviance`
   (`R/tidiers.R:66`). `analyze_dose_response()` carries both tibbles on the
   `qdr_analysis` object.

9. **Output plots — PARTIAL.** `plot_dose_response()` / `autoplot.qdr_analysis`
   plot data points, fitted curve, and bootstrap CI ribbons for exponential and
   approximate beta-Poisson. The exact model breaks in two spots that hardcode the
   approximate-BP `n50` algebra: `bootstrap_prediction_matrix()`
   (`R/plot.R:60-70`) and `autoplot.qdr_bootstrap()` (`R/plot.R:152-166`), both
   using `if (model == "exponential") … else {n50 algebra}`. The exact model's
   bare estimate line works (it routes through `model_probability`), but its CI
   bands and bootstrap-parameter plot do not.

10. **Best-fitting model choice — IMPLEMENTED, two-model scope.** Delivered as a
    sorted per-model tibble with a `recommendation` label
    (`recommended` / `acceptable_alternative` / `preferred_but_inadequate` /
    `not_recommended`) in `build_model_assessment()` (`R/assessment.R:27-67`),
    printed by `print.qdr_analysis` (`R/analysis.R:83-87`). `preferred` is driven
    by the chi-squared deviance test; `consensus_model_decision()` adds an
    experimental three-vote (chi-squared/AIC/BIC) summary. Two-model only — the
    exact model is excluded via `fit_dose_response_models` and the `length(fits)
    != 2L` guard.

### Corrections to the planning docs

- `consensus_model_decision()` is **present and exported** (`R/assessment.R:111-182`,
  `NAMESPACE`). The implementation plan and `HANDOFF.md` claim it was "missing from
  the working tree" — that claim is stale; the function is in the working source.
- The exact beta-Poisson model **exists and is wired into fitting, `effective_dose`,
  and bootstrap** (as of the Phase 3 commit), but is opt-in and absent from the
  "Current implementation notes" below.

### Remaining work, grouped

- **Independent of the three-model decision — DONE (2026-07-20):**
  - Item 1: `validate_dose_response_groups()` enforces the fitting criteria on the
    aggregated dose groups (`R/data.R`) as hard errors with explanatory console
    messages: ≥3 distinct dose groups (the minimum leaving one residual df for a
    two-parameter fit) and >1 positive response in total. The nonzero-response
    rule uses **reading B** ("more than one positive response counted in total",
    `sum(positive) > 1`) rather than reading A ("more than one *responding dose
    group*"). Reading B matches the CAMRA reference, which fits experiments with a
    single responding dose but multiple responders (e.g. Enterovirus 62, E. coli
    39); reading A wrongly rejected them. Initial guard commit `5776110`; reading-B
    correction after validating against the QMRA-wiki target.
  - Item 2: confirmed `dose_trend_test()` correctly implements the directional
    monotonic-increasing pre-screen (Weir `Zca > 1.645` gate). Sharpened the docs
    and added a directional test (increasing passes, decreasing does not). The
    chi-squared reporting columns from the earlier misread of the requirement
    (commit `d8e0fb9`) were reverted after the stakeholder clarified the intent.
  - Item 4: dedicated `test-bootstrap-confint.R` covering percentile bounds,
    level widening, converged-only filtering, and input validation. Commit
    `8f60fb0`.
- **Gated on the three-model / Pinch Point 7 owner decision** (do not auto-decide):
  - Items 5, 9, 10: extend comparison / plotting / best-model selection to include
    the exact beta-Poisson as a third, non-nested, 2-parameter model.

### Reference validation and the multi-start fitting pivot (2026-07-20)

The two-model outputs were validated against the QMRA-wiki Excel target
(`data/raw/QMRA wiki data_v3_Camila.xlsx`, ~120 fittable single-experiment
datasets). The fit / goodness-of-fit / comparison arithmetic reproduces the
reference; the divergences traced to optimizer start values, not formulas.

**Design pivot (documented so end users are not surprised):** the single
geometric-mean-dose `default_start()` was replaced by `candidate_starts()` plus a
multi-start loop in `fit_dose_response()`. A single seed cannot avoid every
start-dependent local optimum — an exponential fit to beta-Poisson-shaped data
(Ward) falls into a spurious saturated optimum, while a low-infectivity pathogen
(e.g. Coxiella, true k ~ 5.7e-11) collapses to k -> 0. Fitting from several seeds
(capped ID50, geometric-mean dose, and the CAMRA fixed `exp(-5)`) and keeping the
highest-likelihood result resolves both. Result on the reference set:
preferred-model match 96% -> 97%, ID50 within 3% 92 -> 114 / 119, beta-Poisson
deviance matches 81 -> 100, and all six tiny-rate fit crashes eliminated.
Multi-start runs only on the initial fit (no `start` supplied); bootstrap
warm-starts pass an explicit `start`, so the 1k-10k refits are unaffected. See
`README.md` (Fitting robustness) for the user-facing note.

With the item-1 nonzero-response screen set to reading B, experiments 62 and 39
(single responding dose, multiple responders) now fit, so **every fittable
reference experiment runs with zero errors** and the preferred-model match is
118/121 (98%). Pooled experiments diverge separately (see item 6).

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
