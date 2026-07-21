# Prototype validation results — exact beta-Poisson (₁F₁) evaluator

Companion to `EXACT_BETA_POISSON_IMPLEMENTATION_PLAN.md`. Summarizes what the
numerical prototyping established. Oracle: `mpmath.hyp1f1` at 50 digits (to be
re-certified against `Rmpfr` inside the package test suite). Quantity measured:
relative error in the **survival** `S = 1 − P = M(α, α+β, −d)`, because the
binomial log-likelihood uses `log(1 − P) = log(S)` for non-responders and `S` is
the hard-to-compute tail at large dose.

## Headline result

The routed evaluator (`exact_beta_poisson.R`) achieves **worst-case relative
error 2.26e-9** over α ∈ [0.05, 2], β ∈ [0.1, 1e5], dose ∈ [1e-3, 1e5]
(420 grid points), and produces a monotone P(dose) ∈ [0, 1]. No new runtime
dependency; pure base R.

## What each candidate method does (measured, not assumed)

| Method | Worst-case rel. err in S | Verdict |
|---|---|---|
| Naïve double-precision series | ~10³¹ (fails at dose ≥ 100) | The GSL/Cephes failure mode, reproduced. |
| Legacy `1 − (1 + d/β)^(−α)` approximation | ~112% at small β | Valid only for large β; **not** an exact-model path. |
| Large-dose asymptotic alone | fails at small dose | Correct only when dose ≫ parameter scale. |
| **Log-space Kummer series** | **2.8e-9 everywhere** | Always correct, never overflows; O(dose) terms (slow only at large dose). |
| **Routed evaluator (final)** | **2.26e-9 everywhere** | Kummer + guarded asymptotic; fast in the realistic regime. |

## Three findings that shaped the design

1. **The legacy approximation is not exact-grade.** `1 − (1 + d/β)^(−α)` only
   reaches 1e-8 accuracy around β ≈ 1e6–1e7; at β = 1e5 it is ~1e-6, and at
   small β its error in the survival tail reaches ~112%. It is excellent as the
   *approximate* beta-Poisson (its intended regime) but must not stand in for the
   exact model. This directly supports Pinch Point 5.

2. **A fixed dose cutoff for the asymptotic is wrong.** Switching to the
   large-dose asymptotic at a fixed dose (e.g. 600) gives ~10⁴ relative error at
   large β, because when β is large a moderate dose is still `dose ≪ β` — not the
   asymptotic regime. The switch must be **regime-aware**. The validated criterion
   accepts the asymptotic only when *all three* hold:
   `dose ≥ α+β` **and** `exp(−dose) < 1e-12` (recessive ~e^(−d) term negligible)
   **and** the 2F0 truncation term `< 1e-12`. Otherwise the log-space Kummer
   series is used.

3. **Integer β creates a false "converged" signal.** For integer β the
   asymptotic's 2F0 series terminates exactly, which can look converged at small
   dose where the asymptotic is invalid (observed error ~10⁶ before the
   `dose ≥ α+β` and `exp(−dose)` guards were added). Both guards are necessary,
   not decorative.

## Performance outlook (relevant to the 10,000-bootstrap target)

The cost driver is the Kummer series, which needs ≈ O(dose) terms. Crucially,
**realistic fits have small/moderate β**, where:
- large gec doses satisfy the asymptotic guards → fast asymptotic (a few terms);
- small doses converge in a short Kummer series.

The slow Kummer path (~10⁵ terms) is reached **only** in the extreme
large-β *and* large-dose corner, which is transient during optimization and rare
in practice. Spot check (noro-like small β, doses 10³–10⁵ gec): every case took
the fast asymptotic path. If benchmarks later show that corner is hit, a
higher-order large-β expansion (corrections in 1/β) would remove it; flagged as a
future optimization.

## Design decisions confirmed by the prototype

- Fit and evaluate natively in **(α, β)** (Pinch Point 1) — matches the
  collaborators' `Probinfbpa(a, b, d)` interface and the literature.
- **Base-R evaluator, `Rmpfr` as a test-only oracle** (Pinch Point 2) — no
  runtime dependency, no compiled code; process-safe for the `mirai` parallel
  backend.
- **Error-driven crossover**, not a magic dose (Pinch Point 5) — and, usefully,
  the realized switch dose can be reported per fit to reproduce a "33k-like"
  number in an auditable way.

## Handoff continuity

Across the asymptotic↔Kummer handoff the survival is continuous to ~1e-6
(relative jump), so the likelihood surface is smooth enough for BFGS. The
fitting-path evaluator (Pinch Point 3) remains the single continuous method to be
finalized once fitting is wired up.

## Artifacts

- `exact_beta_poisson.R` — the validated base-R evaluator + response function.
- `prototype_1f1.py`, `prototype_1f1_v2.py`, `calibrate_routing.py`,
  `verify_routing_v3.py`, `tune_guard.py`, `certify_final.py` — the validation
  harness (oracle comparisons, regime calibration, guard tuning, final
  certification).

## Next steps

1. Restore the working tree locally (git note in the plan, Appendix A) so the
   package builds with `consensus_model_decision` present.
2. Drop `exact_beta_poisson.R` into `R/`, add the `Rmpfr` grid test to
   `tests/testthat/` (skip-if-not-installed), and confirm the R run reproduces
   the ≤ 3e-9 figure natively.
3. Wire the model into `model_probability` / fitting / effective_dose, then
   proceed through the plan's phases.
