# Numerical validation harness — exact beta-Poisson (₁F₁) evaluator

Provenance for the base-R evaluator in `R/` (function `exact_beta_poisson_*`).
These Python scripts were used to design and certify the numerics against a
high-precision oracle **before** the base-R version was written (R was not
available in the prototyping environment). They are kept here as a record; they
are **not** part of the R package build.

## Requirement

```
pip install mpmath        # 50-digit oracle: mpmath.hyp1f1
python3 certify_final.py  # prints the worst-case relative error over the grid
```

## Files (run order / purpose)

| File | Purpose |
|---|---|
| `prototype_1f1.py` | First evaluator drafts vs the mpmath oracle; exposes the naïve-series failure and the legacy-approximation error. |
| `prototype_1f1_v2.py` | Log-space Kummer + asymptotic; adds the handoff/no-gap test. |
| `calibrate_routing.py` | Maps which method is accurate in which (β, dose) regime; calibrates the closed-form approximation error. |
| `verify_r_logic.py` | Mirrors the base-R logic 1:1 to catch transcription bugs; found the fixed-cutoff failure at large β. |
| `verify_routing_v3.py` | Error-driven routing; found the integer-β false-convergence case. |
| `tune_guard.py` | Tunes the asymptotic-acceptance guards (`dose ≥ α+β`, `exp(−dose) < 1e-12`, 2F0 truncation). |
| `certify_final.py` | Final certification: worst-case rel. error 2.26e-9 over 420 points incl. β = 1e5; monotone P(dose) ∈ [0,1]. |
| `PROTOTYPE_VALIDATION_RESULTS.md` | Human-readable summary of the findings. |

## Status

Superseded in the package by an `Rmpfr` skip-if-not-installed test in
`tests/testthat/` that performs the same oracle comparison natively in R. Keep
these for reference / methods-section provenance.
