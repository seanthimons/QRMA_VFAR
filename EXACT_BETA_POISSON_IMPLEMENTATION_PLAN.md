# Implementation Plan: Exact Beta-Poisson (Confluent Hypergeometric) Dose-Response Model in `qrmavfar`

| | |
|---|---|
| **Status** | Draft for review |
| **Audience** | QMRA modeling collaborators + `qrmavfar` maintainers |
| **Target branch** | `feature/package-dose-response-models` (a `dev-hypergeo` topic branch is recommended for the work) |
| **Prepared for** | Review and comment before any code is written |
| **Scope** | Add the exact beta-Poisson dose-response model, evaluated via the confluent hypergeometric function ₁F₁, as a first-class model in the package, on par with the existing `exponential` and (approximate) `beta_poisson` models |

> **How to read this document.** Sections 1–5 establish shared context (the model, the current package, the numerical problem, and the legacy code we are inheriting). Section 6 is the substance for reviewers: each **pinch point** is presented as a decision with explicit options, advantages, disadvantages, and considerations, followed by a provisional recommendation. Sections 7–11 cover validation, testing, phasing, and risks. Nothing here is committed until the pinch-point decisions are signed off.

---

## 1. Purpose and scope

Collaborators have specifically requested the **exact beta-Poisson** dose-response model — the form derived without the large-β approximation, expressed through the confluent hypergeometric function (Kummer's function) ₁F₁. The package today fits only the **exponential** and the **approximate beta-Poisson** models. The known obstacle is numerical: naïve evaluation of ₁F₁ at large dose (or large parameters) fails, which is the origin of the "GSL fails to converge at very large numbers" experience.

We already possess a working, if unpolished, implementation of the exact model in C (a Cephes-based ₁F₁ wrapped for R) together with a hand-tuned guard cascade that sidesteps the numerical failure modes. The goal of this plan is to fold that hard-won numerical knowledge into `qrmavfar` cleanly, as a fully integrated model that participates in fitting, effective-dose calculation, bootstrapping (targeting **10,000 replicates**), tidying, plotting, model comparison, and the consensus workflow — while deciding deliberately about each of the numerical and statistical trade-offs involved.

**In scope:** the exact beta-Poisson response function; its numerically robust evaluation; MLE fitting; `effective_dose`/N50 support; bootstrap support at scale; tidy/glance/augment/predict/plot support; extension of the model-comparison, assessment, and consensus machinery from two models to three; documentation; and tests including a high-precision numerical oracle.

**Out of scope (for now):** other dose-response families (log-probit, log-logistic, Weibull, Gompertz), Bayesian fitting, and correlated-parameter Monte-Carlo uncertainty propagation of the kind seen in the legacy noro scripts (that is a downstream risk-assessment concern, not a model-fitting concern).

---

## 2. Background: the exact beta-Poisson model

### 2.1 Where it sits in the single-hit family

All the mechanistic models in this package descend from **single-hit theory**: each ingested organism independently has some probability of surviving host defenses and initiating infection. If every organism shares a single survival probability `r`, integrating over a Poisson-distributed dose yields the **exponential** model, `P(d) = 1 − exp(−r·d)`. If instead `r` varies from host to host (or organism to organism) according to a **Beta(α, β)** distribution — capturing heterogeneity in susceptibility/infectivity — integrating the single-hit survival over that Beta distribution yields the **beta-Poisson** family.

### 2.2 Exact form

Carrying out that integration exactly gives

```
P(d) = 1 − ₁F₁(α, α + β, −d)
```

where `₁F₁(a, b, z) = M(a, b, z)` is Kummer's confluent hypergeometric function. Here `α` and `β` are precisely the two shape parameters of the Beta(α, β) infectivity distribution. This is the form the collaborators want and the form the legacy C code computes (`1 − hyperg(a, a+b, −dose)`).

### 2.3 Relationship to the approximate model already in the package

The package's current `beta_poisson` model is the **approximate** beta-Poisson,

```
P(d) = 1 − (1 + d/β)^(−α)          [equivalently, in the package's N50 parameterization]
P(d) = 1 − [1 + (d / N50)(2^(1/α) − 1)]^(−α)
```

The approximation is valid when `β ≫ 1` and `β ≫ α`. In that regime the exact and approximate forms coincide to high accuracy. **The exact model earns its keep precisely where the approximation breaks down:** small β (strong heterogeneity), where the ₁F₁ form and the power-law approximation diverge materially. A useful framing for reviewers: keeping both models in the package is not redundant — their agreement (or disagreement) is itself diagnostic of whether the large-β approximation is safe for a given pathogen/dataset.

### 2.4 Note on the parameterization the collaborators use

The legacy interface `Probinfbpa(a, b, d) = .C("dr1f1", a, b, d)` takes `(a, b) = (α, β)` directly and returns `1 − ₁F₁(α, α+β, −d)`. In the noro scripts, `α` and `β` are themselves reconstructed at each Monte-Carlo draw from correlated logit-/log-normal "transformed parameters." The operational consequence for us: **the collaborators' code speaks native (α, β)**, not N50. This weighs heavily on Pinch Point 1.

---

## 3. Current package state (integration surface)

The following is the inventory of everything the new model must touch. Reviewers familiar with the package can skim; it is included so the blast radius is explicit.

| File | Function(s) | What must change |
|---|---|---|
| `R/models.R` | `exponential_response`, `beta_poisson_response`, `model_probability` | Add `exact_beta_poisson_response(dose, alpha, beta)` and a stable ₁F₁ evaluator; add a `switch` branch in `model_probability`. |
| `R/fit.R` | `fit_dose_response` (`match.arg`), `fit_dose_response_models`, `model_parameters`, `default_start`, `model_label`, `effective_dose` | Register the model, its parameter names `(alpha, beta)`, sensible start values, a label, and a root-solving branch in `effective_dose`. |
| `R/fit.R` | `compare_dose_response_models` | Currently **hard-coded to exactly two models** with a single nested chi-square. Must generalize to three. See Pinch Point 7. |
| `R/assessment.R` | `assess_dose_response_models`, `build_model_assessment`, `consensus_model_decision` | Generalize assessment to N models. **`consensus_model_decision` is currently missing from the working tree** (present only in `qrmavfar.Rcheck/00_pkg_src/…`) and must be restored, then broadened. See Pinch Point 10. |
| `R/bootstrap.R` | `bootstrap_dose_response`, `bootstrap_refit`, `bootstrap_dose_response_models` | Parameter-name handling is already generic via `model_parameters()`; verify `ed10`/`ed50` work through the exact `effective_dose`. Performance at 10k replicates is the concern (Pinch Point 8). |
| `R/tidiers.R` | `tidy`, `glance`, `augment`, `predict` | Mostly generic; confirm Wald SEs and `predict` route through the new response function. |
| `R/plot.R` | `prediction_curve`, `bootstrap_prediction_matrix`, `autoplot.qdr_bootstrap` | `bootstrap_prediction_matrix` has a model-specific `if (model == "exponential") … else …` block assuming the approximate BP algebra; needs an exact branch. `autoplot.qdr_bootstrap` similarly special-cases parameters. |
| `NAMESPACE`, `man/` | exports, roxygen | New exports/docs; regenerate. |
| `DESCRIPTION` | `Suggests` | Add `Rmpfr` (test-only oracle) and possibly `gsl`/`hypergeo` depending on Pinch Point 2. |
| `tests/testthat/` | new + existing | New numerical-accuracy and fitting tests; the existing two-model comparison/consensus tests must be preserved or consciously updated. |

---

## 4. The numerical problem

### 4.1 Why direct ₁F₁ evaluation fails

Kummer's function has the series

```
₁F₁(a, b, z) = Σ_{n≥0} (a)_n / (b)_n · z^n / n!
```

For our application `z = −d ≤ 0`. When `d` is large, the terms alternate in sign and grow enormous before they decay, so the partial sums involve **catastrophic cancellation**: the true result is a small number obtained by subtracting huge nearly-equal quantities, and double precision is annihilated. The classic remedy is the **asymptotic expansion** valid for large `|z|`, but it is itself only asymptotic (non-convergent) and accurate in a different region. Between "series still accurate" and "asymptotic already accurate" there is a **gap** where neither representation delivers full precision.

### 4.2 Why GSL specifically fails — it is not a GSL bug

The `gsl_sf_hyperg_1F1` routine implements exactly this two-pronged strategy (power series + asymptotic, pick the better). The Cephes routine in our legacy code (`hyperg.c`, Moshier) does the same: it sums the series (`hy1f1p`), also tries the asymptotic (`hy1f1a`), estimates the cancellation error of each, and returns whichever is better — printing an error when the self-estimated relative error exceeds `1e-12`. Because GSL and Cephes are the **same algorithm family**, they fail in the **same region**. Swapping one library for the other does not solve the problem. The solution is to *avoid calling ₁F₁ in the failure region at all*, which is exactly what the legacy `drfunc` cascade does (Section 5).

### 4.3 The crossover direction (and a discrepancy worth resolving)

A methods note for the collaborators: the numerically correct crossover uses the **exact ₁F₁ at small dose** and a **closed-form asymptotic/approximation at large dose** — because ₁F₁ fails at *large* argument, not small. This is the direction the legacy code implements. If a cited paper describes "Pfaff transformation below 33k gec and ₁F₁ above 33k gec," that is either (a) a different parameterization, (b) describing a ₂F₁ (Gauss) computation rather than ₁F₁ (Pfaff is properly a ₂F₁ identity; the confluent analogue is Kummer's transformation), or (c) an inverted description. **We should reconcile the paper's wording against our code before publication** so the methods section is internally consistent. This is tracked as an open question in Section 10.

---

## 5. Legacy code analysis (`dev/hypergeo`, `dev/hgdr.new`)

### 5.1 The Cephes core (`dev/hgdr.new/`)

`hyperg.c` is Stephen Moshier's Cephes confluent hypergeometric routine (Release 2.1, 1988), supported by the usual Cephes machinery (`gamma.c`, `pow.c`, `powi.c`, `polevl.c`, `floor.c`, `const.c`, `mtherr.c`). It is compiled with the other translation units into `hypg.so`/`hypg.dll`. Its own documentation states the accuracy is validated for `a, b, x ∈ [0, 30]` and warns of larger errors "when b is near a negative integer or zero" and for "combinations of arguments [that] yield serious cancellation." This is the numerical fine print behind our failure modes.

### 5.2 The `drfunc` guard cascade (`dev/hgdr.new/hypg.c`) — the real prize

This function is the accumulated numerical wisdom. It is, in effect, the piecewise evaluation strategy we would otherwise have to derive from scratch, already tuned against real pathogen parameter ranges:

```c
double drfunc(double a, double b, double dose)
{
  if (dose < 1e-4)                             return dose*a/(a+b);        // (0) linear small-dose limit
  if (a>1e3 && b<a/100)                        return 1-exp(-dose);        // (1) degenerates to exponential
  if (a>1e2 && b>1e5)                          return 1-exp(-dose*a/b);    // (2) exponential w/ rate a/b
  if (a>1e1 && b>1e5 && dose*a/b>10.0)         return 1-pow(1+dose/b,-a);  // (3) approx beta-Poisson
  if (a>1.0 && b>20*a && dose>10.0)            return 1-pow(1+dose/b,-a);  // (4)
  if (a>1.0 && b>a     && dose>50.0)           return 1-pow(1+dose/b,-a);  // (5)
  if (a>1.0 && b<a     && dose>20.0)           return 1-pow(1+dose/b,-a);  // (6)
  if (a<1.0 && b>50*a)                         return 1-pow(1+dose/b,-a);  // (7)
  if (a<0.1 && b>20*a)                         return 1-pow(1+dose/b,-a);  // (8)
  if (round(a)-a < 1e-4) a = 1.0001*a;                                     // (9) nudge near-integer a
  if (round(b)-b < 1e-4) b = 1.0001*b;                                     // (10) nudge near-integer b
  return 1 - hyperg(a, a+b, -dose);                                        // (11) exact ₁F₁ fallback
}
```

Reading the cascade:

- **(0)** For essentially zero dose, `P → d·α/(α+β)`, the exact linear response (mean infectivity `α/(α+β)` times dose). Avoids evaluating ₁F₁ near the origin.
- **(1)–(2)** Degenerate regimes where the beta-Poisson collapses to an exponential; use it directly.
- **(3)–(8)** The large-dose / large-β regions where the **approximate beta-Poisson** `1 − (1 + d/β)^(−α)` is accurate. These are exactly the regions where ₁F₁'s series cancels catastrophically. Note the thresholds are *regime-dependent* (dose > 10, > 20, or > 50 depending on the relationship between `a` and `b`), **not a single universal cutoff.**
- **(9)–(10)** A deliberate fudge: when `a` or `b` is within `1e-4` of an integer, nudge it by 0.01% to dodge Cephes' precision collapse near integer `b` (the documented "b near a negative integer or zero" hazard).
- **(11)** Only the well-conditioned residual — small-to-moderate dose, non-integer parameters — reaches the exact ₁F₁ call.

The strategic insight: **by the time control reaches line (11), all the hard cases have been diverted.** The exact evaluator therefore only ever runs in a region where even a modest implementation converges quickly. That is what makes a dependency-free re-implementation feasible (Pinch Point 2).

### 5.3 The R interface (`dev/hypergeo/f1.R`)

```r
dyn.load('hypg.dll')
Probinfbpa <- function(a, b, d)
  .C("dr1f1", as.double(a), as.double(b), as.double(d), pinf = double(1))$pinf
```

A thin `.C` shim over `dr1f1`, which simply calls `drfunc`. Confirms the `(α, β, dose)` signature and the `1 − ₁F₁` convention.

### 5.4 Performance evidence (`dev/hypergeo/predications.R`)

The collaborators already profiled this workload. A naïve loop extrapolated to roughly **16 hours for 10,000 iterations**; after restructuring/parallelization it fell to about **19 minutes** (their note: "~98% decrease"). This is directly relevant to the 10,000-bootstrap target and sets a realistic performance bar (Pinch Point 8).

---

## 6. Pinch points (decisions for reviewers)

Each subsection is a decision. The **recommendation** is provisional and is what we will implement absent reviewer objection.

### Pinch Point 1 — Parameterization: `(α, β)` vs `(α, N50)`

**Options.** (A) Fit natively in `(α, β)`, report N50/ID50 as a derived quantity. (B) Fit in `(α, N50)` to mirror the existing approximate `beta_poisson` API.

**Advantages of (A) `(α, β)`:**
- Matches the collaborators' existing code and the QMRA literature for the *exact* model; parameters are directly comparable to published values.
- `α` and `β` are the mechanistic Beta-distribution shape parameters — interpretable as heterogeneity/mean-infectivity.
- No inner root-solve during optimization (see disadvantage of B).

**Disadvantages of (A):**
- API differs from the approximate model (`n50` there, `beta` here), so downstream code that assumes a shared column name must handle both. (The package already keys off `model_parameters()`, so this is manageable but must be verified in tidy/bootstrap/plot.)

**Advantages of (B) `(α, N50)`:**
- Uniform parameter vocabulary with the approximate model; `effective_dose(…, 0.5)` trivially returns a fitted parameter.

**Disadvantages of (B):**
- There is **no closed form** for N50 under the exact model. Every likelihood evaluation would need to root-solve `₁F₁(α, α+β, −N50) = 0.5` to recover `β`, nesting a solver inside the optimizer inside each of 10,000 bootstraps. Slow and fragile.
- Diverges from the collaborators' parameter conventions, complicating cross-checks.

**Considerations.** The package's machinery is already parameter-name-generic; the main cost of (A) is a handful of `if model == …` branches in `plot.R` and the `autoplot` methods. N50 can still be surfaced everywhere via `effective_dose`.

**Recommendation: (A).** Fit in `(α, β)`; expose N50/ID50 through `effective_dose`. This aligns with collaborators, avoids a nested solver in the hot path, and keeps interpretation literature-standard.

---

### Pinch Point 2 — The ₁F₁ evaluation engine

This is the central numerical decision. **Options:** (A) Port `drfunc` (guards + a Kummer-stabilized series for the residual) to **vectorized base R**, with `Rmpfr` used only in tests as an accuracy oracle. (B) Evaluate ₁F₁ with **`Rmpfr` arbitrary precision** at runtime. (C) **Retain the compiled C** (`hgdr.new`) as `src/` in the package via `.C`/`.Call`. (D) Depend on an existing R package (`gsl`, `hypergeo`).

**Option A — base-R port of `drfunc` (+ Rmpfr in tests).**
- *Advantages:* No runtime dependency, no compiled code, no system libraries; matches the package's existing careful base-R numeric style (`log1pexp`, `log_expm1`, etc.); fully vectorizable over the fixed dose grid; portable across platforms and CRAN-friendly; reproducible under `renv`. Because the guards divert every hard case, the residual exact calls are well-conditioned and cheap. Rmpfr as a *Suggests*-only oracle gives a correctness guarantee in CI without shipping it.
- *Disadvantages:* We re-implement and must re-validate the cascade; risk of transcription error (mitigated by validating against the legacy `.dll` and Rmpfr); the residual base-R ₁F₁ series must be written carefully (Kummer transformation to keep terms positive, capped iterations, cancellation guard).

**Option B — Rmpfr at runtime.**
- *Advantages:* Conceptually simplest correctness story; no crossover discontinuities; robust across the whole domain.
- *Disadvantages:* Rmpfr is ~10²–10³× slower than double precision. At 10,000 bootstraps × ~8 dose groups × dozens of optimizer likelihood evaluations × many pathogens, this is the difference between minutes and many hours (consistent with the collaborators' own 16-hour extrapolation). Adds a GMP/MPFR **system-library** dependency, complicating installation and CI. Overkill given the guards already tame the hard cases.

**Option C — ship the compiled C.**
- *Advantages:* Reuses the exact, already-validated numerics; fastest per-call.
- *Disadvantages:* Introduces a compiled-code toolchain requirement and platform-specific build/CI burden (Windows/macOS/Linux); Cephes licensing/attribution must be handled; the code is C89 with global error handling (`mtherr`) that is awkward to make thread-safe for the `mirai` parallel backend; larger maintenance surface for a package that is currently pure R.

**Option D — `gsl` or `hypergeo` package.**
- *Advantages:* Off-the-shelf; less code to own.
- *Disadvantages:* `gsl` reintroduces exactly the convergence failure we are trying to escape (same algorithm family, Section 4.2), and adds a system dependency. `hypergeo` targets the *Gauss* ₂F₁ and complex arguments; using it for ₁F₁ is indirect and not obviously more robust at large real argument.

**Considerations.** The guards are the key enabler: they make Option A viable by ensuring ₁F₁ is only ever called where it is easy. Whichever engine we choose for the *residual* exact call, the **cascade itself should be ported**, because it encodes which regions to avoid.

**Recommendation: (A),** with the structure: port the `drfunc` guard cascade to vectorized base R; implement the residual exact ₁F₁ via a Kummer-transformed positive-term series with a hard iteration cap and error estimate; keep `Rmpfr` as a `Suggests`-only oracle used by the test suite to certify accuracy (target ≤ 1e-8 relative error) across the parameter/dose grid; and additionally cross-check against the legacy `hypg.dll`/`.so` where a build is available. This keeps the shipped package pure R and fast while giving us a high-precision correctness net.

---

### Pinch Point 3 — Continuity vs. speed: branch-switch kinks under MLE

**The issue.** The `drfunc` branch conditions depend on `α` and `β` (e.g. `b > 20*a`, `b < a`, `a > 1`). During maximum-likelihood optimization the solver moves `(α, β)` around; as it crosses a boundary, the functional form switches (exact ₁F₁ ↔ approximation). The response is *continuous by design* (the branches are constructed to agree at the seams), but its **derivative can have small discontinuities** — kinks in the likelihood surface. Gradient-based optimizers (the package uses `optim(method = "BFGS")` with a numerical Hessian) can stall, misreport the Hessian, or converge slowly on such surfaces.

**Options.** (A) Use one **continuous evaluator** for *fitting* (the residual exact-series path, or Rmpfr), and reserve the fast guarded cascade for *prediction* and *Monte-Carlo*, where doses are fixed data and kinks are harmless. (B) Use the guarded cascade everywhere but switch the optimizer to a **derivative-free** method (Nelder–Mead) for this model. (C) Use the cascade everywhere and **smooth the seams** (blend the two forms over a transition band so first derivatives match).

**Advantages / disadvantages.**
- *(A)* Cleanest likelihood surface and most trustworthy Wald standard errors (which come from the Hessian); costs some fitting speed, but fitting is a tiny fraction of total work compared to 10k bootstrap refits — and each bootstrap refit can *also* use the continuous path since accuracy there matters more than raw speed. Slight architectural complexity (two evaluation modes). **Recommended.**
- *(B)* Trivial to implement; Nelder–Mead tolerates kinks. But it does not produce a Hessian, so we lose the Wald SEs that `tidy()` reports, and it converges more slowly — multiplied over 10k bootstraps, possibly a net loss.
- *(C)* Mathematically nicest (truly smooth), but tuning transition bands for every seam is fiddly and easy to get subtly wrong; hard to validate.

**Considerations.** In practice the kinks are small because the approximation is accurate wherever it is selected. The safe default is (A): correctness and honest uncertainty first, with the fast path retained where it is free of consequences. We should empirically check convergence stability on the real species set either way.

**Recommendation: (A),** with a fitting/prediction split of the evaluator. Keep `optim`'s BFGS for the Hessian-based SEs; feed it the continuous evaluator.

---

### Pinch Point 4 — Near-integer parameters and integer-`β` poles

**The issue.** Cephes ₁F₁ loses precision when `b = α + β` is at or near a non-positive integer; `drfunc` handles this by nudging `a`/`b` by 0.01% (lines 9–10). This is a **deliberate bias** trade — a tiny parameter perturbation exchanged for numerical stability.

**Options.** (A) Preserve the nudge on the fast path; on the fitting path, handle integer/near-integer `b` properly via a limiting form or the Kummer transformation `M(a, b, z) = e^z · M(b−a, b, −z)` (which reindexes the problematic parameter). (B) Preserve the nudge everywhere (simplest, matches legacy exactly). (C) Reparameterize to avoid ever landing on integer `b`.

**Advantages / disadvantages.**
- *(A)* Removes the bias where it matters (parameter estimates and SEs) while keeping the battle-tested behavior for prediction. Slightly more code. **Recommended.**
- *(B)* Bit-for-bit agreement with the collaborators' current outputs — valuable for cross-validation — but bakes a 0.01% bias into fitted parameters, which could matter for tight confidence intervals across 10k bootstraps.
- *(C)* Elegant in principle but there is no natural reparameterization that guarantees non-integer `α + β` without distorting interpretation.

**Considerations.** The magnitude of the nudge bias is almost certainly negligible for point estimates, but bootstrap CIs aggregate many fits and could reveal artifacts. Worth a small sensitivity check.

**Recommendation: (A).** Keep the nudge on the fast/prediction path for continuity with legacy results; use the Kummer-reindexed evaluation (no nudge) on the fitting/bootstrap path. Document both.

---

### Pinch Point 5 — Crossover thresholds: fixed vs adaptive, and matching the literature

**The issue.** The legacy thresholds (dose > 10 / 20 / 50, `b > 20*a`, `b > 1e5`, etc.) are hand-tuned constants. They are dataset-agnostic but were presumably tuned for particular pathogen ranges. A cited paper uses a single "33k gec" cutoff, which does not appear in our code.

**Options.** (A) Port the legacy thresholds verbatim (proven, but opaque). (B) Replace them with an **error-driven switch**: evaluate the residual series, and fall back to the approximation only when the series' self-estimated cancellation error exceeds a tolerance — a principled, dataset-adaptive criterion. (C) Expose the crossover as a documented, overridable argument with the legacy values as defaults.

**Advantages / disadvantages.**
- *(A)* Reproduces collaborators' numbers exactly; zero design risk; but the constants are unexplained and may not generalize to new species in the QMRA-wiki set.
- *(B)* Most defensible scientifically — the switch happens exactly where precision demands it, and it self-documents. Requires implementing a reliable per-term cancellation estimate (Cephes already does this internally; we would surface it). **Recommended as the primary, with (A) as a validated fallback.**
- *(C)* Pragmatic compromise; good for transparency and for the methods section, and lets us reconcile the "33k gec" figure by making the cutoff explicit and reportable.

**Considerations.** Approach (B) also resolves the paper-discrepancy concern: if the switch is error-driven, the methods section can state the tolerance rather than a magic dose, which is both more honest and more general. We should still document, for any given fitted `(α, β)`, the dose at which the switch actually occurs — that reproduces a "33k-like" number in an auditable way.

**Recommendation: (B) primary + (A) as a cross-check.** Implement an error-driven crossover; verify it selects the same regions as the legacy thresholds on representative parameters; report the realized switch dose per fit.

---

### Pinch Point 6 — `effective_dose` / N50 for the exact model

**The issue.** `effective_dose` currently has closed forms for the exponential and approximate BP models. The exact model has no closed-form inverse; `ed10`/`ed50` are needed by the bootstrap (`bootstrap_refit`) and by N50 reporting.

**Options.** (A) Root-solve `P(d) = p` with `uniroot` on `log(d)` (monotonic, so robust). (B) Use the approximate-BP inverse as a fast approximation to N50.

**Advantages / disadvantages.**
- *(A)* Exact and consistent with the fitted model; `uniroot` on a monotone function is cheap and reliable. Adds one solve per effective-dose call — negligible next to a full refit, but multiplied across 10k bootstraps × 2 effective doses it is worth vectorizing/caching. **Recommended.**
- *(B)* Fast but inconsistent: reporting an approximate N50 for an *exact* fit invites confusion and small contradictions with the fitted curve.

**Considerations.** Provide a good bracketing interval from the data's dose range to keep `uniroot` fast and guaranteed to converge. Guard against `p` beyond the achievable range.

**Recommendation: (A).** Monotone `uniroot` in log-dose with data-informed brackets.

---

### Pinch Point 7 — The three-model comparison problem

**The issue.** `compare_dose_response_models()` currently **requires exactly two fits** and performs a single nested chi-square (the simpler model's extra-parameter test against the fuller). With three models — exponential (1 parameter), approximate BP (2), exact BP (2) — this breaks: the two 2-parameter models are **not nested in each other**, and the function's `length(fits) != 2` guard and single simpler/fuller split no longer apply.

**Options.**
- **(A) Exponential-as-baseline nested tests + information criteria across all.** Keep the chi-square as a nested test of exponential (the common sub-model) against each 2-parameter model; rank all three by AIC/BIC. The two 2-parameter models are compared *only* by AIC/BIC (legitimate, since they are non-nested).
- **(B) Keep chi-square strictly pairwise; report the exact model on fit + AIC/BIC only.** Least disruption to existing tests; the exact model is scored absolutely (goodness-of-fit) and by IC, but excluded from the nested comparison.
- **(C) Full generalization to N models** with a documented policy for nested vs non-nested pairs.

**Advantages / disadvantages.**
- *(A)* Statistically coherent, keeps a nested test where one is valid, and gives a clean multi-model ranking. Requires rewriting the two-model logic and **updating the existing tests** (which assert two-row outputs and specific `chi_square_df`, `preferred`, etc.). **Recommended.**
- *(B)* Minimal risk to the current test suite and workflow; but treating the exact model as a second-class citizen is unsatisfying given it is the collaborators' headline request.
- *(C)* Most future-proof (anticipates log-probit, Weibull, …) but the largest change and the easiest to over-engineer now.

**Considerations.** The existing tests encode important behavior (e.g. "chi-squared and IC preferences remain distinct," the Ward dataset preferring beta-Poisson). Any generalization must keep those semantics for the two-model case as a special case, so legacy analyses reproduce. This is the single biggest *statistical* design decision in the plan and deserves explicit sign-off.

**Recommendation: (A),** implemented so that the two-model call path reproduces today's outputs exactly (regression-locked), with the three-model path adding exponential-baseline nested tests plus AIC/BIC ranking across all models.

---

### Pinch Point 8 — Performance and parallelism at 10,000 bootstraps

**The issue.** The target is 10,000 bootstrap replicates per model, per species, across the QMRA-wiki set. Each replicate is a full refit (an `optim` run with many likelihood evaluations, each evaluating the response over all dose groups). The collaborators' own profiling shows this is the dominant cost.

**Options / levers.**
- **Vectorize the response over the dose grid** (evaluate all dose groups at once per likelihood call) — essential regardless.
- **Use the existing `mirai` backend** (`backend = "mirai"`) for parallel refits; the package already supports it. Distributing 10k refits across cores is the primary win and mirrors the collaborators' "~19 minutes" result.
- **Warm-start each bootstrap** from the point estimate (already done: `bootstrap_refit` passes `start = original_fit$coefficients`), which minimizes optimizer iterations.
- **Keep the residual ₁F₁ cheap** by leaning on the guards (Pinch Point 2) so the expensive path is rarely taken.
- **Cache/vectorize `effective_dose`** solves (Pinch Point 6).

**Advantages / disadvantages.** Parallelism gives near-linear speedup but requires the evaluator to be **thread-/process-safe** — a strong additional argument *against* shipping the global-state Cephes C (Option C in PP2) and *for* the pure-R evaluator, which is trivially safe under `mirai` process daemons. Vectorization complicates the code slightly but is the highest-leverage single change.

**Considerations.** We should benchmark a single species end-to-end at 1k and 10k replicates, sequential and mirai, before committing to the full run, and record timings (as the collaborators did) so the methods section can report them.

**Recommendation:** Pure-R vectorized evaluator + `mirai` parallel backend + warm starts. Benchmark early; publish timings.

---

### Pinch Point 9 — Dependency footprint and reproducibility

**The issue.** The package today has a clean, pure-R `Imports` list and uses `renv`. Each engine option in PP2 has a different dependency cost.

**Considerations / recommendation.** Prefer **zero new runtime dependencies** (consistent with PP2-A). Add `Rmpfr` to **`Suggests`** only (tests skip gracefully if absent, via `testthat::skip_if_not_installed`). Avoid compiled code and system libraries (GMP/MPFR, GSL) in the shipped package. Record any test-time additions in `renv.lock`. This keeps installation frictionless for collaborators and CI green across platforms.

---

### Pinch Point 10 — Restore and broaden `consensus_model_decision`

**The issue.** `consensus_model_decision()` is **exported (`NAMESPACE`), documented (`man/consensus_model_decision.Rd`), and tested**, but the function body is **absent from the working `R/assessment.R`** — it survives only in `qrmavfar.Rcheck/00_pkg_src/qrmavfar/R/assessment.R`. The working tree is mid-edit and currently would fail `R CMD check` and its own tests. It treats the chi-square, AIC, and BIC selections as three votes and reports unanimity/majority/no-consensus.

**Options.** (A) Restore it verbatim from the `.Rcheck` source, then extend the vote logic to three models. (B) Rewrite from scratch for N models.

**Advantages / disadvantages.** (A) recovers known-good, test-covered behavior with minimal risk, then layers the three-model extension on top; the voting logic (`criterion_choice`, tie handling, `criteria_available`) already anticipates ≥2 models and generalizes naturally. (B) is unnecessary churn.

**Recommendation: (A).** Restore first (so the package builds and the suite passes), commit that as a discrete fix, *then* broaden voting to include the exact model — keeping the existing two/three-model consensus tests as regression anchors.

---

## 7. Validation strategy

Correctness is verified against **three independent oracles**, in increasing order of authority:

1. **Legacy `hypg.dll`/`.so`.** Where a build is available, compare the new base-R `drfunc` port against `Probinfbpa` on a dense `(α, β, dose)` grid spanning the failure regions. This confirms we faithfully reproduced the collaborators' numbers, guard-for-guard. `driver.c` and `output.txt` provide additional fixed reference points.
2. **`Rmpfr` arbitrary precision.** Compute `1 − ₁F₁(α, α+β, −d)` at high precision and require the shipped evaluator to match to a set relative tolerance (target ≤ 1e-8) everywhere, including deep into the large-dose region. This is the ground-truth check that the guards + residual series are not just self-consistent but actually correct.
3. **Analytic limits.** Small-dose linear limit `P → d·α/(α+β)`; large-β agreement with the approximate beta-Poisson; exponential degeneracy in the appropriate corner. These catch structural errors the numeric oracles might mask.

Statistical validation: refit the **Ward rotavirus** fixture and confirm the exact and approximate models agree closely (Ward's β is large), and confirm existing exponential/approximate estimates are unchanged (regression lock).

---

## 8. Testing plan

- **Numerical accuracy** (`test-hyperg.R`, new): grid comparison vs Rmpfr (skipped if Rmpfr absent), vs analytic limits, and vs legacy reference points. Explicitly include the large-dose and small-β cases that break naïve ₁F₁.
- **Response function** (`test-models.R`, extend): vectorization, monotonicity in dose, `P(0)=0`, `P→1` as `d→∞`, parameter validation errors.
- **Fitting** (`test-models.R`/`test-fit.R`): recover known `(α, β)` from simulated data; convergence flags; Wald SEs finite; warm-started refit stability.
- **Effective dose**: `effective_dose(fit, 0.5)` matches a direct root solve; round-trips with the fitted curve.
- **Three-model comparison/consensus** (extend existing): two-model path reproduces current outputs exactly (regression), three-model path behaves per Pinch Point 7/10.
- **Bootstrap/plot**: `bootstrap_prediction_matrix` exact branch; `autoplot.qdr_bootstrap` parameter plot for `(α, β)`; a small (e.g. 200-replicate) end-to-end bootstrap smoke test.
- **Build health**: restore `consensus_model_decision` so `R CMD check` is clean before feature work lands.

---

## 9. Proposed implementation phases

1. **Repository hygiene.** Restore `consensus_model_decision` (Pinch Point 10); confirm `R CMD check` and the current suite pass. Commit as an isolated fix. *(Also: clear the stray `.git/index.lock` locally — see Section 11.)*
2. **Numerical core.** Implement the base-R `drfunc` port + residual Kummer-stabilized ₁F₁ + error-driven crossover; add the Rmpfr-oracle tests. No package wiring yet — prove the numerics first.
3. **Model registration.** Add `exact_beta_poisson_response`, wire `model_probability`, `model_parameters`, `default_start`, `model_label`, `fit_dose_response` `match.arg`, and the `effective_dose` root-solver.
4. **Fitting/uncertainty.** Continuous evaluator on the fitting path (Pinch Point 3); verify BFGS convergence and SEs on real species; bootstrap plumbing and `bootstrap_prediction_matrix` exact branch.
5. **Comparison/assessment/consensus.** Generalize to three models (Pinch Points 7, 10) with two-model regression locks.
6. **Docs + benchmarks.** Roxygen/`NAMESPACE`, README, `NEWS.md`; benchmark 1k/10k sequential vs mirai on one species; record timings.
7. **Full run + review.** Fit the QMRA-wiki species set; sanity-check exact-vs-approximate divergence as a heterogeneity diagnostic.

---

## 10. Risks and open questions

- **Paper vs code crossover discrepancy (open).** Reconcile the cited "Pfaff below 33k gec / ₁F₁ above 33k gec" description against our exact-at-small-dose / approximate-at-large-dose implementation. Confirm parameterization and whether the reference concerns ₂F₁. *Owner: collaborators.*
- **Three-model comparison semantics (needs sign-off).** Pinch Point 7 changes a core statistical workflow; confirm the exponential-baseline-nested + IC-ranking policy is acceptable and that reproducing the two-model outputs is sufficient.
- **Nudge bias under 10k bootstraps (to measure).** Quantify whether the near-integer nudge perturbs bootstrap CIs; adopt the fitting-path Kummer reindexing if so (Pinch Point 4).
- **Optimizer stability on real data (to measure).** Confirm BFGS + continuous evaluator converges reliably across all species; have Nelder–Mead as a per-species fallback if not.
- **Performance headroom (to benchmark).** Validate that pure-R + mirai meets a practical wall-clock target at 10k replicates before committing to the full species run.
- **Legacy build availability.** The `.dll`/`.so` oracle is most useful if we can run it; if not, we rely on Rmpfr + reference points, which is sufficient but loses the guard-for-guard cross-check.

---

## 11. Recommendations at a glance

| Pinch point | Recommendation |
|---|---|
| 1. Parameterization | Fit in **(α, β)**; derive N50 via `effective_dose`. |
| 2. ₁F₁ engine | **Base-R port** of `drfunc` + Kummer-stabilized residual series; **Rmpfr in tests only**. |
| 3. Fitting continuity | **Continuous evaluator for fitting/bootstrap**; fast guarded cascade for prediction/Monte-Carlo. |
| 4. Near-integer params | Nudge on fast path; **Kummer reindex (no nudge) on fitting path**. |
| 5. Crossover | **Error-driven switch** (primary) validated against legacy thresholds; report realized switch dose. |
| 6. Effective dose | **Monotone `uniroot`** in log-dose with data-informed brackets. |
| 7. Three-model comparison | **Exponential-baseline nested tests + AIC/BIC ranking**; two-model path regression-locked. |
| 8. Performance | **Vectorized pure-R + `mirai`** parallel + warm starts; benchmark early. |
| 9. Dependencies | **No new runtime deps**; `Rmpfr` in **`Suggests`**; no compiled code/system libs. |
| 10. Consensus | **Restore first**, then broaden voting to three models. |

---

## Appendix A — Git state note (environmental)

During setup, the repository's `HEAD` was found pointing at a then-nonexistent `dev-hypergeo` ref ("your current branch appears to be broken"). `HEAD` has been repointed to `feature/package-dose-response-models`, which now resolves correctly, so the branch is no longer broken. One residual item could **not** be completed from the working environment: a stray zero-byte **`.git/index.lock`** remains and the mounted filesystem refuses its removal ("Operation not permitted"), which leaves the *index* stale and makes `git status` report spurious modifications. **On a normal local checkout, delete `.git/index.lock` and run `git reset` (mixed; it does not touch the working tree) to refresh the index.** All working files are intact on disk. A `dev-hypergeo` branch now exists at a real commit and can be used for this work if desired.

## Appendix B — References

- Haas, Rose, Gerba, *Quantitative Microbial Risk Assessment* (single-hit theory; beta-Poisson).
- Teunis & Havelaar (2000), "The Beta Poisson dose-response model is not a single-hit model" (exact vs approximate; ₁F₁ form).
- Moshier, *Cephes Mathematical Library* (`hyperg.c`, confluent hypergeometric implementation and its accuracy caveats).
- Weir et al. (2017), microbial dose-response modeling tool (`dev/Weir et al 2017 …pdf`).
- Legacy implementation: `dev/hgdr.new/` (C core + `drfunc` cascade), `dev/hypergeo/f1.R` (R interface), `dev/hypergeo/predications.R` (performance profiling).
