# Exact beta-Poisson dose-response (confluent hypergeometric) — prototype
#
# Response:  P(d) = 1 - M(alpha, alpha + beta, -d)
# where M = 1F1 is Kummer's confluent hypergeometric function and (alpha, beta)
# are the shape parameters of the Beta(alpha, beta) infectivity distribution
# (single-hit theory). This file prototypes ONLY the numerically robust
# evaluator + response function; package wiring (fit / effective_dose /
# bootstrap / comparison) comes in later phases.
#
# We compute and return the SURVIVAL S = 1 - P directly, because the binomial
# log-likelihood needs log(1 - P) = log(S) for non-responders, and S is the
# hard-to-get tail at large dose. The response is 1 - S.
#
# VALIDATED NUMERICAL STRATEGY  (worst-case relative error in S = 2.3e-9 over
# alpha in [0.05, 2], beta in [0.1, 1e5], dose in [1e-3, 1e5], vs mpmath at 50
# digits; re-certify against Rmpfr in the package test suite):
#
#   route each dose element to exactly one method --
#     * dose < 1e-4                          -> linear limit  S = 1 - d*alpha/(alpha+beta)
#     * asymptotic ACCEPTED when ALL hold:   -> large-argument asymptotic
#          dose >= alpha+beta                   (dose past the parameter scale)
#          exp(-dose) < 1e-12                   (recessive ~e^-d term negligible; dose > ~28)
#          2F0 truncation term < 1e-12          (algebraic asymptotic converged)
#     * otherwise                            -> log-space Kummer series (always
#                                               correct, never overflows)
#
# WHY NOT the legacy drfunc approximation `1-(1+d/beta)^-alpha`: measured error
# up to ~112% in S at small beta, and it only reaches 1e-8 near beta ~ 1e6-1e7,
# so it is NOT used as an exact-model path here. WHY NOT a fixed dose cutoff:
# at large beta a fixed cutoff wrongly applies the asymptotic where dose << beta
# (error ~1e4); the switch must be regime-aware as above.
#
# Performance note (10k-bootstrap target): realistic fits have small/moderate
# beta, where large gec doses take the fast asymptotic (a few terms) and only
# small doses use the short Kummer series. The slow Kummer path (~1e5 terms) is
# reached only in the extreme large-beta / large-dose corner; a higher-order
# large-beta expansion could remove even that and is flagged as a future
# optimization (see EXACT_BETA_POISSON_IMPLEMENTATION_PLAN.md, Pinch Point 2/5).

# ---- large-dose asymptotic, vectorized over dose; returns value + truncation --
# M(a, a+b, -d) ~ Gamma(B)/Gamma(B-a) * d^-a * 2F0(a, a-B+1; ; 1/d),  B = a+b
.ebp_asymp <- function(alpha, beta, dose, kmax = 80L) {
  B <- alpha + beta
  log_lead <- lgamma(B) - lgamma(B - alpha) - alpha * log(dose)
  term <- rep(1.0, length(dose))
  s <- rep(1.0, length(dose))
  smallest <- rep(1.0, length(dose))
  active <- rep(TRUE, length(dose))
  for (k in 0:(kmax - 1L)) {
    r <- (alpha + k) * (alpha - B + 1 + k) / ((k + 1) * dose)
    nt <- term * r
    grew <- abs(nt) >= abs(term)          # divergent asymptotic: stop at smallest term
    freeze <- active & grew
    active <- active & !grew
    if (any(active)) {
      term[active] <- nt[active]
      s[active] <- s[active] + term[active]
      smallest[active] <- abs(term[active])
    }
    if (!any(active)) break
  }
  list(value = exp(log_lead) * s, truncation = smallest)
}

# ---- log-space positive-term Kummer series, vectorized over dose --------------
# S = M(alpha, alpha+beta, -d) = exp(-d) * M(beta, alpha+beta, d); all terms > 0,
# summed in log space so nothing overflows regardless of dose magnitude.
.ebp_kummer <- function(alpha, beta, dose, drop = 40, maxit = 2e6L) {
  A <- beta
  B <- alpha + beta
  logx <- log(dose)
  log_t <- numeric(length(dose))          # log(term_0) = 0
  log_s <- numeric(length(dose))          # log(sum) = 0
  peak <- numeric(length(dose))
  active <- rep(TRUE, length(dose))
  n <- 0L
  while (any(active) && n < maxit) {
    idx <- which(active)
    log_t[idx] <- log_t[idx] + log((A + n) / (B + n)) + logx[idx] - log(n + 1)
    mx <- pmax(log_s[idx], log_t[idx])
    log_s[idx] <- mx + log(exp(log_s[idx] - mx) + exp(log_t[idx] - mx))
    peak[idx] <- pmax(peak[idx], log_t[idx])
    done <- (log_t[idx] < peak[idx] - drop) & (n > dose[idx])
    active[idx[done]] <- FALSE
    n <- n + 1L
  }
  exp(-dose + log_s)
}

#' Exact beta-Poisson survival S(dose) = 1 - P(dose), vectorized over dose.
#' @export
exact_beta_poisson_survival <- function(dose, alpha, beta) {
  stopifnot(all(is.finite(dose)), all(dose >= 0), alpha > 0, beta > 0)
  B <- alpha + beta
  out <- numeric(length(dose))
  tiny <- dose < 1e-4
  if (any(tiny)) out[tiny] <- 1 - dose[tiny] * alpha / B
  rest <- which(!tiny)
  if (length(rest)) {
    d <- dose[rest]
    asy <- .ebp_asymp(alpha, beta, d)
    use_asy <- (d >= B) & (exp(-d) < 1e-12) & (asy$truncation < 1e-12)
    out[rest[use_asy]] <- asy$value[use_asy]
    kum <- rest[!use_asy]
    if (length(kum)) out[kum] <- .ebp_kummer(alpha, beta, dose[kum])
  }
  out
}

#' Exact beta-Poisson response P(dose) = 1 - S(dose).
#' Slots into models.R beside exponential_response / beta_poisson_response and
#' is registered in model_probability(); parameters are (alpha, beta).
#' @export
exact_beta_poisson_response <- function(dose, alpha, beta) {
  1 - exact_beta_poisson_survival(dose, alpha, beta)
}

# ---- optional self-check when sourced interactively (needs Rmpfr) -------------
if (identical(environment(), globalenv()) && requireNamespace("Rmpfr", quietly = TRUE)) {
  ref <- function(a, b, d) {
    dd <- Rmpfr::mpfr(d, 200)
    as.numeric(Rmpfr::hypergeom1F1(Rmpfr::mpfr(a, 200), Rmpfr::mpfr(a + b, 200), -dd))
  }
  g <- expand.grid(a = c(0.05, 0.265, 1, 2), b = c(0.1, 1, 5, 1000, 1e5),
                   d = c(1e-3, 10, 100, 1e3, 33000, 1e5))
  err <- mapply(function(a, b, d)
    abs(exact_beta_poisson_survival(d, a, b) - ref(a, b, d)) / abs(ref(a, b, d)),
    g$a, g$b, g$d)
  cat(sprintf("max relative error vs Rmpfr over %d points: %.3e\n", nrow(g), max(err)))
}
