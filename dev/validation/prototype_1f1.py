"""
Prototype + validation of a numerically robust evaluator for the exact
beta-Poisson response:   P(d) = 1 - M(alpha, alpha+beta, -d)
where M = 1F1 (Kummer confluent hypergeometric).

The decision-relevant quantity is the SURVIVAL  S = 1 - P = M(a, a+b, -d),
because the binomial log-likelihood uses log(1-P)=log(S) for non-responders.
We therefore measure relative error in S (=M), the hard tail.

Oracle: mpmath.hyp1f1 at 50 digits.
Methods under test:
  (naive)   direct double-precision power series           -> expected to fail at large d
  (kummer)  log-space all-positive Kummer-transformed sum  -> M(a,b,-d)=e^{-d}M(b-a,b,d)
  (asymp)   large-d asymptotic (Gamma-ratio * d^-a * 2F0)  -> fast, accurate for large d
  (legacy)  1-(1+d/b)^-a  (the drfunc approximation branch) -> check where it breaks
  (select)  error-driven switch kummer<->asymp             -> the proposed fitting-path evaluator
"""
import math
import mpmath as mp
mp.mp.dps = 50

# ---------- oracle ----------
def S_oracle(a, b, d):
    # S = M(a, a+b, -d) at high precision
    return mp.hyp1f1(a, a+b, -mp.mpf(d))

# ---------- naive double-precision series (the thing that fails) ----------
def S_naive(a, b, d):
    A, B, z = a, a+b, -d
    term = 1.0; s = 1.0
    for n in range(0, 100000):
        term *= (A+n)/(B+n) * z/(n+1)
        s += term
        if abs(term) < 1e-18*abs(s) and n > abs(z):
            break
    return s

# ---------- log-space positive series for M(A,B,x), x>0, A,B>0 (all terms +) ----------
def logM_pos(A, B, x, maxit=2_000_000):
    # returns log( sum_n t_n ), t_0=1, ratio t_{n+1}/t_n = (A+n)/(B+n)*x/(n+1)
    log_t = 0.0          # log of current term (t_0 = 1)
    log_s = 0.0          # log of running sum
    peak = 0.0
    for n in range(0, maxit):
        log_t += math.log((A+n)/(B+n)) + math.log(x) - math.log(n+1)
        # log_s = logaddexp(log_s, log_t)
        m = max(log_s, log_t)
        log_s = m + math.log(math.exp(log_s-m) + math.exp(log_t-m))
        if log_t > peak: peak = log_t
        if log_t < peak - 40.0 and n > x:   # term negligible and past the peak
            break
    return log_s

def S_kummer(a, b, d):
    # M(a,b,-d) = e^{-d} M(b-a, b, d);  b-a = beta, b = alpha+beta  (both > 0)
    logM = logM_pos(b, a+b, d)      # A=b-a=beta, B=a+b, x=d
    return math.exp(-d + logM)

# ---------- large-d asymptotic:  M(a,b,-d) ~ G(b)/G(b-a) d^-a * 2F0(a,a-b+1;;1/d) ----------
def S_asymp(a, b, d, kmax=60):
    B = a+b
    log_lead = math.lgamma(B) - math.lgamma(B-a) - a*math.log(d)   # log[ G(b)/G(b-a) d^-a ]
    # 2F0 asymptotic (divergent) series: t_0=1, t_{k+1}/t_k = (a+k)(a-B+1+k)/((k+1)) * (1/d)
    term = 1.0; s = 1.0; prev = float('inf')
    for k in range(0, kmax):
        r = (a+k)*(a-B+1+k)/((k+1)*d)
        nt = term*r
        if abs(nt) > abs(term):      # divergent series: stop at smallest term
            break
        term = nt; s += term
    return math.exp(log_lead)*s

def S_legacy(a, b, d):
    # drfunc approximation branch: P = 1-(1+d/b)^-a  ->  S = (1+d/b)^-a
    return (1.0 + d/b)**(-a)

# ---------- proposed selector (fitting-path evaluator) ----------
def S_select(a, b, d):
    if d < 1e-4:
        return 1.0 - d*a/(a+b)                 # linear limit (P small)
    # try asymptotic; if its smallest-term (truncation) error is small, trust it
    B = a+b
    # crude asymptotic reliability: last included term magnitude
    log_lead = math.lgamma(B) - math.lgamma(B-a) - a*math.log(d)
    term = 1.0; s = 1.0; smallest = 1.0
    for k in range(0, 60):
        r = (a+k)*(a-B+1+k)/((k+1)*d)
        nt = term*r
        if abs(nt) > abs(term): break
        term = nt; s += term; smallest = abs(term)
    if smallest < 1e-12:                       # asymptotic converged well
        return math.exp(log_lead)*s
    return S_kummer(a, b, d)                    # else all-positive Kummer series

# ---------------- validation grid ----------------
alphas = [0.05, 0.1, 0.265, 0.5, 1.0, 2.0]
betas  = [0.1, 1.0, 5.0, 50.0, 1000.0, 1e5]
doses  = [1e-3, 1.0, 10.0, 100.0, 1e3, 1e4, 33000.0, 1e5]

def relerr(x, ref):
    ref = float(ref)
    if ref == 0: return abs(x)
    return abs((x-ref)/ref)

import collections
worst = collections.defaultdict(float)
fail_examples = collections.defaultdict(list)

for a in alphas:
    for b in betas:
        for d in doses:
            ref = S_oracle(a,b,d)
            for name, fn in [("naive",S_naive),("kummer",S_kummer),
                             ("asymp",S_asymp),("legacy",S_legacy),("select",S_select)]:
                try:
                    val = fn(a,b,d)
                    e = relerr(val, ref)
                except Exception as ex:
                    e = float('inf')
                if e > worst[name]:
                    worst[name] = e
                if e > 1e-6 and len(fail_examples[name]) < 6:
                    fail_examples[name].append((a,b,d,e))

print("=== worst-case relative error in S=1-P (=M) over the grid ===")
for name in ["naive","kummer","asymp","legacy","select"]:
    print(f"  {name:8s}: {worst[name]:.3e}")

print("\n=== where each method exceeds 1e-6 (a, beta, dose, relerr) ===")
for name in ["naive","kummer","asymp","legacy","select"]:
    print(f"  [{name}]")
    if not fail_examples[name]:
        print("     (none)")
    for (a,b,d,e) in fail_examples[name]:
        print(f"     a={a:<6} beta={b:<8} dose={d:<9} relerr={e:.2e}")
