"""
v2: Kummer log-space series is accurate everywhere (v1 showed 2.8e-9 worst).
So the selector switches to the asymptotic ONLY for large dose, where (a) the
asymptotic is valid and (b) the Kummer series would need ~O(d) terms (slow).
Fixes v1 selector bug (asymptotic wrongly trusted at small dose when beta~1).
"""
import math, time
import mpmath as mp
mp.mp.dps = 50

def S_oracle(a,b,d): return mp.hyp1f1(a, a+b, -mp.mpf(d))

def logM_pos(A,B,x,maxit=2_000_000):
    log_t=0.0; log_s=0.0; peak=0.0; n=0
    for n in range(maxit):
        log_t += math.log((A+n)/(B+n)) + math.log(x) - math.log(n+1)
        m=max(log_s,log_t); log_s=m+math.log(math.exp(log_s-m)+math.exp(log_t-m))
        if log_t>peak: peak=log_t
        if log_t < peak-40.0 and n>x: break
    return log_s, n+1

def S_kummer(a,b,d):
    logM,nterms = logM_pos(b, a+b, d)     # M(a,b,-d)=e^-d M(b-a,b,d); b-a=beta,b=a+b
    return math.exp(-d+logM), nterms

def asymp(a,b,d,kmax=80):
    B=a+b
    log_lead=math.lgamma(B)-math.lgamma(B-a)-a*math.log(d)
    term=1.0; s=1.0; smallest=1.0
    for k in range(kmax):
        r=(a+k)*(a-B+1+k)/((k+1)*d)
        nt=term*r
        if abs(nt)>abs(term): break        # divergent -> stop at smallest term
        term=nt; s+=term; smallest=abs(term)
    return math.exp(log_lead)*s, smallest

D_SWITCH = 2000.0   # above this dose, prefer asymptotic if it is reliable
def S_select(a,b,d):
    if d < 1e-4:
        return 1.0 - d*a/(a+b)
    if d > D_SWITCH:
        val,smallest = asymp(a,b,d)
        if smallest < 1e-12:               # asymptotic truncation error tiny
            return val
    return S_kummer(a,b,d)[0]

alphas=[0.05,0.1,0.265,0.5,1.0,2.0]
betas =[0.1,1.0,5.0,50.0,1000.0,1e5]
doses =[1e-3,1.0,10.0,100.0,1e3,1e4,33000.0,1e5,1e6]

def relerr(x,ref):
    ref=float(ref); return abs(x) if ref==0 else abs((x-ref)/ref)

worst_sel=0.0; worst_kum=0.0; fails=[]; maxterms=0; slow=[]
t0=time.time()
for a in alphas:
    for b in betas:
        for d in doses:
            ref=S_oracle(a,b,d)
            es=relerr(S_select(a,b,d),ref)
            if es>worst_sel: worst_sel=es
            if es>1e-8 and len(fails)<12: fails.append((a,b,d,es))
            # kummer accuracy + cost (skip d=1e6 cost: ~1e6 terms)
            if d<=1e5:
                val,nt=S_kummer(a,b,d)
                ek=relerr(val,ref)
                if ek>worst_kum: worst_kum=ek
                if nt>maxterms: maxterms=nt
                if nt>50000: slow.append((a,b,d,nt))
print(f"grid points: {len(alphas)*len(betas)*len(doses)}   elapsed {time.time()-t0:.1f}s")
print(f"worst rel err  SELECT : {worst_sel:.3e}")
print(f"worst rel err  KUMMER : {worst_kum:.3e}  (max terms used: {maxterms})")
print("select failures >1e-8:", fails if fails else "(none)")
print("kummer term-count >50k (cost driver for the asymptotic switch):")
for (a,b,d,nt) in slow[:8]: print(f"   a={a} beta={b} dose={d:g} -> {nt} terms")

# realized asymptotic switch dose for a representative small-beta fit
print("\nrealized crossover check (a=0.265,beta=5): dose where asymptotic first reliable")
a,b=0.265,5.0
for d in [500,1000,2000,5000,1e4,3.3e4]:
    _,sm=asymp(a,b,d)
    print(f"   dose={d:>7g}  asymptotic smallest-term={sm:.1e}  reliable={sm<1e-12}")

# ---- bridging test: is the asymptotic accurate for ALL d>=600 (non-integer beta)? ----
print("\n=== handoff test: asymptotic rel-err for non-integer beta at moderate/large dose ===")
worst=0.0
for a in [0.05,0.265,0.5,1.0,2.0]:
    for b in [0.1,0.7,2.3,5.5,17.9,123.4]:     # deliberately non-integer beta
        for d in [600,700,800,1000,2000,1e4,1e5]:
            ref=float(S_oracle(a,b,d))
            val,_=asymp(a,b,d)
            e=abs((val-ref)/ref) if ref else abs(val)
            if e>worst: worst=e
            if e>1e-8:
                print(f"   GAP a={a} beta={b} dose={d:g} relerr={e:.2e}")
print(f"worst asymptotic rel-err for d>=600, non-integer beta: {worst:.3e}")
print("=> if <1e-8, base-R hot path can use: direct cumprod Kummer (d<700) + asymptotic (d>=700)")
