# Mirror the base-R exact_beta_poisson.R logic 1:1 to catch transcription bugs.
import math, mpmath as mp
mp.mp.dps = 50
def ref(a,b,d): return float(mp.hyp1f1(a, a+b, -mp.mpf(d)))

def kummer(alpha, beta, dose, reltol=1e-16):   # mirrors .ebp_survival_kummer (scalar d)
    A=beta; B=alpha+beta
    n_max=int(math.ceil(dose)+40*math.ceil(math.sqrt(dose)+1))
    term=1.0; m=1.0
    for n in range(1, n_max+1):
        term *= ((A+n-1)/(B+n-1))*(dose/n)
        m += term
        if term < reltol*m and n > dose: break
    return math.exp(-dose)*m

def asymp(alpha, beta, dose, kmax=80):          # mirrors .ebp_survival_asymp (scalar d)
    B=alpha+beta
    log_lead=math.lgamma(B)-math.lgamma(B-alpha)-alpha*math.log(dose)
    term=1.0; s=1.0; active=True
    for k in range(kmax):
        r=(alpha+k)*(alpha-B+1+k)/((k+1)*dose)
        nt=term*r
        if abs(nt)>abs(term):
            active=False
        if not active: break
        term=nt; s+=term
    return math.exp(log_lead)*s

def survival(dose, alpha, beta):                # mirrors exact_beta_poisson_survival
    if dose < 1e-4:  return 1 - dose*alpha/(alpha+beta)
    if dose >= 600:  return asymp(alpha,beta,dose)
    return kummer(alpha,beta,dose)

worst=0.0; args=None
for a in [0.05,0.1,0.265,0.5,1.0,2.0]:
    for b in [0.1,0.7,1.0,5.0,50.0,1000.0,1e5]:
        for d in [1e-3,1.0,10.0,100.0,599.0,600.0,1e3,1e4,33000.0,1e5]:
            e=abs((survival(d,a,b)-ref(a,b,d))/ref(a,b,d))
            if e>worst: worst=e; args=(a,b,d)
print(f"base-R-logic worst rel err in S: {worst:.3e}  at (alpha,beta,dose)={args}")
# continuity across the 600 handoff (no jump that would kink the likelihood)
for a,b in [(0.265,0.7),(0.5,5.5),(1.0,123.4)]:
    lo=survival(599.999,a,b); hi=survival(600.001,a,b)
    print(f"  handoff a={a} b={b}: S(599.999)={lo:.6e}  S(600.001)={hi:.6e}  jump={abs(lo-hi)/hi:.2e}")
