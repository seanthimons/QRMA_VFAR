import math, mpmath as mp
mp.mp.dps=50
def ref(a,b,d): return float(mp.hyp1f1(a,a+b,-mp.mpf(d)))

def asymp(a,b,dose,kmax=80):            # mirrors .ebp_asymp per element
    B=a+b; log_lead=math.lgamma(B)-math.lgamma(B-a)-a*math.log(dose)
    term=1.0;s=1.0;smallest=1.0;active=True
    for k in range(kmax):
        r=(a+k)*(a-B+1+k)/((k+1)*dose); nt=term*r
        grew=abs(nt)>=abs(term)
        if active and grew: active=False
        elif active: term=nt; s+=term; smallest=abs(term)
        if not active: break
    return math.exp(log_lead)*s, smallest

def kummer(a,b,dose,drop=40,maxit=2_000_000):   # mirrors .ebp_kummer per element
    A=b;B=a+b;logx=math.log(dose);log_t=0.0;log_s=0.0;peak=0.0;n=0
    while n<maxit:
        log_t += math.log((A+n)/(B+n))+logx-math.log(n+1)
        m=max(log_s,log_t); log_s=m+math.log(math.exp(log_s-m)+math.exp(log_t-m))
        if log_t>peak: peak=log_t
        if (log_t<peak-drop) and (n>dose): break
        n+=1
    return math.exp(-dose+log_s)

def survival(a,b,dose):                 # mirrors exact_beta_poisson_survival per element
    B=a+b
    if dose<1e-4: return 1-dose*a/B
    val,trunc=asymp(a,b,dose)
    if (dose>=B) and (math.exp(-dose)<1e-12) and (trunc<1e-12): return val
    return kummer(a,b,dose)

worst=0.0;wargs=None
for a in [0.05,0.1,0.265,0.5,1.0,2.0]:
    for b in [0.1,0.7,1.0,5.0,50.0,1000.0,1e5]:
        for d in [1e-3,1.0,10.0,30.0,100.0,600.0,1e3,1e4,33000.0,1e5]:
            e=abs((survival(a,b,d)-ref(a,b,d))/ref(a,b,d))
            if e>worst: worst=e;wargs=(a,b,d)
print(f"FINAL worst rel err in S (S=1-P) over 420 pts incl beta=1e5: {worst:.3e} at {wargs}")

# monotonicity + range sanity for a representative fit (P must rise 0->1 in dose)
a,b=0.265,5.0
P=[1-survival(a,b,d) for d in [0.01,1,10,100,1000,1e4,1e5]]
print("P(dose) monotone increasing in [0,1]:", all(0<=p<=1 for p in P) and all(x<y for x,y in zip(P,P[1:])))
print("  sample P:", [f"{p:.4f}" for p in P])
