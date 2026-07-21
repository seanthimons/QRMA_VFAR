import math, mpmath as mp
mp.mp.dps=50
def ref(a,b,d): return float(mp.hyp1f1(a,a+b,-mp.mpf(d)))
def logsumexp2(x,y):
    m=max(x,y); return m+math.log(math.exp(x-m)+math.exp(y-m))
def logM_pos(A,B,x,maxit=2_000_000):
    log_t=0.0; log_s=0.0; peak=0.0; n=0
    for n in range(maxit):
        log_t += math.log((A+n)/(B+n))+math.log(x)-math.log(n+1)
        log_s=logsumexp2(log_s,log_t)
        if log_t>peak: peak=log_t
        if log_t<peak-40.0 and n>x: break
    return log_s,n+1
def S_kummer(a,b,d):
    logM,nt=logM_pos(b,a+b,d); return math.exp(-d+logM),nt
def asymp(a,b,d,kmax=80):
    B=a+b; log_lead=math.lgamma(B)-math.lgamma(B-a)-a*math.log(d)
    term=1.0; s=1.0; smallest=1.0
    for k in range(kmax):
        r=(a+k)*(a-B+1+k)/((k+1)*d); nt=term*r
        if abs(nt)>abs(term): break
        term=nt; s+=term; smallest=abs(term)
    return math.exp(log_lead)*s, smallest

def make_survival(dfloor, recessive_tol):
    def survival(a,b,d):
        if d<1e-4: return 1-d*a/(a+b),0
        val,smallest=asymp(a,b,d)
        if d>=(a+b) and d>dfloor and math.exp(-d)<recessive_tol and smallest<1e-12:
            return val,1
        return S_kummer(a,b,d)
    return survival

grid=[(a,b,d) for a in [0.05,0.1,0.265,0.5,1.0,2.0]
              for b in [0.1,0.7,1.0,5.0,50.0,1000.0,1e5]
              for d in [1e-3,1.0,10.0,30.0,100.0,600.0,1e3,1e4,33000.0,1e5]]
refs={(a,b,d):ref(a,b,d) for (a,b,d) in grid}

for dfloor,rtol in [(28,1e-12),(35,1e-14),(40,1e-16)]:
    surv=make_survival(dfloor,rtol)
    worst=0.0;wargs=None;fastbig=0;slowbig=0
    for (a,b,d) in grid:
        v,nt=surv(a,b,d); e=abs((v-refs[(a,b,d)])/refs[(a,b,d)])
        if e>worst: worst=e;wargs=(a,b,d)
        if d>=1e3:                          # large-dose cases: want asymptotic(fast)
            if nt==1: fastbig+=1
            else: slowbig+=1
    print(f"dfloor={dfloor} rtol={rtol:g}: worst={worst:.2e} at {wargs} | large-dose fast={fastbig} slow={slowbig}")
