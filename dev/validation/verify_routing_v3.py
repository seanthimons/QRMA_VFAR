import math, mpmath as mp
mp.mp.dps=50
def ref(a,b,d): return float(mp.hyp1f1(a,a+b,-mp.mpf(d)))

def logsumexp2(x,y):
    m=max(x,y); return m+math.log(math.exp(x-m)+math.exp(y-m))

def logM_pos(A,B,x,maxit=2_000_000):
    log_t=0.0; log_s=0.0; peak=0.0; n=0
    for n in range(maxit):
        log_t += math.log((A+n)/(B+n))+math.log(x)-math.log(n+1)
        log_s = logsumexp2(log_s,log_t)
        if log_t>peak: peak=log_t
        if log_t<peak-40.0 and n>x: break
    return log_s, n+1

def S_kummer(a,b,d):
    logM,nt=logM_pos(b,a+b,d)         # M(a,b,-d)=e^-d M(b-a,b,d)
    return math.exp(-d+logM), nt

def asymp(a,b,d,kmax=80):
    B=a+b; log_lead=math.lgamma(B)-math.lgamma(B-a)-a*math.log(d)
    term=1.0; s=1.0; smallest=1.0
    for k in range(kmax):
        r=(a+k)*(a-B+1+k)/((k+1)*d); nt=term*r
        if abs(nt)>abs(term): break
        term=nt; s+=term; smallest=abs(term)
    return math.exp(log_lead)*s, smallest

def survival(a,b,d):                  # proposed routing
    if d<1e-4: return 1-d*a/(a+b), 0
    val,smallest=asymp(a,b,d)
    if d>=(a+b) and smallest<1e-12: return val, 1  # asymptotic valid AND reliable
    v,nt=S_kummer(a,b,d); return v, nt

worst=0.0; wargs=None; maxterms=0; slowcases=[]
for a in [0.05,0.1,0.265,0.5,1.0,2.0]:
    for b in [0.1,0.7,1.0,5.0,50.0,1000.0,1e5]:
        for d in [1e-3,1.0,10.0,100.0,600.0,1e3,1e4,33000.0,1e5]:
            v,nt=survival(a,b,d); r=ref(a,b,d)
            e=abs((v-r)/r)
            if e>worst: worst=e; wargs=(a,b,d)
            if nt>maxterms: maxterms=nt
            if nt>2000: slowcases.append((a,b,d,nt))
print(f"worst rel err in S over full grid (incl beta=1e5): {worst:.3e} at {wargs}")
print(f"max Kummer terms used anywhere: {maxterms}")
print("cases needing >2000 Kummer terms (the slow corner):")
for (a,b,d,nt) in slowcases: print(f"   a={a} beta={b:g} dose={d:g} -> {nt} terms")
if not slowcases: print("   (none)")

# realistic-fit speed check: small/moderate beta with large doses -> should be FAST (asymptotic)
print("\nrealistic regime (noro-like small beta, large gec doses):")
for a,b in [(0.2,0.9),(0.5,4.0),(0.04,0.055)]:
    for d in [1e3,1e4,33000,1e5]:
        v,nt=survival(a,b,d)
        tag="asymptotic(fast)" if nt==1 else f"kummer:{nt} terms"
        print(f"   a={a} beta={b} dose={d:g}: {tag}")
