import math, mpmath as mp
mp.mp.dps=50
def ref(a,b,d): return float(mp.hyp1f1(a,a+b,-mp.mpf(d)))

# candidate large-b (large-beta) closed forms for S=M(a,a+b,-d)
def approx_beta(a,b,d):   return (1+d/b)**(-a)          # legacy: denom beta
def approx_B(a,b,d):      return (1+d/(a+b))**(-a)      # denom alpha+beta

# large-b asymptotic with first correction (DLMF 13.8.4-style, leading + 1/B term)
def largeB(a,b,d,order=3):
    B=a+b; z=-d
    # M(a,B,z) ~ sum_k C_k ; use ratio form of the (1 - z/B)^-a expansion w/ corrections
    base=(1 - z/B)**(-a)
    # first correction term ~ -a(a+1)/2 * z^2/B^2 * (1-z/B)^(-a-2) ... approximate numerically
    # (kept simple; we mainly test whether base alone suffices at various B)
    return base

print("=== approximation error vs beta (find where each closed form is < 1e-8) ===")
for a in [0.05,0.265,1.0,2.0]:
    for b in [1e3,1e4,1e5,1e6,1e7]:
        # worst over a dose sweep
        wa=wb=wL=0.0
        for d in [1,10,100,1e3,1e4,1e5,1e6]:
            r=ref(a,b,d)
            wa=max(wa,abs(approx_beta(a,b,d)-r)/r)
            wb=max(wb,abs(approx_B(a,b,d)-r)/r)
            wL=max(wL,abs(largeB(a,b,d)-r)/r)
        print(f"a={a:<5} b={b:>7g}  (1+d/beta)^-a:{wa:.1e}  (1+d/(a+b))^-a:{wb:.1e}  largeB:{wL:.1e}")

print("\n=== asymptotic reliability (smallest-term<1e-12) vs (beta,dose): when is large-d asymptotic usable ===")
def asymp_ok(a,b,d):
    B=a+b; term=1.0; smallest=1.0
    for k in range(80):
        r=(a+k)*(a-B+1+k)/((k+1)*d); nt=term*r
        if abs(nt)>abs(term): break
        term=nt; smallest=abs(term)
    return smallest<1e-12
for b in [0.5,5,50,1000,1e5]:
    row=[]
    for d in [100,1e3,1e4,1e5,1e6]:
        row.append(f"d={d:g}:{'Y' if asymp_ok(0.265,b,d) else 'n'}")
    print(f"  beta={b:<7g} "+"  ".join(row))
