/* To launch this program from within R use:
First compile the code with
R CMD SHLIB hypg.c hyperg.c gamma.c pow.c powi.c polevl.c floor.c const.c mtherr.c
and load the library (hypg.so) with
dyn.load("~/stat/dr/hgdr.new/hypg.so")
or whichever path exists to the produced lib object
then define the loaded external function as
 bp <- function(a,b,d)
 .C("dr1f1",
  as.double(a),
  as.double(b),
  as.double(d),
  pinf = double(1))$pinf
and call as
bp(a,b,dose)
 */

#include "mth_wrap.h"
#include <R.h>
#include <Rmath.h>

double drfunc(double a, double b, double dose)
{
  if (dose<1e-4) return dose*a/(a+b);
  if (a>1e3 && b<a/100) return 1-WRAP_exp(-dose);
  if (a>1e2 && b>1e5) return 1-WRAP_exp(-dose*a/b);
  if (a>1e1 && b>1e5 && dose*a/b>10.0) return 1-WRAP_pow(1+dose/b, -a);
  if (a>1.0 && b>20*a && dose>10.0) return 1-WRAP_pow(1+dose/b, -a);
  if (a>1.0 && b>a && dose>50.0) return 1-WRAP_pow(1+dose/b, -a);
  if (a>1.0 && b<a && dose>20.0) return 1-WRAP_pow(1+dose/b, -a);
  if (a<1.0 && b>50*a) return 1-WRAP_pow(1+dose/b, -a);
  if (a<0.1 && b>20*a) return 1-WRAP_pow(1+dose/b, -a);
  if (round(a)-a<1e-4) a=1.0001*a;
  if (round(b)-b<1e-4) b=1.0001*b;
  return 1-WRAP_hyperg(a, a+b, -dose);
}

void dr1f1(double *a, double *b, double *dose, double *pinf)
{
  *pinf = drfunc(*a, *b, *dose);
}
