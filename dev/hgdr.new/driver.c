
#include <math.h>
#include <stdio.h>
#include "mth_wrap.h"

void do_stress(int iter);
double guppie (double a, double b, double x);

/******************************************/
int main()
{

fprintf(stdout, "pow = %s\n"
	, (mth_pow==WRAP_pow) ? "local" : "system-supplied" );
fprintf(stdout, "fabs = %s\n"
	, (mth_fabs==WRAP_fabs) ? "local" : "system-supplied" );
fprintf(stdout, "gamma = %s\n"
	, (mth_gamma==WRAP_gamma) ? "local" : "system-supplied" );
fprintf(stdout, "lgamma = %s\n"
	, (mth_lgamma==WRAP_lgamma) ? "local" : "system-supplied" );

#if 1||WANT_JOSTI_BAND
{
int ii;
double aa, xx, y0, y1, diff;
fprintf(stdout, "\npow : ref | mtm | diff\n" );
aa = 3.0;
for(ii = -10; ii <= 10; ii++ ) {
	xx = ii;
	y0 = pow(aa, xx);
	y1 = mth_pow(aa, xx);
	diff = y1 - y0;
	fprintf(stdout, "%lf : %lf | %lf | %lf\n"
	 , xx ,y0 , y1, diff );
	}

fprintf(stdout, "\npow : pow | powi | diff\n" );
aa = 3.0;
for(ii = -10; ii <= 10; ii++ ) {
	xx = ii;
	y0 = pow(aa, xx);
	y1 = mth_powi(aa, ii);
	diff = y1 - y0;
	fprintf(stdout, "%lf : %lf | %lf | %lf\n"
	 , xx ,y0 , y1, diff );
	}

fprintf(stdout, "\ngamma : ref | mth | diff\n" );
aa = 3.0;
for(ii = -10; ii <= 10; ii++ ) {
	xx = mth_powi(aa, ii);
	y0 = tgamma(xx);
	y1 = mth_gamma(xx);
	diff = y1 - y0;
	fprintf(stdout, "%lf : %lf | %lf | %lf\n"
	 , xx ,y0 , y1, diff );
	}

fprintf(stdout, "\nlgamma : ref | mth | diff\n" );
aa = 3.0;
for(ii = -10; ii <= 10; ii++ ) {
	xx = mth_powi(aa, ii);
	y0 = lgamma(xx);
	y1 = mth_lgamma(xx);
	diff = y1 - y0;
	fprintf(stdout, "%lf : %lf | %lf | %lf\n"
	 , xx ,y0 , y1, diff );
	}


fprintf(stdout, "\nhyperg(: (1,1,x) | (2,2,x)\n" );
aa = 3.0;
for(ii = -20; ii <= 20; ii++ ) {
	xx = mth_powi(aa, ii);
	y0 = guppie(1.0 , 1.0 , xx);
	y1 = guppie(2.0 , 2.0 , xx);
	fprintf(stdout, "%lf : %lf | %lf\n"
	 , xx ,y0 , y1);
	}
}
#endif /* want_JOSTI_BAND */

do_stress(10*1000*1000);

return 0;
}

/******************************************/
void do_stress(int iter)
{
int idx;
int ii;
double xx, yy;

for (idx =0; idx < iter; idx++) {
	ii = rand();
	if (ii < 0) xx = (double) 123/ ii;
	else xx = (double) ii /123 ;
	yy = guppie(1.0 , 1.000001 , xx);
	}
}

double guppie (double a, double b, double x)
{
return 1.0 - mth_hyperg(a, a+b, -x) ;
}
/******************************************/
