/*							hyperg.c
 *	Confluent hypergeometric function
 * SYNOPSIS:
 * double a, b, x, y, hyperg();
 * y = hyperg( a, b, x );
 * DESCRIPTION:
 * Computes the confluent hypergeometric function
 *                          1           2
 *                       a x    a(a+1) x
 *   F ( a,b;x )  =  1 + ---- + --------- + ...
 *  1 1                  b 1!   b(b+1) 2!
 * Many higher transcendental functions are special cases of
 * this power series.
 * As is evident from the formula, b must not be a negative
 * integer or zero unless a is an integer with 0 >= a > b.
 * The routine attempts both a direct summation of the series
 * and an asymptotic expansion.  In each case error due to
 * roundoff, cancellation, and nonconvergence is estimated.
 * The result with smaller estimated error is returned.
 * ACCURACY:
 *
 * Tested at random points (a, b, x), all three variables
 * ranging from 0 to 30.
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    DEC       0,30         2000       1.2e-15     1.3e-16
 *    IEEE      0,30        30000       1.8e-14     1.1e-15
 * Larger errors can be observed when b is near a negative
 * integer or zero.  Certain combinations of arguments yield
 * serious cancellation error in the power series summation
 * and also are not in the region of near convergence of the
 * asymptotic series.  An error message is printed if the
 * self-estimated relative error is greater than 1.0e-12.
 */

/*							hyperg.c
Cephes Math Library Release 2.1:  November, 1988
Copyright 1984, 1987, 1988 by Stephen L. Moshier
Direct inquiries to 30 Frost Street, Cambridge, MA 02140
*/

#include <stdio.h>

#include "mconf.h"
#include "protos.h"
#include "const.h"


static double hy1f1p( double a, double b, double x, double *err );
static double hy1f1a( double a, double b, double x, double *err );
static double hyp2f0( double a, double b, double x, int type, double *err );

double mth_hyperg( double a, double b, double x)
{
double asum, psum, acanc, pcanc;

/* See if a Kummer transformation will help */
/* Note AvK. A Kummer transformation will NOT help.
** In this form, It leads to infinite recursion.
*/
#if (0||WANT_KUMMER)
{
double temp ;
temp = b - a;
if ( WRAP_fabs(temp) < WRAP_fabs(b) )
	return WRAP_exp(x) * mth_hyperg( temp, b, -x ) ;
}
#endif

psum = hy1f1p( a, b, x, &pcanc );
if ( pcanc < 1.0e-15 )
	goto done;

/* try asymptotic series */
asum = hy1f1a( a, b, x, &acanc );

/* Pick the result with less estimated error */
if ( acanc < pcanc )
	{
	pcanc = acanc;
	psum = asum;
	}

done:
#if WANT_DUMP_INCIDENTS
if ( pcanc > 1.0e-12 ) {
   mtherr( "mth_hyperg", MTHE_PLOSS );
   printf("hy1f1(a=%f, b=%f; x=%f)=%f  err %1.12f\n",a,b,x,psum,pcanc);
   }
#endif
return psum ;
}

/* Power series summation for confluent hypergeometric function		*/

static double hy1f1p( double a, double b, double x, double *err )
{
#define MTH_MAXITER 5000
double a0, sum, t, u, temp;
double an, bn, maxt, pcanc;
double nn;

/* set up for power series summation */
an = a;
bn = b;
a0 = 1.0;
sum = 1.0;
pcanc = 1.0;	/* estimate 100% error */

nn = 1.0;
maxt = 0.0;

#if 1
{
int clicks;
double frac, whole;
frac = WRAP_modf(b, &whole);
if (whole < 0 && whole + MTH_MAXITER >= 0) {
	if (frac < MTH_MACHEP) goto blowup; /* Estimate B+nn==close-to-zero */
	}

frac = WRAP_modf(a, &whole);
if (whole < 0 && whole + MTH_MAXITER >= 0) {
	if (frac < MTH_MACHEP) clicks = 0 - whole; /* Idem A */
	else clicks = MTH_MAXITER;
	}
else clicks = MTH_MAXITER;
/* now, clicks can be used as iteration count */

for (t = 1.0; t > MTH_MACHEP; ) {
	if (clicks-- <= 0) goto pdone;

	u = x * ( an / (bn * nn) ); /* !Hotspot! */

	/* check for blowup */
	temp = WRAP_fabs(u); /* !Hotspot1! */
	if ( (temp > 1.0 ) && (maxt > (MTH_MAXNUM/temp)) ) /* !Hotspot1! */
		{
		goto blowup;
		}
	a0 *= u;
	sum += a0;
	t = WRAP_fabs(a0);
	if ( t > maxt )
		maxt = t;
	if ( (maxt/WRAP_fabs(sum)) > 1.0e17 ) { /* !Hotspot! */
		goto blowup;
		}
	an += 1.0;
	bn += 1.0;
	nn += 1.0;
	}
goto pdone;
}
#else

if (an > 0 && bn > 0) {
for (t = 1.0; t > MTH_MACHEP; ) {
	/* check bn first since if both
	** an and bn are zero it is a singularity
	*/
	if (nn > MTH_MAXITER)
		goto pdone;
	u = x * ( an / (bn * nn) ); /* !Hotspot! */

	/* check for blowup */
	temp = WRAP_fabs(u); /* !Hotspot1! */
	if ( (temp > 1.0 ) && (maxt > (MTH_MAXNUM/temp)) ) /* !Hotspot1! */
		{
		goto blowup;
		}
	a0 *= u;
	sum += a0;
	t = WRAP_fabs(a0);
	if ( t > maxt )
		maxt = t;
	if ( (maxt/WRAP_fabs(sum)) > 1.0e17 ) { /* !Hotspot! */
		goto blowup;
		}
	an += 1.0;
	bn += 1.0;
	nn += 1.0;
	}
goto pdone;
} else { /* An <= 0 || bn <= 0.0 */
for (t = 1.0; t > MTH_MACHEP; ) {
	/* check bn first since if both
	** an and bn are zero it is a singularity
	*/
	if (bn == 0) {
		mtherr( "mth_hyperg", MTHE_SING );
		sum = MTH_MAXNUM ;
		goto blowup;
		}
	if (an == 0)
		goto pdone ;
	if (nn > MTH_MAXITER)
		goto pdone;
	u = x * ( an / (bn * nn) ); /* !Hotspot! */

	/* check for blowup */
	temp = WRAP_fabs(u); /* !Hotspot1! */
	if ( (temp > 1.0 ) && (maxt > (MTH_MAXNUM/temp)) ) /* !Hotspot1! */
		{
		goto blowup;
		}
	a0 *= u;
	sum += a0;
	t = WRAP_fabs(a0);
	if ( t > maxt )
		maxt = t;
	if ( (maxt/WRAP_fabs(sum)) > 1.0e17 ) { /* !Hotspot! */
		goto blowup;
		}
	an += 1.0;
	bn += 1.0;
	nn += 1.0;
	}
}
#endif

pdone:

/* estimate error due to roundoff and cancellation */
if ( sum != 0.0 )
	maxt /= WRAP_fabs(sum);
maxt *= MTH_MACHEP; 	/* this way avoids multiply overflow */
pcanc = WRAP_fabs( MTH_MACHEP * nn  +  maxt );

blowup:

*err = pcanc;

return sum ;
}

/*							hy1f1a()	*/
/* asymptotic formula for hypergeometric function:
 *
 *        (    -a                         
 *  --    ( |z|                           
 * |  (b) ( -------- 2f0( a, 1+a-b, -1/x )
 *        (  --                           
 *        ( |  (b-a)                      
 *                                x    a-b                     )
 *                               e  |x|                        )
 *                             + -------- 2f0( b-a, 1-a, 1/x ) )
 *                                --                           )
 *                               |  (a)                        )
 */

static double hy1f1a( double a, double b, double x, double *err )
{
double h1, h2, t, u, temp, acanc, asum, err1, err2;

if ( x == 0 ) {
	acanc = 1.0;
	asum = MTH_MAXNUM;
	goto adone;
	}
temp = WRAP_log( WRAP_fabs(x) );
t = x + temp * (a-b);
u = -temp * a;

if ( b > 0 ) {
	temp = WRAP_lgamma(b);
	t += temp;
	u += temp;
	}

h1 = hyp2f0( a, a-b+1, -1.0/x, 1, &err1 );

temp = WRAP_exp(u) / WRAP_gamma(b-a);
h1 *= temp;
err1 *= temp;

h2 = hyp2f0( b-a, 1.0-a, 1.0/x, 2, &err2 );

if (a < 0)
	temp = WRAP_exp(t) / WRAP_gamma(a);
else
	temp = WRAP_exp( t-WRAP_lgamma(a) );

h2 *= temp;
err2 *= temp;

if ( x < 0.0 )
	asum = h1;
else
	asum = h2;

acanc = WRAP_fabs(err1) + WRAP_fabs(err2);

if ( b < 0 ) {
	temp = WRAP_gamma(b);
	asum *= temp;
	acanc *= WRAP_fabs(temp);
	}


if ( asum != 0.0 )
	acanc /= WRAP_fabs(asum);

	/* fudge factor, since error of asymptotic formula
	 * often seems this much larger than advertised.
	*/
acanc *= 30.0;

adone:

*err = acanc;
return asum ;
}

/*							hyp2f0()	*/

static double hyp2f0( double a, double b, double x, int type, double *err )
{
double a0, alast, t, tlast, maxt;
double n, an, bn, u, sum, temp;

an = a;
bn = b;
a0 = 1.0e0;
alast = 1.0e0;
sum = 0.0;
n = 1.0e0;

tlast = 1.0e9;
maxt = 0.0;
for (t = 1.0; t > MTH_MACHEP; ) {
	if ( an == 0 )
		goto pdone;
	if ( bn == 0 )
		goto pdone;

	u = an * (bn * x / n);

	/* check for blowup */
	temp = WRAP_fabs(u);
	if ( (temp > 1.0 ) && (maxt > (MTH_MAXNUM/temp)) )
		goto error;

	a0 *= u;
	t = WRAP_fabs(a0);

	/* terminating condition for asymptotic series */
	if ( t > tlast )
		goto ndone;

	tlast = t;
	sum += alast;	/* the sum is one term behind */
	alast = a0;

	if ( n > 200 )
		goto ndone;

	an += 1.0e0;
	bn += 1.0e0;
	n += 1.0e0;
	if ( t > maxt )
		maxt = t;
	} 

pdone:	/* series converged! */

/* estimate error due to roundoff and cancellation */
*err = WRAP_fabs( MTH_MACHEP * (n+maxt) );

alast = a0;
goto done;

ndone:	/* series did not converge */

/* The following "Converging factors" are supposed to improve accuracy,
 * but do not actually seem to accomplish very much. */

n -= 1.0;
x = 1.0/x;

switch( type )	/* "type" given as subroutine argument */
{
case 1:
	alast *= ( 0.5 + (0.125 + 0.25*b - 0.5*a + 0.25*x - 0.25*n)/x );
	break;

case 2:
	alast *= 2.0/3.0 - b + 2.0*a + x - n;
	break;

default:
	;
}

	/* estimate error due to roundoff, cancellation, and nonconvergence */
*err = MTH_MACHEP * (n + maxt)  +  WRAP_fabs ( a0 );


done:
sum += alast;
return sum ;

	/* series blew up: */
error:
*err = MTH_MAXNUM;
mtherr( "mth_hyperg", MTHE_TLOSS );
return sum ;
}
