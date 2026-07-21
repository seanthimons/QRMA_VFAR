/*							mth_ceil()
 *							mth_floor()
 *							mth_frexp()
 *							mth_ldexp()
 *							mth_signbit()
 *							mth_isnan()
 *							mth_isfinite()
 *
 *	Floating point numeric utilities
 *
 *
 *
 * SYNOPSIS:
 *
 * double mth_ceil(), mth_floor(), mth_frexp(), mth_ldexp();
 * int mth_signbit(), mth_isnan(), mth_isfinite();
 * double x, y;
 * int expnt, n;
 *
 * y = mth_floor(x);
 * y = mth_ceil(x);
 * y = mth_frexp( x, &expnt );
 * y = mth_ldexp( x, n );
 * n = mth_signbit(x);
 * n = mth_isnan(x);
 * n = mth_isfinite(x);
 *
 * y = mth_fabs(x);
 *
 *
 *
 * DESCRIPTION:
 *
 * All four routines return a double precision floating point
 * result.
 *
 * mth_floor() returns the largest integer less than or equal to x.
 * It truncates toward minus infinity.
 *
 * mth_ceil() returns the smallest integer greater than or equal
 * to x.  It truncates toward plus infinity.
 *
 * mth_frexp() extracts the exponent from x.  It returns an integer
 * power of two to expnt and the significand between 0.5 and 1
 * to y.  Thus  x = y * 2**expn.
 *
 * mth_ldexp() multiplies x by 2**n.
 *
 * mth_signbit(x) returns 1 if the sign bit of x is 1, else 0.
 *
 * These functions are part of the standard C run time library
 * for many but not all C compilers.  The ones supplied are
 * written in C for either DEC or IEEE arithmetic.  They should
 * be used only if your compiler library does not already have
 * them.
 *
 * The IEEE versions assume that denormal numbers are implemented
 * in the arithmetic.  Some modifications will be required if
 * the arithmetic has abrupt rather than gradual underflow.
 */

/*
Cephes Math Library Release 2.3:  March, 1995
Copyright 1984, 1995 by Stephen L. Moshier
*/


#include "mconf.h"
#include "protos.h"
#include "const.h"

#ifdef ARCH_UNKNOWN
/* ceil(), floor(), frexp(), ldexp() may need to be rewritten. */
#undef ARCH_UNKNOWN
#if BIGENDIAN
#define ARCH_MIEEE 1
#else
#define ARCH_INTEL_X86 1
#endif
#endif

#ifdef DEC
#define EXPMSK 0x807f
#define MEXP 255
#define NBITS 56
#endif

#ifdef ARCH_INTEL_X86
#define EXPMSK 0x800f
#define MEXP 0x7ff
#define NBITS 53
#endif

#ifdef ARCH_MIEEE
#define EXPMSK 0x800f
#define MEXP 0x7ff
#define NBITS 53
#endif


double mth_ceil(double x)
{
double y;

#ifdef ARCH_UNKNOWN
mtherr( "mth_ceil", MTHE_DOMAIN );
return 0.0;
#endif

y = WRAP_floor(x);
if (y < x) y += 1.0;

#ifdef WANT_MINUSZERO
if (y == 0.0 && x < 0.0)
	return MTH_NEGZERO;
#endif
return y;
}

/* Bit clearing masks: */

static unsigned short bmask[] = {
0xffff,
0xfffe,
0xfffc,
0xfff8,
0xfff0,
0xffe0,
0xffc0,
0xff80,
0xff00,
0xfe00,
0xfc00,
0xf800,
0xf000,
0xe000,
0xc000,
0x8000,
0x0000,
};


double mth_floor(double x)
{
union
	{
	double y;
	unsigned short sh[4];
	} u;
unsigned short *p;
int e;

#ifdef ARCH_UNKNOWN
mtherr( "mth_floor", MTHE_DOMAIN );
return 0.0;
#endif

u.y = x;
/* find the exponent (power of 2) */
#ifdef ARCH_DEC_VAX
p = (unsigned short *)&u.sh[0];
e = ((*p  >> 7) & 0377) - 0201;
p += 3;
#endif

#ifdef ARCH_INTEL_X86
p = (unsigned short *)&u.sh[3];
e = ((*p >> 4) & 0x7ff) - 0x3ff;
p -= 3;
#endif

#ifdef ARCH_MIEEE
p = (unsigned short *)&u.sh[0];
e = ((*p >> 4) & 0x7ff) - 0x3ff;
p += 3;
#endif

if (e < 0) {
	if (u.y < 0.0)
		return -1.0;
	else
		return 0.0;
	}

e = (NBITS -1) - e;
/* clean out 16 bits at a time */
while (e >= 16) {
#ifdef ARCH_INTEL_X86
	*p++ = 0;
#endif

#ifdef ARCH_DEC_VAX
	*p-- = 0;
#endif

#ifdef ARCH_MIEEE
	*p-- = 0;
#endif
	e -= 16;
	}

/* clear the remaining bits */
if (e > 0)
	*p &= bmask[e];

if ((x < 0) && (u.y != x))
	u.y -= 1.0;

return u.y;
}



double mth_frexp( double x, int *pw2)
{
union
	{
	double y;
	unsigned short sh[4];
	} u;
int i;
#ifdef WANT_DENORMAL
int k;
#endif
short *q;

u.y = x;

#ifdef ARCH_UNKNOWN
mtherr( "mth_frexp", MTHE_DOMAIN );
return 0.0;
#endif

#ifdef ARCH_INTEL_X86
q = (short *)&u.sh[3];
#endif

#ifdef ARCH_DEC_VAX
q = (short *)&u.sh[0];
#endif

#ifdef ARCH_MIEEE
q = (short *)&u.sh[0];
#endif

/* find the exponent (power of 2) */
#ifdef ARCH_DEC_VAX
i = ( *q >> 7) & 0377;
if (i == 0) {
	*pw2 = 0;
	return 0.0;
	}
i -= 0200;
*pw2 = i;
*q &= 0x807f;	/* strip all exponent bits */
*q |= 040000;	/* mantissa between 0.5 and 1 */
return u.y;
#endif

#ifdef ARCH_INTEL_X86
i  = ( *q >> 4) & 0x7ff;
if (i != 0)
	goto ieeedon;
#endif

#ifdef ARCH_MIEEE
i  =  *q >> 4;
i &= 0x7ff;
if (i != 0)
	goto ieeedon;
#ifdef WANT_DENORMAL

#else
*pw2 = 0;
return 0.0;
#endif

#endif


#ifndef ARCH_DEC_VAX
/* Number is denormal or zero */
#ifdef WANT_DENORMAL
if (u.y == 0.0) {
	*pw2 = 0;
	return 0.0;
	}


/* Handle denormal number. */
do	{
	u.y *= 2.0;
	i -= 1;
	k  = ( *q >> 4) & 0x7ff;
	} while (k == 0 );

i = i + k;
#endif /* WANT_DENORMAL */

ieeedon:

i -= 0x3fe;
*pw2 = i;
*q &= 0x800f;
*q |= 0x3fe0;
return u.y;
#endif
}


double mth_ldexp( double x, int pw2)
{
union
	{
	double y;
	unsigned short sh[4];
	} u;
short *q;
int e;

#ifdef ARCH_UNKNOWN
mtherr( "mth_ldexp", MTHE_DOMAIN );
return 0.0;
#endif

u.y = x;
#ifdef ARCH_DEC_VAX
q = (short *)&u.sh[0];
e  = (*q >> 7) & 0377;
if (e == 0) return 0.0;
#else

#ifdef ARCH_INTEL_X86
q = (short *)&u.sh[3];
#endif
#ifdef ARCH_MIEEE
q = (short *)&u.sh[0];
#endif

while ((e = (*q & 0x7ff0) >> 4) == 0 ) {
	if (u.y == 0.0) {
		return 0.0 ;
		}
/* Input is denormal. */
	if (pw2 > 0) {
		u.y *= 2.0;
		pw2 -= 1;
		}
	if (pw2 < 0) {
		if (pw2 < -53) return 0.0;
		u.y /= 2.0;
		pw2 += 1;
		}
	if (pw2 == 0)
		return u.y;
	}
#endif /* not ARCH_DEC_VAX */

e += pw2;

/* Handle overflow */
#ifdef ARCH_DEC_VAX
if (e > MEXP)
	return MTH_MAXNUM;
#else
if (e >= MEXP)
	return 2.0*MTH_MAXNUM;
#endif

/* Handle denormalized results */
if (e < 1) {
#ifdef WANT_DENORMAL
	if (e < -53 )
		return 0.0;
	*q &= 0x800f;
	*q |= 0x10;
	while (e < 1) {
		u.y /= 2.0;
		e += 1;
		}
	return u.y;
#else
	return 0.0;
#endif
	}
else
	{
#ifdef ARCH_DEC_VAX
	*q &= 0x807f;	/* strip all exponent bits */
	*q |= (e & 0xff) << 7;
#else
	*q &= 0x800f;
	*q |= (e & 0x7ff) << 4;
#endif
	return u.y;
	}
}

/* AvK: added this for completeness
** In real life, this could be a macro
*/

double mth_fabs(double x)
{
if (WRAP_signbit(x)) return -x;
else return x;
}

/* Return 1 if the sign bit of x is 1, else 0.  */

int mth_signbit(double x)
{
union
	{
	double d;
	short s[4];
	int i[2];
	} u;

u.d = x;

if (sizeof(int) == 4)
	{
#ifdef ARCH_INTEL_X86
	return u.i[1] < 0 ;
#endif
#ifdef ARCH_DEC_VAX
	return u.s[3] < 0 ;
#endif
#ifdef ARCH_MIEEE
	return u.i[0] < 0 ;
#endif
	}
else
	{
#ifdef ARCH_INTEL_X86
	return u.s[3] < 0 ;
#endif
#ifdef ARCH_DEC_VAX
	return u.s[3] < 0 ;
#endif
#ifdef ARCH_MIEEE
	return u.s[0] < 0 ;
#endif
	}
}


/* Return 1 if x is a number that is Not a Number, else return 0.  */

int mth_isnan( double x)
{
#ifdef WANT_NANS
union
	{
	double d;
	unsigned short s[4];
	unsigned int i[2];
	} u;

u.d = x;

if (sizeof(int) == 4)
	{
#ifdef ARCH_INTEL_X86
	if (((u.i[1] & 0x7ff00000) == 0x7ff00000)
	    && (((u.i[1] & 0x000fffff) != 0) || (u.i[0] != 0)))
		return 1;
#endif
#ifdef ARCH_DEC_VAX
	if ((u.s[1] & 0x7fff) == 0) {
		if ((u.s[2] | u.s[1] | u.s[0]) != 0)
			return 1;
		}
#endif
#ifdef ARCH_MIEEE
	if (((u.i[0] & 0x7ff00000) == 0x7ff00000)
	    && (((u.i[0] & 0x000fffff) != 0) || (u.i[1] != 0)))
		return 1;
#endif
	return 0;
	}
else
	{ /* size int not 4 */
#ifdef ARCH_INTEL_X86
	if ((u.s[3] & 0x7ff0) == 0x7ff0) {
		if (((u.s[3] & 0x000f) | u.s[2] | u.s[1] | u.s[0]) != 0)
			return 1;
		}
#endif
#ifdef ARCH_DEC_VAX
	if ((u.s[3] & 0x7fff) == 0) {
		if ((u.s[2] | u.s[1] | u.s[0]) != 0)
			return 1;
		}
#endif
#ifdef ARCH_MIEEE
	if ((u.s[0] & 0x7ff0) == 0x7ff0) {
		if (((u.s[0] & 0x000f) | u.s[1] | u.s[2] | u.s[3]) != 0)
			return 1;
		}
#endif
	return 0;
	} /* size int not 4 */

#else
/* No WANT_NANS.  */
return 0;
#endif
}


/* Return 1 if x is not infinite and is not a NaN.  */

int mth_isfinite(double x)
{
#ifdef WANT_INFINITIES
union
	{
	double d;
	unsigned short s[4];
	unsigned int i[2];
	} u;

u.d = x;

if (sizeof(int) == 4) {
#ifdef ARCH_INTEL_X86
	if ((u.i[1] & 0x7ff00000) != 0x7ff00000) return 1;
#endif

#ifdef ARCH_DEC_VAX
	if ((u.s[3] & 0x7fff) != 0) return 1;
#endif

#ifdef ARCH_MIEEE
	if ((u.i[0] & 0x7ff00000) != 0x7ff00000) return 1;
#endif
	return 0;
	}
else
	{
#ifdef ARCH_INTEL_X86
	if ((u.s[3] & 0x7ff0) != 0x7ff0) return 1;
#endif

#ifdef ARCH_DEC_VAX
	if ((u.s[3] & 0x7fff) != 0) return 1;
#endif

#ifdef ARCH_MIEEE
	if ((u.s[0] & 0x7ff0) != 0x7ff0) return 1;
#endif
	return 0;
	}
#else /* No INFINITY.  */
return 1;
#endif
}
