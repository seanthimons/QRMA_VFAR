/*							const.c
 *
 *	Globally declared constants
 *
 *
 *
 * SYNOPSIS:
 *
 * extern double nameofconstant;
 *
 *
 *
 *
 * DESCRIPTION:
 *
 * This file contains a number of mathematical constants and
 * also some needed size parameters of the computer arithmetic.
 * The values are supplied as arrays of hexadecimal integers
 * for IEEE arithmetic; arrays of octal constants for ARCH_DEC_VAX
 * arithmetic; and in a normal decimal scientific notation for
 * other machines.  The particular notation used is determined
 * by a symbol (ARCH_DEC_VAX, ARCH_INTEL_X86, or ARCH_UNKNOWN) defined in the include file
 * mconf.h.
 *
 * The default size parameters are as follows.
 *
 * For ARCH_DEC_VAX and ARCH_UNKNOWN modes:
 * MACHEP =  1.38777878078144567553E-17       2**-56
 * MAXLOG =  8.8029691931113054295988E1       log(2**127)
 * MINLOG = -8.872283911167299960540E1        log(2**-128)
 * MAXNUM =  1.701411834604692317316873e38    2**127
 *
 * For IEEE arithmetic (ARCH_INTEL_X86):
 * MACHEP =  1.11022302462515654042E-16       2**-53
 * MAXLOG =  7.09782712893383996843E2         log(2**1024)
 * MINLOG = -7.08396418532264106224E2         log(2**-1022)
 * MAXNUM =  1.7976931348623158E308           2**1024
 *
 * The global symbols for mathematical constants are
 * MTH_PI     =  3.14159265358979323846           pi
 * MTH_PIO2   =  1.57079632679489661923           pi/2
 * MTH_PIO4   =  7.85398163397448309616E-1        pi/4
 * MTH_SQRT2  =  1.41421356237309504880           sqrt(2)
 * MTH_SQRTH  =  7.07106781186547524401E-1        sqrt(2)/2
 * MTH_LOG2E  =  1.4426950408889634073599         1/log(2)
 * MTH_SQ2OPI =  7.9788456080286535587989E-1      sqrt( 2/pi )
 * MTH_LOGE2  =  6.93147180559945309417E-1        log(2)
 * MTH_LOGSQ2 =  3.46573590279972654709E-1        log(2)/2
 * MTH_THPIO4 =  2.35619449019234492885           3*pi/4
 * MTH_TWOOPI =  6.36619772367581343075535E-1     2/pi
 *
 * These lists are subject to change.
 */

/*							const.c */

/*
Cephes Math Library Release 2.3:  March, 1995
Copyright 1984, 1995 by Stephen L. Moshier
*/

	/* AvK, 20060530
	** we need to set CONST_C here , because
	** const.h lies about the type of the
	** constants.
	** Conceptually, they are of course 'double'
	** ,but to make the initialisers easyer, they
	** are actually defined HERE as 'unsigned short []'
	** and declared THERE as 'extern double'.
	** const.h is guarded against CONST_C
	*/
#define CONST_C 1
#include "mconf.h"

#ifdef ARCH_UNKNOWN
#if 1
double MTH_MACHEP =  1.11022302462515654042E-16;   /* 2**-53 */
#else
double MTH_MACHEP =  1.38777878078144567553E-17;   /* 2**-56 */
#endif

double MTH_UFLOWTHRESH =  2.22507385850720138309E-308; /* 2**-1022 */

#ifdef WANT_DENORMAL
double MTH_MAXLOG =  7.09782712893383996732E2;     /* log(MTH_MAXNUM) */
/* double MTH_MINLOG = -7.44440071921381262314E2; */     /* log(2**-1074) */
double MTH_MINLOG = -7.451332191019412076235E2;     /* log(2**-1075) */
#else
double MTH_MAXLOG =  7.08396418532264106224E2;     /* log 2**1022 */
double MTH_MINLOG = -7.08396418532264106224E2;     /* log 2**-1022 */
#endif

double MTH_MAXNUM =  1.79769313486231570815E308;    /* 2**1024*(1-MACHEP) */
double MTH_PI     =  3.14159265358979323846;       /* pi */
double MTH_PIO2   =  1.57079632679489661923;       /* pi/2 */
double MTH_PIO4   =  7.85398163397448309616E-1;    /* pi/4 */
double MTH_SQRT2  =  1.41421356237309504880;       /* sqrt(2) */
double MTH_SQRTH  =  7.07106781186547524401E-1;    /* sqrt(2)/2 */
double MTH_LOG2E  =  1.4426950408889634073599;     /* 1/log(2) */
double MTH_SQ2OPI =  7.9788456080286535587989E-1;  /* sqrt( 2/pi ) */
double MTH_LOGE2  =  6.93147180559945309417E-1;    /* log(2) */
double MTH_LOGSQ2 =  3.46573590279972654709E-1;    /* log(2)/2 */
double MTH_THPIO4 =  2.35619449019234492885;       /* 3*pi/4 */
double MTH_TWOOPI =  6.36619772367581343075535E-1; /* 2/pi */

#ifdef WANT_INFINITIES
double MTH_INFINITY = 99e999;
#else
double MTH_INFINITY =  1.79769313486231570815E308;    /* 2**1024*(1-MTH_MACHEP) */
#endif

#ifdef WANT_NANS
double MTH_NAN = 0.0/0.0;
#else
double MTH_NAN = 0.0;
#endif

#ifdef WANT_MINUSZERO
double MTH_NEGZERO = -0.0;
#else
double MTH_NEGZERO = 0.0;
#endif

#endif /* MACINE_UNKNOWN */

#ifdef ARCH_INTEL_X86
			/* 2**-53 =  1.11022302462515654042E-16 */
unsigned short MTH_MACHEP[4] = {0x0000,0x0000,0x0000,0x3ca0};
unsigned short MTH_UFLOWTHRESH[4] = {0x0000,0x0000,0x0000,0x0010};

#ifdef WANT_DENORMAL
			/* log(MTH_MAXNUM) =  7.09782712893383996732224E2 */
unsigned short MTH_MAXLOG[4] = {0x39ef,0xfefa,0x2e42,0x4086};
			/* log(2**-1074) = - -7.44440071921381262314E2 */
/*unsigned short MTH_MINLOG[4] = {0x71c3,0x446d,0x4385,0xc087};*/
unsigned short MTH_MINLOG[4] = {0x3052,0xd52d,0x4910,0xc087};
#else
			/* log(2**1022) =   7.08396418532264106224E2 */
unsigned short MTH_MAXLOG[4] = {0xbcd2,0xdd7a,0x232b,0x4086};
			/* log(2**-1022) = - 7.08396418532264106224E2 */
unsigned short MTH_MINLOG[4] = {0xbcd2,0xdd7a,0x232b,0xc086};
#endif
			/* 2**1024*(1-MACHEP) =  1.7976931348623158E308 */
unsigned short MTH_MAXNUM[4] = {0xffff,0xffff,0xffff,0x7fef};
unsigned short MTH_PI[4]     = {0x2d18,0x5444,0x21fb,0x4009};
unsigned short MTH_PIO2[4]   = {0x2d18,0x5444,0x21fb,0x3ff9};
unsigned short MTH_PIO4[4]   = {0x2d18,0x5444,0x21fb,0x3fe9};
unsigned short MTH_SQRT2[4]  = {0x3bcd,0x667f,0xa09e,0x3ff6};
unsigned short MTH_SQRTH[4]  = {0x3bcd,0x667f,0xa09e,0x3fe6};
unsigned short MTH_LOG2E[4]  = {0x82fe,0x652b,0x1547,0x3ff7};
unsigned short MTH_SQ2OPI[4] = {0x3651,0x33d4,0x8845,0x3fe9};
unsigned short MTH_LOGE2[4]  = {0x39ef,0xfefa,0x2e42,0x3fe6};
unsigned short MTH_LOGSQ2[4] = {0x39ef,0xfefa,0x2e42,0x3fd6};
unsigned short MTH_THPIO4[4] = {0x21d2,0x7f33,0xd97c,0x4002};
unsigned short MTH_TWOOPI[4] = {0xc883,0x6dc9,0x5f30,0x3fe4};

#ifdef WANT_INFINITIES
unsigned short MTH_INFINITY[4] = {0x0000,0x0000,0x0000,0x7ff0};
#else
unsigned short MTH_INFINITY[4] = {0xffff,0xffff,0xffff,0x7fef};
#endif

#ifdef WANT_NANS
unsigned short MTH_NAN[4] = {0x0000,0x0000,0x0000,0x7ffc};
#else
unsigned short MTH_NAN[4] = {0x0000,0x0000,0x0000,0x0000};
#endif

#ifdef WANT_MINUSZERO
unsigned short MTH_NEGZERO[4] = {0x0000,0x0000,0x0000,0x8000};
#else
unsigned short MTH_NEGZERO[4] = {0x0000,0x0000,0x0000,0x0000};
#endif

#endif

#ifdef ARCH_MIEEE
		/* 2**-53 =  1.11022302462515654042E-16 */
unsigned short MTH_MACHEP[4] = {0x3ca0,0x0000,0x0000,0x0000};
unsigned short MTH_UFLOWTHRESH[4] = {0x0010,0x0000,0x0000,0x0000};
#ifdef WANT_DENORMAL
		/* log(2**1024) =   7.09782712893383996843E2 */
unsigned short MTH_MAXLOG[4] = {0x4086,0x2e42,0xfefa,0x39ef};
		/* log(2**-1074) = - -7.44440071921381262314E2 */
/* unsigned short MTH_MINLOG[4] = {0xc087,0x4385,0x446d,0x71c3}; */
unsigned short MTH_MINLOG[4] = {0xc087,0x4910,0xd52d,0x3052};
#else
		/* log(2**1022) =  7.08396418532264106224E2 */
unsigned short MTH_MAXLOG[4] = {0x4086,0x232b,0xdd7a,0xbcd2};
		/* log(2**-1022) = - 7.08396418532264106224E2 */
unsigned short MTH_MINLOG[4] = {0xc086,0x232b,0xdd7a,0xbcd2};
#endif

		/* 2**1024*(1-MACHEP) =  1.7976931348623158E308 */
unsigned short MTH_MAXNUM[4] = {0x7fef,0xffff,0xffff,0xffff};
unsigned short MTH_PI[4]     = {0x4009,0x21fb,0x5444,0x2d18};
unsigned short MTH_PIO2[4]   = {0x3ff9,0x21fb,0x5444,0x2d18};
unsigned short MTH_PIO4[4]   = {0x3fe9,0x21fb,0x5444,0x2d18};
unsigned short MTH_SQRT2[4]  = {0x3ff6,0xa09e,0x667f,0x3bcd};
unsigned short MTH_SQRTH[4]  = {0x3fe6,0xa09e,0x667f,0x3bcd};
unsigned short MTH_LOG2E[4]  = {0x3ff7,0x1547,0x652b,0x82fe};
unsigned short MTH_SQ2OPI[4] = {0x3fe9,0x8845,0x33d4,0x3651};
unsigned short MTH_LOGE2[4]  = {0x3fe6,0x2e42,0xfefa,0x39ef};
unsigned short MTH_LOGSQ2[4] = {0x3fd6,0x2e42,0xfefa,0x39ef};
unsigned short MTH_THPIO4[4] = {0x4002,0xd97c,0x7f33,0x21d2};
unsigned short MTH_TWOOPI[4] = {0x3fe4,0x5f30,0x6dc9,0xc883};

#ifdef WANT_INFINITIES
unsigned short MTH_INFINITY[4] = {0x7ff0,0x0000,0x0000,0x0000};
#else
unsigned short MTH_INFINITY[4] = {0x7fef,0xffff,0xffff,0xffff};
#endif

#ifdef WANT_NANS
unsigned short MTH_NAN[4] = {0x7ff8,0x0000,0x0000,0x0000};
#else
unsigned short MTH_NAN[4] = {0x0000,0x0000,0x0000,0x0000};
#endif

#ifdef WANT_MINUSZERO
unsigned short MTH_NEGZERO[4] = {0x8000,0x0000,0x0000,0x0000};
#else
unsigned short MTH_NEGZERO[4] = {0x0000,0x0000,0x0000,0x0000};
#endif

#endif

#ifdef ARCH_DEC_VAX
		/* 2**-56 =  1.38777878078144567553E-17 */
unsigned short MTH_MACHEP[4] = {0022200,0000000,0000000,0000000};
unsigned short MTH_UFLOWTHRESH[4] = {0x0080,0x0000,0x0000,0x0000};
		/* log 2**127 = 88.029691931113054295988 */
unsigned short MTH_MAXLOG[4] = {041660,007463,0143742,025733,};
		/* log 2**-128 = -88.72283911167299960540 */
unsigned short MTH_MINLOG[4] = {0141661,071027,0173721,0147572,};
		/* 2**127 = 1.701411834604692317316873e38 */
unsigned short MTH_MAXNUM[4] = {077777,0177777,0177777,0177777,};
unsigned short MTH_PI[4]     = {040511,007732,0121041,064302,};
unsigned short MTH_PIO2[4]   = {040311,007732,0121041,064302,};
unsigned short MTH_PIO4[4]   = {040111,007732,0121041,064302,};
unsigned short MTH_SQRT2[4]  = {040265,002363,031771,0157145,};
unsigned short MTH_SQRTH[4]  = {040065,002363,031771,0157144,};
unsigned short MTH_LOG2E[4]  = {040270,0125073,024534,013761,};
unsigned short MTH_SQ2OPI[4] = {040114,041051,0117241,0131204,};
unsigned short MTH_LOGE2[4]  = {040061,071027,0173721,0147572,};
unsigned short MTH_LOGSQ2[4] = {037661,071027,0173721,0147572,};
unsigned short MTH_THPIO4[4] = {040426,0145743,0174631,007222,};
unsigned short MTH_TWOOPI[4] = {040042,0174603,067116,042025,};

/* Approximate infinity by MTH_MAXNUM.  */
unsigned short MTH_INFINITY[4] = {077777,0177777,0177777,0177777,};
unsigned short MTH_NAN[4] = {0000000,0000000,0000000,0000000};

#ifdef WANT_MINUSZERO
unsigned short MTH_NEGZERO[4] = {0000000,0000000,0000000,0100000};
#else
unsigned short MTH_NEGZERO[4] = {0000000,0000000,0000000,0000000};
#endif

#endif

#ifndef ARCH_UNKNOWN
extern unsigned short MTH_MACHEP[];
extern unsigned short MTH_UFLOWTHRESH[];
extern unsigned short MTH_MAXLOG[];
extern unsigned short MTH_UNDLOG[];
extern unsigned short MTH_MINLOG[];
extern unsigned short MTH_MAXNUM[];
extern unsigned short MTH_PI[];
extern unsigned short MTH_PIO2[];
extern unsigned short MTH_PIO4[];
extern unsigned short MTH_SQRT2[];
extern unsigned short MTH_SQRTH[];
extern unsigned short MTH_LOG2E[];
extern unsigned short MTH_SQ2OPI[];
extern unsigned short MTH_LOGE2[];
extern unsigned short MTH_LOGSQ2[];
extern unsigned short MTH_THPIO4[];
extern unsigned short MTH_TWOOPI[];
extern unsigned short MTH_INFINITY[];
extern unsigned short MTH_NAN[];
extern unsigned short MTH_NEGZERO[];
#endif
