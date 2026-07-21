/*
 *   Created: Tue May 30 2006, AvK
 *  This is intended as a wrapper between mathlab
 *  , user-functions, and (more or less) std math library.
 * The WRAP_xxx Macros provided here function as toggles:
 * they either point to libm-functions
 * , or to the (locally defined) mth_xxx() versions.
 *
 * Note : for historical reasons the gamma() function in UNIX's
 * libm is called tgamma(): True gamma.
 */
#ifndef MTH_WRAP_H
#define MTH_WRAP_H 1

#include "protos.h"

#define WRAP_hyperg mth_hyperg

#define WRAP_gamma tgamma
#define WRAP_lgamma lgamma

#define WRAP_frexp frexp
#define WRAP_ldexp ldexp

#define WRAP_round round
#define WRAP_fabs fabs
#define WRAP_floor floor
#define WRAP_ceil ceil
#define WRAP_modf modf

#define WRAP_isnan isnan
#define WRAP_isfinite isfinite
#define WRAP_signbit signbit


#define WRAP_pow mth_pow
#define WRAP_powi mth_powi

	/* we don't have local versions for these */
#define WRAP_sin sin
#define WRAP_log log
#define WRAP_exp exp

#endif
