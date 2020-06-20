/*
 * Copyright 2016-2018 Leo McCormack
 *
 * Permission to use, copy, modify, and/or distribute this software for any
 * purpose with or without fee is hereby granted, provided that the above
 * copyright notice and this permission notice appear in all copies.
 *
 * THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES WITH
 * REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY
 * AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY SPECIAL, DIRECT,
 * INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM
 * LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR
 * OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR
 * PERFORMANCE OF THIS SOFTWARE.
 */

/**
 * @file saf_utilities.h  
 * @brief Main header for the utilities module (#SAF_UTILITIES_MODULE)
 *
 * Contains a collection of useful memory allocation functions, cross-platform
 * complex number wrappers, and optimised linear algebra routines utilising BLAS
 * and LAPACK.
 *
 * ## Dependencies
 *   A performance library comprising CBLAS and LAPACK routines is required by
 *   the module and, thus, also by the SAF framework as a whole. Add one of the
 *   following FLAGS to your project's preprocessor definitions list, in order
 *   to enable one of these suitable performance libraries, which must also be
 *   linked correctly to your project.
 *   - SAF_USE_INTEL_MKL:
 *       to enable Intel's Math Kernal Library with Fortran LAPACK interface
 *   - SAF_USE_ATLAS:
 *       to enable ATLAS BLAS routines and ATLAS's CLAPACK interface
 *   - SAF_USE_OPENBLAS_WITH_LAPACKE:
 *       to enable OpenBLAS with LAPACKE interface
 *
 * @see More information can be found here:
 *      https://github.com/leomccormack/Spatial_Audio_Framework
 * @note MacOSX users only: saf_utilities will employ Apple's Accelerate library
 *       by default, if none of the above FLAGS are defined.
 *
 * @author Leo McCormack
 * @date 11.07.2016
 */

#ifndef __SAF_UTILITIES_H_INCLUDED__
#define __SAF_UTILITIES_H_INCLUDED__
    
/* ========================================================================== */
/*                        Performance Library to Employ                       */
/* ========================================================================== */

#if defined(SAF_USE_INTEL_MKL)
/*
 * Using Intel's Math Kernal Library (MKL)
 * (Generally the fastest library for x86/amd64 based architectures)
 */
# define VECLIB_USE_LAPACK_FORTRAN_INTERFACE /**< Fortan interface of LAPACK */
# include "mkl.h"

#elif defined(SAF_USE_OPEN_BLAS_AND_LAPACKE)
/*
 * Using OpenBLAS and the LAPACKE interface
 * (A good option for ARM based architectures)
 */
# define VECLIB_USE_LAPACKE_INTERFACE /**< LAPACKE interface */
# include "cblas.h"
# include "lapacke.h"

#elif defined(SAF_USE_ATLAS)
/*
 * Using the Automatically Tuned Linear Algebra Software (ATLAS) library
 * (Not recommended, since some saf_veclib functions do not support ATLAS)
 */
# define VECLIB_USE_CLAPACK_INTERFACE /**< CLAPACK interface */
# include "cblas-atlas.h"
# include "clapack.h" /* NOTE: CLAPACK does not include all LAPACK functions! */

#elif defined(__APPLE__)
/*
 * Using Apple's Accelerate library (vDSP)
 * (Solid choice, but only works under MacOSX)
 */
# define VECLIB_USE_LAPACK_FORTRAN_INTERFACE /**< Fortan interface of LAPACK */
# include "Accelerate/Accelerate.h"
#else
# error "SAF requires a library (or libraries) which supports CBLAS and LAPACK"
#endif


/* ========================================================================== */
/*                        Macros and Global Constants                         */
/* ========================================================================== */

#ifndef NUM_EARS
/** 2 (true for most humans) */
# define NUM_EARS 2
#endif
#ifndef MIN
/** Returns the minimum of the two values */
# define MIN(a,b) (( (a) < (b) ) ? (a) : (b))
#endif
#ifndef MAX
/** Returns the maximum of the two values */
# define MAX(a,b) (( (a) > (b) ) ? (a) : (b))
#endif
#ifndef CLAMP
/** Ensures value "a" is clamped between the "min" and "max" values */
# define CLAMP(a,min,max) (MAX(min, MIN(max, a)))
#endif
#ifndef M_PI
/** pi constant (single precision) */
# define M_PI ( 3.14159265358979323846264338327950288f )
#endif
#ifndef M_PId
/** pi constant (double precision) */
# define M_PId ( 3.14159265358979323846264338327950288 )
#endif
#ifndef SAF_PI
/** pi constant (single precision) */
# define SAF_PI ( 3.14159265358979323846264338327950288f )
#endif
#ifndef SAF_PId
/** pi constant (double precision) */
# define SAF_PId ( 3.14159265358979323846264338327950288 )
#endif
/** Returns 0 if "x" is not a power of 2 */
#define SAF_ISPOW2(x) (((x & ~(x-1))==x) ? x : 0);
#ifndef ISEVEN
/** Returns 1 if "n" is even valued, and 0 if it is not */
# define ISEVEN(n)   ((n%2 == 0) ? 1 : 0)
#endif
#ifndef ISODD
/** Returns 1 if "n" is odd valued, and 0 if it is not */
# define ISODD(n)    ((n%2 != 0) ? 1 : 0)
#endif
/** sqrt(4pi) (single precision) */
#define SQRT4PI ( 3.544907701811032f )
/** Converts elevation to inclincation, (in radians) */
#define ELEV2INCL(E) ( (M_PI/2 - E) )


/* ========================================================================== */
/*                 External Resources and SAF Utility Headers                 */
/* ========================================================================== */

/* (In correct linker order) */
/* for error message handling */
#include "saf_utility_error.h"
/* for allocating multi-dimensional arrays */
#include "../resources/md_malloc/md_malloc.h"
/* default FFT implementation, if no optimised implementation is available */
#include "../resources/kissFFT//kiss_fftr.h"
#include "../resources/kissFFT/kiss_fft.h"
/* for generating 3-D convex hulls */
#include "../resources/convhull_3d/convhull_3d.h"
/* for cross-platform complex numbers wrapper */
#include "saf_utility_complex.h"
/* for sorting vectors */
#include "saf_utility_sort.h"
/* filter coefficients (IIR/FIR) */
#include "saf_utility_filters.h"
/* Many handy linear algebra functions based on CBLAS/LAPACK/IntelMKL/Accelerate
 */
#include "saf_utility_veclib.h"
/* For computing spherical/cylindrical Bessel and Hankel functions */
#include "saf_utility_bessel.h"
/* optimised FFT routines */
#include "saf_utility_fft.h"
/* matrix convolver */
#include "saf_utility_matrixConv.h"
/* pitch shifting algorithms */
#include "saf_utility_pitch.h"
/* for decorrelators */
#include "saf_utility_decor.h"
/* for determining ERBs */
#include "saf_utility_erb.h"
/* for misc. functions */
#include "saf_utility_misc.h"
/* various presets for loudspeaker arrays and uniform distributions of points on
 * spheres. */
#include "saf_utility_loudspeaker_presets.h"
/* various presets for microphone and hydrophone arrays. */
#include "saf_utility_sensorarray_presets.h"
/* modified version of the alias-free STFT implementation by Juha Vilkamo.
 * The original source code can be found here:
 *   https://github.com/jvilkamo/afSTFT
 * and the design is also detailed in chapter 1 of [1]
 *
 * @see [1] Pulkki, V., Delikaris-Manias, S. and Politis, A. 2018. Parametric
 *          time--frequency domain spatial audio. John Wiley & Sons,
 *          Incorporated.
 */
#include "../resources/afSTFT/afSTFTlib.h"


#endif /* __SAF_UTILITIES_H_INCLUDED__ */
