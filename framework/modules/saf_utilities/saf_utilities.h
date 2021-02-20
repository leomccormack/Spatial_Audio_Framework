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
 *@addtogroup Utilities
 *@{
 * @file saf_utilities.h
 * @brief Main header for the utilities module (#SAF_UTILITIES_MODULE)
 *
 * A collection of useful utility functions; including: cross-platform complex
 * number wrappers; optimised linear algebra routines based on CBLAS and LAPACK;
 * FFT wrappers and STFT implementation; IIR/FIR filter coefficients and filter
 * bank designs; lists of common loudspeaker and microphone array coordinates;
 * multi-channel and matrix convolvers; spherical Bessel/Hankel functions
 * (including their derivatives); etc.
 *
 * @author Leo McCormack
 * @date 11.07.2016
 * @license ISC
 */

#ifndef __SAF_UTILITIES_H_INCLUDED__
#define __SAF_UTILITIES_H_INCLUDED__

/* The usual suspects: */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include <assert.h>
#include <limits.h>

/* ========================================================================== */
/*                        Macros and Global Constants                         */
/* ========================================================================== */

/** 2 (true for most humans) */
#define NUM_EARS 2

/** Returns the minimum of the two values */
#define SAF_MIN(a,b) (( (a) < (b) ) ? (a) : (b))

/** Returns the maximum of the two values */
#define SAF_MAX(a,b) (( (a) > (b) ) ? (a) : (b))

/** Ensures value "a" is clamped between the "min" and "max" values */
#define SAF_CLAMP(a,min,max) (SAF_MAX(min, SAF_MIN(max, a)))

/** Boolean true */
#define SAF_TRUE ( 1 )

/** Boolean false */
#define SAF_FALSE ( 0 )

/** pi constant (single precision) */
# define SAF_PI ( 3.14159265358979323846264338327950288f )

/** pi constant (double precision) */
# define SAF_PId ( 3.14159265358979323846264338327950288 )

/** Returns 0 if "x" is not a power of 2 */
#define SAF_ISPOW2(x) (((x & ~(x-1))==x) ? x : 0)

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

/** 4pi (single precision) */
#define FOURPI ( 12.566370614359172f )

/** Converts elevation to inclincation, (in radians) */
#define ELEV2INCL(E) ( (SAF_PI/2.0f - E) )

#ifndef DEG2RAD
/** Converts degrees to radians */
# define DEG2RAD(x) (x * SAF_PI / 180.0f)
#endif

#ifndef RAD2DEG
/** Converts radians to degrees  */
# define RAD2DEG(x) (x * 180.0f / SAF_PI)
#endif

#ifndef MKSTRING_
/** Used to make strings inside of Macros */
# define MKSTRING_(s) #s
#endif

#ifndef MKSTRING
/** Used to make strings inside of Macros */
# define MKSTRING(s) MKSTRING_(s)
#endif

/** Indicates that a particular variable is unused (& squelches any warnings) */
#define SAF_UNUSED(x) (void)(x)

#ifndef NDEBUG /* If debug mode: */
/** Macro to print a warning message along with the filename and line number */
# define saf_print_warning(message) {fprintf(stdout, \
                                    "SAF WARNING: %s [%s LINE %u] \n", message,\
                                    __FILE__, __LINE__);}

/** Macro to print a error message along with the filename and line number */
# define saf_print_error(message) {fprintf(stderr, \
                                  "SAF ERROR: %s [%s LINE %u] \n", message, \
                                  __FILE__, __LINE__); \
                                  exit(EXIT_FAILURE);}

/** Macro to make an assertion, along with a string explaining its purpose */
# define saf_assert(x, message) if (!(x)) \
                    {printf("SAF ASSERTION FAILED: (%s), %s [%s LINE %u].\n", \
                     MKSTRING(x), message, __FILE__, __LINE__); \
                     exit(EXIT_FAILURE); }

#else /* ...otherwise macros do nothing, or just the standard behaviour: */
# define saf_print_warning(message)
# define saf_print_error(message) exit(EXIT_FAILURE)
# define saf_assert(x, message) assert(x)
#endif


/* ========================================================================== */
/*                 External Resources and SAF Utility Headers                 */
/* ========================================================================== */

/* For allocating contiguous multi-dimensional arrays.
 * The original source code can be found here (MIT license):
 *   https://github.com/leomccormack/md_malloc
 */
#include "../../resources/md_malloc/md_malloc.h"

/* The default FFT implementation, for when no optimised implementation is
 * available.
 * The original source code can be found here (BSD-3-Clause license):
 *   https://github.com/mborgerding/kissfft
 */
#include "../../resources/kissFFT/kiss_fftr.h"
#include "../../resources/kissFFT/kiss_fft.h"

/* For computing N-dimensional convex hulls and Delaunay meshes.
 * The original source code can be found here (MIT license):
 *   https://github.com/leomccormack/convhull_3d
 */
#include "../../resources/convhull_3d/convhull_3d.h"

/* For resampling audio and FIR filters.
 * The original source code can be found here (BSD-3-Clause license):
 *   https://github.com/xiph/speexdsp/
 */
#include "../../resources/speex_resampler/speex_resampler.h"

/* For data compression/decompression (Lempel–Ziv 1977 and Huffman coding based)
 * The original source code can be found here (BSD-3-Clause license):
 *   https://github.com/madler/zlib
 */
#include "../../resources/zlib/zlib.h"

/* For cross-platform complex number support */
#include "saf_utility_complex.h"

/* For sorting vectors */
#include "saf_utility_sort.h"

/* Filter coefficients and filterbanks (IIR/FIR) */
#include "saf_utility_filters.h"

/* Many handy linear algebra functions based on CBLAS/LAPACK. Additionally, some
 * optimised proprietary Intel MKL and Apple Accelerate routines are employed
 * (for e.g. vector-vector products, addition etc.) if available; otherwise
 * reverting to default implementations. */
#include "saf_utility_veclib.h"

/* For computing spherical/cylindrical Bessel and Hankel functions and their
 * derivatives */
#include "saf_utility_bessel.h"

/* Wrappers for different FFT implementations */
#include "saf_utility_fft.h"

/* Matrix and multi-channel convolvers */
#include "saf_utility_matrixConv.h"

/* Pitch shifting algorithms */
#include "saf_utility_pitch.h"

/* A collection of signal decorrelators */
#include "saf_utility_decor.h"

/* For misc. functions */
#include "saf_utility_misc.h"

/* For computational geometry functions */
#include "saf_utility_geometry.h"

/* For an implementation of the hybrid complex quadrature mirror filterbank */
#include "saf_utility_qmf.h"

/* Various presets for loudspeaker arrays and uniform distributions of points on
 * spheres. */
#include "saf_utility_loudspeaker_presets.h"

/* Various presets for microphone and hydrophone arrays. */
#include "saf_utility_sensorarray_presets.h"

/* A modified version of the alias-free STFT implementation by Juha Vilkamo.
 * The original source code can be found here (MIT license):
 *   https://github.com/jvilkamo/afSTFT
 * The design is also detailed in chapter 1 of [1]
 *
 * @see [1] Pulkki, V., Delikaris-Manias, S. and Politis, A. 2018. Parametric
 *          time--frequency domain spatial audio. John Wiley & Sons,
 *          Incorporated.
 */
#include "../../resources/afSTFT/afSTFTlib.h"

/* Distance variation function filter coefficient functions. */
#include "saf_utility_dvf.h"


#endif /* __SAF_UTILITIES_H_INCLUDED__ */

/**@} */ /* doxygen addtogroup Utilities */
