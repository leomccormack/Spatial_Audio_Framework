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

/*
 * Filename: saf_utilities.h (include header)
 * ------------------------------------------
 * Contains a collection of useful memory allocation functions and cross-
 * platform complex number wrappers. Optimised linear algebra routines utilising
 * BLAS and LAPACK are also included.
 *
 * Dependencies:
 *     A performance library comprising CBLAS and LAPACK routines is required
 *     by the module and, thus, also by the SAF framework as a whole.
 *     Add one of the following FLAGS to your project's preprocessor definitions
 *     list, in order to enable one of these suitable performance libraries,
 *     which must also be linked correctly to your project.
 *         SAF_USE_INTEL_MKL
 *             to enable Intel's Math Kernal Library
 *         SAF_USE_ATLAS
 *             to enable ATLAS BLAS routines with CLAPACK
 *         SAF_USE_OPEN_BLAS_AND_LAPACKE
 *             to enable OpenBLAS and use the LAPACK(E) interface
 *     More information can be found here:
 *         https://github.com/leomccormack/Spatial_Audio_Framework
 *     Mac users only:
 *         saf_utilities will employ Apple's Accelerate library by default, if
 *         none of the above FLAGS are defined.
 * Author, date created:
 *     Leo McCormack, 11.07.2016
 */

#ifndef __SAF_UTILITIES_H_INCLUDED__
#define __SAF_UTILITIES_H_INCLUDED__

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */
    
/* ========================================================================== */
/*                        Performance Library to Employ                       */
/* ========================================================================== */

#if defined(SAF_USE_INTEL_MKL)
# define VECLIB_USE_LAPACK_FORTRAN_INTERFACE  
# include "mkl.h"
#elif defined(SAF_USE_ATLAS)
# define VECLIB_USE_CLAPACK_INTERFACE
# include "cblas-atlas.h"
/* note: ATLAS's CLAPACK does not include some LAPACK functions */
# include "clapack.h"
#elif defined(SAF_USE_OPEN_BLAS_AND_LAPACKE)
# define VECLIB_USE_LAPACKE_INTERFACE
# include "cblas.h"
# include "lapacke.h"
#elif defined(__APPLE__)
# define VECLIB_USE_LAPACK_FORTRAN_INTERFACE
# include "Accelerate/Accelerate.h"
#else
# error "SAF requires a library (or libraries) which supports CBLAS and LAPACK"
#endif
    
    
/* ========================================================================== */
/*                        Macros and Global Constants                         */
/* ========================================================================== */
    
#ifndef MIN
# define MIN(a,b) (( (a) < (b) ) ? (a) : (b))
#endif
#ifndef MAX
# define MAX(a,b) (( (a) > (b) ) ? (a) : (b))
#endif
#ifndef CLAMP
# define CLAMP(a,min,max) (MAX(min, MIN(max, a)))
#endif
#ifndef PI
# define PI ( 3.14159265358979323846264338327950288f )
#endif
#ifndef M_PI
# define M_PI ( 3.14159265358979323846264338327950288f )
#endif
#define SAF_ISPOW2(x) (((x & ~(x-1))==x) ? x : 0);
 
    
/* ========================================================================== */
/*                       Resources and Utility Headers                        */
/* ========================================================================== */
    
/* for error message handling */
#include "../saf_utilities/saf_error.h"
/* for allocating multi-dimensional arrays */
#include "../resources/md_malloc/md_malloc.h"
/* default FFT implementation, if no optimised implementation is available */
#include "../resources/kissFFT/kiss_fftr.h"
/* for generating 3-D convex hulls */
#include "../resources/convhull_3d/convhull_3d.h"
/* for cross-platform complex numbers wrapper */
#include "../saf_utilities/saf_complex.h"
/* for sorting vectors */
#include "../saf_utilities/saf_sort.h"
/* filter coefficients (IIR/FIR) */
#include "../saf_utilities/saf_filters.h"
/* Many handy linear algebra functions based on CBLAS/LAPACK/IntelMKL/Accelerate */
#include "../saf_utilities/saf_veclib.h"
/* optimised FFT routines */
#include "../saf_utilities/saf_fft.h"
/* matrix convolver */
#include "../saf_utilities/saf_matrixConv.h"
/* for decorrelators */
#include "../saf_utilities/saf_decor.h"
/* for determining ERBs */
#include "../saf_utilities/saf_erb.h"
/* various presets for loudspeaker arrays and uniform distributions of points on
 * spheres. */
#include "../saf_utilities/saf_loudspeaker_presets.h"
/* various presets for microphone and hydrophone arrays. */
#include "../saf_utilities/saf_sensorarray_presets.h"

    
#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */

#endif /* __SAF_UTILITIES_H_INCLUDED__ */
