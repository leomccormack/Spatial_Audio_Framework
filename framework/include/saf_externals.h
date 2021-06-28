/*
 * Copyright 2020 Leo McCormack
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
 * @file saf_externals.h
 * @brief Include header for SAF externals
 *
 * @note Including this header is optional and only needed if you wish to have
 *       access to these external libraries in your own project.
 * @warning Using ATLAS (SAF_USE_ATLAS) as the performance library is not
 *          recommended, since some LAPACK routines are not implemented by the
 *          library! However, if you don't mind losing some SAF functionality
 *          (namely: certain veclib utility functions), then it may still be a
 *          good choice for your particular project.
 *
 * ## Required Dependencies
 *   A performance library comprising CBLAS and LAPACK routines is required by
 *   SAF. Add one of the following FLAGS to your project's preprocessor
 *   definitions in order to enable one of these suitable performance libraries,
 *   (which must be correctly linked when building SAF):
 *   - SAF_USE_INTEL_MKL_LP64:
 *       to enable Intel's Math Kernal Library with the Fortran LAPACK interface
 *   - SAF_USE_INTEL_MKL_ILP64
 *       same as SAF_USE_INTEL_MKL except using int64 and LAPACKE interface
 *   - SAF_USE_OPENBLAS_WITH_LAPACKE:
 *       to enable OpenBLAS with the LAPACKE interface
 *   - SAF_USE_APPLE_ACCELERATE:
 *       to enable the Accelerate framework with the Fortran LAPACK interface
 *   - SAF_USE_ATLAS:
 *       to enable ATLAS BLAS routines and ATLAS's CLAPACK interface
 *
 * ## Optional dependencies
 *   netcdf is required by the optional saf_sofa_reader module, which is enabled
 *   with the following pre-processor flag: SAF_ENABLE_SOFA_READER_MODULE
 *
 *   Intel IPP may be optionally used with the flag: SAF_USE_INTEL_IPP
 *
 *   SSE/SSE2/SSE3 intrinsics support may be enabled with: SAF_ENABLE_SIMD
 *   Alternatively, AVX/AVX2 intrinsics may be employed with the additional
 *   flag: SAF_USE_AVX
 *   (Note that intrinsics require a supported x86_64 processor)
 *
 * @see More information can be found in the docs folder regarding dependencies
 *
 * @author Leo McCormack
 * @date 06.08.2020
 */

#ifndef __SAF_EXTERNALS_H_INCLUDED__
#define __SAF_EXTERNALS_H_INCLUDED__

/* Assert that only one CBLAS/LAPACK performance library has been specified */
#if (defined(SAF_USE_INTEL_MKL_LP64) + \
     defined(SAF_USE_INTEL_MKL_ILP64) + \
     defined(SAF_USE_OPEN_BLAS_AND_LAPACKE) + \
     defined(SAF_USE_ATLAS) + \
     defined(SAF_USE_APPLE_ACCELERATE)) > 1
# error Only one performance library flag can be defined!
#endif

/* If SAF_USE_AVX is defined, then so must SAF_ENABLE_SIMD */
#if defined(SAF_USE_AVX) && !defined(SAF_ENABLE_SIMD)
# error SAF_USE_AVX requies SAF_ENABLE_SIMD to also be defined!
#endif

/* ========================================================================== */
/*                        Performance Library to Employ                       */
/* ========================================================================== */

#if defined(SAF_USE_INTEL_MKL_LP64)
/*
 * Using Intel's Math Kernal Library (MKL) LP64 configuration (32-bit int)
 * (Generally the fastest library for x86 based architectures)
 */
/* Note that Intel MKL LP64 supports Fortran LAPACK and LAPACKE interfaces: */
# define SAF_VECLIB_USE_LAPACK_FORTRAN_INTERFACE /**< LAPACK interface */
# include "mkl.h"

#elif defined(SAF_USE_INTEL_MKL_ILP64)
/*
 * Using Intel's Math Kernal Library (MKL) ILP64 configuration (64-bit int)
 * (Generally the fastest library for x86 based architectures. This 64-bit int
 * version is the one employed by e.g. MATLAB. Therefore it is required if you
 * wish to build MEX objects using SAF code; see e.g. extras/safmex)
 */
/* Note that Intel MKL ILP64 will only work with the LAPACKE interface: */
# define SAF_VECLIB_USE_LAPACKE_INTERFACE /**< LAPACK interface */
# define MKL_ILP64
# include "mkl.h"

#elif defined(SAF_USE_OPEN_BLAS_AND_LAPACKE)
/*
 * Using OpenBLAS and the LAPACKE interface
 * (A good option for both x86 and ARM based architectures)
 */
# define SAF_VECLIB_USE_LAPACKE_INTERFACE /**< LAPACK interface */
# include "cblas.h"
# include "lapacke.h"

#elif defined(SAF_USE_ATLAS)
/*
 * Using the Automatically Tuned Linear Algebra Software (ATLAS) library
 * (Not recommended, since some saf_veclib functions do not work with ATLAS)
 */ 
# define SAF_VECLIB_USE_CLAPACK_INTERFACE /**< LAPACK interface */
# include "cblas-atlas.h"
# include "clapack.h"
# warning Note: CLAPACK does not include all LAPACK routines!

#elif defined(__APPLE__) && defined(SAF_USE_APPLE_ACCELERATE)
/*
 * Using Apple's Accelerate library (vDSP)
 * (Solid choice for both x86 and ARM, but only works under MacOSX and is not as
 * fast as Intel MKL for x86 systems)
 */
# define SAF_VECLIB_USE_LAPACK_FORTRAN_INTERFACE /**< LAPACK interface */
# include "Accelerate/Accelerate.h"

#else
/*
 * If you would like to use some other CBLAS/LAPACK supporting library then
 * please get in touch! :-)
 */
# error SAF requires a library (or libraries) which supports CBLAS and LAPACK
#endif


/* ========================================================================== */
/*                        Optional External Libraries                         */
/* ========================================================================== */

#if defined(SAF_USE_INTEL_IPP)
/*
 * The use of Intel's Integrated Performance Primitives (IPP) is optional, but
 * does lead to some minor performance improvements for saf_fft (compared with
 * Intel MKL).
 */
# include "ipp.h"
#endif

#if defined(SAF_ENABLE_SIMD)
/*
 * SAF heavily favours the use of optimised routines provided by e.g. Intel MKL
 * or Apple Accelerate, which internally employ SIMD intrinsics (SSE/AVX etc.).
 * However, in cases where the employed performance library does not offer an
 * implementation for a particular routine, SAF provides fall-back option(s).
 * SIMD accelerated fall-back options may be enabled with: SAF_ENABLE_SIMD
 *
 * By default SSE-based intrinsics are employed, unless the additional flag
 * SAF_USE_AVX is defined, in which case AVX/AVX2 intrinsics are employed.
 */
# if defined(SAF_USE_AVX)
/*
 * Note that AVX/AVX2 requires the '-mavx' and '-mavx2' compiler flags and a
 * supported CPU (x86_64 architecture)
 */
#  include <immintrin.h> /* for AVX and AVX2 */
# else
/*
 * Note that SSE intrinsics require a supported CPU (x86_64 architecture)
 */
#  include <xmmintrin.h> /* for SSE  */
#  include <emmintrin.h> /* for SSE2 */
#  include <pmmintrin.h> /* for SSE3 */
# endif
#endif


/* ========================================================================== */
/*     External Libraries Required by the Optional saf_sofa_reader Module     */
/* ========================================================================== */

#if defined(SAF_ENABLE_SOFA_READER_MODULE)
/*
 * If your compiler stopped at this point, then please add the path for the
 * netcdf include file to your project's include header paths.
 * Instructions for linking the required "netcdf" library may also be found
 * here: docs/SOFA_READER_MODULE_DEPENDENCIES.md
 */
# include <netcdf.h>
#endif


#endif /* __SAF_EXTERNALS_H_INCLUDED__ */
