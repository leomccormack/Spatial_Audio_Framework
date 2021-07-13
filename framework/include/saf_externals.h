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
 *       to enable Intel's Math Kernel Library with the Fortran LAPACK interface
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
 *   FFTW may be optionally used with the flag: SAF_USE_FFTW
 *
 *   SIMD intrinsics support may be enabled with: SAF_ENABLE_SIMD
 *    - SSE/SSE2/SSE3 intrinsics are used by default
 *    - AVX/AVX2 intrinsics are enabled with compiler flag: -mavx2
 *    - AVX-512  intrinsics are enabled with compiler flag: -mavx512f
 *   (Note that intrinsics require a CPU that supports them)
 *
 * @see More information can be found in the docs folder regarding dependencies
 *
 * @author Leo McCormack
 * @date 06.08.2020
 */

#ifndef __SAF_EXTERNALS_H_INCLUDED__
#define __SAF_EXTERNALS_H_INCLUDED__

/* ========================================================================== */
/*                        Performance Library to Employ                       */
/* ========================================================================== */

/* Assert that only one CBLAS/LAPACK performance library has been specified */
#if (defined(SAF_USE_INTEL_MKL_LP64) + \
     defined(SAF_USE_INTEL_MKL_ILP64) + \
     defined(SAF_USE_OPEN_BLAS_AND_LAPACKE) + \
     defined(SAF_USE_ATLAS) + \
     defined(SAF_USE_APPLE_ACCELERATE)) > 1
# error Only one performance library flag can be defined!
#endif

/*
 * Due to the nature of spatial/multi-channel audio signal processing, SAF is
 * required to employ heavy use of linear algebra operations. Therefore, the
 * framework has been written from the ground up to conform to the CBLAS and
 * LAPACK standards, of which there a number of highly optimised performance
 * libraries available:
 */
#if defined(SAF_USE_INTEL_MKL_LP64)
/*
 * Using Intel's Math Kernel Library (MKL) LP64 configuration (32-bit int)
 * (Generally the fastest library for x86 based architectures)
 *
 * Note that Intel MKL not only supports CBLAS and LAPACK, but also offers:
 *  - a highly optimised discrete/fast Fourier transform (DFT/FFT), which is
 *    used by the saf_utility_fft wrapper.
 *  - a number of additional vector, vector-vector, vector-scalar operations
 *    that are not covered by the CBLAS standard; such as: hadamard products,
 *    element-wise additions/substractions, and the modulus or reciprical of
 *    all vector elements, etc.
 */
# include "mkl.h"

#elif defined(SAF_USE_INTEL_MKL_ILP64)
/*
 * Using Intel's Math Kernel Library (MKL) ILP64 configuration (64-bit int)
 * (Generally the fastest library for x86 based architectures)
 *
 * This 64-bit int version is the one employed by e.g. MATLAB. Therefore, it is
 * required if you wish to build MEX objects using SAF code; see e.g.
 * extras/safmex. In general, the performance of this option is practically the
 * same as "SAF_USE_INTEL_MKL_LP64", but it is slower in some very rare special
 * cases. Therefore, "SAF_USE_INTEL_MKL_LP64" is still the favoured option if
 * you are not planning on building MEX objects using SAF.
 */
# define MKL_ILP64 /**< Indicates to MKL that we've linked the ILP64 variant */
# include "mkl.h"

#elif defined(SAF_USE_OPEN_BLAS_AND_LAPACKE)
/*
 * Using OpenBLAS and the LAPACKE interface
 * (A decent option for both x86 and ARM based architectures)
 *
 * This option provides implementations of the CBLAS/LAPACK functions which have
 * decent performance. However, unlike Intel MKL or Apple Accelerate, it does
 * not offer an optimised DFT/FFT or any other linear algebra functions outside
 * of these standards. Therefore, if you are using this option, consider also
 * using Intel's IPP library for the FFT, along with the SSE/AVX/AVX-512
 * fallback implementations for the other linear algebra options, by also
 * defining: "SAF_USE_INTEL_IPP" and "SAF_ENABLE_SIMD" (see instructions below).
 */
# include "cblas.h"
# include "lapacke.h"

#elif defined(SAF_USE_ATLAS)
/*
 * Using the Automatically Tuned Linear Algebra Software (ATLAS) library
 * (Not recommended, since some saf_utility_veclib functions do not work with
 * ATLAS)
 *
 * Basically, do not use this unless you have to, and if you do, be aware that
 * some linear algebra functions in saf_utility_veclib will exit the program if
 * they are called.
 */
# include "cblas-atlas.h"
# include "clapack.h"
# warning Note: CLAPACK does not include all LAPACK routines!

#elif defined(__APPLE__) && defined(SAF_USE_APPLE_ACCELERATE)
/*
 * Using Apple's Accelerate library
 * (Solid choice for both x86 and ARM, but only works under MacOSX and is not
 * quite as fast as Intel MKL for x86 systems)
 *
 * Note that Apple Accelerate not only supports CBLAS and LAPACK, but also
 * offers:
 *  - an optimised discrete/fast Fourier transform (DFT/FFT), which is used by
 *    the saf_utility_fft wrapper.
 *  - a number of additional vector, vector-vector, vector-scalar operations
 *    that are not covered by the CBLAS standard; such as hadamard products,
 *    element-wise additions/substractions, etc.
 *
 * Unlike Intel MKL, not all even number DFT lengths are supported by vDSP.
 * Therefore, be aware that the default kissFFT library (included in
 * framework/resources) is still used as a fall-back option in such cases.
 */
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
 * does lead to improvements in the following:
 *   - slightly faster DFT/FFT (for saf_utility_fft) compared with the
 *     implementation found in Intel MKL, which are both faster than the DFT/FFT
 *     implementations found in Apple Accelerate vDSP and FFTW.
 *   - this overrides the included resources/speex_resampler with the IPP
 *     resampler, which is marginally faster and more accurate.
 *   - this also overrides certain vector-vector, and vector-scalar operations,
 *     such as element-wise multiplication, addition, scaling etc.
 *
 * Note that the IPP DFT/FFT is overriden by FFTW if SAF_USE_FFTW is defined.
 */
# include "ipp.h"
#endif

#if defined(SAF_USE_FFTW)
/*
 * The use of FFTW is optional, but it is faster than the default kissFFT
 * DFT/FFT implementation. However, if you are on an x86 CPU then the DFT/FFT
 * implementations found in Intel IPP, Intel MKL and Apple Accelerate are
 * usually faster options.
 *
 * Note, SAF uses the single-precision version (fftw3f.a), which is built with:
 *   $ ./configure --enable-float
 *   $ make
 *
 * If SAF_USE_FFTW is defined, then FFTW overrides all of the other available
 * DFT/FFT implementations in the saf_utility_fft wrapper.
 */
# include "fftw3.h"
#endif

#if defined(SAF_ENABLE_SIMD)
/*
 * SAF heavily favours the use of optimised linear algebra routines provided by
 * e.g. Intel MKL or Accelerate, since they optimally employ vectorisation
 * (with SSE/AVX etc.). However, in cases where the employed performance library
 * does not offer an implementation for a particular routine, SAF provides fall-
 * back option(s).
 * SIMD accelerated fall-back options may be enabled with: SAF_ENABLE_SIMD
 *
 * By default SSE, SSE2, and SSE3 intrinsics are employed, unless one of the
 * following compiler flags are defined:
 *    - AVX/AVX2 intrinsics are enabled with: -mavx2
 *    - AVX-512  intrinsics are enabled with: -mavx512f
 *
 * Note that intrinsics require a CPU that supports them (x86_64 architecture)
 * To find out which SIMD intrinsics are supported by your own CPU, use the
 * following terminal command on macOS: $ sysctl -a | grep machdep.cpu.features
 * Or on Linux, use: $ lscpu
 */
# if (defined(__AVX__) && defined(__AVX2__)) || defined(__AVX512F__)
/*
 * Note that AVX/AVX2 requires the '-mavx2' compiler flag
 * Whereas AVX-512 requires the '-mavx512f' compiler flag
 */
#  include <immintrin.h> /* for AVX, AVX2, and/or AVX-512 */
# endif
# if defined(__SSE__) && defined(__SSE2__) && defined(__SSE3__)
#  include <xmmintrin.h>  /* for SSE  */
#  include <emmintrin.h>  /* for SSE2 */
#  include <pmmintrin.h>  /* for SSE3 */
# else
#  error SAF_ENABLE_SIMD requires at least SSE, SSE2 and SSE3 support
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
