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
 * @note More information can be found in the docs folder regarding where to
 *       find and how to link dependencies
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
 *   If the optional saf_sofa_reader module is enabled and SAF_ENABLE_NETCDF is
 *   defined, then the netcdf library must also be linked along with saf.
 *
 *   Intel IPP may be optionally employed with the flag: SAF_USE_INTEL_IPP
 *
 *   FFTW may be optionally employed with the flag: SAF_USE_FFTW
 *
 *   SIMD intrinsics utilisation may be enabled with: SAF_ENABLE_SIMD
 *    - SSE/SSE2/SSE3 intrinsics are used by default
 *    - AVX/AVX2 intrinsics are enabled with compiler flag: -mavx2
 *    - AVX-512  intrinsics are enabled with compiler flag: -mavx512f
 *   (Note that intrinsics require a CPU that supports them)
 *
 * @author Leo McCormack
 * @date 06.08.2020
 * @license ISC
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
     defined(SAF_USE_GSL) + \
     defined(SAF_USE_APPLE_ACCELERATE)) != 1
# error One (and only one) performance library flag should be defined!
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
 *    used by the saf_utility_fft wrapper [unless this is overriden by the Intel
 *    IPP implementation (SAF_USE_INTEL_IPP), or FFTW (SAF_USE_FFTW)].
 *  - a number of additional vector, vector-vector, vector-scalar operations
 *    that are not covered by the CBLAS standard; such as: hadamard products,
 *    element-wise additions/subtractions, and the modulus or reciprical of
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
 * of these standards. Therefore, consider also using Intel's IPP library or
 * FFTW for the DFT/FFT with: "SAF_USE_INTEL_IPP" or "SAF_USE_FFTW"
 *
 * Note that "SAF_USE_INTEL_IPP" also offers support for certain linear algebra
 * operations not covered by the CBLAS/LAPACK standards, which SAF can leverage.
 *
 * Alternatively, SSE/AVX/AVX-512 fallback implementations for certain linear
 * algebra operations may be enabled with: "SAF_ENABLE_SIMD"
 *
 * More information regarding these additional options can be found below.
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
 *    the saf_utility_fft wrapper  [unless this is overriden by the Intel IPP
 *    implementation (SAF_USE_INTEL_IPP), or FFTW (SAF_USE_FFTW)].
 *  - a number of additional vector, vector-vector, vector-scalar operations
 *    that are not covered by the CBLAS standard; such as hadamard products,
 *    element-wise additions/subtractions, etc.
 *
 * Unlike e.g. Intel MKL's DFT, not all even number DFT lengths are supported by
 * vDSP. Therefore, be aware that the default kissFFT library (included in
 * framework/resources) is still used as a fall-back option in such cases.
 */
# include "Accelerate/Accelerate.h"

#elif defined(SAF_USE_GSL)
/*
 * Using the GNU Scientific Library (GSL)
 *
 * Please feel free to try it out and report back. Note also that certain LAPACK
 * functions used in saf_utility_veclib will need to be swapped out for
 * equivalent functions in GSL
 */
# error Using GNU Scientific Library (GSL) is currently unsupported/incomplete
# include "gsl_cblas.h"

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
 * Note, SAF uses the single-precision version (fftw3f), which is built with:
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
 * following compiler flags are given:
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

#if defined(SAF_ENABLE_SOFA_READER_MODULE)
/*
 * The built-in saf_sofa_open() SOFA file reader has two implementations:
 *    - By default, the function wraps around the "libmysofa" library
 *      (BSD-3-Clause license), which depends on only zlib (which is included
 *      in framework/resources/zlib). The downsides of this option, is that zlib
 *      has file size limits for each chunk (<4GB) and it is quite slow at
 *      decompressing large files.
 *    - If SAF_ENABLE_NETCDF is defined, then an alternative SOFA reader may be
 *      used. This version requires netcdf to be linked to SAF, along with its
 *      dependencies. The netcdf loader gets around the file size limits of
 *      the libmysofa loader and is also approximately 3 times faster.
 *      Therefore, if you intend to load many large SOFA files
 *      (especially microphone arrays or Ambisonic IRs), then this alternative
 *      SOFA reader is either required (to get around the file size limit) or
 *      may be preferred due to the shorter loading times. The downsides of
 *      using the netcdf option is that it is NOT thread-safe! and requires
 *      these additional external libraries to be linked to SAF.
 *
 * Note that the "mysofa" interface, e.g. mysofa_load(), may also be called
 * directly, rather than using saf_sofa_open().
 */
# ifdef SAF_ENABLE_NETCDF
#  include <netcdf.h>
# endif
#endif


/* ========================================================================== */
/*                   Configuration and Status Flags/Strings                   */
/* ========================================================================== */

/* Currently employed performance library: */
#if defined(SAF_USE_INTEL_MKL_LP64)
# define SAF_CURRENT_PERFORMANCE_LIBRARY_STRING "Intel MKL (LP64)"
#elif defined(SAF_USE_INTEL_MKL_ILP64)
# define SAF_CURRENT_PERFORMANCE_LIBRARY_STRING "Intel MKL (ILP64)"
#elif defined(SAF_USE_OPEN_BLAS_AND_LAPACKE)
# define SAF_CURRENT_PERFORMANCE_LIBRARY_STRING "OpenBLAS with LAPACKE"
#elif defined(SAF_USE_ATLAS)
# define SAF_CURRENT_PERFORMANCE_LIBRARY_STRING "ATLAS"
#elif defined(__APPLE__) && defined(SAF_USE_APPLE_ACCELERATE)
# define SAF_CURRENT_PERFORMANCE_LIBRARY_STRING "Apple Accelerate"
#else
# define SAF_CURRENT_PERFORMANCE_LIBRARY_STRING "NONE"
#endif

/* Status of Intel IPP */
#if defined(SAF_USE_INTEL_IPP)
# define SAF_INTEL_IPP_STATUS_STRING "Enabled"
#else
# define SAF_INTEL_IPP_STATUS_STRING "Disabled"
#endif

/* Status of FFTW */
#if defined(SAF_USE_FFTW)
# define SAF_FFTW_STATUS_STRING "Enabled"
#else
# define SAF_FFTW_STATUS_STRING "Disabled"
#endif

/* Status of SIMD intrinsics */
#if defined(SAF_ENABLE_SIMD)
# define SAF_SIMD_STATUS_STRING "Enabled"
/* Which SIMD intrinsics are currently enabled? */
# if defined(__AVX512F__)
#  define SAF_ENABLED_SIMD_INTRINSICS_STRING "SSE, SSE2, SSE3, AVX, AVX2, AVX512F"
# elif defined(__AVX__) && defined(__AVX2__)
#  define SAF_ENABLED_SIMD_INTRINSICS_STRING "SSE, SSE2, SSE3, AVX, AVX2"
# elif defined(__SSE__) && defined(__SSE2__) && defined(__SSE3__)
#  define SAF_ENABLED_SIMD_INTRINSICS_STRING "SSE, SSE2, SSE3"
# else
#  define SAF_ENABLED_SIMD_INTRINSICS_STRING "None"
# endif
#else
# define SAF_SIMD_STATUS_STRING "Disabled"
# define SAF_ENABLED_SIMD_INTRINSICS_STRING "None"
#endif

/* Status of netCDF */
#if defined(SAF_ENABLE_NETCDF)
# define SAF_NETCDF_STATUS_STRING "Enabled"
#else
# define SAF_NETCDF_STATUS_STRING "Disabled"
#endif

/** Current configuration information */
#define SAF_EXTERNALS_CONFIGURATION_STRING  \
    "Current SAF externals configuration: "                               "\n" \
    " - Performance library: " SAF_CURRENT_PERFORMANCE_LIBRARY_STRING     "\n" \
    " - Intel IPP status:    " SAF_INTEL_IPP_STATUS_STRING                "\n" \
    " - FFTW status:         " SAF_FFTW_STATUS_STRING                     "\n" \
    " - SIMD status:         " SAF_SIMD_STATUS_STRING                     "\n" \
    " - Enabled intrinsics:  " SAF_ENABLED_SIMD_INTRINSICS_STRING         "\n" \
    " - netCDF status:       " SAF_NETCDF_STATUS_STRING                   "\n"


#endif /* __SAF_EXTERNALS_H_INCLUDED__ */
