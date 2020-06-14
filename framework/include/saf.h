/*
 * Copyright 2018-2020 Leo McCormack
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
 * @file saf.h
 * @brief Main include header for the Spatial_Audio_Framework
 *
 * ## Core modules
 *   saf_utilities, saf_cdf4sap, saf_hoa, saf_hrir, saf_sh, saf_vbap, saf_reverb
 *
 * ## Optional modules
 *   saf_sofa_reader
 *
 * @author Leo McCormack
 * @date 06.04.2018
 */

#ifndef SAF_H_INCLUDED
#define SAF_H_INCLUDED

#define SAF_VERSION_MAJOR 1       /**< Major version */
#define SAF_VERSION_MINOR 1       /**< Minor version */
#define SAF_VERSION_PATCH 0       /**< Patch version */
#define SAF_VERSION_SPECIAL ""    /**< Append text, "alpha", "beta", "" etc. */
#define MKSTRING_(s) #s           /**< Stringify */
#define MKSTRING(s) MKSTRING_(s)  /**< Stringify */

/** The Spatial_Audio_Framework Version as a number */
#define SAF_VERSION ((SAF_VERSION_MAJOR << 16) + (SAF_VERSION_MINOR << 8) + \
                    SAF_VERSION_PATCH)

/** The Spatial_Audio_Framework Version as a string */
#define SAF_VERSION_STRING MKSTRING(SAF_VERSION_MAJOR)"."MKSTRING( \
                           SAF_VERSION_MINOR)"."MKSTRING(SAF_VERSION_PATCH)""\
                           SAF_VERSION_SPECIAL

/** The Spatial_Audio_Framework Version as a banner */
#define SAF_VERSION_BANNER \
    "      _____     _____     ______                                    \n" \
    "    (  ____ ) (  ___  ) (  ____  )                                  \n" \
    "    | (    |/ | (   ) | | (    |/                                   \n" \
    "    | (_____  | (___) | | (__                                       \n" \
    "    (_____  ) |  ___  | |  __)                                      \n" \
    "          ) | | (   ) | | (                                         \n" \
    "    /|____) | | |   | | | |                                         \n" \
    "    (_______) |_|   |_| |_|        (Version "SAF_VERSION_STRING") \n\n"


/* ========================================================================== */
/*                                Core Modules                                */
/* ========================================================================== */

/**
 * SAF Module: Utilities
 *
 * Contains a collection of useful memory allocation functions and cross-
 * platform complex number wrappers. Optimised linear algebra routines utilising
 * BLAS and LAPACK are also included.
 *
 * ## Dependencies
 *   A performance library comprising CBLAS and LAPACK routines is required by
 *   the module and, thus, also by the SAF framework as a whole. Add one of the
 *   following FLAGS to your project's preprocessor definitions list in order
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
 */
#define SAF_MODULE_UTILITIES
#include "saf_utilities.h"

/**
 * SAF Module: CDF4SAP
 *
 * Covarience Domain Framework for Spatial Audio Processing (CDF4SAP). A C
 * implementation of the framework described in [1].
 *
 * ## Dependencies
 *   saf_utilities
 *
 * @see [1] Vilkamo, J., Backstrom, T., & Kuntz, A. (2013). Optimized covariance
 *          domain framework for time--frequency processing of spatial audio.
 *          Journal of the Audio Engineering Society, 61(6), 403-411.
 */
#define SAF_MODULE_CDF4SAP
#include "saf_cdf4sap.h"

/**
 * SAF Module: HOA
 *
 * A collection of higher-order Ambisonics related functions. Some of which are
 * derived from the Matlab library by Archontis Politis, found here:
 *  - https://github.com/polarch/Higher-Order-Ambisonics
 *
 * ## Dependencies
 *   saf_utilities, saf_vbap, saf_sh
 */
#define SAF_MODULE_HOA
#include "saf_hoa.h"

/**
 * SAF Module: HRIR
 *
 * A collection of head-related impulse-response (HRIR) functions. Including
 * estimation of the interaural time differences (ITDs), conversion of HRIRs to
 * HRTF filterbank coefficients, and HRTF interpolation utilising amplitude-
 * normalised VBAP gains.
 *
 * ## Dependencies
 *   saf_utilities
 */
#define SAF_MODULE_HRIR
#include "saf_hrir.h"

/**
 * SAF Module: REVERB
 *
 * A collection of reverb and room simulation algorithms.
 *
 * ## Dependencies
 *   saf_utilities, saf_sh
 */
#define SAF_MODULE_REVERB
#include "saf_reverb.h"

/**
 * SAF Module: SH
 *
 * A collection of spherical harmonic related functions. Many of which have been
 * derived from Matlab libraries by Archontis Politis, found here:
 *   - https://github.com/polarch/Spherical-Harmonic-Transform
 *   - https://github.com/polarch/Array-Response-Simulator
 *   - https://github.com/polarch/Spherical-Array-Processing
 *
 * ## Dependencies
 *   saf_utilities
 */
#define SAF_MODULE_SH
#include "saf_sh.h"

/**
 * SAF Module: VBAP
 *
 * VBAP functions largely derived from the MATLAB library by Archontis Politis,
 * found here:
 *  - https://github.com/polarch/Vector-Base-Amplitude-Panning
 *
 * ## Dependencies
 *  saf_utilities
 */
#define SAF_MODULE_VBAP
#include "saf_vbap.h"


/* ========================================================================== */
/*                              Optional Modules                              */
/* ========================================================================== */

/**
 * SAF Module: SOFA Reader
 *
 * A simple SOFA file reader that returns only the bare minimum needed to
 * load HRIR data.
 *
 * ## Enable instructions
 *   Add this pre-processor definition to your project: SAF_ENABLE_SOFA_READER
 *   and ensure that netcdf is also linked to your project
 * ## Dependencies
 *   saf_utilities, saf_hrir, netcdf
 */
#define SAF_MODULE_SOFA_READER
#ifdef  SAF_ENABLE_SOFA_READER
# include "saf_sofa_reader.h"
#endif /* SAF_ENABLE_SOFA_READER */


#endif /* SAF_H_INCLUDED */
