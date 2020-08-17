/*
 * Copyright 2018-2020 Leo McCormack
 *
 * This software is dual-licensed.
 *
 * ISC License ( https://choosealicense.com/licenses/isc/ ):
 *   By default, this software is provided permissively under the terms of the
 *   ISC license; since all of the core (non-optional) modules are licensed as
 *   such. Also note that all third-party code, which has been adapted or
 *   adopted by these core modules, is also provided under similar permissive
 *   licensing terms (MIT, BSD, WOL, etc.). The license employed by third-party
 *   code is stated clearly in the code comments at the point of inclusion.
 *
 * GNU GPLv3 License ( https://choosealicense.com/licenses/gpl-3.0/ ):
 *   Including and enabling certain optional modules, which are instead provided
 *   under the copy-left GNU GPLv3 license, will mean that the use of this
 *   software is instead governed by the GNU GPLv3 licencing terms.
 *
 * Note that the license employed by each module is stated clearly below. For
 * more information, refer to the license terms found at the top of each source
 * and header file.
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
 * @brief Main include header for the Spatial_Audio_Framework (SAF)
 *
 * ## Core modules
 *   saf_utilities, saf_cdf4sap, saf_hoa, saf_hrir, saf_sh, saf_vbap, saf_reverb
 *
 * @note All core modules are provided permissively under the terms of the ISC
 *       license. Therefore, by default, the use of this framework is governed
 *       by the ISC licencing terms.
 *
 * ## Optional modules
 *   saf_sofa_reader (ISC License), example (GPLv3 License)
 *
 * @note Some optional modules are provided under the copy-left GPLv3 license.
 *       Therefore, including and enabling these specific optional modules will
 *       mean that the use of this framework is instead governed by the GPLv3
 *       licencing terms.
 *
 * @author Leo McCormack
 * @date 06.04.2018
 */

#ifndef SAF_H_INCLUDED
#define SAF_H_INCLUDED

/*  doxygen module definitions: */
/** @defgroup Utilities    saf_utilities
 *  @brief    Utilities module */
/** @defgroup CDF4SAP      saf_cdf4sap
 *  @brief    Covariance domain framework processing module */
/** @defgroup HOA          saf_hoa
 *  @brief    Higher-order Ambisonics module */
/** @defgroup SH           saf_sh
 *  @brief    Spherical harmonic domain processing module */
/** @defgroup HRIR         saf_hrir
 *  @brief    HRIR/HRTF processing module */
/** @defgroup Reverb       saf_reverb
 *  @brief    Reverb and room simulation module */
/** @defgroup VBAP         saf_vbap
 *  @brief    Vector-base Amplitude Panning (VBAP) module */
/** @defgroup SOFA_Reader  saf_sofa_reader
 *  @brief    SOFA file reading module */


/* ========================================================================== */
/*                                SAF Version                                 */
/* ========================================================================== */
/* Note that this version info is also scanned and employed by CMakeLists.txt */

#define SAF_VERSION_MAJOR 1       /**< Major version */
#define SAF_VERSION_MINOR 1       /**< Minor version */
#define SAF_VERSION_PATCH 3       /**< Patch version */
#define SAF_VERSION_SPECIAL ""    /**< Append text ("alpha", "beta", "" etc.) */
#define MKSTRING_(s) #s           /**< Stringify */
#define MKSTRING(s) MKSTRING_(s)  /**< Stringify */

/** The Spatial_Audio_Framework Version as a number */
#define SAF_VERSION ((SAF_VERSION_MAJOR << 16) + (SAF_VERSION_MINOR << 8) + \
                    SAF_VERSION_PATCH)

/** The Spatial_Audio_Framework Version as a string */
#define SAF_VERSION_STRING MKSTRING(SAF_VERSION_MAJOR) "." MKSTRING( \
                           SAF_VERSION_MINOR) "." MKSTRING(SAF_VERSION_PATCH) \
                           "" SAF_VERSION_SPECIAL

/** The Spatial_Audio_Framework Version as a banner */
#define SAF_VERSION_BANNER \
    "      _____     _____     ______                                    \n" \
    "    (  ____ ) (  ___  ) (  ____  )                                  \n" \
    "    | (    |/ | (   ) | | (    |/                                   \n" \
    "    | (_____  | (___) | | (__                                       \n" \
    "    (_____  ) |  ___  | |  __)                                      \n" \
    "          ) | | (   ) | | (                                         \n" \
    "    /|____) | | |   | | | |                                         \n" \
    "    (_______) |_|   |_| |_|        (Version " SAF_VERSION_STRING ") \n\n"


/* ========================================================================== */
/*                                Core Modules                                */
/* ========================================================================== */

/**
 * SAF Module: Utilities
 *
 * A collection of useful utility functions; including: cross-platform complex
 * number wrappers; optimised linear algebra routines based on CBLAS and LAPACK;
 * FFT wrappers and STFT implementation; IIR/FIR filter coefficients and filter
 * bank designs; lists of common loudspeaker and microphone array coordinates;
 * multi-channel and matrix convolvers; spherical Bessel/Hankel functions
 * (including their derivatives); etc.
 *
 * License: ISC
 */
#define SAF_UTILITIES_MODULE
#include "../modules/saf_utilities/saf_utilities.h"

/**
 * SAF Module: CDf4SAP
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
 *
 * License: ISC
 */
#define SAF_CDF4SAP_MODULE
#include "../modules/saf_cdf4sap/saf_cdf4sap.h"

/**
 * SAF Module: HOA
 *
 * A collection of higher-order Ambisonics related functions; many of which are
 * derived from the Matlab library found in [1] (BSD-3-Clause License).
 *
 * ## Dependencies
 *   saf_utilities, saf_vbap, saf_sh
 *
 * @see [1] https://github.com/polarch/Higher-Order-Ambisonics
 *
 * License: ISC
 */
#define SAF_HOA_MODULE
#include "../modules/saf_hoa/saf_hoa.h"

/**
 * SAF Module: HRIR
 *
 * A collection of head-related impulse-response (HRIR) functions; including:
 * estimation of the interaural time differences (ITDs), conversion of HRIRs to
 * HRTF filterbank coefficients, and HRTF interpolation.
 *
 * ## Dependencies
 *   saf_utilities
 *
 * License: ISC
 */
#define SAF_HRIR_MODULE
#include "../modules/saf_hrir/saf_hrir.h"

/**
 * SAF Module: Reverb
 *
 * A collection of reverb and room simulation algorithms.
 *
 * ## Dependencies
 *   saf_utilities, saf_sh
 *
 * License: ISC
 */
#define SAF_REVERB_MODULE
#include "../modules/saf_reverb/saf_reverb.h"

/**
 * SAF Module: SH
 *
 * A collection of spherical harmonic related functions; many of which have been
 * derived from the Matlab libraries found in [1-3] (BSD-3-Clause License).
 *
 * ## Dependencies
 *   saf_utilities
 *
 * @see [1] https://github.com/polarch/Spherical-Harmonic-Transform
 * @see [2] https://github.com/polarch/Array-Response-Simulator
 * @see [3] https://github.com/polarch/Spherical-Array-Processing
 *
 * License: ISC
 */
#define SAF_SH_MODULE
#include "../modules/saf_sh/saf_sh.h"

/**
 * SAF Module: VBAP
 *
 * VBAP functions largely derived from the Matlab library found in [1]
 * (BSD-3-Clause License).
 *
 * ## Dependencies
 *  saf_utilities
 *
 * @see [1] https://github.com/polarch/Vector-Base-Amplitude-Panning
 *
 * License: ISC
 */
#define SAF_VBAP_MODULE
#include "../modules/saf_vbap/saf_vbap.h"


/* ========================================================================== */
/*                              Optional Modules                              */
/* ========================================================================== */

/**
 * SAF Module: SOFA_Reader
 *
 * A simple SOFA file reader that returns only the bare minimum needed to
 * load HRIR data.
 *
 * ## Enable instructions
 *   Add this pre-processor definition to your project:
 *       SAF_ENABLE_SOFA_READER_MODULE
 *   and ensure that the netcdf library is also linked to your project
 * ## Dependencies
 *   saf_utilities, saf_hrir, netcdf
 *
 * License: ISC
 */
#define SAF_SOFA_READER_MODULE
#ifdef  SAF_ENABLE_SOFA_READER_MODULE
# include "../modules/saf_sofa_reader/saf_sofa_reader.h"
#endif /* SAF_ENABLE_SOFA_READER_MODULE */

/**
* SAF Module: EXAMPLE
*
* Description of module
*
* ## Enable instructions
*   Add this pre-processor definition to your project:
*       SAF_ENABLE_X
* ## Dependencies
*   saf_utilities ...
*
* License: GNU GPLv3
*/
#define SAF_X_MODULE
#ifdef  SAF_ENABLE_X_MODULE
# include "../modules/saf_X/saf_X.h"
#endif /* SAF_ENABLE_X_MODULE */


#endif /* SAF_H_INCLUDED */
