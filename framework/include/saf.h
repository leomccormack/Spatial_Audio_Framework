/*
 * Copyright 2018-2019 Leo McCormack
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
 * Filename: saf.h (include header)
 * --------------------------------
 * Main include header for the spatial audio framework. Instructions on how
 * to enable the framework modules is provided below.
 *
 * Included modules:
 *     afSTFTlib, saf_cdf4sap, saf_hoa, saf_hrir, saf_sh, saf_sofa_reader,
 *     saf_utilities, saf_vbap
 * Author, date created:
 *     Leo McCormack, 06.04.2018
 */

#ifndef SAF_H_INCLUDED
#define SAF_H_INCLUDED

/* ========================================================================== */
/*                                Core Modules                                */
/* ========================================================================== */

/*
 * SAF Module: Utilities
 * ---------------------
 * Contains a collection of useful memory allocation functions and cross-
 * platform complex number wrappers. Optimised linear algebra routines utilising
 * BLAS and LAPACK are also included.
 *
 * Enable instructions:
 *     Cannot be disabled.
 * Dependencies:
 *     A performance library comprising CBLAS and LAPACK routines is required
 *     by the module and, thus, also by the SAF framework as a whole.
 *     Add one of the following FLAGS to your project's preprocessor definitions
 *     list, in order to enable one of these suitable performance libraries,
 *     which must also be linked correctly to your project.
 *         SAF_USE_INTEL_MKL
 *             to enable Intel's Math Kernal Library
 *         SAF_USE_ATLAS_WITH_LAPACK
 *             to enable ATLAS BLAS routines with netlib's LAPACK
 *             see: http://math-atlas.sourceforge.net/atlas_install/node8.html
 *         SAF_USE_OPENBLAS_WITH_LAPACK
 *             to enable OpenBLAS and use netlib's LAPACK
 *     More information can be found here:
 *         https://github.com/leomccormack/Spatial_Audio_Framework
 *     Mac users only:
 *         saf_utilities will employ Apple's Accelerate library by default, if
 *         none of the above FLAGS are defined.
 */
#include "../modules/saf_utilities/saf_utilities.h"


/* ========================================================================== */
/*                              Optional Modules                              */
/* ========================================================================== */

/*
 * SAF Module: afSTFT
 * ------------------
 * The Alias-free STFT implementation by Juha Vilkamo, with some minor changes.
 * The Original source code can be found here:
 *     https://github.com/jvilkamo/afSTFT
 *
 * Enable instructions:
 *     Add this pre-processor definition to your project: SAF_ENABLE_AFSTFT
 * Dependencies:
 *     saf_utilities
 */
#ifdef SAF_ENABLE_AFSTFT
# include "../resources/afSTFT/afSTFTlib.h"
#endif

/*
 * SAF Module: CDF4SAP
 * -------------------
 * Covarience Domain Framework for Spatial Audio Processing (CDF4SAP). A C
 * implementation of the framework described in [1].
 *
 * Enable instructions:
 *     Add this pre-processor definition to your project: SAF_ENABLE_CDF4SAP
 * Dependencies:
 *     saf_utilities
 *
 * [1] Vilkamo, J., Bäckström, T., & Kuntz, A. (2013). Optimized covariance
 *     domain framework for time–frequency processing of spatial audio. Journal
 *     of the Audio Engineering Society, 61(6), 403-411.
 */
#ifdef SAF_ENABLE_CDF4SAP
# include "../modules/saf_cdf4sap/saf_cdf4sap.h"
#endif

/*
 * SAF Module: HOA
 * ---------------
 * A collection of higher-order Ambisonics related functions. Some of which are
 * derived from the Matlab library by Archontis Politis, found here:
 *     https://github.com/polarch/Higher-Order-Ambisonics
 *
 * Enable instructions:
 *     Add this pre-processor definition to your project: SAF_ENABLE_HOA
 * Dependencies:
 *     saf_utilities, saf_vbap, saf_sh
 */
#ifdef SAF_ENABLE_HOA
# include "../modules/saf_hoa/saf_hoa.h"
#endif

/*
 * SAF Module: HRIR
 * ----------------
 * A collection of head-related impulse-response (HRIR) functions. Including
 * estimation of the interaural time differences (ITDs), conversion of HRIRs to
 * HRTF filterbank coefficients, and HRTF interpolation utilising amplitude-
 * normalised VBAP gains.
 *
 * Enable instructions:
 *     Add this pre-processor definition to your project: SAF_ENABLE_HRIR
 * Dependencies:
 *     saf_utilities, afSTFTlib
 */
#ifdef SAF_ENABLE_HRIR
# include "../modules/saf_hrir/saf_hrir.h"
#endif

/*
 * SAF Module: SOFA Reader
 * -----------------------
 * A simple SOFA file reader that returns only the bare minimum.
 *
 * Enable instructions:
 *     Add this pre-processor definition to your project: SAF_ENABLE_SOFA_READER
 * Dependencies:
 *     netcdf library
 */
#ifdef SAF_ENABLE_SOFA_READER
# include "../modules/saf_hrir/saf_sofa_reader.h"
#endif

/*
 * SAF Module: SH
 * --------------
 * A collection of spherical harmonic related functions. Many of which have been
 * derived from Matlab libraries by Archontis Politis; found here:
 *     https://github.com/polarch/Spherical-Harmonic-Transform
 *     https://github.com/polarch/Array-Response-Simulator
 *     https://github.com/polarch/Spherical-Array-Processing
 *
 * Enable instructions:
 *     Add this pre-processor definition to your project: SAF_ENABLE_SH
 * Dependencies:
 *     saf_utilities
 */
#ifdef SAF_ENABLE_SH
# include "../modules/saf_sh/saf_sh.h"
#endif

/*
 * SAF Module: VBAP
 * ----------------
 * VBAP functions largely derived from the MATLAB library by Archontis Politis,
 * found here: https://github.com/polarch/Vector-Base-Amplitude-Panning
 *
 * Enable instructions:
 *     Add this pre-processor definition to your project: SAF_ENABLE_VBAP 
 * Dependencies:
 *     saf_utilities
 */
#ifdef SAF_ENABLE_VBAP
# include "../modules/saf_vbap/saf_vbap.h"
#endif


#endif /* SAF_H_INCLUDED */
