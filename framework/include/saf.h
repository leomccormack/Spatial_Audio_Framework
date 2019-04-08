/*
 Copyright 2018 Leo McCormack
 
 Permission to use, copy, modify, and/or distribute this software for any purpose with or
 without fee is hereby granted, provided that the above copyright notice and this permission
 notice appear in all copies.
 
 THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES WITH REGARD TO
 THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS. IN NO EVENT
 SHALL THE AUTHOR BE LIABLE FOR ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR
 ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF
 CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE
 OR PERFORMANCE OF THIS SOFTWARE.
*/
/*
 * Filename:
 *     saf.h (include header)
 * Description:
 *     Main include header for the spatial audio framework. Instructions on how to enable
 *     the framework features is provided below.
 * Dependencies:
 *     afSTFTlib, saf_cdf4sap, saf_hoa, saf_hrir, saf_sh, saf_sofa_reader, saf_utilities,
 *     saf_vbap
 * Author, date created:
 *     Leo McCormack, 06.04.2018
 */

#ifndef SAF_H_INCLUDED
#define SAF_H_INCLUDED


/* SAF Utilities:
 *     Contains a collection of useful memory allocation functions and cross-platform
 *     complex number wrappers. Optimised linear algebra routines utilising BLAS and LAPACK
 *     are also included.
 * Enable instructions:
 *     Cannot be disabled.
 * Dependencies:
 *     Windows users only: custom Intel MKL '.lib/.dll' files are required.
 *     Mac users only: saf_utilities will utilise Apple's Accelerate library by default.
 *     However, Mac users may elect to use a custom Intel MKL '.dylib' instead.
 *     Further instructions for both Windows/Mac users can be found here:
 *     https://github.com/leomccormack/Spatial_Audio_Framework 
 */
#include "saf_utilities.h"


/* afSTFT:
 *     The Alias-free STFT implementation by Juha Vilkamo, with some minor changes. The
 *     Orginal source code can be found here: https://github.com/jvilkamo/afSTFT
 * Enable instructions:
 *     Place: #define SAF_ENABLE_AFSTFT, before: #include "saf.h"
 * Dependencies:
 *     saf_utilities
 */
#ifdef SAF_ENABLE_AFSTFT
  #include "afSTFTlib.h"
#endif


/* SAF CDF4SAP:
 *     Covarience Domain Framework for Spatial Audio Processing (CDF4SAP). A C
 *     implementation of the framework described in:
 *         Vilkamo, J., Bäckström, T., & Kuntz, A. (2013). Optimized covariance domain
 *         framework for time–frequency processing of spatial audio. Journal of the Audio
 *         Engineering Society, 61(6), 403-411.
 * Enable instructions:
 *     Place: #define SAF_ENABLE_CDF4SAP, before: #include "saf.h"
 * Dependencies:
 *     saf_utilities
 */
#ifdef SAF_ENABLE_CDF4SAP
  #include "saf_cdf4sap.h"
#endif


/* SAF HOA:
 *     A collection of higher-order Ambisonics related functions. Largely derived from the
 *     Matlab library by Archontis Politis; found here:
 *     https://github.com/polarch/Higher-Order-Ambisonics
 * Enable instructions:
 *     Place: #define SAF_ENABLE_HOA, before: #include "saf.h"
 * Dependencies:
 *     saf_utilities, saf_vbap, saf_sh
 */
#ifdef SAF_ENABLE_HOA
  #include "saf_hoa.h"
#endif


/* SAF HRIR:
 *     A collection of head-related impulse-response (HRIR) functions. Including estimation
 *     of the interaural time differences (ITDs), conversion of HRIRs to HRTF filterbank
 *     coefficients, and HRTF interpolation utilising amplitude-normalised VBAP gains.
 * Enable instructions:
 *     Place: #define SAF_ENABLE_HRIR, before: #include "saf.h"
 * Dependencies:
 *     saf_utilities, afSTFTlib
 */
#ifdef SAF_ENABLE_HRIR
  #include "saf_hrir.h"
#endif


/* SAF SOFA Reader:
 *     A simple header only SOFA file reader that returns only the bare minimum.
 * Enable instructions:
 *     Place: #define SAF_ENABLE_SOFA_READER, before: #include "saf.h"
 * Dependencies:
 *     netcdf library
 */
#ifdef SAF_ENABLE_SOFA_READER
  #include "saf_sofa_reader.h"
#endif


/* SAF SH:
 *     A collection of spherical harmonic related functions. Largely derived from the
 *     Matlab library by Archontis Politis; found here:
 *     https://github.com/polarch/Spherical-Harmonic-Transform
 *     and MATLAB code by Symeon Delikaris-Manias.
 * Enable instructions:
 *     Place: #define SAF_ENABLE_SH, before: #include "saf.h"
 * Dependencies:
 *     saf_utilities
 */
#ifdef SAF_ENABLE_SH
  #include "saf_sh.h"
#endif


/* SAF VBAP:
 *     A collection of vector-based amplitude panning (VBAP) functions. Largely derived from
 *     the Matlab library by Archontis Politis; found here:
 *     https://github.com/polarch/Vector-Base-Amplitude-Panning
 * Enable instructions:
 *     Place: #define SAF_ENABLE_VBAP, before: #include "saf.h"
 * Dependencies:
 *     saf_utilities
 */
#ifdef SAF_ENABLE_VBAP 
  #include "saf_vbap.h"
#endif


#endif /* SAF_H_INCLUDED */




