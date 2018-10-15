/*
 Copyright 2016-2018 Leo McCormack
 
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
 *     saf_sensorarray_presets.c
 * Description:
 *     Contains a collection of useful memory allocation functions and cross-platform
 *     complex number wrappers. Optimised linear algebra routines utilising BLAS and LAPACK
 *     are also included.
 * Dependencies:
 *     Windows users only: Intel's MKL must be installed, which can be freely aquired via:
 *     https://software.intel.com/en-us/articles/free-ipsxe-tools-and-libraries
 *     Mac users only: saf_utilities will utilise Apple's Accelerate library.
 * Author, date created:
 *     Leo McCormack, 11.07.2016
 */

#ifndef __SAF_SENSORARRAY_PRESETS_INCLUDED__
#define __SAF_SENSORARRAY_PRESETS_INCLUDED__

#ifdef __cplusplus
extern "C" {
#endif
    
#ifndef UTIL_DEFAULT_SH_ORDER
  #define UTIL_DEFAULT_SH_ORDER ( 7 )
#endif
    
/* sensor array coordinates */
extern const float __Aalto_Hydrophone_coords_rad[4][2];
extern const float __Sennheiser_Ambeo_coords_rad[4][2];
extern const float __Core_Sound_TetraMic_coords_rad[4][2];
extern const float __Sound_field_SPS200_coords_rad[4][2];
extern const float __Zylia1D_coords_rad[19][2];
extern const float __Eigenmike32_coords_rad[32][2];
extern const float __DTU_mic_coords_rad[52][2];
extern const float __default_coords_rad[(UTIL_DEFAULT_SH_ORDER+1)*(UTIL_DEFAULT_SH_ORDER+1)][2];
extern const float __default_SENSORcoords64_rad[64][2];

/* sensor array maximum order */
extern const int __Aalto_Hydrophone_maxOrder;
extern const int __Sennheiser_Ambeo_maxOrder;
extern const int __Core_Sound_TetraMic_maxOrder;
extern const int __Sound_field_SPS200_maxOrder;
extern const int __Zylia_maxOrder;
extern const int __Eigenmike32_maxOrder;
extern const int __DTU_mic_maxOrder;
  
/* sensor array frequency ranges for each SH order - should only be used as a rough estimate.
 * The upper frequency limits were selected as the point where the spatial correlation went <0.9.
 * The lower frequency limits were selected as the point where the level difference exceeded 6dB
 * (assuming a 15dB maximum amplification with the Tikhonov regularisation method for all mics).
 *
 * For more information on determining the usable frequency range per spherical harmonic order, for
 * a given microphone array, the reader is directed to:
 * Moreau, S., Daniel, J., & Bertet, S. (2006, May). 3D sound field recording with higher order
 * ambisonics–Objective measurements and validation of a 4th order spherical microphone.
 * In 120th Convention of the AES (pp. 20-23).
 */
extern const float __Zylia_freqRange[4];
extern const float __Eigenmike32_freqRange[6];
extern const float __DTU_mic_freqRange[10];

#ifdef __cplusplus
}
#endif


#endif /* __SAF_SENSORARRAY_PRESETS_INCLUDED__ */


