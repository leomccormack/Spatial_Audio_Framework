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
 * @file saf_utility_sensorarray_presets.h
 * @brief A collection of microphone array sensor directions
 * 
 * @author Leo McCormack
 * @date 11.07.2016
 * @license ISC
 */

#ifndef __SAF_SENSORARRAY_PRESETS_INCLUDED__
#define __SAF_SENSORARRAY_PRESETS_INCLUDED__

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */
    
/* ========================================================================== */
/*                 Microphone/Hydrophone Array Configurations                 */
/* ========================================================================== */
/*
 * NOTE: All microphone array sensor directions are given in radians, and in
 * the [azimuth, elevation] convention. [0 0] is looking directly in-front
 * (positive x-axis), with positive elevations looking upwards (positive z-axis)
 * and positive azimuth angles looking leftwards (positive y-axis).
 * In other words: the convention used by SAF follows the "right-hand-rule".
 */
    
/**
 * Sensor array coordinates for the custom hydrophone array made at Aalto
 * University [1]
 *
 * @see [1] Delikaris-Manias, S., McCormack, L., Huhtakallio, I., & Pulkki, V.
 *          (2018, May). Real-time underwater spatial audio: a feasibility
 *          study. In Audio Engineering Society Convention 144. Audio
 *          Engineering Society.
 */
extern const float __Aalto_Hydrophone_coords_rad[4][2];
/**
 * Sensor array coordinates for the Sennheiser Ambeo */
extern const float __Sennheiser_Ambeo_coords_rad[4][2];
/**
 * Sensor array coordinates for the Core Sound TetraMic */
extern const float __Core_Sound_TetraMic_coords_rad[4][2];
/**
 * Sensor array coordinates for the Sound-field SPS200 */
extern const float __Sound_field_SPS200_coords_rad[4][2];
/**
 * Sensor array coordinates for the Zoom H3VR */
extern const float __Zoom_H3VR_coords_rad[4][2];
/**
 * Sensor array coordinates for the Zylia mic */
extern const float __Zylia1D_coords_rad[19][2];
/**
 * Sensor array coordinates for the Eigenmike32 */
extern const float __Eigenmike32_coords_rad[32][2];
/**
 * Sensor array coordinates for the custom 52-sensor array built at the
 * Technical University of Denmark (DTU) [1]
 *
 * @see [1] Marschall, M., 2014. Capturing and reproducing realistic acoustic
 *          scenes for hearing research (Doctoral dissertation, Technical
 *          University of Denmark, Department of Electrical Engineering). */
extern const float __DTU_mic_coords_rad[52][2];
/**
 * Default sensor array coordinates */
extern const float __default_SENSORcoords64_rad[64][2];
/**
 * Max spherical harmonic order for the custom hydrophone array made at Aalto
 * University */
extern const int __Aalto_Hydrophone_maxOrder;
/**
 * Max spherical harmonic order for the Sennheiser Ambeo */
extern const int __Sennheiser_Ambeo_maxOrder;
/**
 * Max spherical harmonic order for the Core Sound TetraMic */
extern const int __Core_Sound_TetraMic_maxOrder;
/**
 * Max spherical harmonic order for the Sound-field SPS200 */
extern const int __Sound_field_SPS200_maxOrder;
/**
 * Max spherical harmonic order for the Zylia mic */
extern const int __Zylia_maxOrder;
/**
 * Max spherical harmonic order for the Eigenmike32 */
extern const int __Eigenmike32_maxOrder;
/**
 * Max spherical harmonic order for the custom 52-sensor array built at the
 * Technical University of Denmark (DTU) */
extern const int __DTU_mic_maxOrder;
/**
 * Sensor array frequency ranges for each SH order, for the Zylia array (should
 * only be used as a rough estimate).
 *
 * The upper frequency limits were selected as the point where
 * the spatial correlation went <0.9. The lower frequency limits were selected
 * as the point where the level difference exceeded 6dB (assuming a 15dB maximum
 * amplification with the Tikhonov regularisation method for all mics).
 *
 * For more information on determining the usable frequency range per spherical
 * harmonic order, for a given microphone array, the reader is directed to [1].
 *
 * @see [1] Moreau, S., Daniel, J., & Bertet, S. (2006, May). 3D sound field
 *          recording with higher order ambisonics--objective measurements and
 *          validation of a 4th order spherical microphone. In 120th Convention
 *          of the AES (pp. 20-23).
 */
extern const float __Zylia_freqRange[4];
/**
 * Sensor array frequency ranges for each SH order, for the Eigenmike32 (should
 * only be used as a rough estimate) */
extern const float __Eigenmike32_freqRange[6];
/**
 * Sensor array frequency ranges for each SH order, for the DTU mic (should
 * only be used as a rough estimate) */
extern const float __DTU_mic_freqRange[10];


#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */

#endif /* __SAF_SENSORARRAY_PRESETS_INCLUDED__ */

/**@} */ /* doxygen addtogroup Utilities */
