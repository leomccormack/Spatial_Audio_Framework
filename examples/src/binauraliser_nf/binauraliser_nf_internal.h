/*
 * Copyright 2022 Michael McCrea, Leo McCormack
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
 * @file: binauraliser_nf_internal.h
 * @brief Convolves input audio (up to 64 channels) with interpolated HRTFs in
 *        the time-frequency domain, and applies optional near-field binaural
 *        filtering, as described in [1].
 *
 * The HRTFs are interpolated by applying amplitude-preserving VBAP gains to the
 * HRTF magnitude responses and inter-aural time differences (ITDs)
 * individually, before being re-combined. The example also allows the user to
 * specify an external SOFA file for the convolution, and rotations of the
 * source directions to accomodate head-tracking.
 *
 * @see [1] S. Spagnol, E. Tavazzi, and F. Avanzini, “Distance rendering and
 *          perception of nearby virtual sound sources with a near-field filter
 *          model,” Applied Acoustics, vol. 115, pp. 61–73, Jan. 2017,
 *          doi: 10.1016/j.apacoust.2016.08.015.
 *
 * @author Michael McCrea, Leo McCormack
 * @date 22.02.2022
 * @license ISC
 */

#ifndef __BINAURALISER_NF_INTERNAL_H_INCLUDED__
#define __BINAURALISER_NF_INTERNAL_H_INCLUDED__

#include "../binauraliser/binauraliser_internal.h"
#include "binauraliser_nf.h"    /* Include header for this example */
#include "saf.h"                /* Main include header for SAF */
#include "saf_externals.h"      /* To also include SAF dependencies (cblas etc.) */

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/* ========================================================================== */
/*                                 Structures                                 */
/* ========================================================================== */

/**
 * Main structure for binauraliserNF. Contains all variables from binauraliser
 * (audio buffers, afSTFT, HRTFs, internal variables, flags, user parameters) plus
 * those specific to the near field variant.
 * FREQUENCY DOMAIN implementation.
 */
typedef struct _binauraliserNF {
    //struct _binauraliser; /**< "inherit" member vars of binauraliser struct - requires -fms-extensions compile flag */

    /* The following variables MUST match those of the _binauraliser struct */

    /* audio buffers */
    float** inputFrameTD;            /**< time-domain input frame; #MAX_NUM_INPUTS x #BINAURALISER_FRAME_SIZE */
    float** outframeTD;              /**< time-domain output frame; #NUM_EARS x #BINAURALISER_FRAME_SIZE */
    float_complex*** inputframeTF;   /**< time-frequency domain input frame; #HYBRID_BANDS x #MAX_NUM_INPUTS x #TIME_SLOTS */
    float_complex*** outputframeTF;  /**< time-frequency domain input frame; #HYBRID_BANDS x #NUM_EARS x #TIME_SLOTS */
    _Atomic_FLOAT32 fs;              /**< Host sampling rate, in Hz */
    int firstInit;                   /**< flag, 1: `_init()` function has never been called, 0: `_init()` function has been called */
    float freqVector[HYBRID_BANDS];  /**< Frequency vector (filterbank centre frequencies) */
    void* hSTFT;                     /**< afSTFT handle */
    
    /* sofa file info */
    char* sofa_filepath;             /**< absolute/relevative file path for a sofa file */
    float* hrirs;                    /**< time domain HRIRs; FLAT: N_hrir_dirs x #NUM_EARS x hrir_len */
    float* hrir_dirs_deg;            /**< directions of the HRIRs in degrees [azi elev]; FLAT: N_hrir_dirs x 2 */
    _Atomic_INT32 N_hrir_dirs;       /**< number of HRIR directions in the current sofa file */
    _Atomic_INT32 hrir_len;          /**< length of the loaded HRIRs, in samples */
    _Atomic_INT32 hrir_orig_fs;      /**< Can be different from hrir_fs if HRIRs were resampled */
    _Atomic_INT32 hrir_fs;           /**< sampling rate of the HRIRs being used for processing (after any resampling) */
    float* weights;                  /**< Integration weights for the HRIR measurement grid; N_hrir_dirs x 1 */
    
    /* vbap gain table */
    int hrtf_vbapTableRes[2];        /**< [0] azimuth, and [1] elevation grid resolution, in degrees */
    int N_hrtf_vbap_gtable;          /**< Number of interpolation weights/directions */
    int* hrtf_vbap_gtableIdx;        /**< N_hrtf_vbap_gtable x 3 */
    float* hrtf_vbap_gtableComp;     /**< N_hrtf_vbap_gtable x 3 */
    
    /* hrir filterbank coefficients */
    float* itds_s;                   /**< interaural-time differences for each HRIR (in seconds); nBands x 1 */
    float_complex* hrtf_fb;          /**< hrtf filterbank coefficients; nBands x nCH x N_hrirs */
    float* hrtf_fb_mag;              /**< magnitudes of the hrtf filterbank coefficients; nBands x nCH x N_hrirs */
    float_complex hrtf_interp[MAX_NUM_INPUTS][HYBRID_BANDS][NUM_EARS]; /**< Interpolated HRTFs */
    
    /* flags/status */
    _Atomic_CODEC_STATUS codecStatus;          /**< see #CODEC_STATUS */
    _Atomic_FLOAT32 progressBar0_1;            /**< Current (re)initialisation progress, between [0..1] */
    char* progressBarText;                     /**< Current (re)initialisation step, string */
    _Atomic_PROC_STATUS procStatus;            /**< see #PROC_STATUS */
    _Atomic_INT32 recalc_hrtf_interpFLAG[MAX_NUM_INPUTS]; /**< 1: re-calculate/interpolate the HRTF, 0: do not */
    _Atomic_INT32 reInitHRTFsAndGainTables;    /**< 1: reinitialise the HRTFs and interpolation tables, 0: do not */
    _Atomic_INT32 recalc_M_rotFLAG;            /**< 1: re-calculate the rotation matrix, 0: do not */
    
    /* misc. */
    float src_dirs_rot_deg[MAX_NUM_INPUTS][2]; /**< Intermediate rotated source directions, in degrees */
    float src_dirs_rot_xyz[MAX_NUM_INPUTS][3]; /**< Intermediate rotated source directions, as unit-length Cartesian coordinates */
    float src_dirs_xyz[MAX_NUM_INPUTS][3];     /**< Intermediate source directions, as unit-length Cartesian coordinates  */
    int nTriangles;                            /**< Number of triangles in the convex hull of the spherical arrangement of HRIR directions/points */
    _Atomic_INT32 new_nSources;                /**< New number of input/source signals (current value will be replaced by this after next re-init) */

    /* user parameters */
    _Atomic_INT32 nSources;                            /**< Current number of input/source signals */
    _Atomic_FLOAT32 src_dirs_deg[MAX_NUM_INPUTS][2];   /**< Current source/panning directions, in degrees */
    _Atomic_INTERP_MODES interpMode;                   /**< see #INTERP_MODES */
    _Atomic_INT32 useDefaultHRIRsFLAG;                 /**< 1: use default HRIRs in database, 0: use those from SOFA file */
    _Atomic_INT32 enableHRIRsDiffuseEQ;                /**< flag to diffuse-field equalisation to the currently loaded HRTFs */
    _Atomic_INT32 enableRotation;                      /**< 1: enable rotation, 0: disable */
    _Atomic_FLOAT32 yaw;                               /**< yaw (Euler) rotation angle, in degrees */
    _Atomic_FLOAT32 roll;                              /**< roll (Euler) rotation angle, in degrees */
    _Atomic_FLOAT32 pitch;                             /**< pitch (Euler) rotation angle, in degrees */
    _Atomic_INT32 bFlipYaw;                            /**< flag to flip the sign of the yaw rotation angle */
    _Atomic_INT32 bFlipPitch;                          /**< flag to flip the sign of the pitch rotation angle */
    _Atomic_INT32 bFlipRoll;                           /**< flag to flip the sign of the roll rotation angle */
    _Atomic_INT32 useRollPitchYawFlag;                 /**< rotation order flag, 1: r-p-y, 0: y-p-r */
    _Atomic_FLOAT32 src_gains[MAX_NUM_INPUTS];         /**< Gains applied per source */


    /* End copied _binauraliser struct members. The following are unique to the _binauraliserNF struct */

    float b_dvf[MAX_NUM_INPUTS][NUM_EARS][2];                   /**< shelf IIR numerator coefficients for each input, left and right. */
    float a_dvf[MAX_NUM_INPUTS][NUM_EARS][2];                   /**< shelf IIR denominator coefficients for each input, left and right. */
    float dvfmags[MAX_NUM_INPUTS][NUM_EARS][HYBRID_BANDS];      /**< DVF filter frequency band magnitudes. */
    float dvfphases[MAX_NUM_INPUTS][NUM_EARS][HYBRID_BANDS];    /**< DVF filter frequency band phases. */

    /* misc. */
    _Atomic_FLOAT32 src_dists_m[MAX_NUM_INPUTS];  /**< Source distance,  meters. */
    float farfield_thresh_m;            /**< Distance considered to be far field (no near field filtering),  meters. */
    float farfield_headroom;            /**< Scale factor applied to farfield_thresh_m when resetting to the far field, and for UI range, meters. */
    float nearfield_limit_m;            /**< Minimum distance allowed for near-field filtering, from head _center_, meters, def. 0.15. */
    float head_radius;                  /**< Head radius, used calculate normalized source distance meters, def. 0.09096. */
    float head_radius_recip;            /**< Reciprocal of head radius. */
    float src_dirs_cur[MAX_NUM_INPUTS][2];    /**< Pointer to assign to the current HRTF directions being operated on (non/rotated directions switch). */

    /* flags/status */
    _Atomic_INT32 recalc_dvfCoeffFLAG[MAX_NUM_INPUTS]; /**< 1: re-calculate the DVF coefficients on change in distance, 0: do not. */

} binauraliserNF_data;


/* ========================================================================== */
/*                             Internal Functions                             */
/* ========================================================================== */

/**
 * Initialise the filterbank used by binauraliserNF.
 *
 * @note Call this function before binauraliser_initHRTFsAndGainTables()
 */
void binauraliserNF_initTFT(void* const hBin);

/**
 * Resets the source distances to the default far field distance.
 */
void binauraliserNF_resetSourceDistances(void* const hBin);


#ifdef __cplusplus
} /* extern "C" { */
#endif /* __cplusplus */

#endif /* __BINAURALISER_INTERNAL_NF_H_INCLUDED__ */
