/*
 * Copyright 2017-2018 Leo McCormack
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
 * @file: binauraliser_internal.h
 * @brief Convolves input audio (up to 64 channels) with interpolated HRTFs in
 *        the time-frequency domain.
 *
 * The HRTFs are interpolated by applying amplitude-preserving VBAP gains to the
 * HRTF magnitude responses and inter-aural time differences (ITDs)
 * individually, before being re-combined. The example also allows the user to
 * specify an external SOFA file for the convolution, and rotations of the
 * source directions to accomodate head-tracking.
 *
 * @author Leo McCormack
 * @date 25.09.2017
 * @license ISC
 */

#ifndef __BINAURALISER_INTERNAL_H_INCLUDED__
#define __BINAURALISER_INTERNAL_H_INCLUDED__

#include "binauraliser.h"  /* Include header for this example */
#include "saf.h"           /* Main include header for SAF */
#include "saf_externals.h" /* To also include SAF dependencies (cblas etc.) */

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/* ========================================================================== */
/*                            Internal Parameters                             */
/* ========================================================================== */

#if !defined(BINAURALISER_FRAME_SIZE)
# if defined(FRAME_SIZE) /* Use the global framesize if it is specified: */
#  define BINAURALISER_FRAME_SIZE ( FRAME_SIZE )          /**< Framesize, in time-domain samples */
# else /* Otherwise, the default framesize for this example is: */
#  define BINAURALISER_FRAME_SIZE ( 128 )                 /**< Framesize, in time-domain samples */
# endif
#endif
#define HOP_SIZE ( 128 )                                  /**< STFT hop size */
#define HYBRID_BANDS ( HOP_SIZE + 5 )                     /**< Number of frequency bands */
#define TIME_SLOTS ( BINAURALISER_FRAME_SIZE / HOP_SIZE ) /**< Number of STFT timeslots */

/* Checks: */
#if (BINAURALISER_FRAME_SIZE % HOP_SIZE != 0)
# error "BINAURALISER_FRAME_SIZE must be an integer multiple of HOP_SIZE"
#endif

/* ========================================================================== */
/*                                 Structures                                 */
/* ========================================================================== */

/**
 * Main structure for binauraliser. Contains variables for audio buffers,
 * afSTFT, HRTFs, internal variables, flags, user parameters.
 * Note: if this is modified, identically modify _binauraliserNF struct.
 */
typedef struct _binauraliser
{
    /* audio buffers */
    float** inputFrameTD;            /**< time-domain input frame; #MAX_NUM_INPUTS x #BINAURALISER_FRAME_SIZE */
    float** outframeTD;              /**< time-domain output frame; #NUM_EARS x #BINAURALISER_FRAME_SIZE */
    float_complex*** inputframeTF;   /**< time-frequency domain input frame; #HYBRID_BANDS x #MAX_NUM_INPUTS x #TIME_SLOTS */
    float_complex*** outputframeTF;  /**< time-frequency domain input frame; #HYBRID_BANDS x #NUM_EARS x #TIME_SLOTS */
    int fs;                          /**< Host sampling rate, in Hz */
    float freqVector[HYBRID_BANDS];  /**< Frequency vector (filterbank centre frequencies) */
    void* hSTFT;                     /**< afSTFT handle */
    
    /* sofa file info */
    char* sofa_filepath;             /**< absolute/relevative file path for a sofa file */
    float* hrirs;                    /**< time domain HRIRs; FLAT: N_hrir_dirs x #NUM_EARS x hrir_len */
    float* hrir_dirs_deg;            /**< directions of the HRIRs in degrees [azi elev]; FLAT: N_hrir_dirs x 2 */
    int N_hrir_dirs;                 /**< number of HRIR directions in the current sofa file */
    int hrir_loaded_len;             /**< length of the loaded HRIRs, in samples */
    int hrir_runtime_len;            /**< length of the HRIRs being used for processing (after any resampling), in samples */
    int hrir_loaded_fs;              /**< sampling rate of the loaded HRIRs  */
    int hrir_runtime_fs;             /**< sampling rate of the HRIRs being used for processing (after any resampling) */
    float* weights;                  /**< Integration weights for the HRIR measurement grid */
    
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
    CODEC_STATUS codecStatus;        /**< see #CODEC_STATUS */
    float progressBar0_1;            /**< Current (re)initialisation progress, between [0..1] */
    char* progressBarText;           /**< Current (re)initialisation step, string */
    PROC_STATUS procStatus;          /**< see #PROC_STATUS */
    int recalc_hrtf_interpFLAG[MAX_NUM_INPUTS]; /**< 1: re-calculate/interpolate the HRTF, 0: do not */
    int reInitHRTFsAndGainTables;    /**< 1: reinitialise the HRTFs and interpolation tables, 0: do not */
    int recalc_M_rotFLAG;            /**< 1: re-calculate the rotation matrix, 0: do not */
    
    /* misc. */
    float src_dirs_rot_deg[MAX_NUM_INPUTS][2]; /**< Intermediate rotated source directions, in degrees */
    float src_dirs_rot_xyz[MAX_NUM_INPUTS][3]; /**< Intermediate rotated source directions, as unit-length Cartesian coordinates */
    float src_dirs_xyz[MAX_NUM_INPUTS][3];     /**< Intermediate source directions, as unit-length Cartesian coordinates  */
    int nTriangles;                            /**< Number of triangles in the convex hull of the spherical arrangement of HRIR directions/points */
    int new_nSources;                          /**< New number of input/source signals (current value will be replaced by this after next re-init) */

    /* user parameters */
    int nSources;                            /**< Current number of input/source signals */
    float src_dirs_deg[MAX_NUM_INPUTS][2];   /**< Current source/panning directions, in degrees */
    INTERP_MODES interpMode;                 /**< see #INTERP_MODES */
    int useDefaultHRIRsFLAG;                 /**< 1: use default HRIRs in database, 0: use those from SOFA file */
    int enableHRIRsDiffuseEQ;                /**< flag to diffuse-field equalisation to the currently loaded HRTFs */
    int enableRotation;                      /**< 1: enable rotation, 0: disable */
    float yaw;                               /**< yaw (Euler) rotation angle, in degrees */
    float roll;                              /**< roll (Euler) rotation angle, in degrees */
    float pitch;                             /**< pitch (Euler) rotation angle, in degrees */
    int bFlipYaw;                            /**< flag to flip the sign of the yaw rotation angle */
    int bFlipPitch;                          /**< flag to flip the sign of the pitch rotation angle */
    int bFlipRoll;                           /**< flag to flip the sign of the roll rotation angle */
    int useRollPitchYawFlag;                 /**< rotation order flag, 1: r-p-y, 0: y-p-r */
    float src_gains[MAX_NUM_INPUTS];         /**< Gains applied per source */

} binauraliser_data;


/* ========================================================================== */
/*                             Internal Functions                             */
/* ========================================================================== */

/** Sets codec status (see #CODEC_STATUS enum) */
void binauraliser_setCodecStatus(void* const hBin,
                                 CODEC_STATUS newStatus);

/**
 * Interpolates between (up to) 3 HRTFs via amplitude-normalised VBAP gains.
 *
 * The HRTF magnitude responses and HRIR ITDs are interpolated seperately before
 * re-introducing the phase.
 *
 * @param[in]  hBin          binauraliser handle
 * @param[in]  mode          see #INTERP_MODES 
 * @param[in]  azimuth_deg   Source azimuth in DEGREES
 * @param[in]  elevation_deg Source elevation in DEGREES
 * @param[out] h_intrp       Interpolated HRTF
 */
void binauraliser_interpHRTFs(void* const hBin,
                              INTERP_MODES mode,
                              float azimuth_deg,
                              float elevation_deg,
                              float_complex h_intrp[HYBRID_BANDS][NUM_EARS]);

/**
 * Initialise the HRTFs: either loading the default set or loading from a SOFA
 * file; and then generate a VBAP gain table for interpolation.
 *
 * @note Call binauraliser_initTFT() (if needed) before calling this function
 */
void binauraliser_initHRTFsAndGainTables(void* const hBin);

/**
 * Initialise the filterbank used by binauraliser.
 *
 * @note Call this function before binauraliser_initHRTFsAndGainTables()
 */
void binauraliser_initTFT(void* const hBin);

/**
 * Returns the source directions for a specified source config preset.
 *
 * The function also returns the number of source in the configuration
 * Note: default uniformly distributed points are used to pad the
 * dirs_deg matrix up to the #MAX_NUM_INPUTS, if nCH is less than
 * this. This can help avoid scenarios of many sources being panned in the same
 * direction, or triangulations errors.
 *
 * @param[in]  preset   See #SOURCE_CONFIG_PRESETS enum.
 * @param[out] dirs_deg Source directions, [azimuth elevation] convention, in
 *                      DEGREES;
 * @param[out] newNCH   (&) new number of channels
 * @param[out] nDims    (&) estimate of the number of dimensions (2 or 3)
 */
void binauraliser_loadPreset(SOURCE_CONFIG_PRESETS preset,
                             float dirs_deg[MAX_NUM_INPUTS][2],
                             int* newNCH,
                             int* nDims);


#ifdef __cplusplus
} /* extern "C" { */
#endif /* __cplusplus */

#endif /* __BINAURALISER_INTERNAL_H_INCLUDED__ */
