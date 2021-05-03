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
 */

#ifndef __BINAURALISER_INTERNAL_H_INCLUDED__
#define __BINAURALISER_INTERNAL_H_INCLUDED__

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "binauraliser.h"
#include "saf.h"
#include "saf_externals.h" /* to also include saf dependencies (cblas etc.) */

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/* ========================================================================== */
/*                            Internal Parameters                             */
/* ========================================================================== */

#ifndef FRAME_SIZE
# define FRAME_SIZE ( 128 )
#endif
#define HOP_SIZE ( 128 )                                    /* STFT hop size = nBands */
#define HYBRID_BANDS ( HOP_SIZE + 5 )                       /* hybrid mode incurs an additional 5 bands  */
#define TIME_SLOTS ( FRAME_SIZE / HOP_SIZE )                /* 4/8/16 */
#ifndef DEG2RAD
# define DEG2RAD(x) (x * M_PI / 180.0f)
#endif
#ifndef RAD2DEG
# define RAD2DEG(x) (x * 180.0f / M_PI)
#endif
#if (FRAME_SIZE % HOP_SIZE != 0)
# error "FRAME_SIZE must be an integer multiple of HOP_SIZE"
#endif

/* ========================================================================== */
/*                                 Structures                                 */
/* ========================================================================== */

/**
 * Main structure for binauraliser. Contains variables for audio buffers,
 * afSTFT, HRTFs, internal variables, flags, user parameters
 */
typedef struct _binauraliser
{
    /* audio buffers */
    float** inputFrameTD;
    float** outframeTD;
    float_complex*** inputframeTF;
    float_complex*** outputframeTF;
    int fs;
    float freqVector[HYBRID_BANDS]; 
    void* hSTFT;
    
    /* sofa file info */
    char* sofa_filepath; 
    float* hrirs;
    float* hrir_dirs_deg;
    int N_hrir_dirs;
    int hrir_len;
    int hrir_fs;
    float* weights;
    
    /* vbap gain table */
    int hrtf_vbapTableRes[2];
    int N_hrtf_vbap_gtable;
    int* hrtf_vbap_gtableIdx;        /**< N_hrtf_vbap_gtable x 3 */
    float* hrtf_vbap_gtableComp;     /**< N_hrtf_vbap_gtable x 3 */
    
    /* hrir filterbank coefficients */
    float* itds_s;                   /**< interaural-time differences for each HRIR (in seconds); nBands x 1 */
    float_complex* hrtf_fb;          /**< hrtf filterbank coefficients; nBands x nCH x N_hrirs */
    float* hrtf_fb_mag;              /**< magnitudes of the hrtf filterbank coefficients; nBands x nCH x N_hrirs */
    float_complex hrtf_interp[MAX_NUM_INPUTS][HYBRID_BANDS][NUM_EARS];
    
    /* flags/status */
    CODEC_STATUS codecStatus;
    float progressBar0_1;
    char* progressBarText;
    PROC_STATUS procStatus;
    int recalc_hrtf_interpFLAG[MAX_NUM_INPUTS];
    int reInitHRTFsAndGainTables;
    int recalc_M_rotFLAG;
    
    /* misc. */
    float src_dirs_rot_deg[MAX_NUM_INPUTS][2];
    float src_dirs_rot_xyz[MAX_NUM_INPUTS][3];
    float src_dirs_xyz[MAX_NUM_INPUTS][3]; 
    int nTriangles;
    int input_nDims;  
    int output_nDims;
    
    /* user parameters */
    int nSources;
    int new_nSources;
    float src_dirs_deg[MAX_NUM_INPUTS][2];
    INTERP_MODES interpMode;
    int useDefaultHRIRsFLAG;                 /**< 1: use default HRIRs in database, 0: use those from SOFA file */
    int enableHRIRsPreProc;                  /**< flag to apply pre-processing to the currently loaded HRTFs */
    int enableRotation;
    float yaw, roll, pitch;                  /**< rotation angles in degrees */
    int bFlipYaw, bFlipPitch, bFlipRoll;     /**< flag to flip the sign of the individual rotation angles */
    int useRollPitchYawFlag;                 /**< rotation order flag, 1: r-p-y, 0: y-p-r */
    
} binauraliser_data;


/* ========================================================================== */
/*                             Internal Functions                             */
/* ========================================================================== */

/**
 * Sets codec status (see #CODEC_STATUS enum)
 */
void binauraliser_setCodecStatus(void* const hBin,
                                 CODEC_STATUS newStatus);

/**
 * Interpolates between (up to) 3 HRTFs via amplitude-normalised VBAP gains.
 *
 * The HRTF magnitude responses and HRIR ITDs are interpolated seperately before
 * re-introducing the phase.
 *
 * @param[in]  hBin          binauraliser handle
 * @param[in]  azimuth_deg   Source azimuth in DEGREES
 * @param[in]  elevation_deg Source elevation in DEGREES
 * @param[out] h_intrp       Interpolated HRTF
 */
void binauraliser_interpHRTFs(void* const hBin,
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
