/*
 * Copyright 2021 Leo McCormack
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
 * @file: spreader_internal.h
 * @brief An arbitrary array panner (HRIRs, microphone array IRs, etc.) with
 *        coherent and incoherent spreading modes, as described in [1].
 *
 * @see [1] McCormack, L. Politis, A., and Pulkki, V., 2021, October. Rendering
 *          of source spread for arbitrary playback setups based on spatial
 *          covariance matching. In 2021 IEEE Workshop on Applications of Signal
 *          Processing to Audio and Acoustics (WASPAA). IEEE
 *
 * @author Leo McCormack
 * @date 07.04.2021
 * @license ISC
 */

#ifndef __SPREADER_INTERNAL_H_INCLUDED__
#define __SPREADER_INTERNAL_H_INCLUDED__

#include "spreader.h"      /* Include header for this example */
#include "saf.h"           /* Main include header for SAF */
#include "saf_externals.h" /* To also include SAF dependencies (cblas etc.) */

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/* ========================================================================== */
/*                            Internal Parameters                             */
/* ========================================================================== */

#if !defined(SPREADER_FRAME_SIZE)
# if defined(FRAME_SIZE) /* Use the global framesize if it is specified: */
#  define SPREADER_FRAME_SIZE ( FRAME_SIZE )          /**< Framesize, in time-domain samples */
# else /* Otherwise, the default framesize for this example is: */
#  define SPREADER_FRAME_SIZE ( 512 )                 /**< Framesize, in time-domain samples */
# endif
#endif
#define MAX_SPREAD_FREQ ( 16e3f )                     /**< Maximum spread frequency, above which no spreading occurs */
#define HOP_SIZE ( 128 )                              /**< STFT hop size */
#define HYBRID_BANDS ( HOP_SIZE + 5 )                 /**< Number of frequency bands */
#define TIME_SLOTS ( SPREADER_FRAME_SIZE / HOP_SIZE ) /**< Number of STFT timeslots */

/* Checks: */
#if (SPREADER_FRAME_SIZE % HOP_SIZE != 0)
# error "SPREADER_FRAME_SIZE must be an integer multiple of HOP_SIZE"
#endif

/* ========================================================================== */
/*                                 Structures                                 */
/* ========================================================================== */

/**
 * Main structure for spreader. Contains variables for audio buffers,
 * afSTFT, HRTFs, internal variables, flags, user parameters
 */
typedef struct _spreader
{
    /* audio buffers and time-frequency transform */
    float** inputFrameTD;              /**< time-domain input frame; #MAX_NUM_INPUTS x #SPREADER_FRAME_SIZE */
    float** outframeTD;                /**< time-domain output frame; #MAX_NUM_OUTPUTS x #SPREADER_FRAME_SIZE */
    float_complex*** inputframeTF;     /**< time-frequency domain input frame; #HYBRID_BANDS x #MAX_NUM_INPUTS x #TIME_SLOTS */
    float_complex*** protoframeTF;     /**< time-frequency domain prototype frame; #HYBRID_BANDS x #MAX_NUM_OUTPUTS x #TIME_SLOTS */
    float_complex*** decorframeTF;     /**< time-frequency domain decorrelated frame; #HYBRID_BANDS x #MAX_NUM_OUTPUTS x #TIME_SLOTS */
    float_complex*** spreadframeTF;    /**< time-frequency domain spread frame; #HYBRID_BANDS x #MAX_NUM_OUTPUTS x #TIME_SLOTS */
    float_complex*** outputframeTF;    /**< time-frequency domain output frame; #HYBRID_BANDS x #MAX_NUM_OUTPUTS x #TIME_SLOTS */
    int fs;                            /**< Host sampling rate, in Hz */
    float freqVector[HYBRID_BANDS];    /**< Frequency vector (filterbank centre frequencies) */
    void* hSTFT;                       /**< afSTFT handle */

    /* Internal */
    int Q;                             /**< Number of channels in the target playback setup; for example: 2 for binaural */
    int nGrid;                         /**< Number of directions/measurements/HRTFs etc. */
    int h_len;                         /**< Length of time-domain filters, in samples */
    float h_fs;                        /**< Sample rate used to measure the filters */
    float* h_grid;                     /**< FLAT: nGrid x Q x h_len */
    float_complex* H_grid;             /**< FLAT: HYBRID_BANDS x Q x nGrid */
    float_complex** HHH[HYBRID_BANDS]; /**< Pre-computed array outer-products; HYBRID_BANDS x nGrid x FLAT: (Q x Q) */
    float* grid_dirs_deg;              /**< Grid directions, in degrees; FLAT: nGrid x 2 */
    float* grid_dirs_xyz;              /**< Grid directions as unit-length Cartesian coordinates; FLAT: nGrid x 3 */
    float* weights;                    /**< Integration weights; nGrid x 1 */
    void* hDecor[SPREADER_MAX_NUM_SOURCES]; /**< handles for decorrelators */
    float* angles;                     /**< angles; nGrid x 1 */
    float_complex** Cproto[SPREADER_MAX_NUM_SOURCES]; /**< Current prototype covariance matrices; HYBRID_BANDS x FLAT:(Q x Q) */
    float_complex** Cy[SPREADER_MAX_NUM_SOURCES];     /**< Target covariance matrices; HYBRID_BANDS x FLAT:(Q x Q) */
    float_complex** prev_M[SPREADER_MAX_NUM_SOURCES]; /**< previous mixing matrices; HYBRID_BANDS x FLAT:(Q x Q) */
    float** prev_Mr[SPREADER_MAX_NUM_SOURCES];        /**< previous residual mixing matrices; HYBRID_BANDS x FLAT:(Q x Q) */
    float_complex** new_M;             /**< mixing matrices; HYBRID_BANDS x FLAT:(Q x Q) */
    float** new_Mr;                    /**< residual mixing matrices; HYBRID_BANDS x FLAT:(Q x Q) */
    float_complex* interp_M;           /**< Interpolated mixing matrix; FLAT:(Q x Q) */
    float* interp_Mr;                  /**< Interpolated residual mixing matrix; FLAT:(Q x Q) */
    float_complex* interp_Mr_cmplx;    /**< Complex variant of interp_Mr */ 
    float interpolatorFadeIn[TIME_SLOTS];  /**< Linear Interpolator - Fade in */
    float interpolatorFadeOut[TIME_SLOTS]; /**< Linear Interpolator - Fade out */

    /* For visualisation */
    int* dirActive[SPREADER_MAX_NUM_SOURCES]; /**< 1: IR direction currently used for spreading, 0: not */

    /* Optimal mixing solution */
    void* hCdf;                        /**< covariance domain framework handle */
    void* hCdf_res;                    /**< covariance domain framework handle for the residual */
    float* Qmix;                       /**< Identity; FLAT: Q x Q */
    float_complex* Qmix_cmplx;         /**< Identity; FLAT: Q x Q */
    float* Cr;                         /**< Residual covariance; FLAT: Q x Q */
    float_complex* Cr_cmplx;           /**< Residual covariance; FLAT: Q x Q */
 
    /* flags/status */
    CODEC_STATUS codecStatus;          /**< see #CODEC_STATUS */
    float progressBar0_1;              /**< Current (re)initialisation progress, between [0..1] */
    char* progressBarText;             /**< Current (re)initialisation step, string */
    PROC_STATUS procStatus;            /**< see #PROC_STATUS */
    int new_nSources;                  /**< New number of input signals (current value will be replaced by this after next re-init) */
    SPREADER_PROC_MODES new_procMode;  /**< See #SPREADER_PROC_MODES (current value will be replaced by this after next re-init) */

    /* user parameters */
    SPREADER_PROC_MODES procMode;      /**< See #SPREADER_PROC_MODES */
    char* sofa_filepath;               /**< SOFA file path */
    int nSources;                      /**< Current number of input signals */
    float src_spread[SPREADER_MAX_NUM_SOURCES];      /**< Source spreading, in degrees */
    float src_dirs_deg[SPREADER_MAX_NUM_SOURCES][2]; /**< Source directions, in degrees */
    int useDefaultHRIRsFLAG;           /**< 1: use default HRIRs in database, 0: use the measurements from SOFA file (can be anything, not just HRTFs) */
    float covAvgCoeff;                 /**< Covariance matrix averaging coefficient, [0..1] */

} spreader_data;


/* ========================================================================== */
/*                             Internal Functions                             */
/* ========================================================================== */

/** Sets codec status (see #CODEC_STATUS enum) */
void spreader_setCodecStatus(void* const hSpr,
                             CODEC_STATUS newStatus);


#ifdef __cplusplus
} /* extern "C" { */
#endif /* __cplusplus */

#endif /* __SPREADER_INTERNAL_H_INCLUDED__ */
