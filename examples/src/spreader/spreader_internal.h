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
 * @file: spreader.h
 * @brief spreader
 * @author Leo McCormack
 * @date 07.04.2021
 */

#ifndef __SPREADER_INTERNAL_H_INCLUDED__
#define __SPREADER_INTERNAL_H_INCLUDED__

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "spreader.h"
#include "saf.h"
#include "saf_externals.h" /* to also include saf dependencies (cblas etc.) */

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/* ========================================================================== */
/*                            Internal Parameters                             */
/* ========================================================================== */

#ifndef FRAME_SIZE
# define FRAME_SIZE ( 512 )
#endif
#define MAX_SPREAD_FREQ ( 16e3f )
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
 * Main structure for spreader. Contains variables for audio buffers,
 * afSTFT, HRTFs, internal variables, flags, user parameters
 */
typedef struct _spreader
{
    /* audio buffers and time-frequency transform */
    float** inputFrameTD;
    float** outframeTD;
    float_complex*** inputframeTF;
    float_complex*** protoframeTF;
    float_complex*** decorframeTF;
    float_complex*** spreadframeTF;
    float_complex*** outputframeTF;
    int fs;
    float freqVector[HYBRID_BANDS]; 
    void* hSTFT;

    /* Internal */
    int Q;                             /**< Number of channels in the target; e.g. 2 for binaural */
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
    float_complex* inFrame_t;          /**< Temporary input frame; Q x 1 */
    float interpolatorFadeIn[TIME_SLOTS];  /**< Interpolator */
    float interpolatorFadeOut[TIME_SLOTS]; /**< Interpolator */

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
    CODEC_STATUS codecStatus;
    float progressBar0_1;
    char* progressBarText;
    PROC_STATUS procStatus;
    int new_nSources;
    SPREADER_PROC_MODES new_procMode;

    /* user parameters */
    SPREADER_PROC_MODES procMode;     /**< See #SPREADER_PROC_MODES */
    char* sofa_filepath;
    int nSources;
    float src_spread[SPREADER_MAX_NUM_SOURCES];
    float src_dirs_deg[SPREADER_MAX_NUM_SOURCES][2]; 
    int useDefaultHRIRsFLAG;          /**< 1: use default HRIRs in database, 0: use the measurements from SOFA file (can be anything, not just HRTFs) */
    float covAvgCoeff;

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
