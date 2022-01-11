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
 * @file ambi_enc_internal.h
 * @brief A basic Ambisonic encoder
 *
 * @author Leo McCormack
 * @date 07.10.2016
 * @license ISC
 */

#ifndef __AMBI_ENC_INTERNAL_H_INCLUDED__
#define __AMBI_ENC_INTERNAL_H_INCLUDED__

#include "ambi_enc.h"      /* Include header for this example */
#include "saf.h"           /* Main include header for SAF */
#include "saf_externals.h" /* To also include SAF dependencies (cblas etc.) */

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/* ========================================================================== */
/*                            Internal Parameters                             */
/* ========================================================================== */

#if !defined(AMBI_ENC_FRAME_SIZE)
# if defined(FRAME_SIZE) /* Use the global framesize if it is specified: */
#  define AMBI_ENC_FRAME_SIZE ( FRAME_SIZE )   /**< Framesize, in time-domain samples */
# else /* Otherwise, the default framesize for this example is: */
#  define AMBI_ENC_FRAME_SIZE ( 64 )           /**< Framesize, in time-domain samples */
# endif
#endif

/* ========================================================================== */
/*                                 Structures                                 */
/* ========================================================================== */

/**
 * Main structure for ambi_enc. Contains variables for audio buffers, internal
 * variables, user parameters
 */
typedef struct _ambi_enc
{
    /* Internal audio buffers */
    float inputFrameTD[MAX_NUM_INPUTS][AMBI_ENC_FRAME_SIZE];              /**< Input frame of signals */
    float prev_inputFrameTD[MAX_NUM_INPUTS][AMBI_ENC_FRAME_SIZE];         /**< Previous frame of signals */
    float tempFrame_fadeOut[MAX_NUM_SH_SIGNALS][AMBI_ENC_FRAME_SIZE];     /**< Temporary frame with linear interpolation (fade-out) applied */
    float tempFrame[MAX_NUM_SH_SIGNALS][AMBI_ENC_FRAME_SIZE];             /**< Temporary frame */
    float outputFrameTD_fadeIn[MAX_NUM_SH_SIGNALS][AMBI_ENC_FRAME_SIZE];  /**< Output frame of SH signals with linear interpolation (fade-in) applied */
    float outputFrameTD[MAX_NUM_SH_SIGNALS][AMBI_ENC_FRAME_SIZE];         /**< Output frame of SH signals */

    /* Internal variables */
    float fs;                                                    /**< Host sampling rate */
    int recalc_SH_FLAG[MAX_NUM_INPUTS];                          /**< Flags, 1: recalc SH weights, 0: do not */
    float Y[MAX_NUM_SH_SIGNALS][MAX_NUM_INPUTS];                 /**< SH weights */
    float prev_Y[MAX_NUM_SH_SIGNALS][MAX_NUM_INPUTS];            /**< Previous SH weights */
    float interpolator_fadeIn[AMBI_ENC_FRAME_SIZE];              /**< Linear Interpolator (fade-in) */
    float interpolator_fadeOut[AMBI_ENC_FRAME_SIZE];             /**< Linear Interpolator (fade-out) */
    int new_nSources;                                            /**< New number of input signals (current value will be replaced by this after next re-init) */
    
    /* user parameters */
    int nSources;                                                /**< Current number of input signals */
    float src_dirs_deg[MAX_NUM_INPUTS][2];                       /**< Source directions, in degrees */
    CH_ORDER chOrdering;                                         /**< Ambisonic channel order convention (see #CH_ORDER) */
    NORM_TYPES norm;                                             /**< Ambisonic normalisation convention (see #NORM_TYPES) */
    SH_ORDERS order;                                             /**< Current SH encoding order */
    int enablePostScaling;                                       /**< Flag 1: output signals scaled by 1/sqrt(nSources), 0: disabled */
    float src_gains[MAX_NUM_INPUTS];                             /**< Gains applied per source */

} ambi_enc_data;
    

/* ========================================================================== */
/*                             Internal Functions                             */
/* ========================================================================== */

/**
 * Returns the source directions for a specified source config preset.
 *
 * The function also returns the number of source in the configuration
 * Note: default uniformly distributed points are used to pad the
 * dirs_deg matrix up to the #MAX_NUM_INPUTS, if nCH is less than this. This can
 * help avoid scenarios of many sources being panned in the same direction, or
 * triangulations errors.
 *
 * @param[in]  preset   See #SOURCE_CONFIG_PRESETS enum
 * @param[out] dirs_deg Source directions, [azimuth elevation] convention, in
 *                      DEGREES;
 * @param[out] nCH      (&) number of source directions in the configuration
 */
void loadSourceConfigPreset(SOURCE_CONFIG_PRESETS preset,
                            float dirs_deg[MAX_NUM_INPUTS][2],
                            int* nCH);


#ifdef __cplusplus
} /* extern "C" { */
#endif /* __cplusplus */

#endif /* __AMBI_ENC_INTERNAL_H_INCLUDED__ */
