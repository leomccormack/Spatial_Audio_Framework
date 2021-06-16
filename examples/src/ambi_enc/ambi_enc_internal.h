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

#ifndef FRAME_SIZE
# define FRAME_SIZE ( 64 ) 
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
    float inputFrameTD[MAX_NUM_INPUTS][FRAME_SIZE];
    float prev_inputFrameTD[MAX_NUM_INPUTS][FRAME_SIZE];
    float tempFrame_fadeOut[MAX_NUM_SH_SIGNALS][FRAME_SIZE];
    float tempFrame[MAX_NUM_SH_SIGNALS][FRAME_SIZE];
    float outputFrameTD_fadeIn[MAX_NUM_SH_SIGNALS][FRAME_SIZE];
    float outputFrameTD[MAX_NUM_SH_SIGNALS][FRAME_SIZE];

    /* Internal variables */
    float fs;
    int recalc_SH_FLAG[MAX_NUM_INPUTS];
    float Y[MAX_NUM_SH_SIGNALS][MAX_NUM_INPUTS];
    float prev_Y[MAX_NUM_SH_SIGNALS][MAX_NUM_INPUTS];
    float interpolator_fadeIn[FRAME_SIZE];
    float interpolator_fadeOut[FRAME_SIZE];
    int new_nSources;
    
    /* user parameters */
    int nSources;
    float src_dirs_deg[MAX_NUM_INPUTS][2];
    CH_ORDER chOrdering;                 /**< Ambisonic channel order convention (see #CH_ORDER) */
    NORM_TYPES norm;                     /**< Ambisonic normalisation convention (see #NORM_TYPES) */
    SH_ORDERS order;
    int enablePostScaling;
    
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
