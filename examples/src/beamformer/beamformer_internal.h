/*
 * Copyright 2019 Leo McCormack
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
 * @file beamformer_internal.h
 * @brief Generates beamformers/virtual microphones in arbitrary directions
 *        with several different beam patterns to choose from
 *
 * @author Leo McCormack
 * @date 17.05.2019
 */

#ifndef __BEAMFORMER_INTERNAL_H_INCLUDED__
#define __BEAMFORMER_INTERNAL_H_INCLUDED__

#include "beamformer.h"    /* Include header for this example */
#include "saf.h"           /* Main include header for SAF */
#include "saf_externals.h" /* To also include SAF dependencies (cblas etc.) */

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/* ========================================================================== */
/*                            Internal Parameters                             */
/* ========================================================================== */

#ifndef FRAME_SIZE
# define FRAME_SIZE ( 128 )               /**< Framesize, in time-domain samples */
#endif
#define MAX_NUM_BEAMS ( MAX_NUM_OUTPUTS ) /**< Maximum permitted number of beams/output channels */

/* ========================================================================== */
/*                                 Structures                                 */
/* ========================================================================== */

/**
 * Main structure for beamformer. Contains variables for audio buffers,
 * beamforming weights, internal variables, flags, user parameters
 */
typedef struct _beamformer
{
    /* Internal audio buffers */
    float SHFrameTD[MAX_NUM_SH_SIGNALS][FRAME_SIZE];
    float prev_SHFrameTD[MAX_NUM_SH_SIGNALS][FRAME_SIZE];
    float tempFrame[MAX_NUM_BEAMS][FRAME_SIZE];
    float tempFrame_fadeOut[MAX_NUM_SH_SIGNALS][FRAME_SIZE];
    float outputFrameTD[MAX_NUM_BEAMS][FRAME_SIZE];
    float outputFrameTD_fadeIn[MAX_NUM_SH_SIGNALS][FRAME_SIZE];

    /* internal variables */
    int fs;
    float beamWeights[MAX_NUM_BEAMS][MAX_NUM_SH_SIGNALS];
    float prev_beamWeights[MAX_NUM_BEAMS][MAX_NUM_SH_SIGNALS];
    float interpolator_fadeIn[FRAME_SIZE];
    float interpolator_fadeOut[FRAME_SIZE];
    int recalc_beamWeights[MAX_NUM_BEAMS];   /**< 0: no init required, 1: init required */ 
    
    /* user parameters */
    int beamOrder;                           /**< beam order */
    int nBeams;                              /**< number of loudspeakers/virtual loudspeakers */
    float beam_dirs_deg[MAX_NUM_BEAMS][2];   /**< beam directions in degrees [azi, elev] */
    STATIC_BEAM_TYPES beamType;              /**< see #STATIC_BEAM_TYPES enum */
    CH_ORDER chOrdering;                     /**< Ambisonic channel order convention (see #CH_ORDER) */
    NORM_TYPES norm;                         /**< Ambisonic normalisation convention (see #NORM_TYPES) */
    
} beamformer_data;


/* ========================================================================== */
/*                             Internal Functions                             */
/* ========================================================================== */
 

#ifdef __cplusplus
} /* extern "C" { */
#endif /* __cplusplus */

#endif /* __BEAMFORMER_INTERNAL_H_INCLUDED__ */
