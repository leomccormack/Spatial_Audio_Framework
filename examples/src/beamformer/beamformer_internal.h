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

#include <stdio.h>
#include <math.h>
#include <string.h>
#include "beamformer.h"
#include "saf.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/* ========================================================================== */
/*                            Internal Parameters                             */
/* ========================================================================== */

#ifndef FRAME_SIZE
# define FRAME_SIZE ( 128 ) 
#endif
#define HOP_SIZE ( 128 )                         /* STFT hop size = nBands */
#define HYBRID_BANDS ( HOP_SIZE + 5 )            /* hybrid mode incurs an additional 5 bands  */
#define TIME_SLOTS ( FRAME_SIZE / HOP_SIZE )     /* 4/8/16 */
#define MAX_NUM_BEAMS ( MAX_NUM_OUTPUTS ) /* Maximum permitted channels for the VST standard */


/* ========================================================================== */
/*                                 Structures                                 */
/* ========================================================================== */

/**
 * Main structure for beamformer. Contains variables for audio buffers, afSTFT,
 * beamforming weights, internal variables, flags, user parameters
 */
typedef struct _beamformer
{
    /* FIFO buffers */
    int FIFO_idx;
    float inFIFO[MAX_NUM_SH_SIGNALS][FRAME_SIZE];
    float outFIFO[MAX_NUM_BEAMS][FRAME_SIZE];

    /* audio buffers + afSTFT time-frequency transform handle */
    float SHFrameTD[MAX_NUM_SH_SIGNALS][FRAME_SIZE];
    float prev_SHFrameTD[MAX_NUM_SH_SIGNALS][FRAME_SIZE];
    float tempFrame[MAX_NUM_BEAMS][FRAME_SIZE];
    float outputFrameTD[MAX_NUM_BEAMS][FRAME_SIZE];
    int fs;
    
    /* internal variables */ 
    float beamWeights[MAX_NUM_BEAMS][MAX_NUM_SH_SIGNALS];
    float prev_beamWeights[MAX_NUM_BEAMS][MAX_NUM_SH_SIGNALS];
    float interpolator[FRAME_SIZE];
    
    /* flags */
    int recalc_beamWeights[MAX_NUM_BEAMS];   /**< 0: no init required, 1: init required */ 
    
    /* user parameters */
    int beamOrder;                           /**< beam order */
    int nBeams;                              /**< number of loudspeakers/virtual loudspeakers */
    float beam_dirs_deg[MAX_NUM_BEAMS][2];   /**< beam directions in degrees [azi, elev] */
    STATIC_BEAM_TYPES beamType;              /**< see #_STATIC_BEAM_TYPES enum */
    CH_ORDER chOrdering;                     /**< only ACN is supported */
    NORM_TYPES norm;                         /**< N3D or SN3D */
    
} beamformer_data;


/* ========================================================================== */
/*                             Internal Functions                             */
/* ========================================================================== */
 

#ifdef __cplusplus
} /* extern "C" { */
#endif /* __cplusplus */

#endif /* __BEAMFORMER_INTERNAL_H_INCLUDED__ */
