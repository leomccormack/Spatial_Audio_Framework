/*
 * Copyright 2020 Leo McCormack
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
 * @file ambi_roomsim_internal.h
 * @brief A simple shoebox room Ambisonic encoder.
 *
 * @author Leo McCormack
 * @date 10.08.2020
 */

#ifndef __AMBI_ROOMSIM_INTERNAL_H_INCLUDED__
#define __AMBI_ROOMSIM_INTERNAL_H_INCLUDED__

#include <stdio.h>
#include <math.h>
#include <string.h>
#include "ambi_roomsim.h"
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
#define ROOM_SIM_MAX_NUM_SOURCES ( 8 )
#define MAX_NUM_SH_SIGNALS ( (MAX_SH_ORDER + 1)*(MAX_SH_ORDER + 1) ) /* (L+1)^2 */


/* ========================================================================== */
/*                                 Structures                                 */
/* ========================================================================== */

/**
 * Main structure for ambi_roomsim. Contains variables for audio buffers, internal
 * variables, user parameters
 */
typedef struct _ambi_roomsim
{
    /* Internals */
    float inputFrameTD[MAX_NUM_INPUTS][FRAME_SIZE];
    float prev_inputFrameTD[MAX_NUM_INPUTS][FRAME_SIZE];
    float tempFrame[MAX_NUM_SH_SIGNALS][FRAME_SIZE];
    float outputFrameTD[MAX_NUM_SH_SIGNALS][FRAME_SIZE];
    float fs;
    int recalc_SH_FLAG[MAX_NUM_INPUTS];
    float Y[MAX_NUM_SH_SIGNALS][MAX_NUM_INPUTS];
    float prev_Y[MAX_NUM_SH_SIGNALS][MAX_NUM_INPUTS];
    float interpolator[FRAME_SIZE];

    /* sd */
    void* hIms;
    int signalLength;
    int sh_order;
    int nBands;
    float abs_wall[5][6];
    float src_pos[3];
    float src2_pos[3];
    float src3_pos[3];
    float src4_pos[3];
    float rec_pos[3];

    long sourceIDs[4], receiverIDs[1];
    float** src_sigs, **rec_sh_outsigs;

    int reinit_room;
    
    /* user parameters */
    int nSources;
    int new_nSources;
    float src_positions[ROOM_SIM_MAX_NUM_SOURCES][3];
    CH_ORDER chOrdering;
    NORM_TYPES norm;
    SH_ORDERS order;
    int enablePostScaling;
    
} ambi_roomsim_data;
    

/* ========================================================================== */
/*                             Internal Functions                             */
/* ========================================================================== */


#ifdef __cplusplus
} /* extern "C" { */
#endif /* __cplusplus */

#endif /* __AMBI_ROOMSIM_INTERNAL_H_INCLUDED__ */
