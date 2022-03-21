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
 * @file tvconv_internal.h
 * @brief A time-varying multi-channel convolver
 * @author Rapolas Daugintis
 * @date 13.07.2021
 */

#ifndef __TVCONV_INTERNAL_H_INCLUDED__
#define __TVCONV_INTERNAL_H_INCLUDED__

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include "tvconv.h"
#include "saf.h"
#include "saf_externals.h" /* to also include saf dependencies (cblas etc.) */

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/* ========================================================================== */
/*                            Internal Parameters                             */
/* ========================================================================== */

#define MIN_FRAME_SIZE ( 512 )
#define MAX_FRAME_SIZE ( 8192 )
#define NUM_DIMENSIONS ( 3 )

/* ========================================================================== */
/*                                 Structures                                 */
/* ========================================================================== */

/** Structure for a vector */
typedef float vectorND[NUM_DIMENSIONS];

/** Main structure for tvconv  */
typedef struct _tvconv
{
    /* FIFO buffers */
    int FIFO_idx;
    float inFIFO[MAX_NUM_CHANNELS][MAX_FRAME_SIZE];
    float outFIFO[MAX_NUM_CHANNELS][MAX_FRAME_SIZE];

    /* Internal buffers */
    float** inputFrameTD;
    float** outputFrameTD;
    
    /* internal */
    void* hTVConv;     /**< saf_TVConv handle */
    int hostBlockSize;     /**< current host block size */
    int hostBlockSize_clamped; /**< Clamped between MIN and #MAX_FRAME_SIZE */
    int host_fs;           /**< current samplerate of the host */
    int reInitFilters;     /**< FLAG: 0: do not reinit, 1: reinit, 2: reinit in progress */
    int nOutputChannels;   /**< number of output channels (same as the number of channels in the loaded wav) */
    
    int ir_fs;
    float** irs;   /**< npositionsx x (FLAT: nfilters x filter_length) */
    int nIrChannels; /**< number of filters per position */
    int ir_length;
    
    /* positions */
    vectorND* listenerPositions;       /**< The listener positions; nListenerPositions x 3 */
    int nListenerPositions;
    vectorND minDimensions;            /**< Minimum values across all dimensions */
    vectorND maxDimensions;            /**< Maximum values across all dimensions */
    int position_idx;
    vectorND sourcePosition;
    
    /* flags/status */
    CODEC_STATUS codecStatus;
    float progressBar0_1;
    char* progressBarText;
    PROC_STATUS procStatus;
    
    /* user parameters */
    int nInputChannels;        /**< number of input channels */
    vectorND targetPosition;    
    char* sofa_filepath;

} tvconv_data;

/* ========================================================================== */
/*                             Internal Functions                             */
/* ========================================================================== */

/** Sets codec status (see #CODEC_STATUS enum) */
void tvconv_setCodecStatus(void* const hTVCnv,
                                 CODEC_STATUS newStatus);
/** Finds the index holding the nearest neigbour to the selected position */
void tvconv_findNearestNeigbour(void* const hTVCnv);

/**
 * Sets the smallest and the highest position of each dimension from the list of
 * positions
 */
void tvconv_setMinMaxDimensions(void* const hTVCnv);


#ifdef __cplusplus
} /* extern "C" { */
#endif /* __cplusplus */

#endif /* __TVCONV_INTERNAL_H_INCLUDED__ */
