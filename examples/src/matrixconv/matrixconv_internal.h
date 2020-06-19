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
 * @file matrixconv_internal.h
 * @brief A standard matrix convolver
 * @author Leo McCormack
 * @date 30.09.2019
 */

#ifndef __MATRIXCONV_INTERNAL_H_INCLUDED__
#define __MATRIXCONV_INTERNAL_H_INCLUDED__

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include "matrixconv.h"
#include "saf.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */
    
/* ========================================================================== */
/*                            Internal Parameters                             */
/* ========================================================================== */

#define MIN_FRAME_SIZE ( 512 )
#define MAX_FRAME_SIZE ( 8192 ) 
#define MAX_NUM_CHANNELS_FOR_WAV ( 1024 )
    
    
/* ========================================================================== */
/*                                 Structures                                 */
/* ========================================================================== */

/**
 * Main structure for matrixconv.
 */
typedef struct _matrixconv
{
    /* FIFO buffers */
    int FIFO_idx;
    float inFIFO[MAX_NUM_CHANNELS][MAX_FRAME_SIZE];
    float outFIFO[MAX_NUM_CHANNELS][MAX_FRAME_SIZE];

    /* input/output buffers */
    float** inputFrameTD;
    float** outputFrameTD;
    
    /* internal */
    void* hMatrixConv;     /**< saf_matrixConv handle */
    int hostBlockSize;     /**< current host block size */
    int hostBlockSize_clamped; /**< Clamped between MIN and #MAX_FRAME_SIZE */
    float* filters;        /**< the matrix of filters; FLAT: nOutputChannels x nInputChannels x filter_length */
    int nfilters;          /**< the number of filters (nOutputChannels x nInputChannels) */
    int input_wav_length;  /**< length of the wav files loaded in samples (inputs are concatenated) */
    int filter_length;     /**< length of the filters (i.e. input_wav_length/nInputChannels) */
    int filter_fs;         /**< current samplerate of the filters */
    int host_fs;           /**< current samplerate of the host */
    int reInitFilters;     /**< FLAG: 0: do not reinit, 1: reinit, 2: reinit in progress */
    int nOutputChannels;   /**< number of output channels (same as the number of channels in the loaded wav) */
    
    /* user parameters */
    int nInputChannels;        /**< number of input channels */
    int enablePartitionedConv; /**< 0: disabled, 1: enabled */
    
} matrixconv_data;
    
    
#ifdef __cplusplus
} /* extern "C" { */
#endif /* __cplusplus */

#endif /* __MATRIXCONV_INTERNAL_H_INCLUDED__ */
