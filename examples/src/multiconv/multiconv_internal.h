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
 * @file multiconv_internal.h
 * @brief A multi-channel convolver
 * @author Leo McCormack
 * @date 23.09.2019
 */

#ifndef __MULTICONV_INTERNAL_H_INCLUDED__
#define __MULTICONV_INTERNAL_H_INCLUDED__

#include "multiconv.h"     /* Include header for this example */
#include "saf.h"           /* Main include header for SAF */
#include "saf_externals.h" /* To also include SAF dependencies (cblas etc.) */

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */
    
/* ========================================================================== */
/*                            Internal Parameters                             */
/* ========================================================================== */

#define MIN_FRAME_SIZE ( 512 )   /**< Minimum framesize, in time-domain samples */
#define MAX_FRAME_SIZE ( 8192 )  /**< Maximum framesize, in time-domain samples */
    
/* ========================================================================== */
/*                                 Structures                                 */
/* ========================================================================== */

/**
 * Main structure for multiconv.
 */
typedef struct _multiconv
{
    /* FIFO buffers */
    int FIFO_idx;                /**< FIFO buffer index */
    float inFIFO[MAX_NUM_CHANNELS][MAX_FRAME_SIZE];  /**< Input FIFO buffer */
    float outFIFO[MAX_NUM_CHANNELS][MAX_FRAME_SIZE]; /**< Output FIFO buffer */

    /* Internal buffers */
    float** inputFrameTD;
    float** outputFrameTD;
    
    /* internal */
    void* hMultiConv;
    int hostBlockSize;
    int hostBlockSize_clamped; /**< Clamped between MIN and #MAX_FRAME_SIZE */
    float* filters;   /**< FLAT: nfilters x filter_length */
    int nfilters;
    int filter_length;
    int filter_fs;
    int host_fs;
    int reInitFilters;
    
    /* user parameters */
    int nChannels;
    int enablePartitionedConv;
    
} multiconv_data;


#ifdef __cplusplus
} /* extern "C" { */
#endif /* __cplusplus */

#endif /* __MULTICONV_INTERNAL_H_INCLUDED__ */
