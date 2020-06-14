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

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include "multiconv.h"
#include "saf.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */
    
/* ========================================================================== */
/*                            Internal Parameters                             */
/* ========================================================================== */

#define MIN_FRAME_SIZE ( 512 )
#define MAX_FRAME_SIZE ( 8192 ) 
#ifndef DEG2RAD
# define DEG2RAD(x) (x * PI / 180.0f)
#endif
#ifndef RAD2DEG
# define RAD2DEG(x) (x * 180.0f / PI)
#endif
    
    
/* ========================================================================== */
/*                                 Structures                                 */
/* ========================================================================== */

/**
 * Main structure for multiconv.
 */
typedef struct _multiconv
{
    /* FIFO buffers */
    int FIFO_idx;
    float inFIFO[MAX_NUM_CHANNELS][MAX_FRAME_SIZE];
    float outFIFO[MAX_NUM_CHANNELS][MAX_FRAME_SIZE];

    /* Internal buffers */
    float** inputFrameTD;
    float** outputFrameTD;
    
    /* internal */
    void* hMultiConv;
    int hostBlockSize;
    int hostBlockSize_clamped; /**< Clamped between MIN and MAX_FRAME_SIZE */
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
