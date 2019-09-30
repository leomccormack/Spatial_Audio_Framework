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

/*
 * Filename: matrixconv_internal.h
 * -------------------------------
 * A matrix convolver
 *
 * Dependencies:
 *     saf_utilities
 * Author, date created:
 *     Leo McCormack, 30.09.2019
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
 
#define MAX_NUM_CHANNELS ( MATRIXCONV_MAX_NUM_CHANNELS )
#define MAX_NUM_CHANNELS_FOR_WAV ( 1024 )
#ifndef DEG2RAD
# define DEG2RAD(x) (x * PI / 180.0f)
#endif
#ifndef RAD2DEG
# define RAD2DEG(x) (x * 180.0f / PI)
#endif
    
    
/* ========================================================================== */
/*                                 Structures                                 */
/* ========================================================================== */

typedef struct _matrixconv
{
    float** inputFrameTD;
    float** outputFrameTD;
    
    /* internal */
    void* hMatrixConv;
    int hostBlockSize; 
    float* filters;   /* FLAT: nfilters x filter_length, & (nOutputChannels x nInputChannels) x filter_length */
    int nfilters;
    int filter_length;
    int filter_fs;
    int host_fs;
    int reInitFilters;
    int nOutputChannels; 
    
    /* user parameters */
    int nInputChannels;
    int enablePartitionedConv;
    
} matrixconv_data;
    
    
#ifdef __cplusplus
} /* extern "C" { */
#endif /* __cplusplus */

#endif /* __MATRIXCONV_INTERNAL_H_INCLUDED__ */
