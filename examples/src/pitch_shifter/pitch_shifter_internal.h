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
 * @file pitch_shifter_internal.h
 * @brief A very basic multichannel pitch shifter
 *
 * @author Leo McCormack
 * @date 05.05.2020
 */

#ifndef __PITCH_SHIFTER_INTERNAL_H_INCLUDED__
#define __PITCH_SHIFTER_INTERNAL_H_INCLUDED__

#include <stdio.h>
#include <math.h>
#include <string.h>
#include "pitch_shifter.h"
#include "saf.h"
#include "saf_externals.h" /* to also include saf dependencies (cblas etc.) */

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/* ========================================================================== */
/*                            Internal Parameters                             */
/* ========================================================================== */

#ifndef FRAME_SIZE
# define FRAME_SIZE ( 128 ) 
#endif


/* ========================================================================== */
/*                                 Structures                                 */
/* ========================================================================== */

/**
 * Main struct for the pitch_shifter
 */
typedef struct _pitch_shifter
{
    /* FIFO buffers */
    int FIFO_idx;
    float inFIFO[MAX_NUM_CHANNELS][FRAME_SIZE];
    float outFIFO[MAX_NUM_CHANNELS][FRAME_SIZE];

    /* internal */
    void* hSmb;
    CODEC_STATUS codecStatus;
    float progressBar0_1;
    char* progressBarText;
    PROC_STATUS procStatus;
    float sampleRate;
    float inputFrame[MAX_NUM_CHANNELS][FRAME_SIZE];
    float outputFrame[MAX_NUM_CHANNELS][FRAME_SIZE];
    int new_nChannels;
    int fftFrameSize, stepsize;

    /* user parameters */
    int nChannels;
    float pitchShift_factor;   /**< 1: no shift, 0.5: down one octave, 2: up one octave */
    PITCH_SHIFTER_FFTSIZE_OPTIONS fftsize_option;
    PITCH_SHIFTER_OSAMP_OPTIONS osamp_option; 
    
} pitch_shifter_data;


/* ========================================================================== */
/*                             Internal Functions                             */
/* ========================================================================== */

/**
 * Sets codec status (see #_CODEC_STATUS enum)
 */
void pitch_shifter_setCodecStatus(void* const hPS,
                                  CODEC_STATUS newStatus);
    
    
#ifdef __cplusplus
} /* extern "C" { */
#endif /* __cplusplus */

#endif /* __PITCH_SHIFTER_INTERNAL_H_INCLUDED__ */
