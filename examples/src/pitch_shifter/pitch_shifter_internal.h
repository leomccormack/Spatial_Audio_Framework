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
 * @license ISC
 */

#ifndef __PITCH_SHIFTER_INTERNAL_H_INCLUDED__
#define __PITCH_SHIFTER_INTERNAL_H_INCLUDED__

#include "pitch_shifter.h" /* Include header for this example */
#include "saf.h"           /* Main include header for SAF */
#include "saf_externals.h" /* To also include SAF dependencies (cblas etc.) */

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/* ========================================================================== */
/*                            Internal Parameters                             */
/* ========================================================================== */

#if !defined(PITCH_SHIFTER_FRAME_SIZE)
# if defined(FRAME_SIZE) /* Use the global framesize if it is specified: */
#  define PITCH_SHIFTER_FRAME_SIZE ( FRAME_SIZE )  /**< Framesize, in time-domain samples */
# else /* Otherwise, the default framesize for this example is: */
#  define PITCH_SHIFTER_FRAME_SIZE ( 128 )         /**< Framesize, in time-domain samples */
# endif
#endif

/* ========================================================================== */
/*                                 Structures                                 */
/* ========================================================================== */

/** Main struct for the pitch_shifter */
typedef struct _pitch_shifter
{
    /* FIFO buffers */
    int FIFO_idx;                   /**< FIFO buffer index */
    float inFIFO[MAX_NUM_CHANNELS][PITCH_SHIFTER_FRAME_SIZE];  /**< Input FIFO buffer */
    float outFIFO[MAX_NUM_CHANNELS][PITCH_SHIFTER_FRAME_SIZE]; /**< Output FIFO buffer */

    /* internal */
    void* hSmb;                     /**< pitch-shifter handle */
    CODEC_STATUS codecStatus;       /**< see #CODEC_STATUS */
    float progressBar0_1;           /**< Current (re)initialisation progress, between [0..1] */
    char* progressBarText;          /**< Current (re)initialisation step, string */
    PROC_STATUS procStatus;         /**< see #PROC_STATUS */
    float sampleRate;               /**< Host sampling rate, in Hz */
    float inputFrame[MAX_NUM_CHANNELS][PITCH_SHIFTER_FRAME_SIZE];  /**< Current input frame */
    float outputFrame[MAX_NUM_CHANNELS][PITCH_SHIFTER_FRAME_SIZE]; /**< Current output frame */
    int new_nChannels;              /**< (current value will be replaced by this after next re-init) */
    int fftFrameSize;               /**< FFT size */
    int stepsize;                   /**< Hop size in samples*/

    /* user parameters */
    int nChannels;                  /**< Current number of input/output channels */
    float pitchShift_factor;        /**< 1: no shift, 0.5: down one octave, 2: up one octave */
    PITCH_SHIFTER_FFTSIZE_OPTIONS fftsize_option; /**< see #PITCH_SHIFTER_FFTSIZE_OPTIONS */
    PITCH_SHIFTER_OSAMP_OPTIONS osamp_option;     /**< see #PITCH_SHIFTER_OSAMP_OPTIONS */
    
} pitch_shifter_data;


/* ========================================================================== */
/*                             Internal Functions                             */
/* ========================================================================== */

/** Sets codec status (see #CODEC_STATUS enum) */
void pitch_shifter_setCodecStatus(void* const hPS,
                                  CODEC_STATUS newStatus);
    
    
#ifdef __cplusplus
} /* extern "C" { */
#endif /* __cplusplus */

#endif /* __PITCH_SHIFTER_INTERNAL_H_INCLUDED__ */
