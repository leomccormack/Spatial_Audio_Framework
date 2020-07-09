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
 * @file decorrelator_internal.h
 * @brief A multi-channel decorrelator
 *
 * @author Leo McCormack
 * @date 07.07.2020
 */


#ifndef __DECORRELATOR_INTERNAL_H_INCLUDED__
#define __DECORRELATOR_INTERNAL_H_INCLUDED__

#include <stdio.h>
#include <math.h>
#include <string.h>
#include "decorrelator.h"
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
#define HOP_SIZE ( 128 ) /* STFT hop size */
#define HYBRID_BANDS ( 133 )
#define TIME_SLOTS ( FRAME_SIZE / HOP_SIZE )  

    
/* ========================================================================== */
/*                                 Structures                                 */
/* ========================================================================== */

/**
 * Main structure for decorrelator. Contains variables for audio buffers, afSTFT,
 * rotation matrices, internal variables, flags, user parameters
 */
typedef struct _decorrelator
{
    /* audio buffers + afSTFT time-frequency transform handle */
    int fs;                         /**< host sampling rate */ 
    float InputFrameTD[MAX_NUM_SH_SIGNALS][FRAME_SIZE];
    float_complex*** InputFrameTF;
    float_complex*** OutputFrameTF;
    complexVector* STFTInputFrameTF;
    complexVector* STFTOutputFrameTF;
    void* hSTFT;                    /**< afSTFT handle */
    int afSTFTdelay;                /**< for host delay compensation */
    float** tempHopFrameTD;         /**< temporary multi-channel time-domain buffer of size "HOP_SIZE". */
    float freqVector[HYBRID_BANDS]; /**< frequency vector for time-frequency transform, in Hz */
     
    /* our codec configuration */
    void* hDec;
    void* hDec2;
    void* hDucker;
    CODEC_STATUS codecStatus;
    float progressBar0_1;
    char* progressBarText;
    
    /* internal variables */
    PROC_STATUS procStatus;
    int new_nChannels;

    /* user parameters */
    int nChannels;
    int enableTransientDucker;
    
} decorrelator_data;


/* ========================================================================== */
/*                             Internal Functions                             */
/* ========================================================================== */

/**
 * Sets codec status. 
 */
void decorrelator_setCodecStatus(void* const hDecor,
                                 CODEC_STATUS newStatus);


#ifdef __cplusplus
} /* extern "C" { */
#endif /* __cplusplus */

#endif /* __DECORRELATOR_INTERNAL_H_INCLUDED__ */
