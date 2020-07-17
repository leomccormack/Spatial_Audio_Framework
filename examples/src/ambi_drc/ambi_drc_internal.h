/*
 * Copyright 2017-2018 Leo McCormack
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
 * @file ambi_drc_internal.h
 * @brief A frequency-dependent spherical harmonic domain dynamic range
 *        compressor (DRC).
 *
 * The implementation can also keep track of the frequency-dependent gain
 * factors for the omnidirectional component over time, for optional plotting.
 * The design is based on the algorithm presented in [1].
 *
 * The DRC gain factors are determined based on analysing the omnidirectional
 * component. These gain factors are then applied to the higher-order
 * components, in a such a manner as to retain the spatial information within
 * them.
 *
 * @author Leo McCormack
 * @date 07.01.2017
 *
 * @see [1] McCormack, L., & Välimäki, V. (2017). "FFT-Based Dynamic Range
 *          Compression". in Proceedings of the 14th Sound and Music Computing
 *          Conference, July 5-8, Espoo, Finland.
 */

#ifndef __AMBI_DRC_INTERNAL_H_INCLUDED__
#define __AMBI_DRC_INTERNAL_H_INCLUDED__

#include <stdio.h>
#include <math.h>
#include <string.h>
#include "ambi_drc.h" 
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
#define HOP_SIZE ( 128 )  /* STFT hop size, can be flexible, but only 'hybrid' mode afSTFT is supported (i.e. non uniform) */
#define TIME_SLOTS ( FRAME_SIZE/HOP_SIZE ) /* time-frequency domain frame size */
#define HYBRID_BANDS ( HOP_SIZE + 5 ) /* hybrid mode incurs an additional 5 bands  */
#if (FRAME_SIZE % HOP_SIZE != 0)
# error "FRAME_SIZE must be an integer multiple of HOP_SIZE"
#endif

/* ========================================================================== */
/*                                 Structures                                 */
/* ========================================================================== */
    
/**
 * Main structure for ambi_drc. Contains variables for audio buffers, afSTFT,
 * internal variables, user parameters
 */
typedef struct _ambi_drc
{ 
    /* audio buffers and afSTFT handle */
    float** frameTD;
    float_complex*** inputFrameTF;
    float_complex*** outputFrameTF;
    void* hSTFT;
    float freqVector[HYBRID_BANDS];

    /* internal */
    int nSH, new_nSH;
    float fs;
    float yL_z1[HYBRID_BANDS];
    int reInitTFT; /**< 0: no init required, 1: init required, 2: init in progress */

#ifdef ENABLE_TF_DISPLAY
    int wIdx, rIdx;
    int storeIdx;
    float** gainsTF_bank0;
    float** gainsTF_bank1;
#endif

    /* user parameters */
    float theshold, ratio, knee, inGain, outGain, attack_ms, release_ms;
    CH_ORDER chOrdering;
    NORM_TYPES norm;
    SH_ORDERS currentOrder;
    
} ambi_drc_data;
     

/* ========================================================================== */
/*                             Internal Functions                             */
/* ========================================================================== */

float ambi_drc_gainComputer(float xG, float T, float R, float W);

float ambi_drc_smoothPeakDetector(float xL, float yL_z1, float alpha_a, float alpha_r);
    
/**
 * Initialise the filterbank used by ambi_drc.
 */
void ambi_drc_initTFT(void* const hAmbi);

void ambi_drc_setInputOrder(SH_ORDERS inOrder, int* nSH);

    
#ifdef __cplusplus
} /* extern "C" */
#endif/* __cplusplus */

#endif /* __AMBI_DRC_INTERNAL_H_INCLUDED__ */
