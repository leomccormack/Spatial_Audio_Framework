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
 * @brief A frequency-dependent Ambisonic sound scene dynamic range compressor
 *        (DRC)
 *
 * The implementation can also keep track of the frequency-dependent gain
 * factors for the omnidirectional component over time (for optional plotting).
 * The design is based on the algorithm presented in [1].
 *
 * The DRC gain factors per band are determined based on the omnidirectional
 * component, which are then applied to all of the higher-order components;
 * thus, the spatial information of the Ambisonic sound scene is retained
 * (although, your perception of them may change due to the DRC).
 *
 * @see [1] McCormack, L., & Välimäki, V. (2017). "FFT-Based Dynamic Range
 *          Compression". in Proceedings of the 14th Sound and Music Computing
 *          Conference, July 5-8, Espoo, Finland.
 *
 * @author Leo McCormack
 * @date 07.01.2017
 * @license ISC
 */

#ifndef __AMBI_DRC_INTERNAL_H_INCLUDED__
#define __AMBI_DRC_INTERNAL_H_INCLUDED__

#include "ambi_drc.h"      /* Include header for this example */
#include "saf.h"           /* Main include header for SAF */
#include "saf_externals.h" /* To also include SAF dependencies (cblas etc.) */

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/* ========================================================================== */
/*                            Internal Parameters                             */
/* ========================================================================== */

#if !defined(AMBI_DRC_FRAME_SIZE)
# if defined(FRAME_SIZE) /* Use the global framesize if it is specified: */
#  define AMBI_DRC_FRAME_SIZE ( FRAME_SIZE )          /**< Framesize, in time-domain samples */
# else /* Otherwise, the default framesize for this example is: */
#  define AMBI_DRC_FRAME_SIZE ( 128 )                 /**< Framesize, in time-domain samples */
# endif
#endif
#define HOP_SIZE ( 128 )                              /**< STFT hop size */
#define HYBRID_BANDS ( HOP_SIZE + 5 )                 /**< Number of frequency bands */
#define TIME_SLOTS ( AMBI_DRC_FRAME_SIZE / HOP_SIZE ) /**< Number of STFT timeslots */

/* Checks: */
#if (AMBI_DRC_FRAME_SIZE % HOP_SIZE != 0)
# error "AMBI_DRC_FRAME_SIZE must be an integer multiple of HOP_SIZE"
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
    float** frameTD;                 /**< Input/output SH signals, in the time-domain; #MAX_NUM_SH_SIGNALS x #AMBI_DRC_FRAME_SIZE */
    float_complex*** inputFrameTF;   /**< Input SH signals, in the time-frequency domain; #HYBRID_BANDS x #MAX_NUM_SH_SIGNALS x #TIME_SLOTS */
    float_complex*** outputFrameTF;  /**< Output SH signals, in the time-frequency domain; #HYBRID_BANDS x #MAX_NUM_SH_SIGNALS x #TIME_SLOTS */
    void* hSTFT;                     /**< Time-frequency transform handle */
    float freqVector[HYBRID_BANDS];  /**< Frequency vector */

    /* internal */
    int nSH;                         /**< Current number of SH signals */
    int new_nSH;                     /**< New number of SH signals (current value will be replaced by this after next re-init) */
    float fs;                        /**< Host sampling rate, in Hz */
    float yL_z1[HYBRID_BANDS];       /**< Delay elements */
    int reInitTFT;                   /**< 0: no init required, 1: init required, 2: init in progress */

#ifdef ENABLE_TF_DISPLAY
    int wIdx;                        /**< Display slot write index */
    int rIdx;                        /**< Display slot read index */
    int storeIdx;                    /**< Display slot storage index */
    float** gainsTF_bank0;           /**< Display slot "0" DRC gains */
    float** gainsTF_bank1;           /**< Display slot "1" DRC gains */
#endif

    /* user parameters */
    float theshold;                  /**< Threshold parameter, in dB */
    float ratio;                     /**< Compression ratio */
    float knee;                      /**< Knee width, in dB */
    float inGain;                    /**< Pre-gain, in dB */
    float outGain;                   /**< Post-gain, in dB*/
    float attack_ms;                 /**< Attack time, in ms */
    float release_ms;                /**< Release time, in ms */
    CH_ORDER chOrdering;             /**< Ambisonic channel order convention (see #CH_ORDER) */
    NORM_TYPES norm;                 /**< Ambisonic normalisation convention (see #NORM_TYPES) */
    SH_ORDERS currentOrder;          /**< Current input SH order */
    
} ambi_drc_data;
     

/* ========================================================================== */
/*                             Internal Functions                             */
/* ========================================================================== */

/** The DRC gain computer */
float ambi_drc_gainComputer(float xG, float T, float R, float W);

/** The envelope detector */
float ambi_drc_smoothPeakDetector(float xL, float yL_z1, float alpha_a, float alpha_r);
    
/** Initialise the filterbank used by ambi_drc */
void ambi_drc_initTFT(void* const hAmbi);

/** Sets the internal input order */
void ambi_drc_setInputOrder(SH_ORDERS inOrder, int* nSH);

    
#ifdef __cplusplus
} /* extern "C" */
#endif/* __cplusplus */

#endif /* __AMBI_DRC_INTERNAL_H_INCLUDED__ */
