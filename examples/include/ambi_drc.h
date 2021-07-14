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
 * @example ambi_drc.h
 * @brief A frequency-dependent Ambisonic sound scene dynamic range compressor
 *        (DRC)
 *
 * ### Files
 * ambi_drc.h (include), ambi_drc_internal.h, ambi_drc.c, ambi_drc_internal.c
 * ### Include Header
 */

/**
 * @file ambi_drc.h
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

#ifndef __AMBI_DRC_H_INCLUDED__
#define __AMBI_DRC_H_INCLUDED__

#include "_common.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#define ENABLE_TF_DISPLAY /**< Enable TF data display related function */

/* ========================================================================== */
/*                             Presets + Constants                            */
/* ========================================================================== */
 
#ifdef ENABLE_TF_DISPLAY
/** How many seconds the display will show historic TF data */
# define AMBI_DRC_NUM_DISPLAY_SECONDS ( 8 )
/** Number of time slots of historic TF data */
# define AMBI_DRC_NUM_DISPLAY_TIME_SLOTS ( (int)(AMBI_DRC_NUM_DISPLAY_SECONDS*48000.0f/(float)128) )
/** Number of samples to offset when reading TF data */
# define AMBI_DRC_READ_OFFSET ( 200 )
/** Number of frequency bands used during processing */
# define AMBI_DRC_NUM_BANDS ( 133 )
#endif
/** -16dB, maximum gain reduction for a given frequency band */
#define AMBI_DRC_SPECTRAL_FLOOR (0.1585f)

#define AMBI_DRC_IN_GAIN_MIN_VAL ( -40.0f )   /**< Minimum input gain, dB */
#define AMBI_DRC_IN_GAIN_MAX_VAL ( 20.0f )    /**< Maximum input gain, dB */
#define AMBI_DRC_THRESHOLD_MIN_VAL ( -60.0f ) /**< Minimum threshold, dB */
#define AMBI_DRC_THRESHOLD_MAX_VAL ( 0.0f )   /**< Maximum threshold, dB */
#define AMBI_DRC_RATIO_MIN_VAL ( 1.0f )       /**< Minimum ratio, X:1 */
#define AMBI_DRC_RATIO_MAX_VAL ( 30.0f )      /**< Maximum ratio, X:1 */
#define AMBI_DRC_KNEE_MIN_VAL ( 0.0f )        /**< Minimum knee, dB */
#define AMBI_DRC_KNEE_MAX_VAL ( 10.0f )       /**< Maximum knee, dB */
#define AMBI_DRC_ATTACK_MIN_VAL ( 10.0f )     /**< Minimum attack time, ms */
#define AMBI_DRC_ATTACK_MAX_VAL ( 200.0f )    /**< Maximum attack time, ms */
#define AMBI_DRC_RELEASE_MIN_VAL ( 50.0f )    /**< Minimum release time, ms */
#define AMBI_DRC_RELEASE_MAX_VAL ( 1000.0f )  /**< Maximum release time, ms */
#define AMBI_DRC_OUT_GAIN_MIN_VAL ( -20.0f )  /**< Minimum output gain, dB */
#define AMBI_DRC_OUT_GAIN_MAX_VAL ( 40.0f )   /**< Maximum output gain, dB */


/* ========================================================================== */
/*                               Main Functions                               */
/* ========================================================================== */

/**
 * Creates an instance of the ambi_drc
 *
 * @param[in] phAmbi (&) address of ambi_drc handle
 */
void ambi_drc_create(void** const phAmbi);

/**
 * Destroys an instance of the ambi_drc
 *
 * @param[in] phAmbi (&) address of ambi_drc handle
 */
void ambi_drc_destroy(void** const phAmbi);

/**
 * Initialises an instance of ambi_drc with default settings
 *
 * @param[in] hAmbi      ambi_drc handle
 * @param[in] samplerate Host samplerate.
 */
void ambi_drc_init(void* const hAmbi,
                   int samplerate);

/**
 * Applies the frequency-dependent dynamic range compression to the input
 * spherical harmonic signals.
 *
 * @param[in] hAmbi    ambi_drc handle
 * @param[in] inputs   Input channel buffers; 2-D array: nCH x nSamples
 * @param[in] outputs  Output channel buffers; 2-D array: nCH x nSamples
 * @param[in] nCH      Number of input/output channels
 * @param[in] nSamples Number of samples in 'inputs'/'output' matrices
 */
void ambi_drc_process(void* const hAmbi,
                      const float *const * inputs,
                      float** const outputs,
                      int nCH,
                      int nSamples);


/* ========================================================================== */
/*                                Set Functions                               */
/* ========================================================================== */

/**
 * Sets all intialisation flags to 1; re-initialising all settings/variables
 * as ambi_drc is currently configured, at next available opportunity.
 */
void ambi_drc_refreshSettings(void* const hAmbi);

/** Sets the compressor threshold value in DECIBELS */
void ambi_drc_setThreshold(void* const hAmbi, float newValue);

/** Sets the compression ratio */
void ambi_drc_setRatio(void* const hAmbi, float newValue);

/** Sets the compressor knee value; 0: hard knee, >0: soft knee, in DECIBELS */
void ambi_drc_setKnee(void* const hAmbi, float newValue);

/** Sets the compressor input gain value, in DECIBELS */
void ambi_drc_setInGain(void* const hAmbi, float newValue);
    
/** Sets the compressor output gain value, in DECIBELS */
void ambi_drc_setOutGain(void* const hAmbi, float newValue);
    
/** Sets the compressor envelope attack time, in miliseconds */
void ambi_drc_setAttack(void* const hAmbi, float newValue);

/** Sets the compressor envelope release time, in miliseconds */
void ambi_drc_setRelease(void* const hAmbi, float newValue);
    
/**
 * Sets the Ambisonic channel ordering convention to decode with, in order to
 * match the convention employed by the input signals (see #CH_ORDER enum)
 */
void ambi_drc_setChOrder(void* const hAmbi, int newOrder);

/**
 * Sets the Ambisonic normalisation convention to decode with, in order to match
 * with the convention employed by the input signals (see #NORM_TYPES enum)
 */
void ambi_drc_setNormType(void* const hAmbi, int newType);
    
/**
 * Sets processing order.
 *
 * If input order is set higher than the input signal order, the extra required
 * channels are filled with zeros. If the input order is set lower than the
 * input signal order, the number input signals are truncated accordingly (see
 * #SH_ORDERS enum)
 */
void ambi_drc_setInputPreset(void* const hAmbi, SH_ORDERS newPreset);

    
/* ========================================================================== */
/*                                Get Functions                               */
/* ========================================================================== */

/**
 * Returns the processing framesize (i.e., number of samples processed with
 * every _process() call )
 */
int ambi_drc_getFrameSize(void);

/**
 * Returns pointers to historic time-frequency data, which may be used for
 * plotting purposes
 */
float** ambi_drc_getGainTF(void* const hAmbi);

/** Returns current TF gain data write index */
int ambi_drc_getGainTFwIdx(void* const hAmbi);

/** Returns current TF gain data read index */
int ambi_drc_getGainTFrIdx(void* const hAmbi);

/** Returns the frequency vector used by the processing */
float* ambi_drc_getFreqVector(void* const hAmbi, int* nFreqPoints);
    
/** Returns the compressor threshold value, in DECIBELS */
float ambi_drc_getThreshold(void* const hAmbi);

/** Returns the compression ratio */
float ambi_drc_getRatio(void* const hAmbi);

/** Returns the compressor knee value 0: hard knee, >0: soft knee, in DECIBELS*/
float ambi_drc_getKnee(void* const hAmbi);

/** Returns the compressor input gain value, in DECIBELS */
float ambi_drc_getInGain(void* const hAmbi);

/** Returns the compressor output gain value, in DECIBELS */
float ambi_drc_getOutGain(void* const hAmbi);

/** Returns the compressor envelope attack time, in miliseconds */
float ambi_drc_getAttack(void* const hAmbi);

/** Returns the compressor envelope release time, in miliseconds */
float ambi_drc_getRelease(void* const hAmbi);

/**
 * Returns the Ambisonic channel ordering convention currently being used to
 * decode with, which should match the convention employed by the input signals
 * (see #CH_ORDER enum)
 */
int ambi_drc_getChOrder(void* const hAmbi);

/**
 * Returns the Ambisonic normalisation convention currently being usedto decode
 * with, which should match the convention employed by the input signals
 * (see #NORM_TYPES enum)
 */
int ambi_drc_getNormType(void* const hAmbi);
    
/** Returns the current processing order (see #SH_ORDERS enum) */
SH_ORDERS ambi_drc_getInputPreset(void* const hAmbi);
    
/**
 * Returns the number of spherical harmonic signals required by the current
 * decoding order: (current_order+1)^2
 */
int ambi_drc_getNSHrequired(void* const hAmbi);
    
/** Returns the DAW/Host sample rate */
int ambi_drc_getSamplerate(void* const hAmbi);
    
/**
 * Returns the processing delay in samples; may be used for delay compensation
 * features
 */
int ambi_drc_getProcessingDelay(void);
    
    
#ifdef __cplusplus
} /* extern "C" { */
#endif /* __cplusplus */

#endif /* __AMBI_DRC_H_INCLUDED__ */
