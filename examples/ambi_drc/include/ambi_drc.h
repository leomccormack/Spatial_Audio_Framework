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
 * @file ambi_drc.h
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

#ifndef __AMBI_DRC_H_INCLUDED__
#define __AMBI_DRC_H_INCLUDED__

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#define ENABLE_TF_DISPLAY

/* ========================================================================== */
/*                             Presets + Constants                            */
/* ========================================================================== */

#define FRAME_SIZE ( 128 )
#define HOP_SIZE ( 128 )  /* STFT hop size, can be flexible, but only 'hybrid' mode afSTFT is supported (i.e. non uniform) */
#define TIME_SLOTS ( FRAME_SIZE/HOP_SIZE ) /* time-frequency domain frame size */
#define HYBRID_BANDS ( HOP_SIZE + 5 ) /* hybrid mode incurs an additional 5 bands  */
#define SPECTRAL_FLOOR (0.1585) /* -16dB, maximum gain reduction for a given frequency band */
#define AMBI_DRC_MAX_SH_ORDER ( 7 )
#define MAX_ORDER ( AMBI_DRC_MAX_SH_ORDER )
#define MAX_NUM_SH_SIGNALS ( (MAX_ORDER+1)*(MAX_ORDER+1) )
#ifdef ENABLE_TF_DISPLAY
# define NUM_DISPLAY_SECONDS ( 8 ) /* How many seconds the display will show historic TF data */
# define NUM_DISPLAY_TIME_SLOTS ( (int)(NUM_DISPLAY_SECONDS*48000.0f/(float)HOP_SIZE) )
# define READ_OFFSET ( 200 )
#endif

/**
 * Available Ambisonic channel ordering conventions
 *
 * @note CH_FUMA only supported for 1st order input.
 */
typedef enum _AMBI_DRC_CH_ORDER {
    CH_ACN = 1, /**< Ambisonic Channel Numbering (ACN) */
    CH_FUMA     /**< (Obsolete) Furse-Malham/B-format (WXYZ) */
    
} AMBI_DRC_CH_ORDER;

/** Number of channel ordering options */
#define AMBI_DRC_NUM_CH_ORDERINGS ( 2 )

/**
 * Available Ambisonic normalisation conventions
 *
 * @note NORM_FUMA only supported for 1st order input and does NOT have the
 *       1/sqrt(2) scaling on the omni.
 */
typedef enum _AMBI_DRC_NORM_TYPES {
    NORM_N3D = 1, /**< orthonormalised (N3D) */
    NORM_SN3D,    /**< Schmidt semi-normalisation (SN3D) */
    NORM_FUMA     /**< (Obsolete) Same as NORM_SN3D for 1st order */
    
} AMBI_DRC_NORM_TYPES;

/** Number of normalisation options */
#define AMBI_DRC_NUM_NORM_TYPES ( 3 )

/**
 * Available decoding orders
 */
typedef enum _AMBI_DRC_INPUT_ORDER{
    INPUT_ORDER_1 = 1, /**< First-order decoding (4 channel input) */
    INPUT_ORDER_2,     /**< Second-order decoding (9 channel input) */
    INPUT_ORDER_3,     /**< Third-order decoding (16 channel input) */
    INPUT_ORDER_4,     /**< Fourth-order decoding (25 channel input) */
    INPUT_ORDER_5,     /**< Fifth-order decoding (36 channel input) */
    INPUT_ORDER_6,     /**< Sixth-order decoding (49 channel input) */
    INPUT_ORDER_7      /**< Seventh-order decoding (64 channel input) */
    
}AMBI_DRC_INPUT_ORDER;

#define AMBI_DRC_IN_GAIN_MIN_VAL ( -40.0f )
#define AMBI_DRC_IN_GAIN_MAX_VAL ( 20.0f )
#define AMBI_DRC_THRESHOLD_MIN_VAL ( -60.0f )
#define AMBI_DRC_THRESHOLD_MAX_VAL ( 0.0f )
#define AMBI_DRC_RATIO_MIN_VAL ( 1.0f )
#define AMBI_DRC_RATIO_MAX_VAL ( 30.0f )
#define AMBI_DRC_KNEE_MIN_VAL ( 0.0f )
#define AMBI_DRC_KNEE_MAX_VAL ( 10.0f )
#define AMBI_DRC_ATTACK_MIN_VAL ( 10.0f )
#define AMBI_DRC_ATTACK_MAX_VAL ( 200.0f )
#define AMBI_DRC_RELEASE_MIN_VAL ( 50.0f )
#define AMBI_DRC_RELEASE_MAX_VAL ( 1000.0f )
#define AMBI_DRC_OUT_GAIN_MIN_VAL ( -20.0f )
#define AMBI_DRC_OUT_GAIN_MAX_VAL ( 40.0f )


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
                      float** const inputs,
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

/**
 * Sets the compressor threshold value in DECIBELS
 */
void ambi_drc_setThreshold(void* const hAmbi, float newValue);

/**
 * Sets the compression ratio
 */
void ambi_drc_setRatio(void* const hAmbi, float newValue);

/**
 * Sets the compressor knee value; 0: hard knee, >0: soft knee, in DECIBELS
 */
void ambi_drc_setKnee(void* const hAmbi, float newValue);

/**
 * Sets the compressor input gain value, in DECIBELS
 */
void ambi_drc_setInGain(void* const hAmbi, float newValue);
    
/**
 * Sets the compressor output gain value, in DECIBELS
 */
void ambi_drc_setOutGain(void* const hAmbi, float newValue);
    
/**
 * Sets the compressor envelope attack time, in miliseconds
 */
void ambi_drc_setAttack(void* const hAmbi, float newValue);

/**
 * Sets the compressor envelope release time, in miliseconds
 */
void ambi_drc_setRelease(void* const hAmbi, float newValue);
    
/**
 * Sets the Ambisonic channel ordering convention to decode with, in order to
 * match the convention employed by the input signals (see
 * 'AMBI_DRC_CH_ORDER' enum)
 */
void ambi_drc_setChOrder(void* const hAmbi, int newOrder);

/**
 * Sets the Ambisonic normalisation convention to decode with, in order to match
 * with the convention employed by the input signals (see
 * 'AMBI_DRC_NORM_TYPE' enum)
 */
void ambi_drc_setNormType(void* const hAmbi, int newType);
    
/**
 * Sets processing order.
 *
 * If input order is set higher than the input signal order, the extra required
 * channels are filled with zeros. If the input order is set lower than the
 * input signal order, the number input signals are truncated accordingly (see
 * 'AMBI_DRC_INPUT_ORDERS' enum)
 */
void ambi_drc_setInputPreset(void* const hAmbi, AMBI_DRC_INPUT_ORDER newPreset);

    
/* ========================================================================== */
/*                                Get Functions                               */
/* ========================================================================== */
    
#ifdef ENABLE_TF_DISPLAY
float** ambi_drc_getGainTF(void* const hAmbi);

int ambi_drc_getGainTFwIdx(void* const hAmbi);
    
int ambi_drc_getGainTFrIdx(void* const hAmbi);
    
float* ambi_drc_getFreqVector(void* const hAmbi, int* nFreqPoints);
#endif
    
/**
 * Returns the compressor threshold value, in DECIBELS
 */
float ambi_drc_getThreshold(void* const hAmbi);

/**
 * Returns the compression ratio
 */
float ambi_drc_getRatio(void* const hAmbi);

/**
 * Returns the compressor knee value 0: hard knee, >0: soft knee, in DECIBELS
 */
float ambi_drc_getKnee(void* const hAmbi);

/**
 * Returns the compressor input gain value, in DECIBELS
 */
float ambi_drc_getInGain(void* const hAmbi);

/**
 * Returns the compressor output gain value, in DECIBELS
 */
float ambi_drc_getOutGain(void* const hAmbi);

/**
 * Returns the compressor envelope attack time, in miliseconds
 */
float ambi_drc_getAttack(void* const hAmbi);

/**
 * Returns the compressor envelope release time, in miliseconds
 */
float ambi_drc_getRelease(void* const hAmbi);

/**
 * Returns the Ambisonic channel ordering convention currently being used to
 * decode with, which should match the convention employed by the input signals
 * (see 'AMBI_DRC_CH_ORDER' enum)
 */
int ambi_drc_getChOrder(void* const hAmbi);

/**
 * Returns the Ambisonic normalisation convention currently being usedto decode
 * with, which should match the convention employed by the input signals
 * (see 'AMBI_DRC_NORM_TYPE' enum)
 */
int ambi_drc_getNormType(void* const hAmbi);
    
/**
 * Returns the current processing order (see 'AMBI_DRC_INPUT_ORDERS' enum)
 */
AMBI_DRC_INPUT_ORDER ambi_drc_getInputPreset(void* const hAmbi);
    
/**
 * Returns the number of spherical harmonic signals required by the current
 * decoding order: (current_order+1)^2
 */
int ambi_drc_getNSHrequired(void* const hAmbi);
    
/**
 * Returns the DAW/Host sample rate
 */
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
