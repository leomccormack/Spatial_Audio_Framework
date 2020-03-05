/*
 * Copyright 2016-2018 Leo McCormack
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
 * @file ambi_enc.h
 * @brief A simple, but flexible, Ambisonic encoder.
 * 
 * @author Leo McCormack
 * @date 07.10.2016
 */

#ifndef __AMBI_ENC_H_INCLUDED__
#define __AMBI_ENC_H_INCLUDED__

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/* ========================================================================== */
/*                             Presets + Constants                            */
/* ========================================================================== */

/**
 * Available source configurations presets to use for encoding
 */
typedef enum _AMBI_ENC_SOURCE_CONFIG_PRESETS{
    SOURCE_CONFIG_PRESET_DEFAULT = 1,
    SOURCE_CONFIG_PRESET_MONO,
    SOURCE_CONFIG_PRESET_STEREO,
    SOURCE_CONFIG_PRESET_5PX,
    SOURCE_CONFIG_PRESET_7PX,
    SOURCE_CONFIG_PRESET_8PX,
    SOURCE_CONFIG_PRESET_9PX,
    SOURCE_CONFIG_PRESET_10PX,
    SOURCE_CONFIG_PRESET_11PX,
    SOURCE_CONFIG_PRESET_11PX_7_4,
    SOURCE_CONFIG_PRESET_13PX,
    SOURCE_CONFIG_PRESET_22PX,
    SOURCE_CONFIG_PRESET_AALTO_MCC,
    SOURCE_CONFIG_PRESET_AALTO_MCC_SUBSET,
    SOURCE_CONFIG_PRESET_AALTO_APAJA,
    SOURCE_CONFIG_PRESET_AALTO_LR,
    SOURCE_CONFIG_PRESET_DTU_AVIL,
    SOURCE_CONFIG_PRESET_ZYLIA_LAB,
    SOURCE_CONFIG_PRESET_T_DESIGN_4,
    SOURCE_CONFIG_PRESET_T_DESIGN_12,
    SOURCE_CONFIG_PRESET_T_DESIGN_24,
    SOURCE_CONFIG_PRESET_T_DESIGN_36,
    SOURCE_CONFIG_PRESET_T_DESIGN_48,
    SOURCE_CONFIG_PRESET_T_DESIGN_60
    
}AMBI_ENC_SOURCE_CONFIG_PRESETS;

/**
 * Available encoding orders
 */
typedef enum _AMBI_ENC_OUTPUT_ORDERS{
    OUTPUT_ORDER_FIRST = 1, /**< First-order encoding (4 channel output) */
    OUTPUT_ORDER_SECOND,    /**< Second-order encoding (9 channel output) */
    OUTPUT_ORDER_THIRD,     /**< Third-order encoding (16 channel output) */
    OUTPUT_ORDER_FOURTH,    /**< Fourth-order encoding (25 channel output) */
    OUTPUT_ORDER_FIFTH,     /**< Fifth-order encoding (36 channel output) */
    OUTPUT_ORDER_SIXTH,     /**< Sixth-order encoding (49 channel output) */
    OUTPUT_ORDER_SEVENTH    /**< Seventh-order encoding (64 channel output) */
    
}AMBI_ENC_OUTPUT_ORDERS;

/** Maximum supported Ambisonic order */
#define AMBI_ENC_MAX_SH_ORDER ( 7 )

/**
 * Available Ambisonic channel ordering conventions
 *
 * @note CH_FUMA only supported for 1st order input.
 */
typedef enum _AMBI_ENC_CH_ORDER {
    CH_ACN = 1, /**< Ambisonic Channel Numbering (ACN) */
    CH_FUMA     /**< (Obsolete) Furse-Malham/B-format (WXYZ) */
    
} AMBI_ENC_CH_ORDER;

/** Number of channel ordering options */
#define AMBI_ENC_NUM_CH_ORDERINGS ( 2 )

/**
 * Available Ambisonic normalisation conventions
 *
 * @note NORM_FUMA only supported for 1st order input and does NOT have the
 *       1/sqrt(2) scaling on the omni.
 */
typedef enum _AMBI_ENC_NORM_TYPES {
    NORM_N3D = 1, /**< orthonormalised (N3D) */
    NORM_SN3D,    /**< Schmidt semi-normalisation (SN3D) */
    NORM_FUMA     /**< (Obsolete) Same as NORM_SN3D for 1st order */
    
} AMBI_ENC_NORM_TYPES;

/** Number of normalisation options */
#define AMBI_ENC_NUM_NORM_TYPES ( 3 )
    
/** Maximum number of inputs */
#define AMBI_ENC_MAX_NUM_INPUTS ( 64 )


/* ========================================================================== */
/*                               Main Functions                               */
/* ========================================================================== */

/**
 * Creates an instance of ambi_enc
 *
 * @param[in] phAmbi (&) address of ambi_enc handle
 */
void ambi_enc_create(void** const phAmbi);

/**
 * Destroys an instance of ambi_enc
 *
 * @param[in] phAmbi (&) address of ambi_enc handle
 */
void ambi_enc_destroy(void** const phAmbi);

/**
 * Initialises an instance of ambi_enc with default settings
 *
 * @param[in] hAmbi      ambi_enc handle
 * @param[in] samplerate Host samplerate.
 */
void ambi_enc_init(void* const hAmbi,
                     int samplerate);

/**
 * Encodes input signals into spherical harmonic signals, at the specified
 * encoding directions.
 *
 * @param[in] hAmbi    Ambi_enc handle
 * @param[in] inputs   Input channel buffers; 2-D array: nInputs x nSamples
 * @param[in] outputs  Output channel buffers; 2-D array: nOutputs x nSamples
 * @param[in] nInputs  Number of input channels
 * @param[in] nOutputs Number of output channels
 * @param[in] nSamples Number of samples in 'inputs'/'output' matrices
 */
void ambi_enc_process(void* const hAmbi,
                      float** const inputs,
                      float** const outputs,
                      int nInputs,
                      int nOutputs,
                      int nSamples);

    
/* ========================================================================== */
/*                                Set Functions                               */
/* ========================================================================== */

/**
 * Sets all intialisation flags to 1; re-initialising all settings/variables
 * as ambi_enc is currently configured, at next available opportunity.
 */
void ambi_enc_refreshParams(void* const hAmbi);
    
/**
 * Sets the encoding order (see 'AMBI_ENC_OUTPUT_ORDERS' enum)
 */
void ambi_enc_setOutputOrder(void* const hAmbi, int newValue);

/**
 * Sets the azimuth for a specific source index
 *
 * @param[in] hAmbi      ambi_enc handle
 * @param[in] index      Source index
 * @param[in] newAzi_deg New azimuth, in DEGREES
 */
void ambi_enc_setSourceAzi_deg(void* const hAmbi, int index, float newAzi_deg);

/**
 * Sets the elevation for a specific source index
 *
 * @param[in] hAmbi       ambi_enc handle
 * @param[in] index       Source index
 * @param[in] newElev_deg New elevation, in DEGREES
 */
void ambi_enc_setSourceElev_deg(void* const hAmbi, int index, float newElev_deg);

/**
 * Sets the number of input signals/sources to encode.
 */
void ambi_enc_setNumSources(void* const hAmbi, int new_nSources);

/**
 * Sets the input configuration preset (see '_AMBI_ENC_SOURCE_CONFIG_PRESETS'
 * enum)
 */
void ambi_enc_setInputConfigPreset(void* const hAmbi, int newPresetID);

/**
 * Sets the Ambisonic channel ordering convention to encode with (see
 * 'AMBI_ENC_CH_ORDER' enum)
 */
void ambi_enc_setChOrder(void* const hAmbi, int newOrder);

/**
 * Sets the Ambisonic normalisation convention to encode with (see
 * 'AMBI_ENC_NORM_TYPE' enum)
 */
void ambi_enc_setNormType(void* const hAmbi, int newType);

    
/* ========================================================================== */
/*                                Get Functions                               */
/* ========================================================================== */

/**
 * Returns the decoding order (see 'AMBI_ENC_INPUT_ORDERS' enum)
 *
 * If decoding order is higher than the input signal order, the extra required
 * channels are filled with zeros. If the decoding order is lower than the input
 * signal order, the number input signals is truncated accordingly.
 */
int ambi_enc_getOutputOrder(void* const hAmbi);

/**
 * Returns the azimuth for a specific source, in DEGREES
 */
float ambi_enc_getSourceAzi_deg(void* const hAmbi, int index);

/**
 * Returns the elevation for a specific source, in DEGREES
 */
float ambi_enc_getSourceElev_deg(void* const hAmbi, int index);

/**
 * Returns the number of input signals/sources to encode.
 */
int ambi_enc_getNumSources(void* const hAmbi);

/**
 * Returns the maximum number of input signals/sources that can be encoded.
 */
int ambi_enc_getMaxNumSources(void);

/**
 * Returns the number of spherical harmonic signals required by the current
 * decoding order: (current_order+1)^2
 */
int ambi_enc_getNSHrequired(void* const hAmbi);

/**
 * Returns the Ambisonic channel ordering convention currently being used to
 * encode with (see 'AMBI_ENC_CH_ORDER' enum)
 */
int ambi_enc_getChOrder(void* const hAmbi);
    
/**
 * Returns the Ambisonic normalisation convention currently being used to encode
 * with (see 'AMBI_ENC_NORM_TYPE' enum)
 */
int ambi_enc_getNormType(void* const hAmbi);


#ifdef __cplusplus
} /* extern "C" { */
#endif /* __cplusplus */

#endif /* __AMBI_ENC_H_INCLUDED__ */
