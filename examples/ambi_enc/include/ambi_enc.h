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

/*
 * Filename: ambi_enc.h (include header)
 * -------------------------------------
 * A simple, but flexible, Ambisonic encoder.
 *
 * Dependencies:
 *     saf_utilities, saf_sh
 * Author, date created:
 *     Leo McCormack, 07.10.2016
 */

#ifndef __AMBI_ENC_H_INCLUDED__
#define __AMBI_ENC_H_INCLUDED__

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */
    
/* ========================================================================== */
/*                             Presets + Constants                            */
/* ========================================================================== */

/*
 * Enum: SOURCE_CONFIG_PRESETS
 * ---------------------------
 * Available source configurations presets to use for encoding
 */
typedef enum _SOURCE_CONFIG_PRESETS{
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
    
}SOURCE_CONFIG_PRESETS;
    
/*
 * Enum: OUTPUT_ORDERS
 * ------------------
 * Available encoding orders
 *
 * Options:
 *     OUTPUT_ORDER_FIRST   - First-order encoding (4 channel output)
 *     OUTPUT_ORDER_SECOND  - Second-order encoding (9 channel output)
 *     OUTPUT_ORDER_THIRD   - Third-order encoding (16 channel output)
 *     OUTPUT_ORDER_FOURTH  - Fourth-order encoding (25 channel output)
 *     OUTPUT_ORDER_FIFTH   - Fifth-order encoding (36 channel output)
 *     OUTPUT_ORDER_SIXTH   - Sixth-order encoding (49 channel output)
 *     OUTPUT_ORDER_SEVENTH - Seventh-order encoding (64 channel output)
 */
#define AMBI_ENC_MAX_SH_ORDER ( 7 )
typedef enum _OUTPUT_ORDERS{ 
    OUTPUT_ORDER_FIRST = 1,
    OUTPUT_ORDER_SECOND,
    OUTPUT_ORDER_THIRD,
    OUTPUT_ORDER_FOURTH,
    OUTPUT_ORDER_FIFTH,
    OUTPUT_ORDER_SIXTH,
    OUTPUT_ORDER_SEVENTH
    
}OUTPUT_ORDERS;
    
/*
 * Enum: _CH_ORDER
 * ---------------
 * Available Ambisonic channel ordering conventions
 * Note: CH_FUMA only supported for 1st order output.
 * Further note: FuMa: CH_FUMA+NORM_FUMA, AmbiX: CH_ACN+NORM_SN3D
 *
 * Options:
 *     CH_ACN  - Ambisonic Channel Numbering (ACN)
 *     CH_FUMA - (Obsolete) Furse-Malham/B-format (WXYZ)
 */
#define AMBI_ENC_NUM_CH_ORDERINGS ( 2 )
typedef enum _CH_ORDER{
    CH_ACN = 1,
    CH_FUMA     /* first-order only */
    
}CH_ORDER;

/*
 * Enum: NORM_TYPES
 * ---------------
 * Available Ambisonic normalisation conventions
 * Note: NORM_FUMA only supported for 1st order output and does NOT have the
 * 1/sqrt(2) scaling on the omni.
 * Further note: FuMa: CH_FUMA+NORM_FUMA, AmbiX: CH_ACN+NORM_SN3D
 *
 * Options:
 *     NORM_N3D  - orthonormalised (N3D)
 *     NORM_SN3D - Schmidt semi-normalisation (SN3D)
 *     NORM_FUMA - (Obsolete) Same as NORM_SN3D for 1st order
 */
#define AMBI_ENC_NUM_NORM_TYPES ( 3 )
typedef enum _NORM_TYPES{
    NORM_N3D = 1,
    NORM_SN3D,
    NORM_FUMA   /* first-order only */
    
}NORM_TYPES;
    
#define AMBI_ENC_MAX_NUM_INPUTS ( 64 )
    
    
/* ========================================================================== */
/*                               Main Functions                               */
/* ========================================================================== */

/*
 * Function: ambi_enc_create
 * -------------------------
 * Creates an instance of ambi_enc
 *
 * Input Arguments:
 *     phAmbi - & address of ambi_enc handle
 */
void ambi_enc_create(void** const phAmbi);


/*
 * Function: ambi_enc_destroy
 * --------------------------
 * Destroys an instance of ambi_enc
 *
 * Input Arguments:
 *     phAmbi - & address of ambi_enc handle
 */
void ambi_enc_destroy(void** const phAmbi);


/*
 * Function: ambi_enc_init
 * -----------------------
 * Initialises an instance of ambi_enc with default settings
 *
 * Input Arguments:
 *     hAmbi      - ambi_enc handle
 *     samplerate - host samplerate.
 */
void ambi_enc_init(void* const hAmbi,
                     int samplerate);
    
    
/*
 * Function: ambi_enc_process
 * --------------------------
 * Encodes input signals into spherical harmonic signals, at the specified
 * encoding directions.
 *
 * Input Arguments:
 *     hAmbi     - ambi_enc handle
 *     inputs    - input channel buffers; 2-D array: nInputs x nSamples
 *     outputs   - Output channel buffers; 2-D array: nOutputs x nSamples
 *     nInputs   - number of input channels
 *     nOutputs  - number of output channels
 *     nSamples  - number of samples in 'inputs'/'output' matrices
 *     isPlaying - flag to say if there is audio in the input buffers, 0: no
 *                 audio, reduced processing, 1: audio, full processing
 */
void ambi_enc_process(void* const hAmbi,
                        float** const inputs,
                        float** const outputs,
                        int nInputs,
                        int nOutputs,
                        int nSamples,
                        int isPlaying);

    
/* ========================================================================== */
/*                                Set Functions                               */
/* ========================================================================== */

/*
 * Function: ambi_enc_refreshParams
 * --------------------------------
 * Sets all intialisation flags to 1. i.e. re-initialise all settings/variables
 * as ambi_enc is currently configured, at next available opportunity.
 *
 * Input Arguments:
 *     hAmbi - ambi_enc handle
 */
void ambi_enc_refreshParams(void* const hAmbi);
    
/*
 * Function: ambi_enc_setOutputOrder
 * ---------------------------------
 * Sets the encoding order.
 *
 * Input Arguments:
 *     hAmbi    - ambi_enc handle
 *     newValue - new encoding order (see 'OUTPUT_ORDERS' enum)
 */
void ambi_enc_setOutputOrder(void* const hAmbi, int newValue);

/*
 * Function: ambi_enc_setSourceAzi_deg
 * -----------------------------------
 * Sets the azimuth for a specific source index
 *
 * Input Arguments:
 *     hAmbi      - ambi_enc handle
 *     index      - source index
 *     newAzi_deg - new azimuth, in DEGREES
 */
void ambi_enc_setSourceAzi_deg(void* const hAmbi, int index, float newAzi_deg);

/*
 * Function: ambi_enc_setSourceElev_deg
 * ------------------------------------
 * Sets the elevation for a specific source index
 *
 * Input Arguments:
 *     hAmbi       - ambi_enc handle
 *     index       - source index
 *     newElev_deg - new elevation, in DEGREES
 */
void ambi_enc_setSourceElev_deg(void* const hAmbi, int index, float newElev_deg);

/*
 * Function: ambi_enc_setNumSources
 * --------------------------------
 * Sets the number of input signals/sources to encode.
 *
 * Input Arguments:
 *     hAmbi        - ambi_enc handle
 *     new_nSources - new number of input signals/sources
 */
void ambi_enc_setNumSources(void* const hAmbi, int new_nSources);

/*
 * Function: ambi_enc_setInputConfigPreset
 * ---------------------------------------
 * Sets the input configuration preset (such as: 5.x, 9.1x, 22.x etc).
 *
 * Input Arguments:
 *     hAmbi       - ambi_enc handle
 *     newPresetID - the configuration preset (see 'PRESETS' enum)
 */
void ambi_enc_setInputConfigPreset(void* const hAmbi, int newPresetID);

/*
 * Function: ambi_enc_setChOrder
 * -----------------------------
 * Sets the Ambisonic channel ordering convention to encode with.
 *
 * Input Arguments:
 *     hAmbi    - ambi_enc handle
 *     newOrder - convention to use (see 'CH_ORDER' enum)
 */
void ambi_enc_setChOrder(void* const hAmbi, int newOrder);

/*
 * Function: ambi_enc_setNormType
 * ------------------------------
 * Sets the Ambisonic normalisation convention to encode with.
 *
 * Input Arguments:
 *     hAmbi   - ambi_enc handle
 *     newType - convention to use (see 'NORM_TYPE' enum)
 */
void ambi_enc_setNormType(void* const hAmbi, int newType);

    
/* ========================================================================== */
/*                                Get Functions                               */
/* ========================================================================== */

/*
 * Function: ambi_enc_getInputOrderPreset
 * --------------------------------------
 * Returns the decoding order. If decoding order is higher than the input signal
 * order, the extra required channels are filled with zeros. If the decoding
 * order is lower than the input signal order, the number input signals is
 * truncated accordingly.
 *
 * Input Arguments:
 *     hAmbi - ambi_enc handle
 * Returns:
 *     current decoding order (see 'INPUT_ORDERS' enum)
 */
int ambi_enc_getOutputOrder(void* const hAmbi);

/*
 * Function: ambi_enc_getSourceAzi_deg
 * -----------------------------------
 * Returns the azimuth for a specific source
 *
 * Input Arguments:
 *     hAmbi - ambi_enc handle
 *     index - source index
 * Returns:
 *     source azimuth, in DEGREES
 */
float ambi_enc_getSourceAzi_deg(void* const hAmbi, int index);

/*
 * Function: ambi_enc_getSourceElev_deg
 * ------------------------------------
 * Returns the elevation for a specific source
 *
 * Input Arguments:
 *     hAmbi - ambi_enc handle
 *     index - source index
 * Returns:
 *     source elevation, in DEGREES
 */
float ambi_enc_getSourceElev_deg(void* const hAmbi, int index);

/*
 * Function: ambi_enc_getNumSources
 * --------------------------------
 * Returns the number of input signals/sources to encode.
 *
 * Input Arguments:
 *     hAmbi        - ambi_enc handle
 * Returns:
 *     number of input signals/sources
 */
int ambi_enc_getNumSources(void* const hAmbi);

int ambi_enc_getMaxNumSources(void);

/*
 * Function: ambi_enc_getNSHrequired
 * ---------------------------------
 * Returns the number of spherical harmonic signals required by the current
 * decoding order i.e. (current_order+1)^2
 *
 * Input Arguments:
 *     hAmbi - ambi_enc handle
 * Returns:
 *     number of required spherical harmonic signals required by current
 *     decoding order
 */
int ambi_enc_getNSHrequired(void* const hAmbi);

/*
 * Function: ambi_enc_getChOrder
 * -----------------------------
 * Returns the Ambisonic channel ordering convention currently being used to
 * encode with.
 *
 * Input Arguments:
 *     hAmbi - ambi_enc handle
 * Returns:
 *     convention currently being used (see 'CH_ORDER' enum)
 */
int ambi_enc_getChOrder(void* const hAmbi);
    
/*
 * Function: ambi_enc_getNormType
 * ------------------------------
 * Returns the Ambisonic normalisation convention currently being used to encode
 * with.
 *
 * Input Arguments:
 *     hAmbi - ambi_enc handle
 * Returns:
 *     convention currently being used (see 'NORM_TYPE' enum)
 */
int ambi_enc_getNormType(void* const hAmbi);

    
#ifdef __cplusplus
} /* extern "C" { */
#endif /* __cplusplus */

#endif /* __AMBI_ENC_H_INCLUDED__ */
