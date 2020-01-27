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

/*
 * Filename: sldoa.h (include header)
 * ----------------------------------
 * A spatially-localised active-intensity based direction-of-arrival estimator
 * (SLDoA). VBAP gain patterns are imposed on the spherical harmonic signals,
 * such that the DoA can be estimated in a spatially-constrained region; thus
 * mitigating the effect of interferes and reflections arriving from other
 * directions. The DoA is estimated per sector for each frequency band.
 * The algorithms within sldoa were developed in collaboration with Symeon
 * Delikaris-Manias and Angelo Farina, and are explained in more detail in [1].
 * Dependencies:
 *     saf_utilities, afSTFTlib, saf_vbap, saf_sh
 * Author, date created:
 *     Leo McCormack, 18.10.2017
 *
 * [1] McCormack, L., Delikaris-Manias, S., Farina, A., Pinardi, D., and Pulkki,
 *     V., “Real-time conversion of sensor array signals into spherical harmonic
 *     signals with applications to spatially localised sub-band sound-field
 *     analysis,” in Audio Engineering Society Convention 144, Audio Engineering
 *     Society, 2018.
 */

#ifndef __SLDOA_H_INCLUDED__
#define __SLDOA_H_INCLUDED__

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */
    
/* ========================================================================== */
/*                             Presets + Constants                            */
/* ========================================================================== */

#define SLDOA_MAX_NUM_INPUT_CHANNELS ( 64 )
    
/* Microphone/Hydrophone array options */
#define ENABLE_ZYLIA_MIC_PRESET
#define ENABLE_EIGENMIKE32_MIC_PRESET
#define ENABLE_DTU_MIC_MIC_PRESET

/*
 * Enum: MASTER_ORDERS
 * -------------------
 * "Master order" relates to the current maximum order to expect. However, the
 * analysis order can be lower for a given frequency, due to the
 * "analysisOrderPerBand" vector, which can contain lower values than the
 * master order, but not higher.
 *
 * Options:
 *     MASTER_ORDER_FIRST   - First-order analysis (4 channel input)
 *     MASTER_ORDER_SECOND  - Second-order analysis (9 channel input)
 *     MASTER_ORDER_THIRD   - Third-order analysis (16 channel input)
 *     MASTER_ORDER_FOURTH  - Fourth-order analysis (25 channel input)
 *     MASTER_ORDER_FIFTH   - Fifth-order analysis (36 channel input)
 *     MASTER_ORDER_SIXTH   - Sixth-order analysis (49 channel input)
 *     MASTER_ORDER_SEVENTH - Seventh-order analysis (64 channel input)
 */
typedef enum _MASTER_ORDERS{
    MASTER_ORDER_FIRST = 1,
    MASTER_ORDER_SECOND,
    MASTER_ORDER_THIRD,
    MASTER_ORDER_FOURTH,
    MASTER_ORDER_FIFTH,
    MASTER_ORDER_SIXTH,
    MASTER_ORDER_SEVENTH
    
}MASTER_ORDERS;
    
typedef enum _MIC_PRESETS{
    MIC_PRESET_IDEAL = 1
#ifdef ENABLE_ZYLIA_MIC_PRESET
    , MIC_PRESET_ZYLIA
#endif
#ifdef ENABLE_EIGENMIKE32_MIC_PRESET
    , MIC_PRESET_EIGENMIKE32
#endif
#ifdef ENABLE_DTU_MIC_MIC_PRESET
    , MIC_PRESET_DTU_MIC
#endif
}MIC_PRESETS;
    
/*
 * Enum: _CH_ORDER
 * ---------------
 * Available Ambisonic channel ordering conventions
 * Note: CH_FUMA only supported for 1st order input.
 * Further note: FuMa: CH_FUMA+NORM_FUMA, AmbiX: CH_ACN+NORM_SN3D
 *
 * Options:
 *     CH_ACN  - Ambisonic Channel Numbering (ACN)
 *     CH_FUMA - (Obsolete) Furse-Malham/B-format (WXYZ)
 */
typedef enum _CH_ORDER{
    CH_ACN = 1,
    CH_FUMA     /* first-order only */
    
}CH_ORDER;

/*
 * Enum: NORM_TYPES
 * ---------------
 * Available Ambisonic normalisation conventions
 * Note: NORM_FUMA only supported for 1st order input and does NOT have the
 * 1/sqrt(2) scaling on the omni.
 * Further note: FuMa: CH_FUMA+NORM_FUMA, AmbiX: CH_ACN+NORM_SN3D
 *
 * Options:
 *     NORM_N3D  - orthonormalised (N3D)
 *     NORM_SN3D - Schmidt semi-normalisation (SN3D)
 *     NORM_FUMA - (Obsolete) Same as NORM_SN3D for 1st order
 */
typedef enum _NORM_TYPES{
    NORM_N3D = 1,
    NORM_SN3D,
    NORM_FUMA   /* first-order only */
}NORM_TYPES;
    
/*
 * Enum: CODEC_STATUS
 * ------------------
 * Current status of the codec.
 *
 * Options:
 *     CODEC_STATUS_INITIALISED     - Codec is initialised and ready to process
 *                                    input audio.
 *     CODEC_STATUS_NOT_INITIALISED - Codec has not yet been initialised, or
 *                                    the codec configuration has changed. Input
 *                                    audio should not be processed.
 *     CODEC_STATUS_INITIALISING    - Codec is currently being initialised,
 *                                    input audio should not be processed.
 */
typedef enum _CODEC_STATUS{
    CODEC_STATUS_INITIALISED = 0,
    CODEC_STATUS_NOT_INITIALISED,
    CODEC_STATUS_INITIALISING
}CODEC_STATUS;

#define SLDOA_PROGRESSBARTEXT_CHAR_LENGTH ( 256 )
    

/* ========================================================================== */
/*                               Main Functions                               */
/* ========================================================================== */

/*
 * Function: sldoa_create
 * ----------------------
 * Creates an instance of the sldoa
 *
 * Input Arguments:
 *     phSld - & address of sldoa handle
 */
void sldoa_create(void** const phSld);

/*
 * Function: sldoa_destroy
 * -----------------------
 * Destroys an instance of the sldoa
 *
 * Input Arguments:
 *     phSld - & address of sldoa handle
 */
void sldoa_destroy(void** const phSld);

/*
 * Function: sldoa_init
 * --------------------
 * Initialises an instance of sldoa with default settings
 *
 * Input Arguments:
 *     hSld       - sldoa handle
 *     samplerate - host samplerate.
 */
void sldoa_init(void* const hSld,
                float samplerate);
    
/*
 * Function: sldoa_initCodec
 * -------------------------
 * Intialises the codec variables, based on current global/user parameters
 *
 * Input Arguments:
 *     hSld - sldoa handle
 */
void sldoa_initCodec(void* const hSld);
    
/*
 * Function: sldoa_process
 * -----------------------
 * Applies the spatially-localised active-intensity based direction-of-arrival
 * estimator (SLDoA) onto the input signals [1].
 *
 * Input Arguments:
 *     hSld      - sldoa handle
 *     inputs    - input channel buffers; 2-D array: nInputs x nSamples
 *     nInputs   - number of input channels
 *     nSamples  - number of samples in 'inputs'/'output' matrices
 *     isPlaying - flag to say if there is audio in the input buffers, 0: no
 *                 audio, reduced processing, 1: audio, full processing
 *
 * [1] McCormack, L., Delikaris-Manias, S., Farina, A., Pinardi, D., and Pulkki,
 *     V., “Real-time conversion of sensor array signals into spherical harmonic
 *     signals with applications to spatially localised sub-band sound-field
 *     analysis,” in Audio Engineering Society Convention 144, Audio Engineering
 *     Society, 2018.
 */
void sldoa_analysis(void* const hSld,
                    float** const inputs,
                    int nInputs,
                    int nSamples,
                    int isPlaying);


/* ========================================================================== */
/*                                Set Functions                               */
/* ========================================================================== */

/*
 * Function: sldoa_setMasterOrder
 * ------------------------------
 * Sets the maximum input/analysis order
 *
 * Input Arguments:
 *     hSld     - sldoa handle
 *     newValue - new maximum analysis order (see 'MASTER_ORDERS' enum)
 */
void sldoa_setMasterOrder(void* const hSld,  int newValue);

/*
 * Function: sldoa_refreshSettings
 * -------------------------------
 * Sets all intialisation flags to 1. i.e. re-initialise all settings/variables
 * as sldoa is currently configured, at next available opportunity.
 *
 * Input Arguments:
 *     hSld - sldoa handle
 */
void sldoa_refreshSettings(void* const hSld);
   
/*
 * Function: sldoa_setMaxFreq
 * --------------------------
 * Sets the maximum analysis frequency
 *
 * Input Arguments:
 *     hSld    - sldoa handle
 *     newFreq - new maximum analysis frequency, in Hz
 */
void sldoa_setMaxFreq(void* const hSld, float newFreq);

/*
 * Function: sldoa_setMinFreq
 * --------------------------
 * Sets the minimum analysis frequency
 *
 * Input Arguments:
 *     hSld    - sldoa handle
 *     newFreq - new minimum analysis frequency, in Hz
 */
void sldoa_setMinFreq(void* const hSld, float newFreq);
    
/*
 * Function: sldoa_setAvg
 * ----------------------
 * Sets the DoA averaging coefficient
 *
 * Input Arguments:
 *     hSld   - sldoa handle
 *     newAvg - new averaging coefficient
 */
void sldoa_setAvg(void* const hSld, float newAvg);

/*
 * Function: sldoa_setAnaOrder
 * ---------------------------
 * Sets the input/analysis order for one specific frequency band.
 *
 * Input Arguments:
 *     hSld     - sldoa handle
 *     newValue - new analysis order
 *     bandIdx  - frequency band index
 */
void sldoa_setAnaOrder(void* const hSld,  int newValue, int bandIdx);

/*
 * Function: sldoa_setAnaOrderAllBands
 * -----------------------------------
 * Sets the input/analysis order for all frequency bands.
 *
 * Input Arguments:
 *     hSld     - sldoa handle
 *     newValue - new analysis order
 */
void sldoa_setAnaOrderAllBands(void* const hSld,  int newValue);

/*
 * Function: sldoa_setChOrder
 * --------------------------
 * Sets the Ambisonic channel ordering convention to decode with, in order to
 * match the convention employed by the input signals
 *
 * Input Arguments:
 *     hSld     - sldoa handle
 *     newOrder - convention to use (see 'CH_ORDER' enum)
 */
void sldoa_setChOrder(void* const hSld, int newOrder);

/*
 * Function: sldoa_setNormType
 * ---------------------------
 * Sets the Ambisonic normalisation convention to decode with, in order to match
 * with the convention employed by the input signals.
 *
 * Input Arguments:
 *     hSld    - sldoa handle
 *     newType - convention to use (see 'NORM_TYPE' enum)
 */
void sldoa_setNormType(void* const hSld, int newType);

/*
 * Function: sldoa_setSourcePreset
 * -------------------------------
 * Sets an input preset, i.e. the microphone/hyrophone array used to capture
 * the input signals.
 *
 * Input Arguments:
 *     hSld        - sldoa handle
 *     newPresetID - preset to use (see 'MIC_PRESETS' enum)
 */
void sldoa_setSourcePreset(void* const hSld, int newPresetID);


/* ========================================================================== */
/*                                Get Functions                               */
/* ========================================================================== */

/*
 * Function: sldoa_getCodecStatus
 * ------------------------------
 * Returns current codec status.
 *
 * Input Arguments:
 *     hSld - sldoa handle
 * Returns:
 *     codec status (see 'CODEC_STATUS' enum)
 */
CODEC_STATUS sldoa_getCodecStatus(void* const hSld);

/*
 * Function: sldoa_getProgressBar0_1
 * ---------------------------------
 * (Optional) Returns current intialisation/processing progress, between 0..1
 * 0: intialisation/processing has started
 * 1: intialisation/processing has ended
 *
 * Input Arguments:
 *     hSld - sldoa handle
 * Returns:
 *     current progress, 0..1
 */
float sldoa_getProgressBar0_1(void* const hSld);

/*
 * Function: sldoa_getProgressBarText
 * ----------------------------------
 * (Optional) Returns current intialisation/processing progress text
 * Note: "text" string should be (at least) of length:
 *     SLDOA_PROGRESSBARTEXT_CHAR_LENGTH
 *
 * Input Arguments:
 *     hSld - sldoa handle
 * Output Arguments:
 *     text - process bar text; SLDOA_PROGRESSBARTEXT_CHAR_LENGTH x 1
 */
void sldoa_getProgressBarText(void* const hSld, char* text);

/*
 * Function: sldoa_getMasterOrder
 * ------------------------------
 * Returns the current maximum analysis/input order.
 *
 * Input Arguments:
 *     hSld - sldoa handle
 * Returns:
 *     current maximum analysis order (see 'MASTER_ORDERS' enum)
 */
int sldoa_getMasterOrder(void* const hSld);

/*
 * Function: sldoa_getSamplingRate
 * -------------------------------
 * Returns the current sampling rate
 *
 * Input Arguments:
 *     hSld - sldoa handle
 *  Returns:
 *     current sampling rate, in Hz
 */
int sldoa_getSamplingRate(void* const hSld);

/*
 * Function: sldoa_getMaxFreq
 * --------------------------
 * Returns the maximum analysis frequency
 *
 * Input Arguments:
 *     hSld - sldoa handle
 *  Returns:
 *     current maximum analysis frequency, in Hz
 */
float sldoa_getMaxFreq(void* const hSld);

/*
 * Function: sldoa_getMinFreq
 * --------------------------
 * Returns the minimum analysis frequency
 *
 * Input Arguments:
 *     hSld - sldoa handle
 *  Returns:
 *     current minimum analysis frequency, in Hz
 */
float sldoa_getMinFreq(void* const hSld);

/*
 * Function: sldoa_getAvg
 * ----------------------
 * Returns the current DoA averaging coefficient value
 *
 * Input Arguments:
 *     hSld - sldoa handle
 *  Returns:
 *     current DoA averaging coefficient.
 */
float sldoa_getAvg(void* const hSld);

/*
 * Function: sldoa_getNumberOfBands
 * --------------------------------
 * Returns the number frequency bands employed by sldoa
 *
 *  Returns:
 *     the number frequency bands employed by sldoa
 */
int sldoa_getNumberOfBands(void);

/*
 * Function: sldoa_getNSHrequired
 * ------------------------------
 * Returns the number of spherical harmonic signals required by the current
 * analysis order i.e. (current_order + 1)^2
 *
 * Input Arguments:
 *     hSld - sldoa handle
 * Returns:
 *     number of required spherical harmonic signals required by current
 *     analysis order
 */
int sldoa_getNSHrequired(void* const hSld);
    
/*
 * Function: sldoa_getDisplayData
 * ------------------------------
 * Returns the analysis output data. Including the DoAs per frequency, and per
 * sector, accompanied by colour coefficients (red: high frequencies, blue:
 * low frequencies), and alpha coefficients (more opaque: higher energy, more
 * transpararent: less energy).
 *
 * Input Arguments:
 *     hSld - sldoa handle
 * Output Arguments:
 *     pAzi_deg         - & azimuth of estimated DoAs;
 *                        FLAT: pNsectorsPerBand x sldoa_getNumberOfBands()
 *     pElev_deg        - & elevation of estimated DoAs;
 *                        FLAT: pNsectorsPerBand x sldoa_getNumberOfBands()
 *     pColourScale     - & colour scale, 0..1, 1:red, 0: blue
 *                        FLAT: pNsectorsPerBand x sldoa_getNumberOfBands()
 *     pAlphaScale      - & alpha scale, 0..1, 1: opaque, 0: transparent;
 *                        FLAT: pNsectorsPerBand x sldoa_getNumberOfBands()
 *     pNsectorsPerBand - & number of sectors per frequency;
 *                        pNsectorsPerBand x 1
 *     pMaxNumSectors   - & maximum number of sectors
 *     pStartBand       - & band index corresponding to lowest frequency
 *     pEndBand         - & band index corresponding to highest frequency
 */
void sldoa_getDisplayData(void *  const hSld,
                          float** pAzi_deg,
                          float** pElev_deg,
                          float** pColourScale,
                          float** pAlphaScale,
                          int** pNsectorsPerBand,
                          int* pMaxNumSectors,
                          int* pStartBand,
                          int* pEndBand);
    
/*
 * Function: sldoa_getAnaOrder
 * ---------------------------
 * Returns the input/analysis order for one specific frequency band.
 *
 * Input Arguments:
 *     hSld    - sldoa handle
 *     bandIdx - frequency band index
 * Returns:
 *     current analysis order at this band index
 */
int sldoa_getAnaOrder(void* const hSld, int bandIdx);

/*
 * Function: sldoa_getAnaOrderAllBands
 * -----------------------------------
 * Returns the input/analysis order for the first frequency band
 *
 * Input Arguments:
 *     hSld - sldoa handle
 * Returns:
 *     current analysis order
 */
int sldoa_getAnaOrderAllBands(void* const hSld);

/*
 * Function: sldoa_getAnaOrderHandle
 * ---------------------------------
 * Returns the input/analysis order for all frequency bands
 *
 * Input Arguments:
 *     hSld      - sldoa handle
 * Output Arguments:
 *     pX_vector - & frequency vector; pNpoints x 1
 *     pY_values - & input/analysis orders; pNpoints x 1
 *     pNpoints  - & number of frequency bands
 */
void sldoa_getAnaOrderHandle(void* const hSld,
                             float** pX_vector,
                             int** pY_values,
                             int* pNpoints);

/*
 * Function: sldoa_getChOrder
 * --------------------------
 * Returns the Ambisonic channel ordering convention currently being used to
 * decode with, which should match the convention employed by the input signals
 *
 * Input Arguments:
 *     hSld - sldoa handle
 * Returns:
 *     convention currently being used (see 'CH_ORDER' enum)
 */
int sldoa_getChOrder(void* const hSld);

/*
 * Function: sldoa_getNormType
 * ---------------------------
 * Returns the Ambisonic normalisation convention currently being usedto decode
 * with, which should match the convention employed by the input signals.
 *
 * Input Arguments:
 *     hSld - sldoa handle
 * Returns:
 *     convention currently being used (see 'NORM_TYPE' enum)
 */
int sldoa_getNormType(void* const hSld);

    
#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */

#endif /* __SLDOA_H_INCLUDED__ */
