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
 * Filename: ambi_dec.h (include header)
 * -------------------------------------
 * A frequency-dependent Ambisonic decoder for loudspeakers or headphones.
 * Different decoder settings can be specified for the low and high frequencies.
 * When utilising spherical harmonic signals derived from real microphone
 * arrays, this implementation also allows the decoding order per frequency band
 * to be specified. Optionally, a SOFA file may be loaded for personalised
 * headphone listening.
 * The algorithms utilised in this Ambisonic decoder were pieced together and
 * developed in collaboration with Archontis Politis.
 *
 * Dependencies:
 *     saf_utilities, afSTFTlib, saf_hoa, saf_vbap, saf_hrir, saf_sh,
 *     saf_sofa_reader
 * Author, date created:
 *     Leo McCormack, 07.12.2017
 */

#ifndef __AMBI_DEC_H_INCLUDED__
#define __AMBI_DEC_H_INCLUDED__

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */
    
/* ========================================================================== */
/*                             Presets + Constants                            */
/* ========================================================================== */
    
/*
 * Enum: MASTER_ORDERS
 * -------------------
 * Available decoding (upper) orders. However, the decoding order may be lower
 * than the "master" order for any given frequency band.
 *
 * Options:
 *     MASTER_ORDER_FIRST   - First-order decoding (4 channel input)
 *     MASTER_ORDER_SECOND  - Second-order decoding (9 channel input)
 *     MASTER_ORDER_THIRD   - Third-order decoding (16 channel input)
 *     MASTER_ORDER_FOURTH  - Fourth-order decoding (25 channel input)
 *     MASTER_ORDER_FIFTH   - Fifth-order decoding (36 channel input)
 *     MASTER_ORDER_SIXTH   - Sixth-order decoding (49 channel input)
 *     MASTER_ORDER_SEVENTH - Seventh-order decoding (64 channel input)
 */
#define AMBI_DEC_MAX_SH_ORDER ( 7 )
typedef enum _MASTER_ORDERS{
    MASTER_ORDER_FIRST = 1,
    MASTER_ORDER_SECOND,
    MASTER_ORDER_THIRD,
    MASTER_ORDER_FOURTH,
    MASTER_ORDER_FIFTH,
    MASTER_ORDER_SIXTH,
    MASTER_ORDER_SEVENTH
    
}MASTER_ORDERS;
    
/*
 * Enum: DECODING_METHODS
 * ----------------------
 * Available decoding methods. See "saf_hoa_internal.h" for a more indepth
 * description of each approach.
 *
 * Options:
 *     DECODING_METHOD_SAD    - Sampling Ambisonic Decoder (SAD)
 *     DECODING_METHOD_MMD    - Mode-Matching Decoder (MMD)
 *     DECODING_METHOD_EPAD   - Energy-Preserving Ambisonic Decoder (EPAD)
 *     DECODING_METHOD_ALLRAD - All-Round Ambisonic Decoder (AllRAD)
 */
#define AMBI_DEC_NUM_DECODING_METHODS ( 4 )
typedef enum _DECODING_METHODS {
    DECODING_METHOD_SAD = 1,
    DECODING_METHOD_MMD,
    DECODING_METHOD_EPAD,
    DECODING_METHOD_ALLRAD
    
} DECODING_METHODS;
    
/*
 * Enum: MIC_PRESETS
 * -----------------
 * Available microphone array presets. These determine the frequency ranges
 * where the microphone array provides usable spherical harmonic components
 * at each order.
 */
typedef enum _MIC_PRESETS{
    MIC_PRESET_IDEAL = 1,
    MIC_PRESET_ZYLIA,
    MIC_PRESET_EIGENMIKE32,
    MIC_PRESET_DTU_MIC
    
}MIC_PRESETS;
   
/*
 * Enum: LOUDSPEAKER_ARRAY_PRESETS
 * -------------------------------
 * Available loudspeaker array presets
 */
typedef enum _LOUDSPEAKER_ARRAY_PRESETS{
    LOUDSPEAKER_ARRAY_PRESET_DEFAULT = 1,
    LOUDSPEAKER_ARRAY_PRESET_5PX,
    LOUDSPEAKER_ARRAY_PRESET_7PX,
    LOUDSPEAKER_ARRAY_PRESET_8PX,
    LOUDSPEAKER_ARRAY_PRESET_9PX,
    LOUDSPEAKER_ARRAY_PRESET_10PX,
    LOUDSPEAKER_ARRAY_PRESET_11PX,
    LOUDSPEAKER_ARRAY_PRESET_11PX_7_4,
    LOUDSPEAKER_ARRAY_PRESET_13PX,
    LOUDSPEAKER_ARRAY_PRESET_22PX,
    LOUDSPEAKER_ARRAY_PRESET_AALTO_MCC,
    LOUDSPEAKER_ARRAY_PRESET_AALTO_APAJA,
    LOUDSPEAKER_ARRAY_PRESET_AALTO_LR,
    LOUDSPEAKER_ARRAY_PRESET_DTU_AVIL,
    LOUDSPEAKER_ARRAY_PRESET_ZYLIA_LAB,
    LOUDSPEAKER_ARRAY_PRESET_T_DESIGN_4,
    LOUDSPEAKER_ARRAY_PRESET_T_DESIGN_12,
    LOUDSPEAKER_ARRAY_PRESET_T_DESIGN_24,
    LOUDSPEAKER_ARRAY_PRESET_T_DESIGN_36,
    LOUDSPEAKER_ARRAY_PRESET_T_DESIGN_48,
    LOUDSPEAKER_ARRAY_PRESET_T_DESIGN_60
    
}LOUDSPEAKER_ARRAY_PRESETS;

/*
 * Enum: DIFFUSE_FIELD_EQ_APPROACH
 * -------------------------------
 * When using mixed order decoding (i.e. different decoding orders for
 * different frequencies), this equalisation helps maintain equal perceived
 * "loudness". At low frequencies, preserving amplitude is more favourable,
 * whereas for high frequencies, preserving energy is better.
 *
 * Options:
 *     AMPLITUDE_PRESERVING - preserve omni amplitude
 *     ENERGY_PRESERVING    - preserve omni energy
 */
typedef enum _DIFFUSE_FIELD_EQ_APPROACH {
    AMPLITUDE_PRESERVING=1,
    ENERGY_PRESERVING
    
} DIFFUSE_FIELD_EQ_APPROACH;

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
#define AMBI_DEC_NUM_CH_ORDERINGS ( 2 )
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
#define AMBI_DEC_NUM_NORM_TYPES ( 3 )
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
    
#define AMBI_DEC_MAX_NUM_OUTPUTS ( 64 )
#define AMBI_DEC_TRANSITION_MIN_VALUE ( 500.0f )
#define AMBI_DEC_TRANSITION_MAX_VALUE ( 2000.0f )
#define AMBI_DEC_PROGRESSBARTEXT_CHAR_LENGTH 256
    

/* ========================================================================== */
/*                               Main Functions                               */
/* ========================================================================== */

/*
 * Function: ambi_dec_create
 * -------------------------
 * Creates an instance of the ambi_dec
 *
 * Input Arguments:
 *     phAmbi - & address of ambi_dec handle
 */
void ambi_dec_create(void** const phAmbi);
    
/*
 * Function: ambi_dec_destroy
 * --------------------------
 * Destroys an instance of the ambi_dec
 *
 * Input Arguments:
 *     phAmbi - & address of ambi_dec handle
 */
void ambi_dec_destroy(void** const phAmbi);

/*
 * Function: ambi_dec_init
 * -----------------------
 * Initialises an instance of ambi_dec with default settings
 *
 * Input Arguments:
 *     hAmbi      - ambi_dec handle
 *     samplerate - host samplerate.
 */
void ambi_dec_init(void* const hAmbi,
                   int samplerate);

/*
 * Function: ambi_dec_initCodec
 * ----------------------------
 * Intialises the codec variables, based on current global/user parameters
 *
 * Input Arguments:
 *     hAmbi - ambi_dec handle
 */
void ambi_dec_initCodec(void* const hAmbi);

/*
 * Function: ambi_dec_process
 * --------------------------
 * Decodes input spherical harmonic signals to the loudspeaker channels.
 *
 * Input Arguments:
 *     hAmbi     - ambi_dec handle
 *     inputs    - input channel buffers; 2-D array: nInputs x nSamples
 *     outputs   - Output channel buffers; 2-D array: nOutputs x nSamples
 *     nInputs   - number of input channels
 *     nOutputs  - number of output channels
 *     nSamples  - number of samples in 'inputs'/'output' matrices
 */
void ambi_dec_process(void* const hAmbi,
                      float** const inputs,
                      float** const outputs,
                      int nInputs,
                      int nOutputs,
                      int nSamples);

    
/* ========================================================================== */
/*                                Set Functions                               */
/* ========================================================================== */

/*
 * Function: ambi_dec_refreshSettings
 * ----------------------------------
 * Sets all intialisation flags to 1. i.e. re-initialise all settings/variables
 * as ambi_dec is currently configured, at next available opportunity.
 *
 * Input Arguments:
 *     hAmbi - ambi_dec handle
 */
void ambi_dec_refreshSettings(void* const hAmbi);
    
/*
 * Function: ambi_dec_setMasterDecOrder
 * ------------------------------------
 * Sets the master decoding order. However, the decoding order may be lower than
 * this for any given frequency, this is just the maximum.
 * If decoding order is higher than the input signal order, the extra required
 * channels are filled with zeros. If the decoding order is lower than the input
 * signal order, the number input signals is truncated accordingly.
 *
 * Input Arguments:
 *     hAmbi    - ambi_dec handle
 *     newValue - new decoding order (see 'MASTER_ORDERS' enum)
 */
void ambi_dec_setMasterDecOrder(void* const hAmbi,  int newValue);
    
/*
 * Function: ambi_dec_setDecOrder
 * ------------------------------
 * Sets the decoding order for a given frequency band.
 * Note: the maximum order is dictated by "ambi_dec_setMasterDecOrder"
 *
 * Input Arguments:
 *     hAmbi    - ambi_dec handle
 *     newValue - new decoding order
 *     bandIdx  - frequency band index
 */
void ambi_dec_setDecOrder(void* const hAmbi,  int newValue, int bandIdx);

/*
 * Function: ambi_dec_setDecOrderAllBands
 * --------------------------------------
 * Sets the decoding order for all frequency bands
 * Note: the maximum order is dictated by "ambi_dec_setMasterDecOrder"
 *
 * Input Arguments:
 *     hAmbi    - ambi_dec handle
 *     newValue - new decoding order
 */
void ambi_dec_setDecOrderAllBands(void* const hAmbi,  int newValue);

/*
 * Function: ambi_dec_setLoudspeakerAzi_deg
 * ----------------------------------------
 * Sets the azimuth of a specific loudspeaker
 *
 * Input Arguments:
 *     hAmbi      - ambi_dec handle
 *     index      - loudspeaker index
 *     newAzi_deg - new azimuth in DEGREES
 */
void ambi_dec_setLoudspeakerAzi_deg(void* const hAmbi, int index, float newAzi_deg);

/*
 * Function: ambi_dec_setLoudspeakerElev_deg
 * -----------------------------------------
 * Sets the elevation of a specific loudspeaker
 *
 * Input Arguments:
 *     hAmbi       - ambi_dec handle
 *     index       - loudspeaker index
 *     newElev_deg - new elevation in DEGREES
 */
void ambi_dec_setLoudspeakerElev_deg(void* const hAmbi, int index, float newElev_deg);

/*
 * Function: ambi_dec_setNumLoudspeakers
 * -------------------------------------
 * Sets the number of loudspeakers to decode to.
 *
 * Input Arguments:
 *     hAmbi             - ambi_dec handle
 *     new_nLoudspeakers - new number of loudspeakers
 */
void ambi_dec_setNumLoudspeakers(void* const hAmbi, int new_nLoudspeakers);

/*
 * Function: ambi_dec_setUseDefaultHRIRsflag
 * -----------------------------------------
 * Sets flag to dictate whether the default HRIRs in the Spatial_Audio_Framework
 * should be used, or a custom HRIR set loaded via a SOFA file.
 * Note: if the custom set failes to load correctly, ambi_dec will revert to the
 * defualt set. Use 'ambi_dec_getUseDefaultHRIRsflag' to check if loading was
 * successful.
 *
 * Input Arguments:
 *     hAmbi     - ambi_dec handle
 *     newState  - 0: use custom HRIR set, 1: use default HRIR set
 */
void ambi_dec_setUseDefaultHRIRsflag(void* const hAmbi, int newState);

/*
 * Function: ambi_dec_setBinauraliseLSflag
 * -----------------------------------------
 * Sets flag to dictate whether the output loudspeaker signals should be
 * binauralised
 *
 * Input Arguments:
 *     hAmbi     - ambi_dec handle
 *     newState  - 0: output loudspeaker signals, 1: output binaural signals
 */
void ambi_dec_setBinauraliseLSflag(void* const hAmbi, int newState);

/*
 * Function: ambi_dec_setSofaFilePath
 * ----------------------------------
 * Sets the file path for a .sofa file, in order to employ a custom HRIR set for
 * the decoding.
 * Note: if the custom set failes to load correctly, hcompass will revert to the
 * defualt set. Use 'ambi_dec_getUseDefaultHRIRsflag' to check if loading was
 * successful.
 *
 * Input Arguments:
 *     hAmbi - ambi_dec handle
 *     path  - file path to .sofa file (WITH file extension)
 */
void ambi_dec_setSofaFilePath(void* const hAmbi, const char* path);

/*
 * Function: ambi_dec_setSourcePreset
 * ----------------------------------
 * Sets the source preset. By default the decoder will decode at the maximum
 * order for all frequencies. However, in the case of spherical harmonic
 * input derived from microphone arrays, the available order is frequency
 * dependent, therefore, different bands require different decoding orders.
 * For conveinience, presets for a handful of comerically available microphone
 * array are included (see "MIC_PRESETS" enum).
 *
 * Input Arguments:
 *     hAmbi       - ambi_dec handle
 *     newPresetID - new preset (see "MIC_PRESETS" enum)
 */
void ambi_dec_setSourcePreset(void* const hAmbi, int newPresetID);

/*
 * Function: ambi_dec_setOutputConfigPreset
 * ----------------------------------------
 * Sets the output loudspeaker preset.
 * For conveinience, presets for several popular arrangements are included (see
 * "PRESETS" enum).
 *
 * Input Arguments:
 *     hAmbi       - ambi_dec handle
 *     newPresetID - new preset (see "PRESETS" enum)
 */
void ambi_dec_setOutputConfigPreset(void* const hAmbi, int newPresetID);

/*
 * Function: ambi_dec_setChOrder
 * -----------------------------
 * Sets the Ambisonic channel ordering convention to decode with, in order to
 * match the convention employed by the input signals
 *
 * Input Arguments:
 *     hAmbi    - ambi_dec handle
 *     newOrder - convention to use (see 'CH_ORDER' enum)
 */
void ambi_dec_setChOrder(void* const hAmbi, int newOrder);

/*
 * Function: ambi_dec_setNormType
 * ------------------------------
 * Sets the Ambisonic normalisation convention to decode with, in order to match
 * with the convention employed by the input signals.
 *
 * Input Arguments:
 *     hAmbi   - ambi_dec handle
 *     newType - convention to use (see 'NORM_TYPE' enum)
 */
void ambi_dec_setNormType(void* const hAmbi, int newType);

/*
 * Function: ambi_dec_setDecMethod
 * -------------------------------
 * Sets the decoding method for a specific decoder. ambi_dec employs two
 * decoders, one for low frequencies and one for high frequencies.
 * (use "ambi_dec_setTransitionFreq" to dictate the transition frequency)
 *
 * Input Arguments:
 *     hAmbi - ambi_dec handle
 *     index - index for low (0) or high (1) frequency decoder
 *     newID - new decoding method (see 'DECODING_METHODS' enum)
 */
void ambi_dec_setDecMethod(void* const hAmbi, int index, int newID);

/*
 * Function: ambi_dec_setEnableMaxRE
 * ---------------------------------
 * Sets a flag to enable/disable the max_rE weighting for one of the decoders.
 *
 * Input Arguments:
 *     hAmbi - ambi_dec handle
 *     index - index for low (0) or high (1) frequency decoder
 *     newID - 0: disabled, 1: enabled
 */
void ambi_dec_setDecEnableMaxrE(void* const hAmbi, int index, int newID);

/*
 * Function: ambi_dec_setDecNormType
 * ---------------------------------
 * Sets the equalisation approach for one of the decoders. This is used to help
 * keep the perceived loudness consistent, when using mixed decoding orders
 * (i.e. different decoding orders for different frequency bands)
 * ambi_dec either to preserves amplitude or energy for each order.
 * Note: it is suggested to preserve amplitude at low-frequencies and energy
 * at high-frequencies.
 *
 * Input Arguments:
 *     hAmbi - ambi_dec handle
 *     index - index for low (0) or high (1) frequency decoder
 *     newID - see "DIFFUSE_FIELD_EQ_APPROACH" enum
 */
void ambi_dec_setDecNormType(void* const hAmbi, int index, int newID);

/*
 * Function: ambi_dec_setTransitionFreq
 * ------------------------------------
 * Sets the frequeny at which to transition from the low frequency decoder to
 * the high frequency decoder.
 *
 * Input Arguments:
 *     hAmbi    - ambi_dec handle
 *     newValue - new transition frequency in Hz
 */
void ambi_dec_setTransitionFreq(void* const hAmbi, float newValue);

    
/* ========================================================================== */
/*                                Get Functions                               */
/* ========================================================================== */
    
/*
 * Function: ambi_dec_getCodecStatus
 * ---------------------------------
 * Returns current codec status.
 *
 * Input Arguments:
 *     hAmbi - ambi_dec handle
 * Returns:
 *     codec status (see 'CODEC_STATUS' enum)
 */
CODEC_STATUS ambi_dec_getCodecStatus(void* const hAmbi);

/*
 * Function: ambi_dec_getProgressBar0_1
 * ------------------------------------
 * (Optional) Returns current intialisation/processing progress, between 0..1
 * 0: intialisation/processing has started
 * 1: intialisation/processing has ended
 *
 * Input Arguments:
 *     hAmbi - ambi_dec handle
 * Returns:
 *     current progress, 0..1
 */
float ambi_dec_getProgressBar0_1(void* const hAmbi);

/*
 * Function: ambi_dec_getProgressBarText
 * -------------------------------------
 * (Optional) Returns current intialisation/processing progress text
 * Note: "text" string should be (at least) of length:
 *     AMBI_DEC_PROGRESSBARTEXT_CHAR_LENGTH
 *
 * Input Arguments:
 *     hAmbi - ambi_dec handle
 * Output Arguments:
 *     text  - process bar text; AMBI_DEC_PROGRESSBARTEXT_CHAR_LENGTH x 1
 */
void ambi_dec_getProgressBarText(void* const hAmbi, char* text);

/*
 * Function: ambi_dec_getMasterDecOrder
 * ------------------------------------
 * Returns the master/maximum decoding order.
 *
 * Input Arguments:
 *     hAmbi - ambi_dec handle
 * Returns:
 *     current maximum decoding order (see 'MASTER_ORDERS' enum)
 */
int ambi_dec_getMasterDecOrder(void* const hAmbi);
    
/*
 * Function: ambi_dec_getDecOrder
 * ------------------------------
 * Returns the decoding order for a given frequency band
 *
 * Input Arguments:
 *     hAmbi   - ambi_dec handle
 *     bandIdx - frequency band index
 * Returns:
 *     current maximum decoding order (see 'MASTER_ORDERS' enum)
 */
int ambi_dec_getDecOrder(void* const hAmbi, int bandIdx);
  
/*
 * Function: ambi_dec_getDecOrderAllBands
 * --------------------------------------
 * Returns the decoding order for the first band.
 *
 * Input Arguments:
 *     hAmbi   - ambi_dec handle
 * Returns:
 *     current maximum decoding order (see 'MASTER_ORDERS' enum)
 */
int ambi_dec_getDecOrderAllBands(void* const hAmbi);

/*
 * Function: ambi_dec_getDecOrderHandle
 * ------------------------------------
 * Returns handles for the decoding orders and frequency vector.
 *
 * Input Arguments:
 *     hAmbi     - ambi_dec handle
 *     pX_vector - & frequency vector; pNpoints x 1
 *     pY_values - & decoding order per frequency; pNpoints x 1
 *     pNpoints  - & number of grid points.
 */
void ambi_dec_getDecOrderHandle(void* const hAmbi,
                                float** pX_vector,
                                int** pY_values,
                                int* pNpoints);

/*
 * Function: ambi_dec_getNumberOfBands
 * -----------------------------------
 * Returns the number of frequency bands employed by ambi_dec.
 *
 * Returns:
 *     The number of frequency bands
 */
int ambi_dec_getNumberOfBands(void);

/*
 * Function: ambi_dec_getLoudspeakerAzi_deg
 * ----------------------------------------
 * Returns the loudspeaker azimuth for a given index.
 *
 * Input Arguments:
 *     hAmbi - ambi_dec handle
 *     index - loudspeaker index
 * Returns:
 *     loudspeaker azimuth in DEGREES
 */
float ambi_dec_getLoudspeakerAzi_deg(void* const hAmbi, int index);

/*
 * Function: ambi_dec_getLoudspeakerElev_deg
 * -----------------------------------------
 * Returns the loudspeaker elevation for a given index.
 *
 * Input Arguments:
 *     hAmbi - ambi_dec handle
 *     index - loudspeaker index
 * Returns:
 *     loudspeaker elevation in DEGREES
 */
float ambi_dec_getLoudspeakerElev_deg(void* const hAmbi, int index);

/*
 * Function: ambi_dec_getNumLoudspeakers
 * -------------------------------------
 * Returns the number of loudspeakers in the current layout
 *
 * Input Arguments:
 *     hAmbi - ambi_dec handle
 * Returns:
 *     number of loudspeakers
 */
int ambi_dec_getNumLoudspeakers(void* const hAmbi);

/*
 * Function: ambi_dec_getMaxNumLoudspeakers
 * ----------------------------------------
 * Returns the maximum number of loudspeakers supported by ambi_dec
 *
 * Returns:
 *     maximum number of loudspeakers
 */
int ambi_dec_getMaxNumLoudspeakers(void);

/*
 * Function: ambi_dec_getNSHrequired
 * ---------------------------------
 * Returns the number of spherical harmonic signals required by the current
 * decoding order i.e. (current_order + 1)^2
 *
 * Input Arguments:
 *     hAmbi - ambi_dec handle
 * Returns:
 *     number of required spherical harmonic signals required by current
 *     decoding order
 */
int  ambi_dec_getNSHrequired(void* const hAmbi); 

/*
 * Function: ambi_dec_getUseDefaultHRIRsflag
 * -----------------------------------------
 * Returns the value of a flag used to dictate whether the default HRIRs in the
 * Spatial_Audio_Framework should be used, or a custom HRIR set loaded via a
 * SOFA file.
 * Note: if the custom set failes to load correctly, ambi_dec will revert to the
 * defualt set.
 *
 * Input Arguments:
 *     hAmbi - ambi_dec handle
 * Returns:
 *     0: use custom HRIR set, 1: use default HRIR set
 */
int ambi_dec_getUseDefaultHRIRsflag(void* const hAmbi);

/*
 * Function: ambi_dec_getBinauraliseLSflag
 * ---------------------------------------
 * Returns the value of a flag used to dictate whether the loudspeaker signals
 * should be binauralised.
 *
 * Input Arguments:
 *     hAmbi - ambi_dec handle
 * Returns:
 *     0: output loudspeaker signals, 1: output binaural signals.
 */
int ambi_dec_getBinauraliseLSflag(void* const hAmbi);

/*
 * Function: ambi_dec_getSofaFilePath
 * ----------------------------------
 * Returns the file path for a .sofa file.
 * Note: if the custom set failes to load correctly, hcompass will revert to the
 * defualt set. Use 'ambi_dec_getUseDefaultHRIRsflag' to check if loading was
 * successful.
 *
 * Input Arguments:
 *     hAmbi - ambi_dec handle
 * Returns:
 *      file path to .sofa file (WITH file extension)
 */
char* ambi_dec_getSofaFilePath(void* const hAmbi);

/*
 * Function: ambi_dec_getChOrder
 * -----------------------------
 * Returns the Ambisonic channel ordering convention currently being used to
 * decode with, which should match the convention employed by the input signals
 *
 * Input Arguments:
 *     hAmbi - ambi_dec handle
 * Returns:
 *     convention currently being used (see 'CH_ORDER' enum)
 */
int ambi_dec_getChOrder(void* const hAmbi);

/*
 * Function: ambi_dec_getNormType
 * ------------------------------
 * Returns the Ambisonic normalisation convention currently being usedto decode
 * with, which should match the convention employed by the input signals.
 *
 * Input Arguments:
 *     hAmbi - ambi_dec handle
 * Returns:
 *     convention currently being used (see 'NORM_TYPE' enum)
 */
int ambi_dec_getNormType(void* const hAmbi);

/*
 * Function: ambi_dec_getDecodingMethod
 * ------------------------------------
 * Returns the currently selected decoding method
 *
 * Input Arguments:
 *     hAmbi - ambi_dec handle
 * Returns:
 *     current decoding order (see 'DECODING_METHODS' enum)
 */
int ambi_dec_getDecMethod(void* const hAmbi, int index);
    
/*
 * Function: ambi_dec_getDecEnableMaxrE
 * ------------------------------------
 * Returns the value of a flag used to dictate whether the max_rE weighting is
 * being applied by a given decoder
 *
 * Input Arguments:
 *     hAmbi - ambi_dec handle
 *     index - index for low (0) or high (1) frequency decoder
 * Returns
 */
int ambi_dec_getDecEnableMaxrE(void* const hAmbi, int index);

/*
 * Function: ambi_dec_getDecNormType
 * ---------------------------------
 * Returns the current equalisation approach for one of the decoders.
 * Note: it is suggested to preserve amplitude at low-frequencies and energy
 * at high-frequencies.
 *
 * Input Arguments:
 *     hAmbi - ambi_dec handle
 *     index - index for low (0) or high (1) frequency decoder
 * Returns
 *     see "DIFFUSE_FIELD_EQ_APPROACH" enum
 */
int ambi_dec_getDecNormType(void* const hAmbi, int index);

/*
 * Function: ambi_dec_getTransitionFreq
 * ------------------------------------
 * Returns the frequency at which to transition from the low frequency decoder
 * to the high frequency decoder.
 *
 * Input Arguments:
 *     hAmbi    - ambi_dec handle
 * Returns
 *     transition frequency in Hz
 */
float ambi_dec_getTransitionFreq(void* const hAmbi);
    
/*
 * Function: ambi_dec_getHRIRsamplerate
 * ------------------------------------
 * Returns the HRIR sample rate
 *
 * Input Arguments:
 *     hAmbi - ambi_dec handle
 * Returns:
 *     HRIR sampling rate
 */
int ambi_dec_getHRIRsamplerate(void* const hAmbi);

/*
 * Function: ambi_dec_getDAWsamplerate
 * -----------------------------------
 * Returns the DAW/Host sample rate
 *
 * Input Arguments:
 *     hAmbi - ambi_dec handle
 * Returns:
 *     DAW/Host sampling rate
 */
int ambi_dec_getDAWsamplerate(void* const hAmbi);
    
/*
 * Function: ambi_dec_getProcessingDelay
 * -------------------------------------
 * Returns the processing delay in samples. May be used for delay compensation
 * features
 *
 * Returns:
 *     processing delay in samples
 */
int ambi_dec_getProcessingDelay(void);
    
    
#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */

#endif /* __AMBIDEC_H_INCLUDED__ */
