/*
 * Copyright 2018 Leo McCormack
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
 * Filename: ambi_bin.h (include header)
 * -------------------------------------
 * A binaural Ambisonic decoder for reproducing ambisonic signals over
 * headphones. The decoder supports sound-field rotation for head-tracking and
 * may also accomodate custom HRIR sets via the SOFA standard.
 *
 * Dependencies:
 *     saf_utilities, afSTFTlib, saf_hrir, saf_vbap, saf_sh, saf_sofa_reader
 * Author, date created:
 *     Leo McCormack, 14.04.2018
 */

#ifndef __AMBI_BIN_H_INCLUDED__
#define __AMBI_BIN_H_INCLUDED__

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */
 
/* ========================================================================== */
/*                             Presets + Constants                            */
/* ========================================================================== */

/*
 * Enum: INPUT_ORDERS
 * ------------------
 * Available decoding orders
 *
 * Options:
 *     INPUT_ORDER_FIRST   - First-order decoding (4 channel input)
 *     INPUT_ORDER_SECOND  - Second-order decoding (9 channel input)
 *     INPUT_ORDER_THIRD   - Third-order decoding (16 channel input)
 *     INPUT_ORDER_FOURTH  - Fourth-order decoding (25 channel input)
 *     INPUT_ORDER_FIFTH   - Fifth-order decoding (36 channel input)
 *     INPUT_ORDER_SIXTH   - Sixth-order decoding (49 channel input)
 *     INPUT_ORDER_SEVENTH - Seventh-order decoding (64 channel input)
 */
#define AMBI_BIN_MAX_SH_ORDER ( 7 )
typedef enum _INPUT_ORDERS{
    INPUT_ORDER_FIRST = 1,
    INPUT_ORDER_SECOND,
    INPUT_ORDER_THIRD,
    INPUT_ORDER_FOURTH,
    INPUT_ORDER_FIFTH,
    INPUT_ORDER_SIXTH,
    INPUT_ORDER_SEVENTH
    
}INPUT_ORDERS;
    
/*
 * Enum: DECODING_METHODS
 * ----------------------
 * Available decoding methods. See "saf_hoa_internal.h" for a more indepth
 * description of each approach.
 *
 * Options:
 *     DECODING_METHOD_LS       - Least-squares (LS) decoder
 *     DECODING_METHOD_LSDIFFEQ - Least-squares (LS) decoder with diffuse-field
 *                                spectral equalisation
 *     DECODING_METHOD_SPR      - Spatial resampling decoder (on the same lines
 *                                as the virtual loudspeaker approach)
 *     DECODING_METHOD_TA       - Time-alignment
 *     DECODING_METHOD_MAGLS    - Magnitude least-squares decoder
 */
#define AMBI_BIN_NUM_DECODING_METHODS ( 5 )
typedef enum _DECODING_METHODS{
    DECODING_METHOD_LS = 1,
    DECODING_METHOD_LSDIFFEQ,
    DECODING_METHOD_SPR,
    DECODING_METHOD_TA,
    DECODING_METHOD_MAGLS,
    
}DECODING_METHODS;

/*
 * Enum: CH_ORDER
 * --------------
 * Available Ambisonic channel ordering conventions
 * Note: CH_FUMA only supported for 1st order input.
 * Further note: FuMa: CH_FUMA+NORM_FUMA, AmbiX: CH_ACN+NORM_SN3D
 *
 * Options:
 *     CH_ACN  - Ambisonic Channel Numbering (ACN)
 *     CH_FUMA - (Obsolete) Furse-Malham/B-format (WXYZ)
 */
#define AMBI_BIN_NUM_CH_ORDERINGS ( 2 )
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
#define AMBI_BIN_NUM_NORM_TYPES ( 3 )
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

#define AMBI_BIN_PROGRESSBARTEXT_CHAR_LENGTH ( 256 )


/* ========================================================================== */
/*                               Main Functions                               */
/* ========================================================================== */

/*
 * Function: ambi_bin_create
 * -------------------------
 * Creates an instance of ambi_bin
 *
 * Input Arguments:
 *     phAmbi - & address of ambi_bin handle
 */
void ambi_bin_create(void** const phAmbi);

/*
 * Function: ambi_bin_destroy
 * --------------------------
 * Destroys an instance of ambi_bin
 *
 * Input Arguments:
 *     phAmbi - & address of ambi_bin handle
 */
void ambi_bin_destroy(void** const phAmbi);

/*
 * Function: ambi_bin_init
 * -----------------------
 * Initialises ambi_bin with default settings, and samplerate.
 *
 * Input Arguments:
 *     hAmbi      - ambi_bin handle
 *     samplerate - host samplerate.
 */
void ambi_bin_init(void* const hAmbi,
                   int samplerate);

/*
 * Function: ambi_bin_initCodec
 * ----------------------------
 * Intialises the codec variables, based on current global/user parameters
 *
 * Input Arguments:
 *     hAmbi - ambi_bin handle
 */
void ambi_bin_initCodec(void* const hAmbi);

/*
 * Function: ambi_bin_process
 * --------------------------
 * Decodes input spherical harmonic signals to the binaural channels.
 *
 * Input Arguments:
 *     hAmbi     - ambi_bin handle
 *     inputs    - input channel buffers; 2-D array: nInputs x nSamples
 *     outputs   - Output channel buffers; 2-D array: nOutputs x nSamples
 *     nInputs   - number of input channels
 *     nOutputs  - number of output channels
 *     nSamples  - number of samples in 'inputs'/'output' matrices
 *     isPlaying - flag to say if there is audio in the input buffers, 0: no
 *                 audio, reduced processing, 1: audio, full processing
 */
void ambi_bin_process(void* const hAmbi,
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
 * Function: ambi_bin_refreshParams
 * --------------------------------
 * Sets intialisation flags to 1. i.e. re-initialise all settings/variables
 * as ambi_bin is currently configured, at next available opportunity.
 *
 * Input Arguments:
 *     hAmbi - ambi_bin handle
 */
void ambi_bin_refreshParams(void* const hAmbi);

/*
 * Function: ambi_bin_setUseDefaultHRIRsflag
 * -----------------------------------------
 * Sets flag to dictate whether the default HRIRs in the Spatial_Audio_Framework
 * should be used, or a custom HRIR set loaded via a SOFA file.
 * Note: if the custom set failes to load correctly, ambi_bin will revert to the
 * defualt set. Use 'ambi_bin_getUseDefaultHRIRsflag' to check if loading was
 * successful.
 *
 * Input Arguments:
 *     hAmbi     - ambi_bin handle
 *     newState  - 0: use custom HRIR set, 1: use default HRIR set
 */
void ambi_bin_setUseDefaultHRIRsflag(void* const hAmbi, int newState);

/*
 * Function: ambi_bin_setSofaFilePath
 * ----------------------------------
 * Sets the file path for a .sofa file, in order to employ a custom HRIR set for
 * the decoding.
 * Note: if the custom set failes to load correctly, hcompass will revert to the
 * defualt set. Use 'ambi_bin_getUseDefaultHRIRsflag' to check if loading was
 * successful.
 *
 * Input Arguments:
 *     hAmbi - ambi_bin handle
 *     path  - file path to .sofa file (WITH file extension)
 */
void ambi_bin_setSofaFilePath(void* const hAmbi, const char* path);

/*
 * Function: ambi_bin_setInputOrderPreset
 * --------------------------------------
 * Sets the decoding order. If decoding order is higher than the input signal
 * order, the extra required channels are filled with zeros. If the decoding
 * order is lower than the input signal order, the number input signals is
 * truncated accordingly.
 *
 * Input Arguments:
 *     hAmbi    - ambi_bin handle
 *     newValue - new decoding order (see 'INPUT_ORDERS' enum)
 */
void ambi_bin_setInputOrderPreset(void* const hAmbi, INPUT_ORDERS newPreset);

/*
 * Function: ambi_bin_setDecodingMethod
 * ------------------------------------
 * Sets the decoding method
 *
 * Input Arguments:
 *     hAmbi     - ambi_bin handle
 *     newMethod - new decoding method (see 'DECODING_METHODS' enum)
 */
void ambi_bin_setDecodingMethod(void* const hAmbi, DECODING_METHODS newMethod);

/*
 * Function: ambi_bin_setChOrder
 * -----------------------------
 * Sets the Ambisonic channel ordering convention to decode with, in order to
 * match the convention employed by the input signals
 *
 * Input Arguments:
 *     hAmbi    - ambi_bin handle
 *     newOrder - convention to use (see 'CH_ORDER' enum)
 */
void ambi_bin_setChOrder(void* const hAmbi, int newOrder);

/*
 * Function: ambi_bin_setNormType
 * ------------------------------
 * Sets the Ambisonic normalisation convention to decode with, in order to match
 * with the convention employed by the input signals.
 *
 * Input Arguments:
 *     hAmbi   - ambi_bin handle
 *     newType - convention to use (see 'NORM_TYPE' enum)
 */
void ambi_bin_setNormType(void* const hAmbi, int newType);

/*
 * Function: ambi_bin_setEnableMaxRE
 * ---------------------------------
 * Sets a flag to enable/disable the max_rE weighting
 *
 * Input Arguments:
 *     hAmbi    - ambi_bin handle
 *     newState - 0: disabled, 1: enabled
 */
void ambi_bin_setEnableMaxRE(void* const hAmbi, int newState);

/*
 * Function: ambi_bin_setEnableDiffuseMatching
 * -------------------------------------------
 * Sets a flag to enable/disable the diffuseness covariance constraint.
 *
 * Input Arguments:
 *     hAmbi    - ambi_bin handle
 *     newState - 0: disabled, 1: enabled
 */
void ambi_bin_setEnableDiffuseMatching(void* const hAmbi, int newState);

/*
 * Function: ambi_bin_setEnablePhaseWarping
 * ----------------------------------------
 * Not implemented yet! Sets a flag to enable/disable phase warping.
 *
 * Input Arguments:
 *     hAmbi    - ambi_bin handle
 *     newState - 0: disabled, 1: enabled
 */
void ambi_bin_setEnablePhaseWarping(void* const hAmbi, int newState);

/*
 * Function: ambi_bin_setEnableRotation
 * ------------------------------------
 * Sets the flag to enable/disable sound-field rotation.
 *
 * Input Arguments:
 *     hAmbi     - ambi_bin handle
 *     newState - 0: disable, 1: enable
 */
void ambi_bin_setEnableRotation(void* const hAmbi, int newState);

/*
 * Function: ambi_bin_setYaw
 * -------------------------
 * Sets the 'yaw' rotation angle
 *
 * Input Arguments:
 *     hAmbi   - ambi_bin handle
 *     newYaw - the 'yaw' rotation angle, in DEGREES
 */
void ambi_bin_setYaw(void* const hAmbi, float newYaw);

/*
 * Function: ambi_bin_setPitch
 * ---------------------------
 * Sets the 'pitch' rotation angle
 *
 * Input Arguments:
 *     hAmbi     - ambi_bin handle
 *     newPitch - the 'pitch' rotation angle, in DEGREES
 */
void ambi_bin_setPitch(void* const hAmbi, float newPitch);

/*
 * Function: ambi_bin_setRoll
 * --------------------------
 * Sets the 'roll' rotation angle
 *
 * Input Arguments:
 *     hAmbi    - ambi_bin handle
 *     newRoll - the 'roll' rotation angle, in DEGREES
 */
void ambi_bin_setRoll(void* const hAmbi, float newRoll);

/*
 * Function: ambi_bin_setFlipYaw
 * -----------------------------
 * Sets a flag as to whether to "flip" the sign of the current 'yaw' angle
 *
 * Input Arguments:
 *     hAmbi     - ambi_bin handle
 *     newState - 0: do not flip sign, 1: flip the sign
 */
void ambi_bin_setFlipYaw(void* const hAmbi, int newState);

/*
 * Function: ambi_bin_setFlipPitch
 * -------------------------------
 * Sets a flag as to whether to "flip" the sign of the current 'pitch' angle
 *
 * Input Arguments:
 *     hAmbi     - ambi_bin handle
 *     newState - 0: do not flip sign, 1: flip the sign
 */
void ambi_bin_setFlipPitch(void* const hAmbi, int newState);

/*
 * Function: ambi_bin_setFlipRoll
 * ------------------------------
 * Sets a flag as to whether to "flip" the sign of the current 'roll' angle
 *
 * Input Arguments:
 *     hAmbi     - ambi_bin handle
 *     newState - 0: do not flip sign, 1: flip the sign
 */
void ambi_bin_setFlipRoll(void* const hAmbi, int newState);

/*
 * Function: ambi_bin_setRPYflag
 * -----------------------------
 * Sets a flag as to whether to use "yaw-pitch-roll" (0) or "roll-pitch-yaw" (1)
 * rotation order.
 *
 * Input Arguments:
 *     hAmbi     - ambi_bin handle
 *     newState - 0: use "yaw-pitch-roll", 1: use "roll-pitch-yaw"
 */
void ambi_bin_setRPYflag(void* const hAmbi, int newState);


/* ========================================================================== */
/*                                Get Functions                               */
/* ========================================================================== */

/*
 * Function: ambi_bin_getCodecStatus
 * ---------------------------------
 * Returns current codec status.
 *
 * Input Arguments:
 *     hAmbi - ambi_bin handle
 * Returns:
 *     codec status (see 'CODEC_STATUS' enum)
 */
CODEC_STATUS ambi_bin_getCodecStatus(void* const hAmbi);
    
/*
 * Function: ambi_bin_getProgressBar0_1
 * ------------------------------------
 * (Optional) Returns current intialisation/processing progress, between 0..1
 * 0: intialisation/processing has started
 * 1: intialisation/processing has ended
 *
 * Input Arguments:
 *     hAmbi - ambi_bin handle
 * Returns:
 *     current progress, 0..1
 */
float ambi_bin_getProgressBar0_1(void* const hAmbi);
    
/*
 * Function: ambi_bin_getProgressBarText
 * -------------------------------------
 * (Optional) Returns current intialisation/processing progress text
 * Note: "text" string should be (at least) of length:
 *     AMBI_BIN_PROGRESSBARTEXT_CHAR_LENGTH
 *
 * Input Arguments:
 *     hAmbi - ambi_bin handle
 * Output Arguments:
 *     text  - process bar text; AMBI_BIN_PROGRESSBARTEXT_CHAR_LENGTH x 1
 */
void ambi_bin_getProgressBarText(void* const hAmbi, char* text);
    
/*
 * Function: ambi_bin_getUseDefaultHRIRsflag
 * -----------------------------------------
 * Returns the value of a flag used to dictate whether the default HRIRs in the
 * Spatial_Audio_Framework should be used, or a custom HRIR set loaded via a
 * SOFA file.
 * Note: if the custom set failes to load correctly, ambi_bin will revert to the
 * defualt set.
 *
 * Input Arguments:
 *     hAmbi - ambi_bin handle
 * Returns:
 *     0: use custom HRIR set, 1: use default HRIR set
 */
int ambi_bin_getUseDefaultHRIRsflag(void* const hAmbi);

/*
 * Function: ambi_bin_getInputOrderPreset
 * --------------------------------------
 * Returns the decoding order. If decoding order is higher than the input signal
 * order, the extra required channels are filled with zeros. If the decoding
 * order is lower than the input signal order, the number input signals is
 * truncated accordingly.
 *
 * Input Arguments:
 *     hAmbi - ambi_bin handle
 * Returns:
 *     current decoding order (see 'INPUT_ORDERS' enum)
 */
int ambi_bin_getInputOrderPreset(void* const hAmbi);

/*
 * Function: ambi_bin_getDecodingMethod
 * ------------------------------------
 * Returns the currently selected decoding method
 *
 * Input Arguments:
 *     hAmbi - ambi_bin handle
 * Returns:
 *     current decoding order (see 'DECODING_METHODS' enum)
 */
int ambi_bin_getDecodingMethod(void* const hAmbi);

/*
 * Function: ambi_bin_getSofaFilePath
 * ----------------------------------
 * Returns the file path for a .sofa file.
 * Note: if the custom set failes to load correctly, hcompass will revert to the
 * defualt set. Use 'ambi_bin_getUseDefaultHRIRsflag' to check if loading was
 * successful.
 *
 * Input Arguments:
 *     hAmbi - ambi_bin handle
 * Returns:
 *      file path to .sofa file (WITH file extension)
 */
char* ambi_bin_getSofaFilePath(void* const hAmbi);

/*
 * Function: ambi_bin_getChOrder
 * -----------------------------
 * Returns the Ambisonic channel ordering convention currently being used to
 * decode with, which should match the convention employed by the input signals
 *
 * Input Arguments:
 *     hAmbi - ambi_bin handle
 * Returns:
 *     convention currently being used (see 'CH_ORDER' enum)
 */
int ambi_bin_getChOrder(void* const hAmbi);

/*
 * Function: ambi_bin_getNormType
 * ------------------------------
 * Returns the Ambisonic normalisation convention currently being usedto decode
 * with, which should match the convention employed by the input signals.
 *
 * Input Arguments:
 *     hAmbi - ambi_bin handle
 * Returns:
 *     convention currently being used (see 'NORM_TYPE' enum)
 */
int ambi_bin_getNormType(void* const hAmbi); 

/*
 * Function: ambi_bin_getNumEars
 * -----------------------------
 * Returns the number of ears possessed by the average homo sapien
 *
 * Returns:
 *     2
 */
int ambi_bin_getNumEars(void);

/*
 * Function: ambi_bin_getNSHrequired
 * ---------------------------------
 * Returns the number of spherical harmonic signals required by the current
 * decoding order i.e. (current_order+1)^2
 *
 * Input Arguments:
 *     hAmbi - ambi_bin handle
 * Returns:
 *     number of required spherical harmonic signals required by current
 *     decoding order
 */
int ambi_bin_getNSHrequired(void* const hAmbi);

/*
 * Function: ambi_bin_getEnableRotation
 * ------------------------------------
 * Returns the flag value which dictates whether to enable/disable sound-field
 * rotation.
 *
 * Input Arguments:
 *     hAmbi - ambi_bin handle
 * Returns:
 *     0: disabled, 1: enabled
 */
int ambi_bin_getEnableMaxRE(void* const hAmbi);

/*
 * Function: ambi_bin_getEnableDiffuseMatching
 * -------------------------------------------
 * Returns the flag value which dictates whether the diffuse covariance
 * contraint is currently enabled.
 *
 * Input Arguments:
 *     hAmbi - ambi_bin handle
 * Returns:
 *     0: disabled, 1: enabled
 */
int ambi_bin_getEnableDiffuseMatching(void* const hAmbi);

/*
 * Function: ambi_bin_getEnablePhaseWarping
 * ----------------------------------------
 * Returns the flag value which dictates whether the phase warping is currently
 * enabled.
 *
 * Input Arguments:
 *     hAmbi - ambi_bin handle
 * Returns:
 *     0: disabled, 1: enabled
 */
int ambi_bin_getEnablePhaseWarping(void* const hAmbi);

/*
 * Function: ambi_bin_getEnableRotation
 * ------------------------------------
 * Returns the flag value which dictates whether to enable/disable sound-field
 * rotation.
 *
 * Input Arguments:
 *     hAmbi - ambi_bin handle
 * Returns:
 *     0: disabled, 1: enabled
 */
int ambi_bin_getEnableRotation(void* const hAmbi);

/*
 * Function: ambi_bin_getYaw
 * -------------------------
 * Returns the 'yaw' rotation angle
 *
 * Input Arguments:
 *     hAmbi - ambi_bin handle
 * Returns:
 *     the 'yaw' rotation angle, in DEGREES
 */
float ambi_bin_getYaw(void* const hAmbi);

/*
 * Function: ambi_bin_getPitch
 * ---------------------------
 * Returns the 'pitch' rotation angle
 *
 * Input Arguments:
 *     hAmbi - ambi_bin handle
 * Returns:
 *     the 'pitch' rotation angle, in DEGREES
 */
float ambi_bin_getPitch(void* const hAmbi);

/*
 * Function: ambi_bin_getRoll
 * --------------------------
 * Returns the 'roll' rotation angle
 *
 * Input Arguments:
 *     hAmbi - ambi_bin handle
 * Returns:
 *     the 'roll' rotation angle, in DEGREES
 */
float ambi_bin_getRoll(void* const hAmbi);

/*
 * Function: ambi_bin_getFlipYaw
 * -----------------------------
 * Returns a flag as to whether to "flip" the sign of the current 'yaw' angle
 *
 * Input Arguments:
 *     hAmbi - ambi_bin handle
 * Returns:
 *     0: do not flip sign, 1: flip the sign
 */
int ambi_bin_getFlipYaw(void* const hAmbi);

/*
 * Function: ambi_bin_getFlipPitch
 * -------------------------------
 * Returns a flag as to whether to "flip" the sign of the current 'pitch' angle
 *
 * Input Arguments:
 *     hAmbi - ambi_bin handle
 * Returns:
 *     0: do not flip sign, 1: flip the sign
 */
int ambi_bin_getFlipPitch(void* const hAmbi);

/*
 * Function: ambi_bin_getFlipRoll
 * ------------------------------
 * Returns a flag as to whether to "flip" the sign of the current 'roll' angle
 *
 * Input Arguments:
 *     hAmbi - ambi_bin handle
 * Returns:
 *     0: do not flip sign, 1: flip the sign
 */
int ambi_bin_getFlipRoll(void* const hAmbi);

/*
 * Function: ambi_bin_getRPYflag
 * -----------------------------
 * Returns a flag as to whether to use "yaw-pitch-roll" (0) or "roll-pitch-yaw"
 * (1) rotation order.
 *
 * Input Arguments:
 *     hAmbi - ambi_bin handle
 * Returns:
 *     0: use "yaw-pitch-roll", 1: use "roll-pitch-yaw"
 */
int ambi_bin_getRPYflag(void* const hAmbi);

/*
 * Function: ambi_bin_getNDirs
 * ---------------------------
 * Returns the number of directions in the currently used HRIR set
 *
 * Input Arguments:
 *     hAmbi - ambi_bin handle
 * Returns:
 *     number of HRIR directions
 */
int ambi_bin_getNDirs(void* const hAmbi);

/*
 * Function: ambi_bin_getHRIRlength
 * --------------------------------
 * Returns the length of HRIRs in time-domain samples
 *
 * Input Arguments:
 *     hAmbi - ambi_bin handle
 * Returns:
 *     HRIR length in samples
 */
int ambi_bin_getHRIRlength(void* const hAmbi);

/*
 * Function: ambi_bin_getHRIRsamplerate
 * ------------------------------------
 * Returns the HRIR sample rate
 *
 * Input Arguments:
 *     hAmbi - ambi_bin handle
 * Returns:
 *     HRIR sampling rate
 */
int ambi_bin_getHRIRsamplerate(void* const hAmbi);

/*
 * Function: ambi_bin_getDAWsamplerate
 * -----------------------------------
 * Returns the DAW/Host sample rate
 *
 * Input Arguments:
 *     hAmbi - ambi_bin handle
 * Returns:
 *     DAW/Host sampling rate
 */
int ambi_bin_getDAWsamplerate(void* const hAmbi);

/*
 * Function: ambi_bin_getProcessingDelay
 * -------------------------------------
 * Returns the processing delay in samples. May be used for delay compensation
 * features
 *
 * Returns:
 *     processing delay in samples
 */
int ambi_bin_getProcessingDelay(void);

    
#ifdef __cplusplus
} /* extern "C" { */
#endif /* __cplusplus */

#endif /* __AMBI_BIN_H_INCLUDED__ */
