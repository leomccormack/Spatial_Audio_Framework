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

/**
 * @file ambi_bin.h
 * @brief A binaural Ambisonic decoder for reproducing ambisonic signals over
 *        headphones
 *
 * The decoder includes many historic and current state-of-the-art decoding
 * approaches. It also supports sound-field rotation for head-tracking and may
 * also accomodate custom HRIR sets via the SOFA standard.
 *
 * @author Leo McCormack
 * @date 14.04.2018
 */

#ifndef __AMBI_BIN_H_INCLUDED__
#define __AMBI_BIN_H_INCLUDED__

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */
 
/* ========================================================================== */
/*                             Presets + Constants                            */
/* ========================================================================== */

/**
 * Available decoding orders
 */
typedef enum _AMBI_BIN_INPUT_ORDERS{
    INPUT_ORDER_FIRST = 1, /**< First-order decoding (4 channel input) */
    INPUT_ORDER_SECOND,    /**< Second-order decoding (9 channel input) */
    INPUT_ORDER_THIRD,     /**< Third-order decoding (16 channel input) */
    INPUT_ORDER_FOURTH,    /**< Fourth-order decoding (25 channel input) */
    INPUT_ORDER_FIFTH,     /**< Fifth-order decoding (36 channel input) */
    INPUT_ORDER_SIXTH,     /**< Sixth-order decoding (49 channel input) */
    INPUT_ORDER_SEVENTH    /**< Seventh-order decoding (64 channel input) */
    
}AMBI_BIN_INPUT_ORDERS;
    
/** Maximum supported Ambisonic order */
#define AMBI_BIN_MAX_SH_ORDER ( 7 )
    
/**
 * Available decoding methods for the ambi_bin example. See saf_hoa_internal.h
 * for a more indepth descriptions of the approaches.
 */
typedef enum _AMBI_BIN_DECODING_METHODS{
    DECODING_METHOD_LS = 1,   /**< Least-squares (LS) decoder */
    DECODING_METHOD_LSDIFFEQ, /**< Least-squares (LS) decoder with diffuse-field
                               *   spectral equalisation */
    DECODING_METHOD_SPR,      /**< Spatial resampling decoder (on the same lines
                               *   as the virtual loudspeaker approach) */
    DECODING_METHOD_TA,       /**< Time-alignment (TA) */
    DECODING_METHOD_MAGLS,    /**< Magnitude least-squares decoder (MagLS) */
    
}AMBI_BIN_DECODING_METHODS;
    
/** Number of decoding method options */
#define AMBI_BIN_NUM_DECODING_METHODS ( 5 )

/**
 * Available Ambisonic channel ordering conventions
 *
 * @note CH_FUMA only supported for 1st order input.
 */
typedef enum _AMBI_BIN_CH_ORDER {
    CH_ACN = 1, /**< Ambisonic Channel Numbering (ACN) */
    CH_FUMA     /**< (Obsolete) Furse-Malham/B-format (WXYZ) */
    
} AMBI_BIN_CH_ORDER;
    
/** Number of channel ordering options */
#define AMBI_BIN_NUM_CH_ORDERINGS ( 2 )

/**
 * Available Ambisonic normalisation conventions
 *
 * @note NORM_FUMA only supported for 1st order input and does NOT have the
 *       1/sqrt(2) scaling on the omni.
 */
typedef enum _AMBI_BIN_NORM_TYPES {
    NORM_N3D = 1,   /**< orthonormalised (N3D) */
    NORM_SN3D,      /**< Schmidt semi-normalisation (SN3D) */
    NORM_FUMA       /**< (Obsolete) Same as NORM_SN3D for 1st order */
    
} AMBI_BIN_NORM_TYPES;
 
/** Number of normalisation options */
#define AMBI_BIN_NUM_NORM_TYPES ( 3 )
    
/**
 * Current status of the codec.
 */
typedef enum _AMBI_BIN_CODEC_STATUS {
    CODEC_STATUS_INITIALISED = 0, /**< Codec is initialised and ready to process
                                   *   input audio. */
    CODEC_STATUS_NOT_INITIALISED, /**< Codec has not yet been initialised, or
                                   *   the codec configuration has changed.
                                   *   Input audio should not be processed. */
    CODEC_STATUS_INITIALISING     /**< Codec is currently being initialised,
                                   *   input audio should not be processed. */
} AMBI_BIN_CODEC_STATUS;

/** Length of string, returned by ambi_bin_getProgressBarText() */
#define AMBI_BIN_PROGRESSBARTEXT_CHAR_LENGTH ( 256 )


/* ========================================================================== */
/*                               Main Functions                               */
/* ========================================================================== */

/**
 * Creates an instance of ambi_bin
 *
 * @param[in] phAmbi (&) address of ambi_bin handle
 */
void ambi_bin_create(void** const phAmbi);

/**
 * Destroys an instance of ambi_bin
 *
 * @param[in] phAmbi (&) address of ambi_bin handle
 */
void ambi_bin_destroy(void** const phAmbi);

/**
 * Initialises ambi_bin with default settings, and samplerate.
 *
 * @param[in] hAmbi      ambi_bin handle
 * @param[in] samplerate host samplerate.
 */
void ambi_bin_init(void* const hAmbi,
                   int samplerate);

/**
 * Intialises the codec variables, based on current global/user parameters
 *
 * @param[in] hAmbi ambi_bin handle
 */
void ambi_bin_initCodec(void* const hAmbi);

/**
 * Decodes input spherical harmonic signals to the binaural channels.
 *
 * @param[in] hAmbi    ambi_bin handle
 * @param[in] inputs   Input channel buffers; 2-D array: nInputs x nSamples
 * @param[in] outputs  Output channel buffers; 2-D array: nOutputs x nSamples
 * @param[in] nInputs  Number of input channels
 * @param[in] nOutputs Number of output channels
 * @param[in] nSamples Number of samples in 'inputs'/'output' matrices
 */
void ambi_bin_process(void* const hAmbi,
                      float** const inputs,
                      float** const outputs,
                      int nInputs,
                      int nOutputs,
                      int nSamples);


/* ========================================================================== */
/*                                Set Functions                               */
/* ========================================================================== */

/**
 * Sets intialisation flags to 1, so as to re-initialise all settings/variables
 * (as ambi_bin is currently configured), at next available opportunity.
 */
void ambi_bin_refreshParams(void* const hAmbi);

/**
 * Sets flag to dictate whether the default HRIRs in the Spatial_Audio_Framework
 * should be used, or a custom HRIR set loaded via a SOFA file.
 *
 * @note If the custom set failes to load correctly, ambi_bin will revert to the
 *       default set. Use ambi_bin_getUseDefaultHRIRsflag() to check if loading
 *       was successful.
 *
 * @param[in] hAmbi     ambi_bin handle
 * @param[in] newState  '0' use custom HRIR set, '1' use default HRIR set
 */
void ambi_bin_setUseDefaultHRIRsflag(void* const hAmbi, int newState);

/**
 * Sets the file path for a .sofa file, in order to employ a custom HRIR set for
 * the decoding.
 *
 * @note If the custom set failes to load correctly, hcompass will revert to the
 *       defualt set. Use ambi_bin_getUseDefaultHRIRsflag() to check if loading
 *       was successful.
 *
 * @param[in] hAmbi ambi_bin handle
 * @param[in] path  File path to .sofa file (WITH file extension)
 */
void ambi_bin_setSofaFilePath(void* const hAmbi, const char* path);

/**
 * Sets the decoding order (see AMBI_BIN_INPUT_ORDERS enum).
 *
 * @note If decoding order is higher than the input signal order, the extra
 *       required channels are filled with zeros. If the decoding order is lower
 *       than the input signal order, the number input signals is truncated
 *       accordingly.
 */
void ambi_bin_setInputOrderPreset(void* const hAmbi,
                                  AMBI_BIN_INPUT_ORDERS newPreset);

/**
 * Sets the decoding method (see AMBI_BIN_DECODING_METHODS enum)
 */
void ambi_bin_setDecodingMethod(void* const hAmbi,
                                AMBI_BIN_DECODING_METHODS newMethod);

/**
 * Sets the Ambisonic channel ordering convention to decode with, in order to
 * match the convention employed by the input signals
 */
void ambi_bin_setChOrder(void* const hAmbi, int newOrder);

/**
 * Sets the Ambisonic normalisation convention to decode with, in order to match
 * with the convention employed by the input signals.
 */
void ambi_bin_setNormType(void* const hAmbi, int newType);

/**
 * Sets a flag to enable/disable the max_rE weighting
 */
void ambi_bin_setEnableMaxRE(void* const hAmbi, int newState);

/**
 * Sets a flag to enable/disable (1 or 0) the diffuseness covariance constraint
 */
void ambi_bin_setEnableDiffuseMatching(void* const hAmbi, int newState);

/**
 * Not implemented yet! Sets a flag to enable/disable (1 or 0) phase warping
 */
void ambi_bin_setEnablePhaseWarping(void* const hAmbi, int newState);

/**
 * Sets the flag to enable/disable (1 or 0) sound-field rotation
 */
void ambi_bin_setEnableRotation(void* const hAmbi, int newState);

/**
 * Sets the 'yaw' rotation angle, in degrees
 */
void ambi_bin_setYaw(void* const hAmbi, float newYaw_deg);

/**
 * Sets the 'pitch' rotation angle, in degrees
 */
void ambi_bin_setPitch(void* const hAmbi, float newPitch);

/**
 * Sets the 'roll' rotation angle, in degrees
 */
void ambi_bin_setRoll(void* const hAmbi, float newRoll);

/**
 * Sets a flag as to whether to "flip" the sign of the current 'yaw' angle
 */
void ambi_bin_setFlipYaw(void* const hAmbi, int newState);

/**
 * Sets a flag as to whether to "flip" the sign of the current 'pitch' angle
 */
void ambi_bin_setFlipPitch(void* const hAmbi, int newState);

/**
 * Sets a flag as to whether to "flip" the sign of the current 'roll' angle
 */
void ambi_bin_setFlipRoll(void* const hAmbi, int newState);

/**
 * Sets a flag as to whether to use "yaw-pitch-roll" (0) or "roll-pitch-yaw" (1)
 * rotation order.
 */
void ambi_bin_setRPYflag(void* const hAmbi, int newState);


/* ========================================================================== */
/*                                Get Functions                               */
/* ========================================================================== */

/**
 * Returns current codec status.
 */
AMBI_BIN_CODEC_STATUS ambi_bin_getCodecStatus(void* const hAmbi);
    
/**
 * (Optional) Returns current intialisation/processing progress, between 0..1
 */
float ambi_bin_getProgressBar0_1(void* const hAmbi);
    
/**
 * (Optional) Returns current intialisation/processing progress text
 *
 * @note "text" string should be (at least) of length:
 *       AMBI_BIN_PROGRESSBARTEXT_CHAR_LENGTH  
 */
void ambi_bin_getProgressBarText(void* const hAmbi, char* text);
    
/**
 * Returns the value of a flag used to dictate whether the default HRIRs in the
 * Spatial_Audio_Framework should be used (1), or a custom HRIR set loaded via a
 * SOFA file (0).
 *
 * @note If the custom set fails to load correctly, ambi_bin will revert to the
 *       default set and this function will return '1'.
 */
int ambi_bin_getUseDefaultHRIRsflag(void* const hAmbi);

/**
 * Returns the decoding order.
 *
 * If decoding order is higher than the input signal order, the extra required
 * channels are filled with zeros. If the decoding order is lower than the input
 * signal order, the number input signals is truncated accordingly.
 */
int ambi_bin_getInputOrderPreset(void* const hAmbi);

/**
 * Returns the currently selected decoding method (see
 * 'AMBI_BIN_DECODING_METHODS' enum)
 */
int ambi_bin_getDecodingMethod(void* const hAmbi);

/**
 * Returns the file path for a .sofa file.
 *
 * @note If the custom set fails to load correctly, ambi_bin will revert to the
 *       default set. Use ambi_bin_getUseDefaultHRIRsflag() to check if loading
 *       was successful. Also, note the .sofa file extension is included in the
 *       returned string.
 */
char* ambi_bin_getSofaFilePath(void* const hAmbi);

/**
 * Returns the Ambisonic channel ordering convention currently being used to
 * decode with, which should match the convention employed by the input signals
 * (see AMBI_BIN_CH_ORDER enum)
 */
int ambi_bin_getChOrder(void* const hAmbi);

/**
 * Returns the Ambisonic normalisation convention currently being usedto decode
 * with, which should match the convention employed by the input signals (see
 * 'AMBI_BIN_NORM_TYPE' enum).
 */
int ambi_bin_getNormType(void* const hAmbi); 

/**
 * Returns the number of ears possessed by the average homo sapien (2)
 */
int ambi_bin_getNumEars(void);

/**
 * Returns the number of spherical harmonic signals required by the current
 * decoding order: (current_order+1)^2
 */
int ambi_bin_getNSHrequired(void* const hAmbi);

/**
 * Returns the flag value which dictates whether to enable/disable maxrE
 * weighting ('0' disabled, '1' enabled).
 */
int ambi_bin_getEnableMaxRE(void* const hAmbi);

/**
 * Returns the flag value which dictates whether the diffuse covariance
 * contraint is currently enabled ('0' disabled, '1' enabled).
 */
int ambi_bin_getEnableDiffuseMatching(void* const hAmbi);

/**
 * Returns the flag value which dictates whether the phase warping is currently
 * enabled ('0' disabled, '1' enabled).
 */
int ambi_bin_getEnablePhaseWarping(void* const hAmbi);

/**
 * Returns the flag value which dictates whether to enable/disable sound-field
 * rotation ('0' disabled, '1' enabled).
 */
int ambi_bin_getEnableRotation(void* const hAmbi);

/**
 * Returns the 'yaw' rotation angle, in degree
 */
float ambi_bin_getYaw(void* const hAmbi);

/**
 * Returns the 'pitch' rotation angle, in degrees
 */
float ambi_bin_getPitch(void* const hAmbi);

/**
 * Returns the 'roll' rotation angle, in degrees
 */
float ambi_bin_getRoll(void* const hAmbi);

/**
 * Returns a flag as to whether to "flip" the sign of the current 'yaw' angle
 * ('0' do not flip sign, '1' flip the sign)
 */
int ambi_bin_getFlipYaw(void* const hAmbi);

/**
 * Returns a flag as to whether to "flip" the sign of the current 'pitch' angle
 * ('0' do not flip sign, '1' flip the sign)
 */
int ambi_bin_getFlipPitch(void* const hAmbi);

/**
 * Returns a flag as to whether to "flip" the sign of the current 'roll' angle
 * ('0' do not flip sign, '1' flip the sign)
 */
int ambi_bin_getFlipRoll(void* const hAmbi);

/**
 * Returns a flag as to whether to use "yaw-pitch-roll" (0) or "roll-pitch-yaw"
 * (1) rotation order.
 */
int ambi_bin_getRPYflag(void* const hAmbi);

/**
 * Returns the number of directions in the currently used HRIR set
 */
int ambi_bin_getNDirs(void* const hAmbi);

/**
 * Returns the length of HRIRs in time-domain samples
 */
int ambi_bin_getHRIRlength(void* const hAmbi);

/**
 * Returns the HRIR sample rate
 */
int ambi_bin_getHRIRsamplerate(void* const hAmbi);

/**
 * Returns the DAW/Host sample rate
 */
int ambi_bin_getDAWsamplerate(void* const hAmbi);

/**
 * Returns the processing delay in samples (may be used for delay compensation
 * features)
 */
int ambi_bin_getProcessingDelay(void);

    
#ifdef __cplusplus
} /* extern "C" { */
#endif /* __cplusplus */

#endif /* __AMBI_BIN_H_INCLUDED__ */
