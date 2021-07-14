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
 * @example ambi_bin.h
 * @brief A binaural Ambisonic decoder for reproducing Ambisonic sound scenes
 *        over headphones
 *
 * ### Files
 * ambi_bin.h (include), ambi_bin_internal.h, ambi_bin.c, ambi_bin_internal.c
 * ### Example Usage
 * \code{.c}
 * int main(void) {
 *     void* hAmbi;
 *     int frameSize;
 *
 *     // Create an instance of ambi_bin
 *     ambi_bin_create(&hAmbi);
 *
 *     // Call any set functions, e.g.:
 *     ambi_bin_setNormType(hAmbi, NORM_N3D);
 *     ambi_bin_setInputOrderPreset(hAmbi, SH_ORDER_FIRST);
 *     ambi_bin_setEnableRotation(hAmbi, SAF_TRUE);
 *     ambi_bin_setYaw(hAmbi, 180.0f); // turn the listener around
 *
 *     // Note that many set functions, for example ambi_bin_setYaw(), will
 *     // update their value immediately. Whereas, others, for example
 *     // ambi_bin_setInputOrderPreset(), which could cause clicks/hangs with
 *     // the main processing loop, will only trigger a flag that indicates that
 *     // a re-initisation is required. Therefore, ambi_bin_initCodec() should
 *     // be called after calling these particular set functions in order to
 *     // update the run-time settings accordingly.
 *     //
 *     // This function is fully thread-safe, and actually calling this on a
 *     // separate thread is actively encouraged, in order to avoid the
 *     // aforementioned run-time clicks/hangs. ambi_bin_process() is muted
 *     // if the initialisations are still on-going, and initialisations are
 *     // paused until the current ambi_bin_process() call has completed.
 *     ambi_bin_initCodec(hAmbi);
 *
 *     // ambi_bin_init() should be called once before calling
 *     // ambi_bin_process() in order to inform the example of the host sample
 *     // rate and flush run-time buffers with zeros. It is not safe to call
 *     // this during ambi_bin_process()!
 *     ambi_bin_init(hAmbi, hostSamplingRate);
 *
 *     // The framesize of this example is fixed, and can be found with
 *     frameSize = ambi_bin_getFrameSize();
 *
 *     // Processing frame-by-frame
 *     ...
 *     // Load signals into inputSignalBuffer (numberOfInputs x frameSize)
 *     ambi_bin_process(hAmbi, inputSignalBuffer, outputSignalBuffer,
 *                      numberOfInputs, numberOfOutputs, frameSize);
 *     // Copy signals from outputSignalBuffer (numberOfOutputs x frameSize)
 *     ...
 *
 *     // Destroy this instance of ambi_bin
 *     ambi_bin_destroy(&hAmbi);
 * }
 *\endcode
 * ### Include Header
 */

/**
 * @file ambi_bin.h 
 * @brief A binaural Ambisonic decoder for reproducing Ambisonic sound scenes
 *        over headphones
 *
 * The decoder offers choice over many different binaural decoding options [1-4]
 * It also supports sound-field rotation for head-tracking and can accomodate
 * loading custom HRIR sets via the SOFA standard.
 *
 * @test test__saf_example_ambi_bin()
 *
 * @see [1] Z. Ben-Hur, F. Brinkmann, J. Sheaffer, S. Weinzierl, and B. Rafaely,
 *          "Spectral equalization in binaural signals represented by order-
 *          truncated spherical harmonics" The Journal of the Acoustical
 *          Society of America, vol. 141, no. 6, pp. 4087--4096, 2017.
 * @see [2] B. Bernschutz, A. V. Giner, C. Po"rschmann, and J. Arend, "Binaural
 *          reproduction of plane waves with reduced modal order" Acta Acustica
 *          united with Acustica, vol. 100, no. 5, pp. 972--983, 2014.
 * @see [3] Zaunschirm M, Scho"rkhuber C, Ho"ldrich R. Binaural rendering of
 *          Ambisonic signals by head-related impulse response time alignment
 *          and a diffuseness constraint. The Journal of the Acoustical Society
 *          of America. 2018 Jun 19;143(6):3616-27
 * @see [4] Scho"rkhuber C, Zaunschirm M, Ho"ldrich R. Binaural Rendering of
 *          Ambisonic Signals via Magnitude Least Squares. InProceedings of the
 *          DAGA 2018 (Vol. 44, pp. 339-342).
 *
 * @author Leo McCormack
 * @date 14.04.2018
 * @license ISC
 */

#ifndef __AMBI_BIN_H_INCLUDED__
#define __AMBI_BIN_H_INCLUDED__

#include "_common.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */
 
/* ========================================================================== */
/*                             Presets + Constants                            */
/* ========================================================================== */
    
/**
 * Available decoding methods for the ambi_bin example. See saf_hoa_internal.h
 * for a more indepth descriptions of the approaches.
 */
typedef enum {
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

/** Available HRIR pre-preprocessing options */
typedef enum {
    HRIR_PREPROC_OFF = 1,     /**< No pre-processing active */
    HRIR_PREPROC_EQ,          /**< Diffuse-field EQ (compensates CTF) */
    HRIR_PREPROC_PHASE,       /**< Phase simplification based on ITD */
    HRIR_PREPROC_ALL,         /**< Diffuse-field EQ AND phase-simplification */
}AMBI_BIN_PREPROC;

/** Number of HRIR pre-preprocessing options */
#define AMBI_BIN_NUM_HRIR_PREPROC_OPTIONS ( 4 )


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
 * @warning This should not be called while _process() is on-going!
 *
 * @param[in] hAmbi      ambi_bin handle
 * @param[in] samplerate host samplerate.
 */
void ambi_bin_init(void* const hAmbi,
                   int samplerate);

/**
 * Intialises the codec variables, based on current global/user parameters
 *
 * @note This function is fully threadsafe. It can even be called periodically
 *       via a timer on one thread, while calling _process() on another thread.
 *       Since, if a set function is called (that warrants a re-init), then a
 *       flag is triggered internally and the next time this function is called,
 *       it will wait until the current process() function has completed before
 *       reinitialising the relevant parameters. If the _initCodec() takes
 *       longer than the time it takes for process() to be called again, then
 *       process() is simply bypassed until the codec is ready.
 * @note This function does nothing if no re-initialisations are required.
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
                      const float *const * inputs,
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
 * Sets the decoding order (see #SH_ORDERS enum)
 *
 * @note If decoding order is higher than the input signal order, the extra
 *       required channels are filled with zeros. If the decoding order is lower
 *       than the input signal order, the number input signals is truncated
 *       accordingly.
 */
void ambi_bin_setInputOrderPreset(void* const hAmbi,
                                  SH_ORDERS newPreset);

/**
 * Sets the decoding method (see #AMBI_BIN_DECODING_METHODS enum)
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
 * with the convention employed by the input signals
 */
void ambi_bin_setNormType(void* const hAmbi, int newType);

/** Sets a flag to enable/disable the max_rE weighting */
void ambi_bin_setEnableMaxRE(void* const hAmbi, int newState);

/** Sets a flag to enable/disable (1 or 0) the diffuse-covariance constraint */
void ambi_bin_setEnableDiffuseMatching(void* const hAmbi, int newState);

/** Sets a flag to enable/disable (1 or 0) truncation EQ */
void ambi_bin_setEnableTruncationEQ(void* const hAmbi, int newState);

/** Sets HRIR pre-processing strategy (see #AMBI_BIN_PREPROC enum) */
void ambi_bin_setHRIRsPreProc(void* const hAmbi, AMBI_BIN_PREPROC newType);

/** Sets the flag to enable/disable (1 or 0) sound-field rotation */
void ambi_bin_setEnableRotation(void* const hAmbi, int newState);

/** Sets the 'yaw' rotation angle, in degrees */
void ambi_bin_setYaw(void* const hAmbi, float newYaw_deg);

/** Sets the 'pitch' rotation angle, in degrees */
void ambi_bin_setPitch(void* const hAmbi, float newPitch);

/** Sets the 'roll' rotation angle, in degrees */
void ambi_bin_setRoll(void* const hAmbi, float newRoll);

/** Sets a flag as to whether to "flip" the sign of the current 'yaw' angle */
void ambi_bin_setFlipYaw(void* const hAmbi, int newState);

/** Sets a flag as to whether to "flip" the sign of the current 'pitch' angle */
void ambi_bin_setFlipPitch(void* const hAmbi, int newState);

/** Sets a flag as to whether to "flip" the sign of the current 'roll' angle */
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
 * Returns the processing framesize (i.e., number of samples processed with
 * every _process() call )
 */
int ambi_bin_getFrameSize(void);

/** Returns current codec status, see #CODEC_STATUS enum */
CODEC_STATUS ambi_bin_getCodecStatus(void* const hAmbi);
    
/** (Optional) Returns current intialisation/processing progress, between 0..1*/
float ambi_bin_getProgressBar0_1(void* const hAmbi);
    
/**
 * (Optional) Returns current intialisation/processing progress text
 *
 * @note "text" string should be (at least) of length:
 *       #PROGRESSBARTEXT_CHAR_LENGTH
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
 * #AMBI_BIN_DECODING_METHODS enum)
 */
AMBI_BIN_DECODING_METHODS ambi_bin_getDecodingMethod(void* const hAmbi);

/**
 * Returns the file path for a .sofa file
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
 * (see #CH_ORDER enum)
 */
int ambi_bin_getChOrder(void* const hAmbi);

/**
 * Returns the Ambisonic normalisation convention currently being usedto decode
 * with, which should match the convention employed by the input signals (see
 * #NORM_TYPES enum).
 */
int ambi_bin_getNormType(void* const hAmbi); 

/** Returns the number of ears possessed by the average homo sapien (2) */
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
 * Returns the flag value which dictates whether the truncation EQ is currently
 * enabled ('0' disabled, '1' enabled).
 */
int ambi_bin_getEnableTruncationEQ(void* const hAmbi);

/**
 * Returns HRIR pre-processing strategy.  
 * (see #AMBI_BIN_PREPROC enum)
 */
AMBI_BIN_PREPROC ambi_bin_getHRIRsPreProc(void* const hAmbi);

/**
 * Returns the flag value which dictates whether to enable/disable sound-field
 * rotation ('0' disabled, '1' enabled).
 */
int ambi_bin_getEnableRotation(void* const hAmbi);

/** Returns the 'yaw' rotation angle, in degree */
float ambi_bin_getYaw(void* const hAmbi);

/** Returns the 'pitch' rotation angle, in degrees */
float ambi_bin_getPitch(void* const hAmbi);

/** Returns the 'roll' rotation angle, in degrees */
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

/** Returns the number of directions in the currently used HRIR set */
int ambi_bin_getNDirs(void* const hAmbi);

/** Returns the length of HRIRs in time-domain samples */
int ambi_bin_getHRIRlength(void* const hAmbi);

/** Returns the HRIR sample rate */
int ambi_bin_getHRIRsamplerate(void* const hAmbi);

/** Returns the DAW/Host sample rate */
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
