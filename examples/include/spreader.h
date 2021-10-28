/*
 * Copyright 2021 Leo McCormack
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
 * @example spreader.h
 * @brief An arbitrary array panner (HRIRs, microphone array IRs, etc.) with
 *        coherent and incoherent spreading modes.
 *
 * ### Files
 * spreader.h (include), spreader_internal.h, spreader.c, spreader_internal.c
 * ### Include Header
 */

/**
 * @file: spreader.h
 * @brief An arbitrary array panner (HRIRs, microphone array IRs, etc.) with
 *        coherent and incoherent spreading modes, as described in [1].
 *
 * @see [1] McCormack, L. Politis, A., and Pulkki, V., 2021, October. Rendering
 *          of source spread for arbitrary playback setups based on spatial
 *          covariance matching. In 2021 IEEE Workshop on Applications of Signal
 *          Processing to Audio and Acoustics (WASPAA). IEEE
 *
 * @author Leo McCormack
 * @date 07.04.2021
 * @license ISC
 */

#ifndef __SPREADER_H_INCLUDED__
#define __SPREADER_H_INCLUDED__

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include "_common.h"

/** Maximum number of sources supported by the spreader example */
#define SPREADER_MAX_NUM_SOURCES ( 8 )

/** Available processing modes */
typedef enum {
    SPREADER_MODE_NAIVE = 1, /**< Simple coherent copies of the input signal(s)
                              *   areassigned to the spreading areas */
    SPREADER_MODE_OM,        /**< Optimal mixing solution */
    SPREADER_MODE_EVD        /**< Basic solution based on an Eigenvalue
                              *   decomposition */
} SPREADER_PROC_MODES;

/* ========================================================================== */
/*                               Main Functions                               */
/* ========================================================================== */

/**
 * Creates an instance of the spreader
 *
 * @param[in] phSpr (&) address of spreader handle
 */
void spreader_create(void** const phSpr);

/**
 * Destroys an instance of the spreader
 *
 * @param[in] phSpr (&) address of spreader handle
 */
void spreader_destroy(void** const phSpr);

/**
 * Initialises an instance of spreader with default settings
 *
 * @warning This should not be called while _process() is on-going!
 *
 * @param[in] hSpr       spreader handle
 * @param[in] samplerate Host samplerate.
 */
void spreader_init(void* const hSpr,
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
 * @param[in] hSpr spreader handle
 */
void spreader_initCodec(void* const hSpr);

/**
 * Spatialises and spreads the input signals in the user specified directions
 *
 * @param[in] hSpr      spreader handle
 * @param[in] inputs    Input channel buffers; 2-D array: nInputs x nSamples
 * @param[in] outputs   Output channel buffers; 2-D array: nOutputs x nSamples
 * @param[in] nInputs   Number of input channels
 * @param[in] nOutputs  Number of output channels
 * @param[in] nSamples  Number of samples in 'inputs'/'output' matrices
 */
void spreader_process(void* const hSpr,
                      const float *const * inputs,
                      float** const outputs,
                      int nInputs,
                      int nOutputs,
                      int nSamples);


/* ========================================================================== */
/*                                Set Functions                               */
/* ========================================================================== */

/**
 * Sets all intialisation flags to 1; re-initialising all settings/variables
 * as spreader is currently configured, at next available opportunity.
 */
void spreader_refreshSettings(void* const hSpr);

/** Sets the spreading mode (see #SPREADER_PROC_MODES) */
void spreader_setSpreadingMode(void* const hSpr, int newMode);

/** Sets the averaging coefficient [0..1] */
void spreader_setAveragingCoeff(void* const hSpr, float newValue);
    
/** Sets the panning azimuth for a specific channel index, in DEGREES */
void spreader_setSourceAzi_deg(void* const hSpr,
                               int index,
                               float newAzi_deg);

/** Sets the panning elevation for a specific channel index, in DEGREES */
void spreader_setSourceElev_deg(void* const hSpr,
                                int index,
                                float newElev_deg);

/** Sets the source spread for a specific channel index, in DEGREES */
void spreader_setSourceSpread_deg(void* const hSpr,
                                  int index,
                                  float newSpread_deg);

/** Sets the number of input channels/sources to binauralise. */
void spreader_setNumSources(void* const hSpr, int new_nSources);

/**
 * Sets flag to dictate whether the default HRIRs in the Spatial_Audio_Framework
 * should be used (1), or a custom HRIR set loaded via a SOFA file (0).
 *
 * @note If the custom set fails to load correctly, spreader will revert to
 *       the defualt set. Use spreader_getUseDefaultHRIRsflag() to check if
 *       loading was successful.
 */
void spreader_setUseDefaultHRIRsflag(void* const hSpr, int newState);

/**
 * Sets the file path for a .sofa file, in order to employ a custom HRIR set for
 * the decoding.
 *
 * @note If the custom set fails to load correctly, spreader will revert to
 *       the defualt set. Use spreader_getUseDefaultHRIRsflag() to check if
 *       loading was successful.
 *
 * @param[in] hSpr spreader handle
 * @param[in] path File path to .sofa file (WITH file extension)
 */
void spreader_setSofaFilePath(void* const hSpr, const char* path);


/* ========================================================================== */
/*                                Get Functions                               */
/* ========================================================================== */

/**
 * Returns the processing framesize (i.e., number of samples processed with
 * every _process() call )
 */
int spreader_getFrameSize(void);

/** Returns current codec status codec status (see #CODEC_STATUS enum) */
CODEC_STATUS spreader_getCodecStatus(void* const hSpr);

/**
 * (Optional) Returns current intialisation/processing progress, between 0..1
 * - 0: intialisation/processing has started
 * - 1: intialisation/processing has ended
 */
float spreader_getProgressBar0_1(void* const hSpr);

/**
 * (Optional) Returns current intialisation/processing progress text
 *
 * @note "text" string should be (at least) of length:
 *       #PROGRESSBARTEXT_CHAR_LENGTH
 */
void spreader_getProgressBarText(void* const hSpr, char* text);

/**
 * Returns the pointer to a vector describing which directions are currently
 * being used for the spreading, for a given source index
 */
int* spreader_getDirectionActivePtr(void* const hSpr, int index);

/** Returns the spreading mode (see #SPREADER_PROC_MODES) */
int spreader_getSpreadingMode(void* const hSpr);

/** Returns the averaging coefficient [0..1] */
float spreader_getAveragingCoeff(void* const hSpr);

/** Returns the source azimuth for a given source index, in DEGREES */
float spreader_getSourceAzi_deg(void* const hSpr, int index);

/** Returns the source elevation for a given source index, in DEGREES */
float spreader_getSourceElev_deg(void* const hSpr, int index);

/** Returns the source spread for a given source index, in DEGREES */
float spreader_getSourceSpread_deg(void* const hSpr, int index);

/** Returns the number of inputs/sources in the current config */
int spreader_getNumSources(void* const hSpr);

/** Returns the maximum number of input sources supported by spreader */
int spreader_getMaxNumSources(void);

/** Returns the number of ears possessed by the average homo sapien */
int spreader_getNumOutputs(void* const hSpr);

/** Returns the number of directions in the currently used HRIR set */
int spreader_getNDirs(void* const hSpr);

/** Returns the IR/TF azimuth for a given index, in DEGREES */
float spreader_getIRAzi_deg(void* const hSpr, int index);

/** Returns the IR/TF elevation for a given index, in DEGREES */
float spreader_getIRElev_deg(void* const hSpr, int index);

/** Returns the length of IRs in time-domain samples */
int spreader_getIRlength(void* const hSpr);

/** Returns the IR sample rate */
int spreader_getIRsamplerate(void* const hSpr);

/**
 * Returns the value of a flag used to dictate whether the default HRIRs in the
 * Spatial_Audio_Framework should be used (1), or a custom HRIR set loaded via a
 * SOFA file (0).
 *
 * @note If the custom set fails to load correctly, spreader will revert to
 *       the defualt set, so this will be '1'
 */
int spreader_getUseDefaultHRIRsflag(void* const hSpr);

/**
 * Returns the file path for a .sofa file.
 *
 * @note If the custom set fails to load correctly, spreader will revert to
 *       the defualt set. Use 'spreader_getUseDefaultHRIRsflag()' to check
 *       if loading was successful.
 *
 * @param[in] hSpr spreader handle
 * @returns        File path to .sofa file (WITH file extension)
 */
char* spreader_getSofaFilePath(void* const hSpr);

/** Returns the DAW/Host sample rate */
int spreader_getDAWsamplerate(void* const hSpr);

/**
 * Returns the processing delay in samples (may be used for delay compensation
 * purposes)
 */
int spreader_getProcessingDelay(void);


#ifdef __cplusplus
} /* extern "C" { */
#endif /* __cplusplus */

#endif /* __SPREADER_H_INCLUDED__ */
