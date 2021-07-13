/*
 * Copyright 2020 Leo McCormack
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
 * @example decorrelator.h
 * @brief A multi-channel decorrelator
 *
 * ### Files
 * decorrelator.h (include), decorrelator_internal.h, decorrelator.c,
 * decorrelator_internal.c
 * ### Include Header
 */

/**
 * @file decorrelator.h
 * @brief A multi-channel decorrelator
 *
 * @author Leo McCormack
 * @date 07.07.2020
 * @license ISC
 */

#ifndef __DECORRELATOR_H_INCLUDED__
#define __DECORRELATOR_H_INCLUDED__

#include "_common.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/* ========================================================================== */
/*                               Main Functions                               */
/* ========================================================================== */

/**
 * Creates an instance of decorrelator
 *
 * @param[in] phDecor (&) address of decorrelator handle
 */
void decorrelator_create(void** const phDecor);

/**
 * Destroys an instance of decorrelator
 *
 * @param[in] phDecor (&) address of decorrelator handle
 */
void decorrelator_destroy(void** const phDecor);

/**
 * Initialises decorrelator with default settings, and samplerate.
 *
 * @warning This should not be called while _process() is on-going!
 *
 * @param[in] hDecor     decorrelator handle
 * @param[in] samplerate host samplerate.
 */
void decorrelator_init(void* const hDecor,
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
 * @param[in] hDecor decorrelator handle
 */
void decorrelator_initCodec(void* const hDecor);

/**
 * Decorrelates the input signals
 *
 * @param[in] hDecor   decorrelator handle
 * @param[in] inputs   Input channel buffers; 2-D array: nInputs x nSamples
 * @param[in] outputs  Output channel buffers; 2-D array: nOutputs x nSamples
 * @param[in] nInputs  Number of input channels
 * @param[in] nOutputs Number of output channels
 * @param[in] nSamples Number of samples in 'inputs'/'output' matrices
 */
void decorrelator_process(void* const hDecor,
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
 * (as decorrelator is currently configured), at next available opportunity.
 */
void decorrelator_refreshParams(void* const hDecor);

/** Sets the number of input/output channels */
void decorrelator_setNumberOfChannels(void* const hDecor,
                                      int newValue);

/** Sets the decorrelation amount [0..1] */
void decorrelator_setDecorrelationAmount(void* const hDecor,
                                         float newValue);

/** Sets whether to apply level compensation (0 or 1) */
void decorrelator_setLevelCompensationFlag(void* const hDecor,
                                           int newValue);

/** Sets whether to bypass decorrelating the transients (0 or 1) */
void decorrelator_setTransientBypassFlag(void* const hDecor,
                                         int newValue);


/* ========================================================================== */
/*                                Get Functions                               */
/* ========================================================================== */

/**
 * Returns the processing framesize (i.e., number of samples processed with
 * every _process() call )
 */
int decorrelator_getFrameSize(void);

/** Returns current codec status, see #CODEC_STATUS enum */
CODEC_STATUS decorrelator_getCodecStatus(void* const hDecor);
    
/**
 * (Optional) Returns current intialisation/processing progress, between 0..1
 */
float decorrelator_getProgressBar0_1(void* const hDecor);
    
/**
 * (Optional) Returns current intialisation/processing progress text
 *
 * @note "text" string should be (at least) of length:
 *       #PROGRESSBARTEXT_CHAR_LENGTH
 */
void decorrelator_getProgressBarText(void* const hDecor, char* text);

/** Returns the number of input/output channels */
int decorrelator_getNumberOfChannels(void* const hDecor);

/** Returns the decorrelation amount [0..1] */
float decorrelator_getDecorrelationAmount(void* const hDecor);

/** Returns whether to apply level compensation (0 or 1) */
int decorrelator_getLevelCompensationFlag(void* const hDecor);

/** Returns whether to bypass decorrelating the transients (0 or 1) */
int decorrelator_getTransientBypassFlag(void* const hDecor);

/** Returns the DAW/Host sample rate */
int decorrelator_getDAWsamplerate(void* const hDecor);

/**
 * Returns the processing delay in samples (may be used for delay compensation
 * features)
 */
int decorrelator_getProcessingDelay(void);

    
#ifdef __cplusplus
} /* extern "C" { */
#endif /* __cplusplus */

#endif /* __DECORRELATOR_H_INCLUDED__ */
