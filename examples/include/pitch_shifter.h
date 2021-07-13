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
 * @example pitch_shifter.h
 * @brief A very basic multichannel pitch shifter
 * 
 * ### Files
 * pitch_shifter.h (include), pitch_shifter_internal.h, pitch_shifter.c,
 * pitch_shifter_internal.c
 * ### Include Header
 */

/**
 * @file pitch_shifter.h
 * @brief A very basic multichannel pitch shifter
 *
 * @author Leo McCormack
 * @date 05.05.2020
 * @license ISC
 */

#ifndef __PITCH_SHIFTER_H_INCLUDED__
#define __PITCH_SHIFTER_H_INCLUDED__

#include "_common.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */
    
/* ========================================================================== */
/*                             Presets + Constants                            */
/* ========================================================================== */

/**
 * Available FFT size options. The higher it is, the more drastic the pitch
 * shifting factor can be, at the cost of increased latency and cpu requirements
 */
typedef enum {
    PITCH_SHIFTER_FFTSIZE_512 = 1,
    PITCH_SHIFTER_FFTSIZE_1024,
    PITCH_SHIFTER_FFTSIZE_2048,
    PITCH_SHIFTER_FFTSIZE_4096,
    PITCH_SHIFTER_FFTSIZE_8192,
    PITCH_SHIFTER_FFTSIZE_16384

}PITCH_SHIFTER_FFTSIZE_OPTIONS;

/** Number of FFT size options */
#define PITCH_SHIFTER_NUM_FFTSIZE_OPTIONS ( 6 )

/**
 * Available oversampling options. The higher it is, the better the signal
 * fidelity, but at the cost of increased latency and cpu requirements
 */
typedef enum {
    PITCH_SHIFTER_OSAMP_2 = 1,
    PITCH_SHIFTER_OSAMP_4,
    PITCH_SHIFTER_OSAMP_8,
    PITCH_SHIFTER_OSAMP_16,
    PITCH_SHIFTER_OSAMP_32

}PITCH_SHIFTER_OSAMP_OPTIONS;

/** Number of over-sampling options */
#define PITCH_SHIFTER_NUM_OSAMP_OPTIONS ( 5 )

/** Maximum pitch shifting factor */
#define PITCH_SHIFTER_MAX_SHIFT_FACTOR ( 2.0f )

/** Minimum pitch shifting factor */
#define PITCH_SHIFTER_MIN_SHIFT_FACTOR ( 0.5f )


/* ========================================================================== */
/*                               Main Functions                               */
/* ========================================================================== */

/**
 * Creates an instance of pitch_shifter
 *
 * @param[in] phPS (&) address of pitch_shifter handle
 */
void pitch_shifter_create(void** const phPS);

/**
 * Destroys an instance of pitch_shifter
 *
 * @param[in] phPS (&) address of pitch_shifter handle
 */
void pitch_shifter_destroy(void** const phPS);

/**
 * Initialises an instance of pitch_shifter with default settings
 *
 * @warning This should not be called while _process() is on-going!
 *
 * @param[in] hPS       pitch_shifter handle
 * @param[in] samplerate Host samplerate.
 */
void pitch_shifter_init(void* const hPS,
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
 * @param[in] hPS pitch_shifter handle
 */
void pitch_shifter_initCodec(void* const hPS);

/**
 * Pitch shifts the input signals
 *
 * @param[in] hPS      pitch_shifter handle
 * @param[in] inputs   Input channel buffers; 2-D array: nInputs x nSamples
 * @param[in] outputs  Output channel buffers; 2-D array: nOutputs x nSamples
 * @param[in] nInputs  Number of input channels
 * @param[in] nOutputs Number of output channels
 * @param[in] nSamples Number of samples in 'inputs'/'output' matrices
 */
void pitch_shifter_process(void* const hPS,
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
 * as pitch_shifter is currently configured, at next available opportunity.
 *
 * @param[in] hPS pitch_shifter handle
 */
void pitch_shifter_refreshParams(void* const hPS);

/**
 * Sets the pitch shift factor, 1: no change, 2: up one octave, 0.5: down one
 * octave
 */
void pitch_shifter_setPitchShiftFactor(void* const hPS, float newValue);

/** Sets the number channels to pitch shift */
void pitch_shifter_setNumChannels(void* const hPS, int newValue);

/**
 * Sets the FFT size used by the algorithm (see #PITCH_SHIFTER_FFTSIZE_OPTIONS
 * enum)
 */
void pitch_shifter_setFFTSizeOption(void* const hPS,
                                    PITCH_SHIFTER_FFTSIZE_OPTIONS newOption);

/**
 * Sets the oversampling factor used by the algorithm (see
 * #PITCH_SHIFTER_OSAMP_OPTIONS enum)
 */
void pitch_shifter_setOSampOption(void* const hPS,
                                  PITCH_SHIFTER_OSAMP_OPTIONS newOption);


/* ========================================================================== */
/*                                Get Functions                               */
/* ========================================================================== */

/**
 * Returns the processing framesize (i.e., number of samples processed with
 * every _process() call )
 */
int pitch_shifter_getFrameSize(void);

/** Returns current codec status (see #CODEC_STATUS enum) */
CODEC_STATUS pitch_shifter_getCodecStatus(void* const hPS);

/**
 * (Optional) Returns current intialisation/processing progress, between 0..1
 *  - 0: intialisation/processing has started
 *  - 1: intialisation/processing has ended
 */
float pitch_shifter_getProgressBar0_1(void* const hPS);

/**
 * (Optional) Returns current intialisation/processing progress text
 *
 * @note "text" string should be (at least) of length:
 *       #PROGRESSBARTEXT_CHAR_LENGTH
 */
void pitch_shifter_getProgressBarText(void* const hPS, char* text);

/**
 * Returns the pitch shift factor, 1: no change, 2: up one octave, 0.5: down one
 * octave
 */
float pitch_shifter_getPitchShiftFactor(void* const hPS);

/**
 * Returns the FFT size used by the algorithm (see
 * #PITCH_SHIFTER_FFTSIZE_OPTIONS enum)
 */
PITCH_SHIFTER_FFTSIZE_OPTIONS pitch_shifter_getFFTSizeOption(void* const hPS);

/**
 * Returns the oversampling factor used by the algorithm (see
 * #PITCH_SHIFTER_OSAMP_OPTIONS enum)
 */
PITCH_SHIFTER_OSAMP_OPTIONS pitch_shifter_getOSampOption(void* const hPS);

/** Returns the number of channels required by the current configuration */
int pitch_shifter_getNCHrequired(void* const hPS);

/**
 * Returns the processing delay in samples (may be used for delay compensation
 * features)
 */
int pitch_shifter_getProcessingDelay(void* const hPS);


#ifdef __cplusplus
} /* extern "C" { */
#endif /* __cplusplus */

#endif /* __PITCH_SHIFTER_H_INCLUDED__ */
