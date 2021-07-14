/*
 * Copyright 2019 Leo McCormack
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
 * @example multiconv.h
 * @brief A multi-channel convolver
 *
 * ### Files
 * multiconv.h (include), multiconv_internal.h, multiconv.c,
 * multiconv_internal.c
 * ### Include Header
 */

/**
 * @file multiconv.h
 * @brief A multi-channel convolver 
 * @author Leo McCormack
 * @date 23.09.2019
 * @license ISC
 */

#ifndef __MULTICONV_H_INCLUDED__
#define __MULTICONV_H_INCLUDED__

#include "_common.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/* ========================================================================== */
/*                               Main Functions                               */
/* ========================================================================== */

/**
 * Creates an instance of multiconv
 *
 * @param[in] phMCnv (&) address of multiconv handle
 */
void multiconv_create(void** const phMCnv);

/**
 * Destroys an instance of multiconv
 *
 * @param[in] phMCnv (&) address of multiconv handle
 */
void multiconv_destroy(void** const phMCnv);

/**
 * Initialises an instance of multiconv with default settings
 *
 * @param[in] hMCnv         multiconv handle
 * @param[in] samplerate    Host samplerate
 * @param[in] hostBlockSize Host frame/block size
 */
void multiconv_init(void* const hMCnv,
                    int samplerate,
                    int hostBlockSize);

/**
 * Performs the multuchannel convolution processing
 *
 * @param[in] hMCnv     multiconv handle
 * @param[in] inputs    Input channel buffers; 2-D array: nInputs x nSamples
 * @param[in] outputs   Output channel buffers; 2-D array: nOutputs x nSamples
 * @param[in] nInputs   Number of input channels
 * @param[in] nOutputs  Number of output channels
 * @param[in] nSamples  Number of samples in 'inputs'/'output' matrices
 */
void multiconv_process(void* const hMCnv,
                       const float *const * inputs,
                       float** const outputs,
                       int nInputs,
                       int nOutputs,
                       int nSamples);


/* ========================================================================== */
/*                                Set Functions                               */
/* ========================================================================== */

/**
 * Sets all intialisation flags to 1. Re-initialising all settings/variables,
 * as multiconv is currently configured, at next available opportunity.
 */
void multiconv_refreshParams(void* const hMCnv);

/**
 * Checks whether things have to be reinitialised, and does so if it is needed
 */
void multiconv_checkReInit(void* const hMCnv);

/**
 * Loads the multichannel of filters
 *
 * @param[in] hMCnv       multiconv handle
 * @param[in] H           Input channel buffers; 2-D array:
 *                        numChannels x nSamples
 * @param[in] numChannels Number of channels in loaded data (also the number of
 *                        outputs)
 * @param[in] numSamples  Number of samples (per channel) in the loaded data
 * @param[in] sampleRate  Samplerate of the loaded data
 */
void multiconv_setFilters(void* const hMCnv,
                          const float** H,
                          int numChannels,
                          int numSamples,
                          int sampleRate);
    
/** Enable (1), disable (0), partitioned convolution */
void multiconv_setEnablePart(void* const hMCnv, int newState);
    
/** Sets the number of input/output channels */
void multiconv_setNumChannels(void* const hMCnv, int newValue);
    

/* ========================================================================== */
/*                                Get Functions                               */
/* ========================================================================== */

/**
 * Returns the processing framesize (i.e., number of samples processed with
 * every _process() call )
 */
int multiconv_getFrameSize(void);

/**
 * Returns a flag indicating whether partitioned convolution is enabled (1) or
 * disabled (0)
 */
int multiconv_getEnablePart(void* const hMCnv);

/** Returns the number input/output channels */
int multiconv_getNumChannels(void* const hMCnv);

/** Returns the currect host block size */
int multiconv_getHostBlockSize(void* const hMCnv);

/** Returns the number of filters in the loaded wav file */
int multiconv_getNfilters(void* const hMCnv);

/** Returns the current filter length, in samples */
int multiconv_getFilterLength(void* const hMCnv);

/** Returns the samplerate of the loaded filters */
int multiconv_getFilterFs(void* const hMCnv);

/** Returns the samperate of the host */
int multiconv_getHostFs(void* const hMCnv);

/**
 * Returns the processing delay in samples (may be used for delay compensation
 * features)
 */
int multiconv_getProcessingDelay(void* const hMCnv);


#ifdef __cplusplus
} /* extern "C" { */
#endif /* __cplusplus */

#endif /* __MULTICONV_H_INCLUDED__ */
