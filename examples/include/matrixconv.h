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
 * @example matrixconv.h
 * @brief A standard matrix convolver
 * 
 * ### Files
 * matrixconv.h (include), matrixconv_internal.h, matrixconv.c,
 * matrixconv_internal.c
 * ### Include Header
 */

/**
 * @file matrixconv.h
 * @brief A standard matrix convolver
 * @author Leo McCormack
 * @date 30.09.2019
 * @license ISC
 */

#ifndef __MATRIXCONV_H_INCLUDED__
#define __MATRIXCONV_H_INCLUDED__

#include "_common.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/* ========================================================================== */
/*                               Main Functions                               */
/* ========================================================================== */

/**
 * Creates an instance of matrixconv
 *
 * @param[in] phMCnv (&) address of matrixconv handle
 */
void matrixconv_create(void** const phMCnv);

/**
 * Destroys an instance of matrixconv
 *
 * @param[in] phMCnv (&) address of matrixconv handle
 */
void matrixconv_destroy(void** const phMCnv);

/**
 * Initialises an instance of matrixconv with default settings
 *
 * @param[in] hMCnv         matrixconv handle
 * @param[in] samplerate    Host samplerate.
 * @param[in] hostBlockSize Host frame/block size
 */
void matrixconv_init(void* const hMCnv,
                     int samplerate,
                     int hostBlockSize);

/**
 * Performs the matrix convolution processing
 *
 * @param[in] hMCnv     matrixconv handle
 * @param[in] inputs    Input channel buffers; 2-D array: nInputs x nSamples
 * @param[in] outputs   Output channel buffers; 2-D array: nOutputs x nSamples
 * @param[in] nInputs   Number of input channels
 * @param[in] nOutputs  Number of output channels
 * @param[in] nSamples  Number of samples in 'inputs'/'output' matrices
 */
void matrixconv_process(void* const hMCnv,
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
 * as matrixconv is currently configured, at next available opportunity.
 */
void matrixconv_refreshParams(void* const hMCnv);

/**
 * Checks whether things have to be reinitialised, and does so if it is needed
 */
void matrixconv_checkReInit(void* const hMCnv);

/**
 * Loads the matrix of filters, which should have the input filters
 * concatenated for each output.
 *
 * For example, a matrix: 25 x 32 x 512 (numInputs x numOutputs x filterLength)
 * should be loaded as a 25 x 16384   (note 32x512=16384).
 *
 * This is then divided by the number of inputs, which should be user specified
 * to be 32 in this case.
 *
 * @param[in] hMCnv       matrixconv handle
 * @param[in] H           Input channel buffers; 2-D array:
 *                        numChannels x nSamples
 * @param[in] numChannels Number of channels in loaded data (also the number of
 *                        outputs)
 * @param[in] numSamples  Number of samples (per channel) in the loaded data
 * @param[in] sampleRate  Samplerate of the loaded data
 */
void matrixconv_setFilters(void* const hMCnv,
                           const float** H,
                           int numChannels,
                           int numSamples,
                           int sampleRate);

/** Enable (1), disable (0), partitioned convolution */
void matrixconv_setEnablePart(void* const hMCnv, int newState);
    
/**
 * Sets the number of input channels.
 *
 * @note The loaded wav data channels are divided by the number of channels
 *       (into equal lengths). These are interpreted by matrixconv as the
 *       filters to apply to each input channel to acquire the corresponding
 *       output channel
 */
void matrixconv_setNumInputChannels(void* const hMCnv, int newValue);
    

/* ========================================================================== */
/*                                Get Functions                               */
/* ========================================================================== */

/**
 * Returns the processing framesize (i.e., number of samples processed with
 * every _process() call )
 */
int matrixconv_getFrameSize(void);

/**
 * Returns a flag indicating whether partitioned convolution is enabled (1) or
 * disabled (0)
 */
int matrixconv_getEnablePart(void* const hMCnv);
    
/** Returns the number input channels */
int matrixconv_getNumInputChannels(void* const hMCnv);

/**
 * Returns the number of output channels (the same as the number of channels in
 * the loaded wav file)
 */
int matrixconv_getNumOutputChannels(void* const hMCnv);

/** Returns the currect host block size */
int matrixconv_getHostBlockSize(void* const hMCnv);

/**
 * Returns the number of filters in the loaded wav file (number of outputs
 * multiplied by the number of inputs)
 */
int matrixconv_getNfilters(void* const hMCnv);

/** Returns the current filter length, in samples */
int matrixconv_getFilterLength(void* const hMCnv);

/** Returns the samplerate of the loaded filters */
int matrixconv_getFilterFs(void* const hMCnv);

/** Returns the samperate of the host */
int matrixconv_getHostFs(void* const hMCnv);

/**
 * Returns the processing delay in samples (may be used for delay compensation
 * features)
 */
int matrixconv_getProcessingDelay(void* const hMCnv);
    
    
#ifdef __cplusplus
} /* extern "C" { */
#endif /* __cplusplus */

#endif /* __MATRIXCONV_H_INCLUDED__ */
