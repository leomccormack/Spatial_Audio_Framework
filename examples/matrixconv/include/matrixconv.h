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

/*
 * Filename: matrixconv.h (include header)
 * ---------------------------------------
 * A matrix convolver
 *
 * Dependencies:
 *     saf_utilities
 * Author, date created:
 *     Leo McCormack, 30.09.2019
 */

#ifndef __MATRIXCONV_H_INCLUDED__
#define __MATRIXCONV_H_INCLUDED__

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */
    
/* ========================================================================== */
/*                             Presets + Constants                            */
/* ========================================================================== */

#define MATRIXCONV_MAX_NUM_CHANNELS ( 64 )
    

/* ========================================================================== */
/*                               Main Functions                               */
/* ========================================================================== */

/*
 * Function: matrixconv_create
 * ---------------------------
 * Creates an instance of matrixconv
 *
 * Input Arguments:
 *     phMCnv - & address of matrixconv handle
 */
void matrixconv_create(void** const phMCnv);

/*
 * Function: matrixconv_destroy
 * ----------------------------
 * Destroys an instance of matrixconv
 *
 * Input Arguments:
 *     phMCnv - & address of matrixconv handle
 */
void matrixconv_destroy(void** const phMCnv);

/*
 * Function: matrixconv_init
 * -------------------------
 * Initialises an instance of matrixconv with default settings
 *
 * Input Arguments:
 *     hMCnv      - matrixconv handle
 *     samplerate - host samplerate.
 */
void matrixconv_init(void* const hMCnv,
                    int samplerate,
                    int hostBlockSize);
    
/*
 * Function: matrixconv_process
 * ----------------------------
 *
 * Input Arguments:
 *     hMCnv     - matrixconv handle
 *     inputs    - input channel buffers; 2-D array: nInputs x nSamples
 *     outputs   - Output channel buffers; 2-D array: nOutputs x nSamples
 *     nInputs   - number of input channels
 *     nOutputs  - number of output channels
 *     nSamples  - number of samples in 'inputs'/'output' matrices
 */
void matrixconv_process(void* const hMCnv,
                       float** const inputs,
                       float** const outputs,
                       int nInputs,
                       int nOutputs,
                       int nSamples);

    
/* ========================================================================== */
/*                                Set Functions                               */
/* ========================================================================== */
    
void matrixconv_refreshParams(void* const hMCnv);

void matrixconv_checkReInit(void* const hMCnv);

void matrixconv_setFilters(void* const hMCnv,
                          const float** H,
                          int numChannels,
                          int numSamples,
                          int sampleRate);
    
void matrixconv_setEnablePart(void* const hMCnv, int newState);
    
void matrixconv_setNumInputChannels(void* const hMCnv, int newValue);
    

/* ========================================================================== */
/*                                Get Functions                               */
/* ========================================================================== */
 
int matrixconv_getEnablePart(void* const hMCnv);
    
int matrixconv_getNumInputChannels(void* const hMCnv);
    
int matrixconv_getNumOutputChannels(void* const hMCnv);
    
int matrixconv_getHostBlockSize(void* const hMCnv);
    
int matrixconv_getNfilters(void* const hMCnv);
    
int matrixconv_getFilterLength(void* const hMCnv);
    
int matrixconv_getFilterFs(void* const hMCnv);
    
int matrixconv_getHostFs(void* const hMCnv);
    
    
#ifdef __cplusplus
} /* extern "C" { */
#endif /* __cplusplus */

#endif /* __MATRIXCONV_H_INCLUDED__ */
