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
 * Filename: multiconv.h (include header)
 * --------------------------------------
 * A multi-channel convolver
 *
 * Dependencies:
 *     saf_utilities
 * Author, date created:
 *     Leo McCormack, 23.09.2019
 */

#ifndef __MULTICONV_H_INCLUDED__
#define __MULTICONV_H_INCLUDED__

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */
    
/* ========================================================================== */
/*                             Presets + Constants                            */
/* ========================================================================== */

#define MULTICONV_MAX_NUM_CHANNELS ( 64 )
    

/* ========================================================================== */
/*                               Main Functions                               */
/* ========================================================================== */

/*
 * Function: multiconv_create
 * --------------------------
 * Creates an instance of multiconv
 *
 * Input Arguments:
 *     phMCnv - & address of multiconv handle
 */
void multiconv_create(void** const phMCnv);

/*
 * Function: multiconv_destroy
 * ---------------------------
 * Destroys an instance of multiconv
 *
 * Input Arguments:
 *     phMCnv - & address of multiconv handle
 */
void multiconv_destroy(void** const phMCnv);

/*
 * Function: multiconv_init
 * ------------------------
 * Initialises an instance of multiconv with default settings
 *
 * Input Arguments:
 *     hMCnv      - multiconv handle
 *     samplerate - host samplerate.
 */
void multiconv_init(void* const hMCnv,
                    int samplerate,
                    int hostBlockSize);
    
/*
 * Function: multiconv_process
 * ---------------------------
 *
 * Input Arguments:
 *     hMCnv     - multiconv handle
 *     inputs    - input channel buffers; 2-D array: nInputs x nSamples
 *     outputs   - Output channel buffers; 2-D array: nOutputs x nSamples
 *     nInputs   - number of input channels
 *     nOutputs  - number of output channels
 *     nSamples  - number of samples in 'inputs'/'output' matrices
 *     isPlaying - flag to say if there is audio in the input buffers, 0: no
 *                 audio, reduced processing, 1: audio, full processing
 */
void multiconv_process(void* const hMCnv,
                       float** const inputs,
                       float** const outputs,
                       int nInputs,
                       int nOutputs,
                       int nSamples,
                       int isPlaying);

    
/* ========================================================================== */
/*                                Set Functions                               */
/* ========================================================================== */
    
void multiconv_refreshParams(void* const hMCnv);

void multiconv_checkReInit(void* const hMCnv);

void multiconv_setFilters(void* const hMCnv,
                          const float** H,
                          int numChannels,
                          int numSamples,
                          int sampleRate);
    
void multiconv_setEnablePart(void* const hMCnv, int newState);
    

/* ========================================================================== */
/*                                Get Functions                               */
/* ========================================================================== */
 
int multiconv_getEnablePart(void* const hMCnv);
    
    
    
#ifdef __cplusplus
} /* extern "C" { */
#endif /* __cplusplus */

#endif /* __MULTICONV_H_INCLUDED__ */
