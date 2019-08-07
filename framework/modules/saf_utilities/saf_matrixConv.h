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
 * Filename: saf_matrixConv.h
 * --------------------------
 * Matrix convolver. UNTESTED!
 *
 * Dependencies:
 *     Currently, only Intel MKL
 * Author, date created:
 *     Leo McCormack, 06.04.2019
 */

#ifndef SAF_MATRIXCOLV_H_INCLUDED
#define SAF_MATRIXCOLV_H_INCLUDED

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */
    
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
    
    
/* ========================================================================== */
/*                               Main Functions                               */
/* ========================================================================== */

#if defined(INTEL_MKL_VERSION)
    
/* NOTE: Not tested yet!! */
 
/*
 * Function: matrixConv_create
 * ---------------------------
 * Creates an instance of matrixConv
 *
 * Input Arguments:
 *     phMC     - & address of matrixConv handle
 *     hopSize  - hop size in samples.
 *     H        - time-domain filters; FLAT: nCHout x nCHin x length_h
 *     length_h - length of the filters
 *     nCHin    - number of input channels
 *     nCHout   - number of output channels
 */
void matrixConv_create(/* Input Arguments */
                       void ** const phMC,
                       int hopSize,
                       float* H,
                       int length_h,
                       int nCHin,
                       int nCHout);

/*
 * Function: matrixConv_destroy
 * ----------------------------
 * Destroys an instance of matrixConv
 *
 * Input Arguments:
 *     phMC     - & address of matrixConv handle
 */
void matrixConv_destroy(/*Input Arguments*/
                        void ** const phMC);

/*
 * Function: matrixConv_apply
 * --------------------------
 * Performs the matrix convolution.
 *
 * Input Arguments:
 *     hMC        - matrixConv handle
 *     inputSigs  - input signals; FLAT: nCHin x hopSize
 * Output Arguments:
 *     outputSigs - output signals; FLAT: nCHout x hopSize
 */
/* performs the matrix convolution.  */
void matrixConv_apply(void * const hMC,
                      float* inputSigs,
                      float* outputSigs);
    
#else
    /* PLEASE NOTE: matrixConv only works with Intel's MKL. Apple's vDSP does not support non-power-of-2's FFTs. */
#endif /* defined(INTEL_MKL_VERSION) */

    
#ifdef __cplusplus
}/* extern "C" */
#endif /* __cplusplus */

#endif /* SAF_MATRIXCOLV_H_INCLUDED */
