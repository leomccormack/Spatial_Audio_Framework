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
 * Matrix convolver functions mostly stolen from some Matlab scripts by
 * Archontis Politis ;-)
 * (with permission of course)
 *
 * Dependencies:
 *     saf_fft
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

/* Included:
 *   Matrix Convolver
 *     y = H * x; looped/summed over in/output channels, applied block-by-block
 *       where
 *         y: nOutputChannels x blockSize
 *         x: nInputChannels  x blockSize
 *         H: nOutputChannels x nInputChannels x filterLength
 *
 *   Multi Convolver
 *     y = H * x; looped over channels, applied block-by-block
 *       where
 *         y: nChannels x blockSize
 *         x: nChannels x blockSize
 *         H: nChannels x filterLength
 */

/* ========================================================================== */
/*                              Matrix Convolver                              */
/* ========================================================================== */

/*
 * Function: saf_matrixConv_create
 * -------------------------------
 * Creates an instance of matrixConv
 * This is a matrix convolver intended for block-by-block processing.
 *
 * Input Arguments:
 *     phMC        - & address of matrixConv handle
 *     hopSize     - hop size in samples.
 *     H           - time-domain filters; FLAT: nCHout x nCHin x length_h
 *     length_h    - length of the filters
 *     nCHin       - number of input channels
 *     nCHout      - number of output channels
 *     usePartFLAG - 0: normal fft-based convolution, 1: fft-based partitioned
 *                   convolution
 */
void saf_matrixConv_create(/* Input Arguments */
                           void ** const phMC,
                           int hopSize,
                           float* H,
                           int length_h,
                           int nCHin,
                           int nCHout,
                           int usePartFLAG);

/*
 * Function: saf_matrixConv_destroy
 * --------------------------------
 * Destroys an instance of matrixConv
 *
 * Input Arguments:
 *     phMC     - & address of matrixConv handle
 */
void saf_matrixConv_destroy(/*Input Arguments*/
                            void ** const phMC);

/*
 * Function: saf_matrixConv_apply
 * ------------------------------
 * Performs the matrix convolution.
 * Note: if the number of input+output channels, the filters, or the hopsize
 * change: simply destroy and re-create the matrixConv instance
 *
 * Input Arguments:
 *     hMC        - matrixConv handle
 *     inputSigs  - input signals;  FLAT: nCHin  x hopSize
 * Output Arguments:
 *     outputSigs - output signals; FLAT: nCHout x hopSize
 */
void saf_matrixConv_apply(void * const hMC,
                          float* inputSigs,
                          float* outputSigs);


/* ========================================================================== */
/*                            Multi-Channel Convolver                         */
/* ========================================================================== */

/*
 * Function: saf_multiConv_create
 * ------------------------------
 * Creates an instance of multiConv
 * This is a multi-channel convolver intended for block-by-block processing.
 * Note: nCH can just be 1, in which case this is simply a single-channel
 * convolver.
 *
 * Input Arguments:
 *     phMC        - & address of multiConv handle
 *     hopSize     - hop size in samples.
 *     H           - time-domain filters; FLAT: nCH x length_h
 *     length_h    - length of the filters
 *     nCH         - number of filters & input/output channels
 *     usePartFLAG - 0: normal fft-based convolution, 1: fft-based partitioned
 *                   convolution
 */
void saf_multiConv_create(/* Input Arguments */
                          void ** const phMC,
                          int hopSize,
                          float* H,
                          int length_h,
                          int nCH,
                          int usePartFLAG);

/*
 * Function: saf_multiConv_destroy
 * -------------------------------
 * Destroys an instance of multiConv
 *
 * Input Arguments:
 *     phMC     - & address of multiConv handle
 */
void saf_multiConv_destroy(/*Input Arguments*/
                           void ** const phMC);

/*
 * Function: saf_multiConv_apply
 * -----------------------------
 * Performs the multi-channel convolution
 *
 * Input Arguments:
 *     hMC        - multiConv handle
 *     inputSigs  - input signals;  FLAT: nCH x hopSize
 * Output Arguments:
 *     outputSigs - output signals; FLAT: nCH x hopSize
 */
void saf_multiConv_apply(void * const hMC,
                         float* inputSigs,
                         float* outputSigs);


#ifdef __cplusplus
}/* extern "C" */
#endif /* __cplusplus */

#endif /* SAF_MATRIXCOLV_H_INCLUDED */
