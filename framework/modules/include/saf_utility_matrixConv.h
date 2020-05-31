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
 * @file saf_utility_matrixConv.h
 * @brief Utility: Matrix convolver functions
 *
 * These have been mostly derived from some Matlab scripts by Archontis Politis
 *
 * @author Leo McCormack
 * @date 06.04.2019 
 */

#ifndef SAF_UTILITY_MATRIXCOLV_H_INCLUDED
#define SAF_UTILITY_MATRIXCOLV_H_INCLUDED

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

/**
 * Creates an instance of matrixConv
 *
 * This is a matrix convolver intended for block-by-block processing.
 *
 * @param[in] phMC        (&) address of matrixConv handle
 * @param[in] hopSize     Hop size in samples.
 * @param[in] H           Time-domain filters; FLAT: nCHout x nCHin x length_h
 * @param[in] length_h    Length of the filters
 * @param[in] nCHin       Number of input channels
 * @param[in] nCHout      Number of output channels
 * @param[in] usePartFLAG '0': normal fft-based convolution, '1': fft-based
 *                        partitioned convolution
 */
void saf_matrixConv_create(/* Input Arguments */
                           void ** const phMC,
                           int hopSize,
                           float* H,
                           int length_h,
                           int nCHin,
                           int nCHout,
                           int usePartFLAG);

/**
 * Destroys an instance of matrixConv
 *
 * @param[in] phMC (&) address of matrixConv handle
 */
void saf_matrixConv_destroy(/* Input Arguments */
                            void ** const phMC);

/**
 * Performs the matrix convolution.
 *
 * @note If the number of input+output channels, the filters, or the hopsize
 *       need tochange: simply destroy and re-create the matrixConv instance.
 *
 * @param[in]  hMC        matrixConv handle
 * @param[in]  inputSigs  Input signals;  FLAT: nCHin  x hopSize
 * @param[out] outputSigs Output signals; FLAT: nCHout x hopSize
 */
void saf_matrixConv_apply(/* Input Arguments */
                          void * const hMC,
                          float* inputSigs,
                          /* Output Arguments */
                          float* outputSigs);


/* ========================================================================== */
/*                            Multi-Channel Convolver                         */
/* ========================================================================== */

/**
 * Creates an instance of multiConv
 *
 * This is a multi-channel convolver intended for block-by-block processing.
 *
 * @note nCH can just be 1, in which case this is simply a single-channel
 *       convolver.
 *
 * @param[in] phMC        (&) address of multiConv handle
 * @param[in] hopSize     Hop size in samples.
 * @param[in] H           Time-domain filters; FLAT: nCH x length_h
 * @param[in] length_h    Length of the filters
 * @param[in] nCH         Number of filters & input/output channels
 * @param[in] usePartFLAG '0': normal fft-based convolution, '1': fft-based
 *                        partitioned convolution
 */
void saf_multiConv_create(/* Input Arguments */
                          void ** const phMC,
                          int hopSize,
                          float* H,
                          int length_h,
                          int nCH,
                          int usePartFLAG);

/**
 * Destroys an instance of multiConv
 *
 * @param[in] phMC (&) address of multiConv handle
 */
void saf_multiConv_destroy(/* Input Arguments */
                           void ** const phMC);

/**
 * Performs the multi-channel convolution
 *
 * @param[in] hMC         multiConv handle
 * @param[in] inputSigs   Input signals;  FLAT: nCH x hopSize
 * @param[out] outputSigs Output signals; FLAT: nCH x hopSize
 */
void saf_multiConv_apply(/* Input Arguments */
                         void * const hMC,
                         float* inputSigs,
                         /* Output Arguments */
                         float* outputSigs);


#ifdef __cplusplus
}/* extern "C" */
#endif /* __cplusplus */

#endif /* SAF_UTILITY_MATRIXCOLV_H_INCLUDED */
