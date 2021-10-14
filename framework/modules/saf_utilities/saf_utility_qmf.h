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
 *@addtogroup Utilities
 *@{
 * @file saf_utility_qmf.h
 * @brief An implementation of the complex Quadrature Mirror Filterbank (QMF)
 *        described in [1].
 *
 * @see [1] Herre, J., Purnhagen, H., Breebaart, J., Faller, C., Disch, S.,
 *          Kjörling, K., Schuijers, E., Hilpert, J. and Myburg, F., 2005. The
 *          reference model architecture for MPEG spatial audio coding.
 *
 * @author Leo McCormack
 * @date 14.07.2020
 * @license ISC
 */

#ifndef SAF_QMF_H_INCLUDED
#define SAF_QMF_H_INCLUDED

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include "saf_utility_complex.h"

/** Options for how the frequency domain data is permuted when using qmf */
typedef enum {
    QMF_BANDS_CH_TIME, /**< nBands x nChannels x nTimeHops */
    QMF_TIME_CH_BANDS  /**< nTimeHops x nChannels x nBands */

}QMF_FDDATA_FORMAT;

/**
 * Creates an instance of the qmf filterbank
 *
 * @test test__qmf()
 *
 * @param[in] phQMF      (&) address of qmf handle
 * @param[in] nCHin      Number of input channels
 * @param[in] nCHout     Number of output channels
 * @param[in] hopsize    Hop size, in samples
 * @param[in] hybridmode 0: disabled, 1: hybrid-filtering enabled
 * @param[in] format     frequency-domain frame format, see #QMF_FDDATA_FORMAT
 *                       enum
 */
void qmf_create(/* Input Arguments */
                void ** const phQMF,
                int nCHin,
                int nCHout,
                int hopsize,
                int hybridmode,
                QMF_FDDATA_FORMAT format);

/**
 * Destroys an instance of the qmf filterbank
 *
 * @param[in] phQMF (&) address of qmf handle
 */
void qmf_destroy(void ** const phQMF);

/**
 * Performs QMF analysis of the input time-domain signals
 *
 * @param[in]  hQMF      qmf handle
 * @param[in]  dataTD    Time-domain input; nCHin x framesize
 * @param[in]  framesize Frame size of time-domain data
 * @param[out] dataFD    Frequency-domain output; see #QMF_FDDATA_FORMAT enum
 */
void qmf_analysis(/* Input Arguments */
                  void * const hQMF,
                  float** dataTD,
                  int framesize,
                  /* Output Arguments */
                  float_complex*** dataFD);

/**
 * Performs QMF synthesis of the input frequency-domain signals
 *
 * @param[in]  hQMF      qmf handle
 * @param[in]  dataFD    Frequency-domain input; see #QMF_FDDATA_FORMAT enum
 * @param[in]  framesize Frame size of time-domain data
 * @param[out] dataTD    Time-domain output;  nCHout x framesize
 */
void qmf_synthesis(/* Input Arguments */
                   void * const hQMF,
                   float_complex*** dataFD,
                   int framesize,
                   /* Output Arguments */
                   float** dataTD);

/**
 * Changes the number input and/or output channels
 *
 * @param[in] hQMF       qmf handle
 * @param[in] new_nCHin  New number of input channels
 * @param[in] new_nCHout New number of output channels
 */
void qmf_channelChange(void * const hQMF,
                       int new_nCHin,
                       int new_nCHout);

/** Flushes the analysis and synthesis buffers with zeros. */
void qmf_clearBuffers(void * const hQMF);

/**
 * Returns the processing delay in samples
 *
 * @note The QMF filterbank delay is broken down into the following:
 *          analysis delay:         5*hopsize
 *        + hybrid-filtering delay: 6*hopsize     (or 0, if it is disabled)
 *        + synthesis delay         4*hopsize
 */
int qmf_getProcDelay(void * const hQMF);

/** Returns the number of frequency bands */
int qmf_getNBands(void * const hQMF);

/**
 * Computes the QMF/hybrid-QMF centre frequencies
 *
 * @note 'nBands' can be found with qmf_getNBands()
 *
 * @param[in]  hQMF       qmf handle
 * @param[in]  fs         Sampling rate in Hz
 * @param[in]  nBands     Length of 'centreFreq'
 * @param[out] centreFreq The frequency vector: nBands x 1
 */
void qmf_getCentreFreqs(/* Input Arguments */
                        void * const hQMF,
                        float fs,
                        int nBands,
                        /* Output Arguments */
                        float* centreFreq);

/**
 * Converts FIR filters into Filterbank Coefficients by passing them through
 * the QMF filterbank
 *
 * @param[in]  hIR        Time-domain FIR; FLAT: N_dirs x nCH x ir_len
 * @param[in]  N_dirs     Number of FIR sets
 * @param[in]  nCH        Number of channels per FIR set
 * @param[in]  ir_len     Length of the FIR
 * @param[in]  hopSize    Hop size
 * @param[in]  hybridmode 0: disabled, 1:enabled
 * @param[out] hFB        The FIRs as Filterbank coefficients;
 *                        FLAT: N_bands x nCH x N_dirs
 */
void qmf_FIRtoFilterbankCoeffs(/* Input Arguments */
                               float* hIR,
                               int N_dirs,
                               int nCH,
                               int ir_len,
                               int hopSize,
                               int hybridmode,
                               /* Output Arguments */
                               float_complex* hFB);


#ifdef __cplusplus
}/* extern "C" */
#endif /* __cplusplus */

#endif /* SAF_QMF_H_INCLUDED */

/**@} */ /* doxygen addtogroup Utilities */
