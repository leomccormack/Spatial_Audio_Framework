/*
 Copyright (c) 2015 Juha Vilkamo
 
 Permission is hereby granted, free of charge, to any person obtaining a copy
 of this software and associated documentation files (the "Software"), to deal
 in the Software without restriction, including without limitation the rights
 to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:
 
 The above copyright notice and this permission notice shall be included in
 all copies or substantial portions of the Software.
 
 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 THE SOFTWARE.
 */

/**
 * @file afSTFTlib.h
 * @brief A modified version of afSTFTlib
 *
 * The original afSTFT code (by Juha Vilkamo) can be found here:
 * https://github.com/jvilkamo/afSTFT
 * This version is slightly modified to be more in-line with how the rest of SAF
 * is structured.
 * The files afSTFTlib.h/.c act as the interface to afSTFT, which is then
 * implemented in afSTFT_internal.h/.c.
 *
 * This version also adds functionality to change the number of channels on the
 * fly, flush the run-time buffers with zeros, return the current frequency
 * vector and the current processing delay.
 * It also incorporates SAF utilities (for the vectorisation and FFT).
 *
 * The afSTFT design is also described in more detail in [1]
 *
 * @see [1] Vilkamo, J., & Ba"ckstro"m, T. (2018). Time--Frequency Processing:
 *          Methods and Tools. In Parametric Time--Frequency Domain Spatial
 *          Audio. John Wiley & Sons.
 *
 * @author Juha Vilkamo
 * @date 08.04.2015
 * @license MIT
 */

#ifndef __afSTFTlib_INCLUDED__
#define __afSTFTlib_INCLUDED__

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */
   
/**
 * Remove the "AFSTFT_USE_SAF_UTILITIES" definition, and add the vecTools.h/c
 * and fft4g.h/c to your project, if you want to use the original afSTFT vector
 * code. Note the vecTools.h/c and fft4g.h/c files, may be found here:
 *     https://github.com/jvilkamo/afSTFT
 */
#define AFSTFT_USE_SAF_UTILITIES
#ifdef AFSTFT_USE_SAF_UTILITIES
# include "../../modules/saf_utilities/saf_utilities.h"
#else
# include <stdio.h>
# include <stdlib.h>
# include <string.h>
#endif
    
/** Prototype filter used by afSTFTlib */
extern const float __afSTFT_protoFilter1024[10240];

/** Prototype filter used by afSTFTlib (low-delay mode) */
extern const float __afSTFT_protoFilter1024LD[10240];

/** Options for how the frequency domain data is permuted when using afSTFT */
typedef enum {
    AFSTFT_BANDS_CH_TIME, /**< nBands x nChannels x nTimeHops */
    AFSTFT_TIME_CH_BANDS  /**< nTimeHops x nChannels x nBands */

}AFSTFT_FDDATA_FORMAT;
    
/**
 * Creates an instance of afSTFT
 *
 * @test test__afSTFT()
 *
 * @param[in] phSTFT       (&) address of afSTFT handle
 * @param[in] nCHin        Number of input channels
 * @param[in] nCHout       Number of output channels
 * @param[in] hopsize      Hop size, in samples
 * @param[in] lowDelayMode 0: disabled, 1: low-delay mode enabled
 * @param[in] hybridmode   0: disabled, 1: hybrid-filtering enabled
 * @param[in] format       Frequency-domain frame format, see
 *                         #AFSTFT_FDDATA_FORMAT enum
 */
void afSTFT_create(void ** const phSTFT,
                   int nCHin,
                   int nCHout,
                   int hopsize,
                   int lowDelayMode,
                   int hybridmode,
                   AFSTFT_FDDATA_FORMAT format);

/**
 * Destroys an instance of afSTFT
 *
 * @param[in] phSTFT  (&) address of afSTFT handle
 */
void afSTFT_destroy(void ** const phSTFT);

/**
 * Performs forward afSTFT transform
 *
 * @param[in]  hSTFT     afSTFT handle
 * @param[in]  dataTD    Time-domain input; nCHin x framesize
 * @param[in]  framesize Frame size of time-domain data
 * @param[out] dataFD    Frequency-domain output; #AFSTFT_FDDATA_FORMAT
 */
void afSTFT_forward(void * const hSTFT,
                    float** dataTD,
                    int framesize,
                    float_complex*** dataFD);

/**
 * Performs forward afSTFT transform (dataFD dimensions are known)
 *
 * @note If the dimensions of dataFD are known, then this function will use the
 *       same dataFD format as used by afSTFT_forward(), but with the speed
 *       of afSTFT_forward_flat()
 *
 * @param[in]  hSTFT        afSTFT handle
 * @param[in]  dataTD       Time-domain input; nCHin x framesize
 * @param[in]  framesize    Frame size of time-domain data
 * @param[in]  dataFD_nCH   Number of channels dataFD is allocated (the max)
 * @param[in]  dataFD_nHops Number of timeslots dataFD is allocated (the max)
 * @param[out] dataFD       Frequency-domain output; #AFSTFT_FDDATA_FORMAT
 */
void afSTFT_forward_knownDimensions(void * const hSTFT,
                                    float** dataTD,
                                    int framesize,
                                    int dataFD_nCH,
                                    int dataFD_nHops,
                                    float_complex*** dataFD);

/**
 * Performs forward afSTFT transform (flattened arrays)
 *
 * @param[in]  hSTFT     afSTFT handle
 * @param[in]  dataTD    Time-domain input; FLAT: nCHin x framesize
 * @param[in]  framesize Frame size of time-domain data
 * @param[out] dataFD    Frequency-domain output; FLAT: #AFSTFT_FDDATA_FORMAT
 */
void afSTFT_forward_flat(void * const hSTFT,
                         float* dataTD,
                         int framesize,
                         float_complex* dataFD);

/**
 * Performs backward afSTFT transform
 *
 * @param[in]  hSTFT     afSTFT handle
 * @param[in]  dataFD    Frequency-domain input; #AFSTFT_FDDATA_FORMAT
 * @param[in]  framesize Frame size of time-domain data
 * @param[out] dataTD    Time-domain output;  nCHout x framesize
 */
void afSTFT_backward(void * const hSTFT,
                     float_complex*** dataFD,
                     int framesize,
                     float** dataTD);

/**
 * Performs backward afSTFT transform (dataFD dimensions are known)
 *
 * @note If the dimensions of dataFD are known, then this function will use the
 *       same dataFD format as used by afSTFT_backward(), but with the speed
 *       of afSTFT_backward_flat()
 *
 * @param[in]  hSTFT        afSTFT handle
 * @param[in]  dataFD       Frequency-domain input; #AFSTFT_FDDATA_FORMAT
 * @param[in]  framesize    Frame size of time-domain data
 * @param[in]  dataFD_nCH   Number of channels dataFD is allocated (the max)
 * @param[in]  dataFD_nHops Number of timeslots dataFD is allocated (the max)
 * @param[out] dataTD       Time-domain output;  nCHout x framesize
 */
void afSTFT_backward_knownDimensions(void * const hSTFT,
                                     float_complex*** dataFD,
                                     int framesize,
                                     int dataFD_nCH,
                                     int dataFD_nHops,
                                     float** dataTD);

/**
 * Performs backward afSTFT transform (flattened arrays)
 *
 * @param[in]  hSTFT     afSTFT handle
 * @param[in]  dataFD    Frequency-domain input; FLAT: #AFSTFT_FDDATA_FORMAT
 * @param[in]  framesize Frame size of time-domain data
 * @param[out] dataTD    Time-domain output; FLAT: nCHout x framesize
 */
void afSTFT_backward_flat(void * const hSTFT,
                          float_complex* dataFD,
                          int framesize,
                          float* dataTD);

/**
 * Re-allocates memory to support a change in the number of input/output
 * channels
 *
 * @param[in] hSTFT      afSTFT handle
 * @param[in] new_nCHin  New number of input channels
 * @param[in] new_nCHout New number of output channels
 */
void afSTFT_channelChange(void * const hSTFT,
                          int new_nCHin,
                          int new_nCHout);

/** Flushes time-domain buffers with zeros */
void afSTFT_clearBuffers(void * const hSTFT);

/** Returns number of frequency bands */
int afSTFT_getNBands(void * const hSTFT);

/**
 * Returns current processing delay, in samples
 *
 * @note The afSTFT filterbank delay is broken down into the following:
 *          analysis delay:         5*hopsize
 *        + hybrid-filtering delay: 3*hopsize     (or 0, if it is disabled)
 *        + synthesis delay         4*hopsize
 *
 *       If the low-delay mode is enabled, it is instead:
 *          analysis delay:         2*hopsize
 *        + hybrid-filtering delay: 3*hopsize     (or 0, if it is disabled)
 *        + synthesis delay         2*hopsize
 */
int afSTFT_getProcDelay(void * const hSTFT);

/** Returns current frequency vector */
void afSTFT_getCentreFreqs(void * const hSTFT,
                           float fs,
                           int nBands,
                           float* freqVector);

/**
 * Converts FIR filters into Filterbank Coefficients by passing them through
 * the afSTFT filterbank
 *
 * @param[in]  hIR        Time-domain FIR; FLAT: N_dirs x nCH x ir_len
 * @param[in]  N_dirs     Number of FIR sets
 * @param[in]  nCH        Number of channels per FIR set
 * @param[in]  ir_len     Length of the FIR
 * @param[in]  hopSize    Hop size
 * @param[in]  LDmode     0: disabled, 1:enabled
 * @param[in]  hybridmode 0: disabled, 1:enabled
 * @param[out] hFB        The FIRs as Filterbank coefficients;
 *                        FLAT: N_bands x nCH x N_dirs
 */
void afSTFT_FIRtoFilterbankCoeffs(/* Input Arguments */
                                  float* hIR,
                                  int N_dirs,
                                  int nCH,
                                  int ir_len,
                                  int hopSize,
                                  int LDmode,
                                  int hybridmode,
                                  /* Output Arguments */
                                  float_complex* hFB);


#ifdef __cplusplus
}/* extern "C" */
#endif /* __cplusplus */

#endif /* __afSTFTlib_INCLUDED__ */
