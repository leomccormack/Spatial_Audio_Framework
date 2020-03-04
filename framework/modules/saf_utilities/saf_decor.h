/*
 * Copyright 2018 Leo McCormack
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
 * @file saf_decor.h
 * @brief Collection of signal decorrelators.
 *
 * @author Leo McCormack
 * @date 30.07.2018 
 */

#ifndef SAF_DECOR_H_INCLUDED
#define SAF_DECOR_H_INCLUDED

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

/**
 * Returns delay values for multiple channels per frequency, such that once
 * applied to an input signal (via simple frequency-dependent delay lines), the
 * resulting signal is decorrelated w.r.t the original.
 *
 * @note This is a very basic algorithm and sounds particulary bad for transient
 *       signals. Consider using a transient detector to "duck" the decorrelated
 *       signal during such transients, to improve signal fidelity.
 *
 * @param[in]  nChannels  Number of channels
 * @param[in]  freqs      Centre frequencies; nFreqs x 1
 * @param[in]  nFreqs     Number of elements in frequency vector
 * @param[in]  fs         Sampling rate
 * @param[in]  maxTFdelay Max number of time-slots to delay
 * @param[in]  hopSize    STFT hop size
 * @param[out] delayTF    The resulting time delays per channel and frequency;
 *                        FLAT: nFreq x nChannels
 */
void getDecorrelationDelays(int nChannels,
                            float* freqs,
                            int nFreqs,
                            float fs,
                            int maxTFdelay,
                            int hopSize,
                            int* delayTF);

/**
 * Returns quick and dirty exponentially decaying noise bursts.
 *
 * With long T60 times, it can be used to approximate the late reverberation
 * tail of room impulse responses. With much shorter t60 times, it can be used
 * for decorrelation purposes.
 *
 * @param[in]  nChannels   Number of channels
 * @param[in]  fs          Sampling rate
 * @param[in]  t60         T60 times (in seconds) per octave band; nBands x 1
 * @param[in]  fcen_oct    Octave band centre frequencies; nBands x 1
 * @param[in]  nBands      Number of octave bands
 * @param[in]  flattenFLAG '0' nothing, '1' flattens the magnitude response to
 *                         unity
 * @param[out] rir_filt    (&) the shaped noise bursts;
 *                         FLAT: nChannels x rir_len
 * @param[out] rir_len     (&) length of filters, in samples
 */
void synthesiseNoiseReverb(int nChannels,
                           float fs,
                           float* t60,
                           float* fcen_oct,
                           int nBands,
                           int flattenFLAG,
                           float** rir_filt,
                           int* rir_len);


#ifdef __cplusplus
}/* extern "C" */
#endif /* __cplusplus */

#endif /* SAF_DECOR_H_INCLUDED */
