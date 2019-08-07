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

/*
 * Filename: saf_decor.h
 * ---------------------
 * A collection of signal decorrelators.
 *
 * Dependencies:
 *     none
 * Author, date created:
 *     Leo McCormack, 30.07.2018
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
    
/*
 * Function: getDecorrelationDelays
 * --------------------------------
 * Returns delay values for multiple channels per frequency, such that once
 * applied to an input signal (via simple frequency-dependent delay lines), the
 * resulting signal is decorrelated w.r.t the original.
 * Note: this is a very basic algorithm and sounds particulary bad for transient
 * signals. Consider using a transient detector to "duck" the deccorelated
 * signal during such transients, to improve signal fidelity.
 *
 * Input Arguments:
 *     nChannels  - number of channels
 *     freqs      - centre frequencies; nFreqs x 1
 *     nFreqs     - number of elements in frequency vector
 *     fs         - sampling rate
 *     maxTFdelay - max number of time-slots to delay
 *     hopSize    - STFT hop size
 * Output Arguments:
 *     delayTF    - the resulting time delays per channel and frequency;
 *                  FLAT: nFreq x nChannels
 */
void getDecorrelationDelays(int nChannels,
                            float* freqs,
                            int nFreqs,
                            float fs,
                            int maxTFdelay,
                            int hopSize,
                            int* delayTF);
 
 
#ifdef __cplusplus
}/* extern "C" */
#endif /* __cplusplus */

#endif /* SAF_DECOR_H_INCLUDED */
