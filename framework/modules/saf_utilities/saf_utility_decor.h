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
 *@addtogroup Utilities
 *@{
 * @file saf_utility_decor.h
 * @brief A collection of signal decorrelators
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
#include "saf_utility_complex.h"

/** Lattice all-pass filter coeffs (numerator) for 256 channels, 20th order */
extern const float __lattice_coeffs_o20[256][20];
/** Lattice all-pass filter coeffs (numerator) for 256 channels, 18th order */
extern const float __lattice_coeffs_o18[256][18];
/** Lattice all-pass filter coeffs (numerator) for 256 channels, 16th order */
extern const float __lattice_coeffs_o16[256][16];
/** Lattice all-pass filter coeffs (numerator) for 256 channels, 15th order */
extern const float __lattice_coeffs_o15[256][15];
/** Lattice all-pass filter coeffs (numerator) for 256 channels, 14th order */
extern const float __lattice_coeffs_o14[256][14];
/** Lattice all-pass filter coeffs (numerator) for 256 channels, 12th order */
extern const float __lattice_coeffs_o12[256][12];
/** Lattice all-pass filter coeffs (numerator) for 256 channels, 10th order */
extern const float __lattice_coeffs_o10[256][10];
/** Lattice all-pass filter coeffs (numerator) for 256 channels, 8th order */
extern const float __lattice_coeffs_o8[256][8];
/** Lattice all-pass filter coeffs (numerator) for 256 channels, 6th order */
extern const float __lattice_coeffs_o6[256][6];
/** Lattice all-pass filter coeffs (numerator) for 256 channels, 4th order */
extern const float __lattice_coeffs_o4[256][4];
/** Lattice all-pass filter coeffs (numerator) for 256 channels, 3rd order */
extern const float __lattice_coeffs_o3[256][3];
/** Lattice all-pass filter coeffs (numerator) for 256 channels, 2nd order */
extern const float __lattice_coeffs_o2[256][2];

/**
 * Returns delay values for multiple channels per frequency, such that once
 * applied to an input signal (via simple frequency-dependent delay lines), the
 * resulting signal is decorrelated w.r.t the original
 *
 * @note This is a very basic algorithm and sounds particulary bad for transient
 *       signals. Consider using a transient detector to "duck" the decorrelated
 *       signal during such transients, to improve signal fidelity. See e.g.
 *       transientDucker_create()
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
void getDecorrelationDelays(/* Input Arguments */
                            int nChannels,
                            float* freqs,
                            int nFreqs,
                            float fs,
                            int maxTFdelay,
                            int hopSize,
                            /* Output Arguments */
                            int* delayTF);

/**
 * Returns quick and dirty exponentially decaying noise bursts
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
void synthesiseNoiseReverb(/* Input Arguments */
                           int nChannels,
                           float fs,
                           float* t60,
                           float* fcen_oct,
                           int nBands,
                           int flattenFLAG,
                           /* Output Arguments */
                           float** rir_filt,
                           int* rir_len);

/**
 * Creates an instance of the lattice all-pass-filter-based multi-channel
 * signal decorrelator
 *
 * This decorrelator is intended for decorrelating signals in the time-frequency
 * domain, and is therefore well-suited for audio coding [1] or Direction Audio
 * Coding (DirAC) [2] purposes.
 *
 * @note Contrary to how the decorrelator is defined in the standard [3], it is
 *       recommended to use lower orders. Note that one may cascade 2 instances
 *       (each with lower orders), to attain the desired decorrelation without
 *       introducing instabilities, which accompanies the use of the much higher
 *       filter orders. Furthermore, the coefficients found in
 *       saf_utility_latticeCoeffs.h differ from those found in the standard [3]
 *       so the output will also be numerically different.
 * @test test__latticeDecorrelator()
 *
 * @param[in] phDecor      (&) address of lattice decorrelator handle
 * @param[in] nCH          Number of channels
 * @param[in] orders       Lattice all-pass filter orders (2,3,4,6,8,10,12,14,
 *                         15,16 18, or 20) per band grouping; (nCutoffs+1) x 1
 * @param[in] freqCutoffs  Frequency cut-offs defining the band groupings;
 *                         nCutoffs x 1
 * @param[in] fixedDelays  Fixed time-frequency hop delays; (nCutoffs+1) x 1
 * @param[in] nCutoffs     Number of cutoff frequencies
 * @param[in] freqVector   Frequency vector; nBands x 1
 * @param[in] lookupOffset Optional offset for look-up tables (set to 0 if using
 *                         just one instance)
 * @param[in] nBands       Number of bands
 *
 * @see [1] Herre, J., Kjo"rling, K., Breebaart, J., Faller, C., Disch, S.,
 *          Purnhagen, H., Koppens, J., Hilpert, J., Ro"den, J., Oomen, W. and
 *          Linzmeier, K., 2008. MPEG surround-the ISO/MPEG standard for
 *          efficient and compatible multichannel audio coding. Journal of the
 *          Audio Engineering Society, 56(11), pp.932--955
 * @see [2] Pulkki, V., 2007. Spatial sound reproduction with directional audio
 *          coding. Journal of the Audio Engineering Society, 55(6), pp.503-516
 * @see [3] ISO/IEC FCD 23003-1, "MPEG-D (MPEG Audio Technologies)--Part 1: MPEG
 *          Surround" International Standards Organization, Geneva, Switzerland
 *          (2006)
 */
void latticeDecorrelator_create(/* Input Arguments */
                                void** phDecor,
                                int nCH,
                                int* orders,
                                float* freqCutoffs,
                                int* fixedDelays,
                                int nCutoffs,
                                float* freqVector,
                                int lookupOffset,
                                int nBands);

/**
 * Destroys an instance of the lattice all-pass-filter-based multi-channel
 * signal decorrelator
 *
 * @param[in] phDecor (&) address of lattice decorrelator handle
 */
void latticeDecorrelator_destroy(/* Input Arguments */
                                 void** phDecor);

/**
 * Applies the lattice all-pass-filter-based multi-channel signal decorrelator
 *
 * @param[in]  hDecor     lattice decorrelator handle
 * @param[in]  inFrame    input frame; nBands x nCH x nTimeSlots
 * @param[in]  nTimeSlots Number of time slots per frame
 * @param[out] decorFrame decorrelated frame; nBands x nCH x nTimeSlots
 */
void latticeDecorrelator_apply(/* Input Arguments */
                               void* hDecor,
                               float_complex*** inFrame,
                               int nTimeSlots,
                               /* Output Arguments */
                               float_complex*** decorFrame);

/**
 * Creates an instance of the transient ducker
 *
 * @param[in] phDucker (&) address of ducker handle
 * @param[in] nCH       Number of channels
 * @param[in] nBands    Number of frequency bands
 */
void transientDucker_create(/* Input Arguments */
                            void** phDucker,
                            int nCH,
                            int nBands);

/**
 * Destroys an instance of the transient ducker
 *
 * @param[in] phDucker (&) address of ducker handle
 */
void transientDucker_destroy(/* Input Arguments */
                             void** phDucker);

/**
 * Applies the transient ducker
 *
 * @param[in]  hDucker    ducker handle
 * @param[in]  inFrame    input frame; nBands x nCH x nTimeSlots
 * @param[in]  nTimeSlots Number of time slots per frame
 * @param[in]  alpha      alpha value [0,1];
 * @param[in]  beta       beta value [0,1];
 * @param[out] outFrame   ducked frame; nBands x nCH x nTimeSlots
 */
void transientDucker_apply(/* Input Arguments */
                           void* hDucker,
                           float_complex*** inFrame,
                           int nTimeSlots,
                           float alpha,
                           float beta,
                           /* Output Arguments */
                           float_complex*** outFrame);


#ifdef __cplusplus
}/* extern "C" */
#endif /* __cplusplus */

#endif /* SAF_DECOR_H_INCLUDED */

/**@} */ /* doxygen addtogroup Utilities */
