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
 *@addtogroup Utilities
 *@{
 * @file saf_utility_pitch.h
 * @brief A collection of pitch shifting algorithms
 *
 * @author Leo McCormack
 * @date 04.05.2020
 * @license ISC
 */

#ifndef SAF_PITCH_H_INCLUDED
#define SAF_PITCH_H_INCLUDED

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/* ========================================================================== */
/*                              SMB PitchShifter                              */
/* ========================================================================== */

/**
 * Creates an instance of SMB PitchShifter
 *
 * @note Higher FFT sizes will permit more drastic pitch shifts. Increasing the
 *       Oversampling factor will increase latency, but also improve signal
 *       fidelity.
 *
 * @test test__smb_pitchShifter()
 *
 * @param[in] hSmb         (&) address of smb pitchShifter handle
 * @param[in] nCH          number of channels
 * @param[in] fftFrameSize FFT size
 * @param[in] osamp        Oversampling/overlapping factor
 * @param[in] sampleRate   Sampling rate, Hz
 */
void smb_pitchShift_create(/* Input Arguments */
                           void** hSmb,
                           int nCH,
                           int fftFrameSize,
                           int osamp,
                           float sampleRate);

/**
 * Destroys an instance of SMB PitchShifter
 *
 * @param[in] hSmb (&) address of smb pitchShifter handle
 */
void smb_pitchShift_destroy(/* Input Arguments */
                            void ** const hSmb);

/**
 * Performs pitch shifting of the input signals, while retaining the same time
 * duration as the original using the algorithm detailed in [1]
 *
 * This implementation was orginally written by Stephan M. Bernsee (c) 1999-2015
 * distributed under the WOL license. It has been modified to better work with
 * frame-by-frame processing. It also supports multiple input channels and
 * employs saf_utility_fft.h and saf_utility_veclib.h for additional run-time
 * optimisations.
 *
 * @param[in]  hSmb       (&) smb pitchShifter handle
 * @param[in]  pitchShift Pitch shift factor, 0.5: down 1 octave, 1: no shift,
 *                        2: up 1 octave
 * @param[in]  frameSize  Number of samples in frame
 * @param[in]  inFrame    Input frame;  FLAT: nCH x frameSize
 * @param[out] outFrame   Output frame; FLAT: nCH x frameSize
 *
 * @see [1] http://blogs.zynaptiq.com/bernsee/pitch-shifting-using-the-ft/
 *          Copyright 1999-2015 Stephan M. Bernsee, The Wide Open License (WOL)
 */
void smb_pitchShift_apply(/* Input Arguments */
                          void* hSmb,
                          float pitchShift,
                          int frameSize,
                          float *inFrame,
                          /* Output Arguments */
                          float *outFrame);


#ifdef __cplusplus
}/* extern "C" */
#endif /* __cplusplus */

#endif /* SAF_PITCH_H_INCLUDED */

/**@} */ /* doxygen addtogroup Utilities */
