/*
 * Copyright 2017-2018 Leo McCormack
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
 * @example binauraliser_nf.h
 * @brief Convolves input audio (up to 64 channels) with interpolated HRTFs in
 *        the time-frequency domain, and applies optional near-field binaural
 *        filtering.
 *
 * ### Files
 * binauraliser.h (include), binauraliser_internal.h, binauraliser_nf.h (include),
 * binauraliser_nf_internal.h, binauraliser.c, binauraliser_internal.c,
 * binauraliser_nf_internal.c, binauraliser_nf.c
 * ### Include Header
 */

/**
 * @file: binauraliser_nf.h
 * @brief Convolves input audio (up to 64 channels) with interpolated HRTFs in
 *        the time-frequency domain, and applies optional near-field binaural
 *        filtering.
 *
 * The HRTFs are interpolated by applying amplitude-preserving VBAP gains to the
 * HRTF magnitude responses and inter-aural time differences (ITDs)
 * individually, before being re-combined. The example also allows the user to
 * specify an external SOFA file for the convolution, and rotations of the
 * source directions to accomodate head-tracking.
 *
 * @author Leo McCormack
 * @date 25.09.2017
 * @license ISC
 */

#ifndef __BINAURALISER_H_INCLUDED__
#define __BINAURALISER_H_INCLUDED__

#include "_common.h"
#include "binauraliser.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/* ========================================================================== */
/*                             Presets + Constants                            */
/* ========================================================================== */



/* ========================================================================== */
/*                               Main Functions                               */
/* ========================================================================== */

/**
 * Creates an instance of the binauraliser
 *
 * @param[in] phBin (&) address of binauraliser handle
 */
void binauraliserNF_create(void** const phBin);

/**
 * Destroys an instance of the binauraliser
 *
 * @param[in] phBin (&) address of binauraliser handle
 */
void binauraliserNF_destroy(void** const phBin);

/**
 * Initialises an instance of binauraliser with default settings
 *
 * @warning This should not be called while _process() is on-going!
 *
 * @param[in] hBin       binauraliser handle
 * @param[in] samplerate Host samplerate.
 */
void binauraliserNF_init(void* const hBin,
                       int samplerate);

/**
 * Binauralises the input signals at the user specified directions
 *
 * @param[in] hBin      binauraliser handle
 * @param[in] inputs    Input channel buffers; 2-D array: nInputs x nSamples
 * @param[in] outputs   Output channel buffers; 2-D array: nOutputs x nSamples
 * @param[in] nInputs   Number of input channels
 * @param[in] nOutputs  Number of output channels
 * @param[in] nSamples  Number of samples in 'inputs'/'output' matrices
 */
void binauraliserNF_process(void* const hBin,
                          const float *const * inputs,
                          float** const outputs,
                          int nInputs,
                          int nOutputs,
                          int nSamples);


/* ========================================================================== */
/*                                Set Functions                               */
/* ========================================================================== */

/**
 * Sets the panning distance for a specific channel index, in METERS
 */
void binauraliserNF_setSourceDist_m(void* const hBin,
                                    int index,
                                    float newDist_m);

void binauraliserNF_setInputConfigPreset(void* const hBin,
                                         int newPresetID);

void binauraliserNF_resetSourceDistances(void* const hBin);
/* ========================================================================== */
/*                                Get Functions                               */
/* ========================================================================== */

/**
 * Returns the source elevation for a given index, in METERS
 */
float binauraliser_getSourceDist_m(void* const hBin, int index);

/**
* Returns the distance considered to be the far field (beyond which no near field filtering is applied), in METERS
*/
float binauraliser_getFarfieldThresh_m(void* const hBin);

/**
* Returns the scaling factor to give the far field threshold headroom (useful for UI range limits)
*/
float binauraliser_getFarfieldHeadroom(void* const hBin);

/**
* Returns the minimum distance possible for near field filter, in METERS
*/
float binauraliser_getNearfieldLimit_m(void* const hBin);

#ifdef __cplusplus
} /* extern "C" { */
#endif /* __cplusplus */

#endif /* __BINAURALISER_H_INCLUDED__ */
