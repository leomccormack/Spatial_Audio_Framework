/*
 * Copyright 2022 Michael McCrea, Leo McCormack
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
 *        filtering, as described in [1].
 *
 * @see [1] S. Spagnol, E. Tavazzi, and F. Avanzini, “Distance rendering and
 *          perception of nearby virtual sound sources with a near-field filter
 *          model,” Applied Acoustics, vol. 115, pp. 61–73, Jan. 2017,
 *          doi: 10.1016/j.apacoust.2016.08.015.
 *
 * ### Files
 * binauraliser.h (include), binauraliser_internal.h,
 * binauraliser_nf.h (include), binauraliser_nf_internal.h,
 * binauraliser.c, binauraliser_internal.c,
 * binauraliser_nf_internal.c, binauraliser_nf.c
 * ### Include Header
 */

/**
 * @file: binauraliser_nf.h
 * @brief Convolves input audio (up to 64 channels) with interpolated HRTFs in
 *        the time-frequency domain, and applies optional near-field binaural
 *        filtering, as described in [1].
 *
 * The HRTFs are interpolated by applying amplitude-preserving VBAP gains to the
 * HRTF magnitude responses and inter-aural time differences (ITDs)
 * individually, before being re-combined. The example also allows the user to
 * specify an external SOFA file for the convolution, and rotations of the
 * source directions to accomodate head-tracking.
 *
 * @see [1] S. Spagnol, E. Tavazzi, and F. Avanzini, “Distance rendering and
 *          perception of nearby virtual sound sources with a near-field filter
 *          model,” Applied Acoustics, vol. 115, pp. 61–73, Jan. 2017,
 *          doi: 10.1016/j.apacoust.2016.08.015.
 *
 * @author Michael McCrea, Leo McCormack
 * @date 22.02.2022
 * @license ISC
 */

#ifndef __BINAURALISER_NF_H_INCLUDED__
#define __BINAURALISER_NF_H_INCLUDED__

#include <binauraliser.h>

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/* ========================================================================== */
/*                               Main Functions                               */
/* ========================================================================== */

/**
 * Creates an instance of the binauraliser
 *
 * @param[in] phBin (&) address of binauraliserNF handle
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
 * @param[in] hBin       binauraliserNF handle
 * @param[in] samplerate Host samplerate.
 */
void binauraliserNF_init(void* const hBin,
                         int samplerate);

/**
 * Intialises the codec variables, based on current global/user parameters
 *
 * @note This function is fully threadsafe. It can even be called periodically
 *       via a timer on one thread, while calling _process() on another thread.
 *       Since, if a set function is called (that warrants a re-init), then a
 *       flag is triggered internally and the next time this function is called,
 *       it will wait until the current process() function has completed before
 *       reinitialising the relevant parameters. If the _initCodec() takes
 *       longer than the time it takes for process() to be called again, then
 *       process() is simply bypassed until the codec is ready.
 * @note This function does nothing if no re-initialisations are required.
 *
 * @param[in] hBin binauraliser handle
 */
/* See the source definition for a note on redundancy with binauraliser_initCodec. */
void binauraliserNF_initCodec(void* const hBin);

/**
 * Binauralises the input signals at the user specified directions
 *
 * @param[in] hBin      binauraliserNF handle
 * @param[in] inputs    Input channel buffers; 2-D array: nInputs x nSamples
 * @param[in] outputs   Output channel buffers; 2-D array: nOutputs x nSamples
 * @param[in] nInputs   Number of input channels
 * @param[in] nOutputs  Number of output channels
 * @param[in] nSamples  Number of samples in 'inputs'/'output' matrices
 */
void binauraliserNF_process(void* const hBin,
                            const float *const * inputs,
                            float* const* outputs,
                            int nInputs,
                            int nOutputs,
                            int nSamples);

/**
 * Alternate version of binauraliserNF_process() that performs frequency-domain
 * DVF filtering. Not used but kept for posterity.
 */
void binauraliserNF_processFD(void* const hBin,
                          const float *const * inputs,
                          float* const* outputs,
                          int nInputs,
                          int nOutputs,
                          int nSamples);


/* ========================================================================== */
/*                                Set Functions                               */
/* ========================================================================== */

/**
 * Sets the panning distance for a specific channel index, in METERS
 *
 * @param[in] hBin      binauraliserNF handle
 * @param[in] index     source index
 * @param[in] newDist_m source distance in meter
 */
void binauraliserNF_setSourceDist_m(void* const hBin,
                                    int index,
                                    float newDist_m);

/** Loads an input preset (see #SOURCE_CONFIG_PRESETS enum)
 *
 * @param[in] hBin        binauraliserNF handle
 * @param[in] newPresetID index of the source preset
 */
void binauraliserNF_setInputConfigPreset(void* const hBin,
                                         int newPresetID);


/* ========================================================================== */
/*                                Get Functions                               */
/* ========================================================================== */

/**
 * Returns the source elevation for a given index, in METERS
 *
 * @param[in] hBin  binauraliserNF handle
 * @param[in] index source index
 */
float binauraliserNF_getSourceDist_m(void* const hBin, int index);

/**
 * Returns the distance considered to be the far field (beyond which no near
 * field filtering is applied), in METERS
 */
float binauraliserNF_getFarfieldThresh_m(void* const hBin);

/**
 * Returns the scaling factor to give the far field threshold headroom (useful
 * for UI range limits)
 */
float binauraliserNF_getFarfieldHeadroom(void* const hBin);

/**
 * Returns the minimum distance possible for near field filter, in METERS
 */
float binauraliserNF_getNearfieldLimit_m(void* const hBin);

#ifdef __cplusplus
} /* extern "C" { */
#endif /* __cplusplus */

#endif /* __BINAURALISER_NF_H_INCLUDED__ */
