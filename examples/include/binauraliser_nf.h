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
 * binauraliser.h (include), binauraliser_internal.h, binauraliser_nf.h (include), binauraliser_nf_internal.h, binauraliser.c, binauraliser_internal.c, binauraliser_nf_internal.c, binauraliser_nf.c
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

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/* ========================================================================== */
/*                             Presets + Constants                            */
/* ========================================================================== */

/** Available interpolation modes */
typedef enum {
    INTERP_TRI = 1, /**< Triangular interpolation */
    INTERP_TRI_PS   /**< Triangular interpolation (with phase-simplification) */
}INTERP_MODES;


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


/* ========================================================================== */
/*                                Get Functions                               */
/* ========================================================================== */

/**
 * Returns the processing framesize (i.e., number of samples processed with
 * every _process() call )
 */
int binauraliser_getFrameSize(void);

/** Returns current codec status codec status (see #CODEC_STATUS enum) */
CODEC_STATUS binauraliser_getCodecStatus(void* const hBin);

/**
 * (Optional) Returns current intialisation/processing progress, between 0..1
 * - 0: intialisation/processing has started
 * - 1: intialisation/processing has ended
 */
float binauraliser_getProgressBar0_1(void* const hBin);

/**
 * (Optional) Returns current intialisation/processing progress text
 *
 * @note "text" string should be (at least) of length:
 *       #PROGRESSBARTEXT_CHAR_LENGTH
 */
void binauraliser_getProgressBarText(void* const hBin, char* text);

/** Returns the source azimuth for a given index, in DEGREES */
float binauraliser_getSourceAzi_deg(void* const hBin, int index);

/** Returns the source elevation for a given index, in DEGREES */
float binauraliser_getSourceElev_deg(void* const hBin, int index);

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

/**
 * Returns the number of inputs/sources in the current layout
 */
int binauraliser_getNumSources(void* const hBin);

/** Returns the maximum number of input sources supported by binauraliser */
int binauraliser_getMaxNumSources(void);

/** Returns the number of ears possessed by the average homo sapien */
int binauraliser_getNumEars(void);

/** Returns the number of directions in the currently used HRIR set */
int binauraliser_getNDirs(void* const hBin);

/**
 * Returns the number of triangular groupings (faces) returned by the Convex
 * Hull
 */
int binauraliser_getNTriangles(void* const hBin);

/** Returns the HRIR/HRTF azimuth for a given index, in DEGREES */
float binauraliser_getHRIRAzi_deg(void* const hBin, int index);

/** Returns the HRIR/HRTF elevation for a given index, in DEGREES */
float binauraliser_getHRIRElev_deg(void* const hBin, int index);

/** Returns the length of HRIRs in time-domain samples */
int binauraliser_getHRIRlength(void* const hBin);

/** Returns the HRIR sample rate */
int binauraliser_getHRIRsamplerate(void* const hBin);

/**
 * Returns the value of a flag used to dictate whether the default HRIRs in the
 * Spatial_Audio_Framework should be used (1), or a custom HRIR set loaded via a
 * SOFA file (0).
 *
 * @note If the custom set fails to load correctly, binauraliser will revert to
 *       the defualt set, so this will be '1'
 */
int binauraliser_getUseDefaultHRIRsflag(void* const hBin);

/**
 * Returns the file path for a .sofa file.
 *
 * @note If the custom set fails to load correctly, binauraliser will revert to
 *       the defualt set. Use 'binauraliser_getUseDefaultHRIRsflag()' to check
 *       if loading was successful.
 *
 * @param[in] hBin binauraliser handle
 * @returns        File path to .sofa file (WITH file extension)
 */
char* binauraliser_getSofaFilePath(void* const hBin);

/**
 * Returns the flag indicating whether the diffuse-field EQ applied to the HRTFs
 * is enabled (1) or disabled (0).
 */
int binauraliser_getEnableHRIRsDiffuseEQ(void* const hBin);

/** Returns the DAW/Host sample rate */
int binauraliser_getDAWsamplerate(void* const hBin);

/**
 * Returns the flag value which dictates whether to enable/disable sound-field
 * rotation (0: disabled, 1: enabled)
 */
int binauraliser_getEnableRotation(void* const hBin);

/** Returns the 'yaw' rotation angle, in DEGREES */
float binauraliser_getYaw(void* const hBin);

/** Returns the 'pitch' rotation angle, in DEGREES */
float binauraliser_getPitch(void* const hBin);

/** Returns the 'roll' rotation angle, in DEGREES */
float binauraliser_getRoll(void* const hBin);

/**
 * Returns a flag as to whether to "flip" the sign of the current 'yaw' angle
 * (0: do not flip sign, 1: flip the sign)
 */
int binauraliser_getFlipYaw(void* const hBin);

/**
 * Returns a flag as to whether to "flip" the sign of the current 'pitch' angle
 * (0: do not flip sign, 1: flip the sign)
 */
int binauraliser_getFlipPitch(void* const hBin);

/**
 * Returns a flag as to whether to "flip" the sign of the current 'roll' angle
 * (0: do not flip sign, 1: flip the sign)
 */
int binauraliser_getFlipRoll(void* const hBin);

/**
 * Returns a flag as to whether to use "yaw-pitch-roll" (0) or "roll-pitch-yaw"
 * (1) rotation order.
 */
int binauraliser_getRPYflag(void* const hBin);

/** NOT IMPLEMENTED YET */
int binauraliser_getInterpMode(void* const hBin);

/**
 * Returns the processing delay in samples (may be used for delay compensation
 * purposes)
 */
int binauraliser_getProcessingDelay(void);


#ifdef __cplusplus
} /* extern "C" { */
#endif /* __cplusplus */

#endif /* __BINAURALISER_H_INCLUDED__ */
