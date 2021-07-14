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
 * @example rotator.h
 * @brief A basic spherical harmonic signals rotator
 * 
 * ### Files
 * rotator.h (include), rotator_internal.h, rotator.c, rotator_internal.c
 * ### Include Header
 */

/**
 * @file rotator.h
 * @brief A basic spherical harmonic/ Ambisonic signals rotator, based on the
 *        recursive approach detailed in [1]
 *
 * @test test__saf_example_rotator()
 *
 * @see [1] Ivanic, J., Ruedenberg, K. (1998). Rotation Matrices for Real
 *          Spherical Harmonics. Direct Determination by Recursion Page:
 *          Additions and Corrections. Journal of Physical Chemistry A, 102(45),
 *          9099?9100.
 *
 * @author Leo McCormack
 * @date 02.11.2017
 * @license ISC
 */

#ifndef __ROTATOR_H_INCLUDED__
#define __ROTATOR_H_INCLUDED__

#include "_common.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/* ========================================================================== */
/*                               Main Functions                               */
/* ========================================================================== */

/**
 * Creates an instance of rotator
 *
 * @param[in] phRot (&) address of rotator handle
 */
void rotator_create(void** const phRot);

/**
 * Destroys an instance of rotator
 *
 * @param[in] phRot (&) address of rotator handle
 */
void rotator_destroy(void** const phRot);

/**
 * Initialises an instance of rotator with default settings
 *
 * @param[in] hRot       rotator handle
 * @param[in] samplerate Host samplerate.
 */
void rotator_init(void* const hRot,
                  int samplerate);

/**
 * Rotates the input spherical harmonic signals
 *
 * @param[in] hRot     rotator handle
 * @param[in] inputs   Input channel buffers; 2-D array: nInputs x nSamples
 * @param[in] outputs  Output channel buffers; 2-D array: nOutputs x nSamples
 * @param[in] nInputs  Number of input channels
 * @param[in] nOutputs Number of output channels
 * @param[in] nSamples Number of samples in 'inputs'/'output' matrices
 */
void rotator_process(void* const hRot,
                     const float *const * inputs,
                     float** const outputs,
                     int nInputs,
                     int nOutputs,
                     int nSamples);


/* ========================================================================== */
/*                                Set Functions                               */
/* ========================================================================== */

/**
 * Returns the processing framesize (i.e., number of samples processed with
 * every _process() call )
 */
int rotator_getFrameSize(void);

/** Sets the 'yaw' rotation angle, in DEGREES */
void rotator_setYaw(void* const hRot, float newYaw);

/** Sets the 'pitch' rotation angle, in DEGREES */
void rotator_setPitch(void* const hRot, float newPitch);

/** Sets the 'roll' rotation angle , in DEGREES */
void rotator_setRoll(void* const hRot, float newRoll);

/** Sets the quaternion 'W' value [-1..1] */
void rotator_setQuaternionW(void* const hRot, float newValue);

/** Sets the quaternion 'X' value [-1..1] */
void rotator_setQuaternionX(void* const hRot, float newValue);

/** Sets the quaternion 'Y' value [-1..1] */
void rotator_setQuaternionY(void* const hRot, float newValue);

/** Sets the quaternion 'Z' value [-1..1] */
void rotator_setQuaternionZ(void* const hRot, float newValue);

/**
 * Sets a flag as to whether to "flip" the sign of the current 'yaw' angle
 * (0: do not flip sign, 1: flip the sign)
 */
void rotator_setFlipYaw(void* const hRot, int newState);

/**
 * Sets a flag as to whether to "flip" the sign of the current 'pitch' angle
 * (0: do not flip sign, 1: flip the sign)
 */
void rotator_setFlipPitch(void* const hRot, int newState);

/**
 * Sets a flag as to whether to "flip" the sign of the current 'roll' angle
 * (0: do not flip sign, 1: flip the sign)
 */
void rotator_setFlipRoll(void* const hRot, int newState);

/**
 * Sets a flag as to whether to invert the quaternion used for rotation
 * (0: do not flip sign, 1: flip the sign)
 */
void rotator_setFlipQuaternion(void* const hRot, int newState);

/**
 * Sets the Ambisonic channel ordering convention to decode with, in order to
 * match the convention employed by the input signals (see #CH_ORDER enum)
 */
void rotator_setChOrder(void* const hRot, int newOrder);

/**
 * Sets the Ambisonic normalisation convention to decode with, in order to match
 * with the convention employed by the input signals (see #NORM_TYPES enum)
 */
void rotator_setNormType(void* const hRot, int newType);

/** Sets the input/output order (see #SH_ORDERS enum) */
void rotator_setOrder(void* const hRot, int newOrder);

/**
 * Sets a flag as to whether to use "yaw-pitch-roll" (0) or "roll-pitch-yaw" (1)
 * rotation order.
 */
void rotator_setRPYflag(void* const hRot, int newState);


/* ========================================================================== */
/*                                Get Functions                               */
/* ========================================================================== */

/** Returns the 'yaw' rotation angle, in DEGREES */
float rotator_getYaw(void* const hRot);

/** Returns the 'pitch' rotation angle, in DEGREES */
float rotator_getPitch(void* const hRot);

/** Returns the 'roll' rotation angle, in DEGREES */
float rotator_getRoll(void* const hRot);

/** Returns the quaternion 'W' value [-1..1] */
float rotator_getQuaternionW(void* const hRot);

/** Returns the quaternion 'X' value [-1..1] */
float rotator_getQuaternionX(void* const hRot);

/** Returns the quaternion 'Y' value [-1..1] */
float rotator_getQuaternionY(void* const hRot);

/** Returns the quaternion 'Z' value [-1..1] */
float rotator_getQuaternionZ(void* const hRot);

/**
 * Returns a flag as to whether to "flip" the sign of the current 'yaw' angle
 * (0: do not flip sign, 1: flip the sign)
 */
int rotator_getFlipYaw(void* const hRot);

/**
 * Returns a flag as to whether to "flip" the sign of the current 'pitch' angle
 * (0: do not flip sign, 1: flip the sign)
 */
int rotator_getFlipPitch(void* const hRot);

/** 
 * Returns a flag as to whether to "flip" the sign of the current 'roll' angle
 * (0: do not flip sign, 1: flip the sign)
 */
int rotator_getFlipRoll(void* const hRot);

/**
 * Returns a flag as to whether to invert the quaternion used for rotation
 * (0: do not flip sign, 1: flip the sign)
 */
int rotator_getFlipQuaternion(void* const hRot);

/**
 * Returns a flag as to whether to use "yaw-pitch-roll" (0) or "roll-pitch-yaw"
 * (1) rotation order.
 */
int rotator_getRPYflag(void* const hRot);

/**
 * Returns the Ambisonic channel ordering convention currently being used to
 * decode with, which should match the convention employed by the input signals
 * (see #CH_ORDER enum)
 */
int rotator_getChOrder(void* const hRot);

/**
 * Returns the Ambisonic normalisation convention currently being usedto decode
 * with, which should match the convention employed by the input signals
 * (see #NORM_TYPES enum)
 */
int rotator_getNormType(void* const hRot);

/** Returns the input/output order (see #SH_ORDERS enum) */
int rotator_getOrder(void* const hRot);

/**
 * Returns the number of spherical harmonic signals required by the current
 * input/output order: (current_order+1)^2 
 */
int rotator_getNSHrequired(void* const hRot);

/**
 * Returns the processing delay in samples (may be used for delay compensation
 * features)
 */
int rotator_getProcessingDelay(void);
    
    
#ifdef __cplusplus
} /* extern "C" { */
#endif /* __cplusplus */

#endif /* __ROTATOR_H_INCLUDED__ */
