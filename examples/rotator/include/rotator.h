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
 * @file rotator.h
 * @brief  A simple spherical harmonic domain rotator, based on the recursive
 *         approach detailed in [1].
 *
 * @see [1] Ivanic, J., Ruedenberg, K. (1998). Rotation Matrices for Real
 *          Spherical Harmonics. Direct Determination by Recursion Page:
 *          Additions and Corrections. Journal of Physical Chemistry A, 102(45),
 *          9099?9100.
 *
 * @author Leo McCormack
 * @date 02.11.2017
 */

#ifndef __ROTATOR_H_INCLUDED__
#define __ROTATOR_H_INCLUDED__

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */
    
/* ========================================================================== */
/*                             Presets + Constants                            */
/* ========================================================================== */

/**
 * Available input/output orders.
 */
typedef enum _ROTATOR_INPUT_ORDERS {
    INPUT_ORDER_FIRST=1, /**< First-order rotation (4 channel input/output) */
    INPUT_ORDER_SECOND,  /**< Second-order rotation (9 channel input/output) */
    INPUT_ORDER_THIRD,   /**< Third-order rotation (16 channel input/output) */
    INPUT_ORDER_FOURTH,  /**< Fourth-order rotation (25 channel input/output) */
    INPUT_ORDER_FIFTH,   /**< Fifth-order rotation (36 channel input/output) */
    INPUT_ORDER_SIXTH,   /**< Sixth-order rotation (49 channel input/output) */
    INPUT_ORDER_SEVENTH  /**< Seventh-order rotation (64 channel input/output)*/
    
} ROTATOR_INPUT_ORDERS;
    
/** Maximum supported Ambisonic order */
#define ROTATOR_MAX_SH_ORDER ( 7 )

/**
 * Available Ambisonic channel ordering conventions
 *
 * @note CH_FUMA only supported for 1st order input.
 */
typedef enum ROTATOR_CH_ORDER {
    CH_ACN = 1, /**< Ambisonic Channel Numbering (ACN) */
    CH_FUMA     /**< (Obsolete) Furse-Malham/B-format (WXYZ) */
    
} ROTATOR_CH_ORDER;

/** Number of channel ordering options */
#define ROTATOR_NUM_CH_ORDERINGS ( 2 )

/**
 * Available Ambisonic normalisation conventions
 *
 * @note NORM_FUMA only supported for 1st order input and does NOT have the
 *       1/sqrt(2) scaling on the omni.
 */
typedef enum _ROTATOR_NORM_TYPES {
    NORM_N3D = 1, /**< orthonormalised (N3D) */
    NORM_SN3D,    /**< Schmidt semi-normalisation (SN3D) */
    NORM_FUMA     /**< (Obsolete) Same as NORM_SN3D for 1st order */
    
} ROTATOR_NORM_TYPES;

/** Number of normalisation options */
#define ROTATOR_NUM_NORM_TYPES ( 3 )

#define ROTATOR_MAX_NUM_CHANNELS ( 64 )


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
 * Rotates the input spherical harmonic signals.
 *
 * @param[in] hRot     rotator handle
 * @param[in] inputs   Input channel buffers; 2-D array: nInputs x nSamples
 * @param[in] outputs  Output channel buffers; 2-D array: nOutputs x nSamples
 * @param[in] nInputs  Number of input channels
 * @param[in] nOutputs Number of output channels
 * @param[in] nSamples Number of samples in 'inputs'/'output' matrices
 */
void rotator_process(void* const hRot,
                     float** const inputs,
                     float** const outputs,
                     int nInputs,
                     int nOutputs,
                     int nSamples);


/* ========================================================================== */
/*                                Set Functions                               */
/* ========================================================================== */

/**
 * Sets the 'yaw' rotation angle, in DEGREES
 */
void rotator_setYaw(void* const hRot, float newYaw);

/**
 * Sets the 'pitch' rotation angle, in DEGREES
 */
void rotator_setPitch(void* const hRot, float newPitch);

/**
 * Sets the 'roll' rotation angle , in DEGREES
 */
void rotator_setRoll(void* const hRot, float newRoll);

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
 * Sets the Ambisonic channel ordering convention to decode with, in order to
 * match the convention employed by the input signals (see 'ROTATOR_CH_ORDER'
 * enum)
 */
void rotator_setChOrder(void* const hRot, int newOrder);

/**
 * Sets the Ambisonic normalisation convention to decode with, in order to match
 * with the convention employed by the input signals (see 'ROTATOR_NORM_TYPE'
 * enum)
 */
void rotator_setNormType(void* const hRot, int newType);

/**
 * Sets the input/output order (see 'ROTATOR_INPUT_ORDERS' enum)
 */
void rotator_setOrder(void* const hRot, int newOrder);

/**
 * Sets a flag as to whether to use "yaw-pitch-roll" (0) or "roll-pitch-yaw" (1)
 * rotation order.
 */
void rotator_setRPYflag(void* const hRot, int newState);


/* ========================================================================== */
/*                                Get Functions                               */
/* ========================================================================== */

/**
 * Returns the 'yaw' rotation angle, in DEGREES
 */
float rotator_getYaw(void* const hRot);

/**
 * Returns the 'pitch' rotation angle, in DEGREES
 */
float rotator_getPitch(void* const hRot);

/**
 * Returns the 'roll' rotation angle, in DEGREES
 */
float rotator_getRoll(void* const hRot);

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
 * Returns a flag as to whether to use "yaw-pitch-roll" (0) or "roll-pitch-yaw"
 * (1) rotation order.
 */
int rotator_getRPYflag(void* const hRot);

/**
 * Returns the Ambisonic channel ordering convention currently being used to
 * decode with, which should match the convention employed by the input signals
 * (see 'ROTATOR_CH_ORDER' enum)
 */
int rotator_getChOrder(void* const hRot);

/**
 * Returns the Ambisonic normalisation convention currently being usedto decode
 * with, which should match the convention employed by the input signals
 * (see 'ROTATOR_NORM_TYPE' enum)
 */
int rotator_getNormType(void* const hRot);

/**
 * Returns the input/output order (see 'ROTATOR_INPUT_ORDERS' enum)
 */
int rotator_getOrder(void* const hRot);

/**
 * Returns the number of spherical harmonic signals required by the current
 * input/output order: (current_order+1)^2 
 */
int rotator_getNSHrequired(void* const hRot);
    
    
#ifdef __cplusplus
} /* extern "C" { */
#endif /* __cplusplus */

#endif /* __ROTATOR_H_INCLUDED__ */
