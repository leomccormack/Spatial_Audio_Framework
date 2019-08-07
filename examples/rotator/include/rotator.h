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

/*
 * Filename: rotator.h (include header)
 * ------------------------------------
 * A simple spherical harmonic domain rotator, based on the recursive approach
 * detailed in [1].
 *
 * Dependencies:
 *     saf_utilities, saf_sh
 * Author, date created:
 *     Leo McCormack, 02.11.2017
 *
 * [1] Ivanic, J., Ruedenberg, K. (1998). Rotation Matrices for Real Spherical
 *     Harmonics. Direct Determination by Recursion Page: Additions and
 *     Corrections. Journal of Physical Chemistry A, 102(45), 9099?9100.
 */

#ifndef __ROTATOR_H_INCLUDED__
#define __ROTATOR_H_INCLUDED__

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */
    
/* ========================================================================== */
/*                             Presets + Constants                            */
/* ========================================================================== */

/*
 * Enum: INPUT_ORDERS
 * ------------------
 * Available input/output orders.
 *
 * Options:
 *     INPUT_ORDER_FIRST   - First-order rotation (4 channel input/output)
 *     INPUT_ORDER_SECOND  - Second-order rotation (9 channel input/output)
 *     INPUT_ORDER_THIRD   - Third-order rotation (16 channel input/output)
 *     INPUT_ORDER_FOURTH  - Fourth-order rotation (25 channel input/output)
 *     INPUT_ORDER_FIFTH   - Fifth-order rotation (36 channel input/output)
 *     INPUT_ORDER_SIXTH   - Sixth-order rotation (49 channel input/output)
 *     INPUT_ORDER_SEVENTH - Seventh-order rotation (64 channel input/output)
 */
#define ROTATOR_MAX_SH_ORDER ( 7 )
typedef enum _INPUT_ORDERS{
    INPUT_ORDER_FIRST = 1,
    INPUT_ORDER_SECOND,
    INPUT_ORDER_THIRD,
    INPUT_ORDER_FOURTH,
    INPUT_ORDER_FIFTH,
    INPUT_ORDER_SIXTH,
    INPUT_ORDER_SEVENTH
    
}INPUT_ORDERS;

/*
 * Enum: _CH_ORDER
 * ---------------
 * Available Ambisonic channel ordering conventions
 * Note: CH_FUMA only supported for 1st order input.
 * Further note: FuMa: CH_FUMA+NORM_FUMA, AmbiX: CH_ACN+NORM_SN3D
 *
 * Options:
 *     CH_ACN  - Ambisonic Channel Numbering (ACN)
 *     CH_FUMA - (Obsolete) Furse-Malham/B-format (WXYZ)
 */
#define ROTATOR_NUM_CH_ORDERINGS ( 2 )
typedef enum _CH_ORDER{
    CH_ACN = 1,
    CH_FUMA     /* first-order only */
    
}CH_ORDER;

/*
 * Enum: NORM_TYPES
 * ---------------
 * Available Ambisonic normalisation conventions
 * Note: NORM_FUMA only supported for 1st order input and does NOT have the
 * 1/sqrt(2) scaling on the omni.
 * Further note: FuMa: CH_FUMA+NORM_FUMA, AmbiX: CH_ACN+NORM_SN3D
 *
 * Options:
 *     NORM_N3D  - orthonormalised (N3D)
 *     NORM_SN3D - Schmidt semi-normalisation (SN3D)
 *     NORM_FUMA - (Obsolete) Same as NORM_SN3D for 1st order
 */
#define ROTATOR_NUM_NORM_TYPES ( 3 )
typedef enum _NORM_TYPES{
    NORM_N3D = 1,
    NORM_SN3D,
    NORM_FUMA   /* first-order only */
}NORM_TYPES;
    
#define ROTATOR_MAX_NUM_CHANNELS ( 64 )
    

/* ========================================================================== */
/*                               Main Functions                               */
/* ========================================================================== */

/*
 * Function: rotator_create
 * ------------------------
 * Creates an instance of rotator
 *
 * Input Arguments:
 *     phRot - & address of rotator handle
 */
void rotator_create(void** const phRot);

/*
 * Function: rotator_destroy
 * -------------------------
 * Destroys an instance of rotator
 *
 * Input Arguments:
 *     phRot - & address of rotator handle
 */
void rotator_destroy(void** const phRot);

/*
 * Function: rotator_init
 * ----------------------
 * Initialises an instance of rotator with default settings
 *
 * Input Arguments:
 *     hRot      - rotator handle
 *     samplerate - host samplerate.
 */
void rotator_init(void* const hRot,
                  int samplerate);
    
/*
 * Function: rotator_process
 * -------------------------
 * Rotates the input spherical harmonic signals.
 *
 * Input Arguments:
 *     hRot     - rotator handle
 *     inputs    - input channel buffers; 2-D array: nInputs x nSamples
 *     outputs   - Output channel buffers; 2-D array: nOutputs x nSamples
 *     nInputs   - number of input channels
 *     nOutputs  - number of output channels
 *     nSamples  - number of samples in 'inputs'/'output' matrices
 *     isPlaying - flag to say if there is audio in the input buffers, 0: no
 *                 audio, reduced processing, 1: audio, full processing
 */
void rotator_process(void* const hRot,
                     float** const inputs,
                     float** const outputs,
                     int nInputs,
                     int nOutputs,
                     int nSamples,
                     int isPlaying);

    
/* ========================================================================== */
/*                                Set Functions                               */
/* ========================================================================== */

/*
 * Function: rotator_setYaw
 * ------------------------
 * Sets the 'yaw' rotation angle
 *
 * Input Arguments:
 *     hRot   - rotator handle
 *     newYaw - the 'yaw' rotation angle, in DEGREES
 */
void rotator_setYaw(void* const hRot, float newYaw);

/*
 * Function: rotator_setPitch
 * --------------------------
 * Sets the 'pitch' rotation angle
 *
 * Input Arguments:
 *     hRot     - rotator handle
 *     newPitch - the 'pitch' rotation angle, in DEGREES
 */
void rotator_setPitch(void* const hRot, float newPitch);

/*
 * Function: rotator_setRoll
 * --------------------------
 * Sets the 'roll' rotation angle
 *
 * Input Arguments:
 *     hRot    - rotator handle
 *     newRoll - the 'roll' rotation angle, in DEGREES
 */
void rotator_setRoll(void* const hRot, float newRoll);

/*
 * Function: rotator_setFlipYaw
 * ----------------------------
 * Sets a flag as to whether to "flip" the sign of the current 'yaw' angle
 *
 * Input Arguments:
 *     hRot     - rotator handle
 *     newState - 0: do not flip sign, 1: flip the sign
 */
void rotator_setFlipYaw(void* const hRot, int newState);

/*
 * Function: rotator_setFlipPitch
 * ------------------------------
 * Sets a flag as to whether to "flip" the sign of the current 'pitch' angle
 *
 * Input Arguments:
 *     hRot     - rotator handle
 *     newState - 0: do not flip sign, 1: flip the sign
 */
void rotator_setFlipPitch(void* const hRot, int newState);

/*
 * Function: rotator_setFlipRoll
 * -----------------------------
 * Sets a flag as to whether to "flip" the sign of the current 'roll' angle
 *
 * Input Arguments:
 *     hRot     - rotator handle
 *     newState - 0: do not flip sign, 1: flip the sign
 */
void rotator_setFlipRoll(void* const hRot, int newState);

/*
 * Function: rotator_setChOrder
 * ----------------------------
 * Sets the Ambisonic channel ordering convention to decode with, in order to
 * match the convention employed by the input signals
 *
 * Input Arguments:
 *     hRot    - rotator handle
 *     newOrder - convention to use (see 'CH_ORDER' enum)
 */
void rotator_setChOrder(void* const hRot, int newOrder);

/*
 * Function: rotator_setNormType
 * -----------------------------
 * Sets the Ambisonic normalisation convention to decode with, in order to match
 * with the convention employed by the input signals.
 *
 * Input Arguments:
 *     hRot   - rotator handle
 *     newType - convention to use (see 'NORM_TYPE' enum)
 */
void rotator_setNormType(void* const hRot, int newType);

/*
 * Function: rotator_setInputOrderPreset
 * -------------------------------------
 * Sets the input/output order.
 *
 * Input Arguments:
 *     hRot     - rotator handle
 *     newValue - new input/output order (see 'INPUT_ORDERS' enum)
 */
void rotator_setOrder(void* const hRot, int newOrder);

/*
 * Function: rotator_setRPYflag
 * ----------------------------
 * Sets a flag as to whether to use "yaw-pitch-roll" (0) or "roll-pitch-yaw" (1)
 * rotation order.
 *
 * Input Arguments:
 *     hRot     - rotator handle
 *     newState - 0: use "yaw-pitch-roll", 1: use "roll-pitch-yaw"
 */
void rotator_setRPYflag(void* const hRot, int newState);


/* ========================================================================== */
/*                                Get Functions                               */
/* ========================================================================== */

/*
 * Function: rotator_getYaw
 * ------------------------
 * Returns the 'yaw' rotation angle
 *
 * Input Arguments:
 *     hRot - rotator handle
 * Returns:
 *     the 'yaw' rotation angle, in DEGREES
 */
float rotator_getYaw(void* const hRot);

/*
 * Function: rotator_getPitch
 * --------------------------
 * Returns the 'pitch' rotation angle
 *
 * Input Arguments:
 *     hRot - rotator handle
 * Returns:
 *     the 'pitch' rotation angle, in DEGREES
 */
float rotator_getPitch(void* const hRot);

/*
 * Function: rotator_getRoll
 * -------------------------
 * Returns the 'roll' rotation angle
 *
 * Input Arguments:
 *     hRot - rotator handle
 * Returns:
 *     the 'roll' rotation angle, in DEGREES
 */
float rotator_getRoll(void* const hRot);

/*
 * Function: rotator_getFlipYaw
 * ----------------------------
 * Returns a flag as to whether to "flip" the sign of the current 'yaw' angle
 *
 * Input Arguments:
 *     hRot - rotator handle
 * Returns:
 *     0: do not flip sign, 1: flip the sign
 */
int rotator_getFlipYaw(void* const hRot);

/*
 * Function: rotator_getFlipPitch
 * ------------------------------
 * Returns a flag as to whether to "flip" the sign of the current 'pitch' angle
 *
 * Input Arguments:
 *     hRot - rotator handle
 * Returns:
 *     0: do not flip sign, 1: flip the sign
 */
int rotator_getFlipPitch(void* const hRot);

/*
 * Function: rotator_getFlipRoll
 * -----------------------------
 * Returns a flag as to whether to "flip" the sign of the current 'roll' angle
 *
 * Input Arguments:
 *     hRot - rotator handle
 * Returns:
 *     0: do not flip sign, 1: flip the sign
 */
int rotator_getFlipRoll(void* const hRot);

/*
 * Function: rotator_getRPYflag
 * ----------------------------
 * Returns a flag as to whether to use "yaw-pitch-roll" (0) or "roll-pitch-yaw"
 * (1) rotation order.
 *
 * Input Arguments:
 *     hRot - rotator handle
 * Returns:
 *     0: use "yaw-pitch-roll", 1: use "roll-pitch-yaw"
 */
int rotator_getRPYflag(void* const hRot);

/*
 * Function: rotator_getChOrder
 * ----------------------------
 * Returns the Ambisonic channel ordering convention currently being used to
 * decode with, which should match the convention employed by the input signals
 *
 * Input Arguments:
 *     hRot - rotator handle
 * Returns:
 *     convention currently being used (see 'CH_ORDER' enum)
 */
int rotator_getChOrder(void* const hRot);

/*
 * Function: rotator_getNormType
 * -----------------------------
 * Returns the Ambisonic normalisation convention currently being usedto decode
 * with, which should match the convention employed by the input signals.
 *
 * Input Arguments:
 *     hRot - rotator handle
 * Returns:
 *     convention currently being used (see 'NORM_TYPE' enum)
 */
int rotator_getNormType(void* const hRot);

/*
 * Function: rotator_getInputOrderPreset
 * -------------------------------------
 * Returns the input/output order.
 *
 * Input Arguments:
 *     hRot - rotator handle
 * Returns:
 *     current input/output order (see 'INPUT_ORDERS' enum)
 */
int rotator_getOrder(void* const hRot);

/*
 * Function: rotator_getNSHrequired
 * --------------------------------
 * Returns the number of spherical harmonic signals required by the current
 * input/output order i.e. (current_order+1)^2
 *
 * Input Arguments:
 *     hRot - rotator handle
 * Returns:
 *     number of required spherical harmonic signals required by current
 *     input/output order
 */
int rotator_getNSHrequired(void* const hRot);
    
    
#ifdef __cplusplus
} /* extern "C" { */
#endif /* __cplusplus */

#endif /* __ROTATOR_H_INCLUDED__ */
