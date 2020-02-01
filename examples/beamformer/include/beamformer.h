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

/*
 * Filename: beamformer.h (include header)
 * ---------------------------------------
 * Generates beamformers/virtual microphones in arbitrary directions. Several
 * different beam pattern types are included.
 *
 * Dependencies:
 *     saf_utilities, afSTFTlib, saf_sh
 * Author, date created:
 *     Leo McCormack, 17.05.2019
 */

#ifndef __BEAMFORMER_H_INCLUDED__
#define __BEAMFORMER_H_INCLUDED__

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */
    
/* ========================================================================== */
/*                             Presets + Constants                            */
/* ========================================================================== */

/*
 * Enum: BEAM_ORDERS
 * -----------------
 * Available beamforming orders
 *
 * Options:
 *     BEAM_ORDER_FIRST   - First-order beamforming (4 channel input)
 *     BEAM_ORDER_SECOND  - Second-order beamforming (9 channel input)
 *     BEAM_ORDER_THIRD   - Third-order beamforming (16 channel input)
 *     BEAM_ORDER_FOURTH  - Fourth-order beamforming (25 channel input)
 *     BEAM_ORDER_FIFTH   - Fifth-order beamforming (36 channel input)
 *     BEAM_ORDER_SIXTH   - Sixth-order beamforming (49 channel input)
 *     BEAM_ORDER_SEVENTH - Seventh-order beamforming (64 channel input)
 */
#define BEAMFORMER_MAX_SH_ORDER ( 7 )
typedef enum _BEAM_ORDERS{
    BEAM_ORDER_FIRST = 1,
    BEAM_ORDER_SECOND,
    BEAM_ORDER_THIRD,
    BEAM_ORDER_FOURTH,
    BEAM_ORDER_FIFTH,
    BEAM_ORDER_SIXTH,
    BEAM_ORDER_SEVENTH
    
}BEAM_ORDERS;
    
/*
 * Enum: BEAM_TYPES
 * ----------------
 * Available beamforming approaches
 *
 * Options:
 *     BEAM_TYPE_CARDIOID      - cardioid
 *     BEAM_TYPE_HYPERCARDIOID - hyper-cardioid
 *     BEAM_TYPE_MAX_EV        - hyper-cardioid with max_rE weighting
 */
#define BEAMFORMER_NUM_BEAM_TYPES ( 3 )
typedef enum _BEAM_TYPES {
    BEAM_TYPE_CARDIOID = 1,
    BEAM_TYPE_HYPERCARDIOID,
    BEAM_TYPE_MAX_EV
    
} BEAM_TYPES;
 
/*
 * Enum: CH_ORDER
 * --------------
 * Available Ambisonic channel ordering conventions
 * Note: CH_FUMA only supported for 1st order input.
 * Further note: FuMa: CH_FUMA+NORM_FUMA, AmbiX: CH_ACN+NORM_SN3D
 *
 * Options:
 *     CH_ACN  - Ambisonic Channel Numbering (ACN)
 *     CH_FUMA - (Obsolete) Furse-Malham/B-format (WXYZ)
 */
#define BEAMFORMER_NUM_CH_ORDERINGS ( 2 )
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
#define BEAMFORMER_NUM_NORM_TYPES ( 3 )
typedef enum _NORM_TYPES{
    NORM_N3D = 1,
    NORM_SN3D,
    NORM_FUMA   /* first-order only */
}NORM_TYPES;
    
#define BEAMFORMER_MAX_NUM_BEAMS ( 64 )
    
    
/* ========================================================================== */
/*                               Main Functions                               */
/* ========================================================================== */

/*
 * Function: beamformer_create
 * ---------------------------
 * Creates an instance of beamformer
 *
 * Input Arguments:
 *     phBeam - & address of beamformer handle
 */
void beamformer_create(void** const phBeam);

/*
 * Function: beamformer_destroy
 * ----------------------------
 * Destroys an instance of beamformer
 *
 * Input Arguments:
 *     phBeam - & address of beamformer handle
 */
void beamformer_destroy(void** const phBeam);

/*
 * Function: beamformer_init
 * -------------------------
 * Initialises an instance of beamformer with default settings
 *
 * Input Arguments:
 *     hBeam      - beamformer handle
 *     samplerate - host samplerate.
 */
void beamformer_init(void* const hBeam,
                 int samplerate);

/*
 * Function: beamformer_process
 * ----------------------------
 * Generates beamformers/virtual microphones in the specified directions
 *
 * Input Arguments:
 *     hBeam     - beamformer handle
 *     inputs    - input channel buffers; 2-D array: nInputs x nSamples
 *     outputs   - output channel buffers; 2-D array: nOutputs x nSamples
 *     nInputs   - number of input channels
 *     nOutputs  - number of output channels
 *     nSamples  - number of samples in 'inputs'/'output' matrices
 */
void beamformer_process(void* const hBeam,
                        float** const inputs,
                        float** const outputs,
                        int nInputs,
                        int nOutputs,
                        int nSamples);


/* ========================================================================== */
/*                                Set Functions                               */
/* ========================================================================== */

/*
 * Function: beamformer_refreshSettings
 * ------------------------------------
 * Sets all intialisation flags to 1. i.e. re-initialise all settings/variables
 * as beamformer is currently configured, at next available opportunity.
 *
 * Input Arguments:
 *     hBeam - beamformer handle
 */
void beamformer_refreshSettings(void* const hBeam);

/*
 * Function: beamformer_checkReInit
 * --------------------------------
 * Check if any reInit Flags are active, and reinitialise if they are.
 * Note: Only call when playback has stopped.
 *
 * Input Arguments:
 *     hBeam - beamformer handle
 */
void beamformer_checkReInit(void* const hBeam);
    
/*
 * Function: beamformer_setBeamOrder
 * ---------------------------------
 * Sets the beamforming order. If the beamforming order is higher than the
 * input signal order, the extra required channels are filled with zeros. If the
 * beamforming order is lower than the input signal order, the number input
 * signals is truncated accordingly.
 *
 * Input Arguments:
 *     hBeam    - beamformer handle
 *     newValue - new beamforming order (see 'BEAM_ORDERS' enum)
 */
void beamformer_setBeamOrder(void* const hBeam,  int newValue);

/*
 * Function: beamformer_setBeamAzi_deg
 * -----------------------------------
 * Sets a beamformer azimuth direction for a given beamforming index
 *
 * Input Arguments:
 *     hBeam      - beamformer handle
 *     index      - beamformer index
 *     newAzi_deg - new beamforming azimuth, in DEGREES
 */
void beamformer_setBeamAzi_deg(void* const hBeam, int index, float newAzi_deg);

/*
 * Function: beamformer_setBeamElev_deg
 * ------------------------------------
 * Sets a beamformer elevation direction for a given beamforming index
 *
 * Input Arguments:
 *     hBeam       - beamformer handle
 *     index       - beamformer index
 *     newElev_deg - new beamforming elevation, in DEGREES
 */
void beamformer_setBeamElev_deg(void* const hBeam, int index, float newElev_deg);

/*
 * Function: beamformer_setNumBeams
 * --------------------------------
 * Sets the number of beamformers to generate
 *
 * Input Arguments:
 *     hBeam      - beamformer handle
 *     new_nBeams - new number of beamformers
 */
void beamformer_setNumBeams(void* const hBeam, int new_nBeams);

/*
 * Function: beamformer_setChOrder
 * -------------------------------
 * Sets the Ambisonic channel ordering convention to decode with, in order to
 * match the convention employed by the input signals
 *
 * Input Arguments:
 *     hBeam    - beamformer handle
 *     newOrder - convention to use (see 'CH_ORDER' enum)
 */
void beamformer_setChOrder(void* const hBeam, int newOrder);

/*
 * Function: beamformer_setNormType
 * --------------------------------
 * Sets the Ambisonic normalisation convention to decode with, in order to match
 * with the convention employed by the input signals.
 *
 * Input Arguments:
 *     hBeam   - beamformer handle
 *     newType - convention to use (see 'NORM_TYPE' enum)
 */
void beamformer_setNormType(void* const hBeam, int newType);
    
/*
 * Function: beamformer_setBeamType
 * --------------------------------
 * Sets the beamforming approach to employ
 *
 * Input Arguments:
 *     hBeam - beamformer handle
 *     newID - approach to use (see 'BEAM_TYPE' enum)
 */
void beamformer_setBeamType(void* const hBeam, int newID);

    
/* ========================================================================== */
/*                                Get Functions                               */
/* ========================================================================== */

/*
 * Function: beamformer_setBeamOrder
 * ---------------------------------
 * Sets the beamforming order. If the beamforming order is higher than the
 * input signal order, the extra required channels are filled with zeros. If the
 * beamforming order is lower than the input signal order, the number input
 * signals is truncated accordingly.
 *
 * Input Arguments:
 *     hBeam - beamformer handle
 * Returns:
 *     beamforming order (see 'BEAM_ORDERS' enum)
 */
int beamformer_getBeamOrder(void* const hBeam);
    
/*
 * Function: beamformer_getNumberOfBands
 * -------------------------------------
 * Returns the number of frequency bands used by beamformer (only for adaptive
 * beamformer algorithms).
 *
 * Returns:
 *     number of frequency bands
 */
int beamformer_getNumberOfBands(void);

/*
 * Function: beamformer_getBeamAzi_deg
 * -----------------------------------
 * Returns the beamformer azimuth direction for a given beamforming index
 *
 * Input Arguments:
 *     hBeam      - beamformer handle
 *     index      - beamformer index
 * Returns:
 *     new beamforming azimuth, in DEGREES
 */
float beamformer_getBeamAzi_deg(void* const hBeam, int index);

/*
 * Function: beamformer_getBeamElev_deg
 * ------------------------------------
 * Returns the beamformer elevation direction for a given beamforming index
 *
 * Input Arguments:
 *     hBeam       - beamformer handle
 *     index       - beamformer index
 * Returns:
 *     new beamforming elevation, in DEGREES
 */
float beamformer_getBeamElev_deg(void* const hBeam, int index);

/*
 * Function: beamformer_getNumBeams
 * --------------------------------
 * Returns the number of beamformers to generate
 *
 * Input Arguments:
 *     hBeam      - beamformer handle
 * Returns:
 *     new number of beamformers
 */
int beamformer_getNumBeams(void* const hBeam);
    
/*
 * Function: beamformer_getMaxNumBeams
 * -----------------------------------
 * Returns the maximum number of beamformers permitted
 *
 * Returns:
 *     maximum number of beamformers permitted
 */
int beamformer_getMaxNumBeams(void);
    
/*
 * Function: beamformer_getNSHrequired
 * -----------------------------------
 * Returns the number of spherical harmonic signals required by the current
 * beamforming order i.e. (current_order+1)^2
 *
 * Input Arguments:
 *     hBeam - beamformer handle
 * Returns:
 *     number of required spherical harmonic signals required by current
 *     beamforming order
 */
int  beamformer_getNSHrequired(void* const hBeam);

/*
 * Function: beamformer_getChOrder
 * -------------------------------
 * Returns the Ambisonic channel ordering convention currently being used to
 * decode with, which should match the convention employed by the input signals
 *
 * Input Arguments:
 *     hBeam - beamformer handle
 * Returns:
 *     convention currently being used (see 'CH_ORDER' enum)
 */
int beamformer_getChOrder(void* const hBeam);

/*
 * Function: beamformer_getNormType
 * --------------------------------
 * Returns the Ambisonic normalisation convention currently being usedto decode
 * with, which should match the convention employed by the input signals.
 *
 * Input Arguments:
 *     hBeam - beamformer handle
 * Returns:
 *     convention currently being used (see 'NORM_TYPE' enum)
 */
int beamformer_getNormType(void* const hBeam);
    
/*
 * Function: beamformer_getBeamType
 * --------------------------------
 * Returns the beamforming approach to employ
 *
 * Input Arguments:
 *     hBeam - beamformer handle
 * Returns:
 *    approach in use (see 'BEAM_TYPE' enum)
 */
int beamformer_getBeamType(void* const hBeam); 
    
    
#ifdef __cplusplus
} /* extern "C" { */
#endif /* __cplusplus */

#endif /* __BEAMFORMER_H_INCLUDED__ */
