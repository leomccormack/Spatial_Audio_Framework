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
 * Filename: array2sh.h (include header)
 * -------------------------------------
 * Spatially encodes spherical or cylindrical sensor array signals into
 * spherical harmonic signals utilising theoretical encoding filters.
 * The algorithms within array2sh were pieced together and developed in
 * collaboration with Symeon Delikaris-Manias and Angelo Farina.
 * A detailed explanation of the algorithms within array2sh can be found in [1].
 * Also included, is a diffuse-field equalisation option for frequencies past
 * aliasing, developed in collaboration with Archontis Politis, 8.02.2019
 * Note: since the algorithms are based on theory, only array designs where
 * there are analytical solutions available are supported. i.e. only spherical
 * or cylindrical arrays, which have phase-matched sensors.
 *
 * Dependencies:
 *     saf_utilities, afSTFTlib, saf_sh, saf_hoa, saf_vbap
 * Author, date created:
 *     Leo McCormack, 13.09.2017
 *
 * [1] McCormack, L., Delikaris-Manias, S., Farina, A., Pinardi, D., and Pulkki,
 *     V., “Real-time conversion of sensor array signals into spherical harmonic
 *     signals with applications to spatially localised sub-band sound-field
 *     analysis,” in Audio Engineering Society Convention 144, Audio Engineering
 *     Society, 2018.
 */

#ifndef __ARRAY2SH_H_INCLUDED__
#define __ARRAY2SH_H_INCLUDED__

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */
    
/* ========================================================================== */
/*                             Presets + Constants                            */
/* ========================================================================== */
    
/*
 * Enum: ENCODING_ORDERS
 * ---------------------
 * Available encoding orders
 *
 * Options:
 *     ENCODING_ORDER_FIRST   - First-order encoding (4 channel output)
 *     ENCODING_ORDER_SECOND  - Second-order encoding (9 channel output)
 *     ENCODING_ORDER_THIRD   - Third-order encoding (16 channel output)
 *     ENCODING_ORDER_FOURTH  - Fourth-order encoding (25 channel output)
 *     ENCODING_ORDER_FIFTH   - Fifth-order encoding (36 channel output)
 *     ENCODING_ORDER_SIXTH   - Sixth-order encoding (49 channel output)
 *     ENCODING_ORDER_SEVENTH - Seventh-order encoding (64 channel output)
 */
#define ARRAY2SH_MAX_SH_ORDER ( 7 )
typedef enum _ENCODING_ORDERS{
    ENCODING_ORDER_FIRST = 1,
    ENCODING_ORDER_SECOND,
    ENCODING_ORDER_THIRD,
    ENCODING_ORDER_FOURTH,
    ENCODING_ORDER_FIFTH,
    ENCODING_ORDER_SIXTH,
    ENCODING_ORDER_SEVENTH
    
}ENCODING_ORDERS;
    
/*
 * Enum: MICROPHONE_ARRAY_PRESETS
 * -------------------------------
 * Available microphone array presets
 */
typedef enum _MICROPHONE_ARRAY_PRESETS{
    MICROPHONE_ARRAY_PRESET_DEFAULT = 1,
    MICROPHONE_ARRAY_PRESET_AALTO_HYDROPHONE,
    MICROPHONE_ARRAY_PRESET_SENNHEISER_AMBEO,
    MICROPHONE_ARRAY_PRESET_CORE_SOUND_TETRAMIC,
    MICROPHONE_ARRAY_PRESET_ZOOM_H3VR_PRESET,
    MICROPHONE_ARRAY_PRESET_SOUND_FIELD_SPS200,
    MICROPHONE_ARRAY_PRESET_ZYLIA_1D,
    MICROPHONE_ARRAY_PRESET_EIGENMIKE32,
    MICROPHONE_ARRAY_PRESET_DTU_MIC

}MICROPHONE_ARRAY_PRESETS;

/*
 * Enum: FILTER_TYPES
 * ----------------------
 * Available encoding filter approaches.
 *
 * Options:
 *     FILTER_SOFT_LIM      - encoding filters based on a 'soft-limiting'
 *                            regularised inversion of the modal responses [1]
 *     FILTER_TIKHONOV      - encoding filters based on a 'Tikhonov' regularised
 *                            inversion of the modal responses [2]
 *     FILTER_Z_STYLE       - encoding filters based on a linear-phase filter-
 *                            bank approach [3]
 *     FILTER_Z_STYLE_MAXRE - same as 'FILTER_Z_STYLE', only it also has max_rE
 *                            weights baked in
 *
 * [1] Bernschutz, B., Porschmann, C., Spors, S., Weinzierl, S., Versterkung,
 *     B., 2011. Soft-limiting der modalen amplitudenverst?rkung bei sph?rischen
 *     mikrofonarrays im plane wave decomposition verfahren. Proceedings of the
 *     37. Deutsche Jahrestagung fur Akustik (DAGA 2011)
 * [2] Moreau, S., Daniel, J., Bertet, S., 2006, 3D sound field recording with
 *     higher order ambisonics-objective measurements and validation of
 *     spherical microphone. In Audio Engineering Society Convention 120.
 * [3] Zotter, F. A Linear-Phase Filter-Bank Approach to Process Rigid Spherical
 *     Microphone Array Recordings.
 */
#define ARRAY2SH_NUM_FILTER_TYPES ( 4 )
typedef enum _FILTER_TYPES{
    FILTER_SOFT_LIM = 1,
    FILTER_TIKHONOV,
    FILTER_Z_STYLE,
    FILTER_Z_STYLE_MAXRE
    
}FILTER_TYPES;

/*
 * Enum: CH_ORDER
 * --------------
 * Available Ambisonic channel ordering conventions
 * Note: CH_FUMA only supported for 1st order output.
 * Further note: FuMa: CH_FUMA+NORM_FUMA, AmbiX: CH_ACN+NORM_SN3D
 *
 * Options:
 *     CH_ACN  - Ambisonic Channel Numbering (ACN)
 *     CH_FUMA - (Obsolete) Furse-Malham/B-format (WXYZ)
 */
#define ARRAY2SH_NUM_CH_ORDERINGS ( 2 )
typedef enum _CH_ORDER{
    CH_ACN = 1,
    CH_FUMA     /* first-order only */
    
}CH_ORDER;

/*
 * Enum: NORM_TYPES
 * ---------------
 * Available Ambisonic normalisation conventions
 * Note: NORM_FUMA only supported for 1st order output and does NOT have the
 * 1/sqrt(2) scaling on the omni.
 * Further note: FuMa: CH_FUMA+NORM_FUMA, AmbiX: CH_ACN+NORM_SN3D
 *
 * Options:
 *     NORM_N3D  - orthonormalised (N3D)
 *     NORM_SN3D - Schmidt semi-normalisation (SN3D)
 *     NORM_FUMA - (Obsolete) Same as NORM_SN3D for 1st order
 */
#define ARRAY2SH_NUM_NORM_TYPES ( 3 )
typedef enum _NORM_TYPES{
    NORM_N3D = 1,
    NORM_SN3D,
    NORM_FUMA   /* first-order only */
    
}NORM_TYPES;

/*
 * Enum: ARRAY_TYPES
 * -----------------
 * List of supported array types.
 * Note: although supported, cylindrical arrays have not really been tested as
 * we don't own one
 *
 * Options:
 *     ARRAY_SPHERICAL   - spherical arrangement of sensors (open/rigid)
 *     ARRAY_CYLINDRICAL - cylindrial arrangement of sensors (open/rigid)
 */
#define ARRAY2SH_NUM_ARRAY_TYPES ( 2 )
typedef enum _ARRAY_TYPES{
    ARRAY_SPHERICAL = 1,
    ARRAY_CYLINDRICAL
    
}ARRAY_TYPES;

/*
 * Enum: WEIGHT_TYPES
 * ------------------
 * List of supported sensor directivities and array construction types.
 *
 * Options:
 *     WEIGHT_RIGID_OMNI   - rigid baffle construction with omni sensors
 *     WEIGHT_RIGID_CARD   - rigid baffle construction with cardioid sensors
 *     WEIGHT_RIGID_DIPOLE - rigid baffle construction with dipole sensors
 *     WEIGHT_OPEN_OMNI    - open array construction with omni sensors
 *     WEIGHT_OPEN_CARD    - open array construction with cardioid sensors
 *     WEIGHT_OPEN_DIPOLE  - open array construction with dipole sensors
 */
#define ARRAY2SH_NUM_WEIGHT_TYPES ( 6 )
typedef enum _WEIGHT_TYPES{
    WEIGHT_RIGID_OMNI = 1, 
    WEIGHT_RIGID_CARD,
    WEIGHT_RIGID_DIPOLE,
    WEIGHT_OPEN_OMNI,
    WEIGHT_OPEN_CARD,
    WEIGHT_OPEN_DIPOLE
    
}WEIGHT_TYPES;
    
/*
 * Enum: EVAL_STATUS
 * -----------------
 * Current status of the encoder.
 *
 * Options:
 *     EVAL_STATUS_EVALUATED          - Encoder has been evaluated
 *     EVAL_STATUS_RECENTLY_EVALUATED - Encoder has recently been evaluated
 *     EVAL_STATUS_NOT_EVALUATED      - Encoder has not been evaluated
 *     EVAL_STATUS_EVALUATING         - Encoder is being evaluated
 */
typedef enum _EVAL_STATUS{
    EVAL_STATUS_EVALUATED = 0,
    EVAL_STATUS_RECENTLY_EVALUATED,
    EVAL_STATUS_NOT_EVALUATED,
    EVAL_STATUS_EVALUATING
}EVAL_STATUS;
    
#define ARRAY2SH_MAX_NUM_SENSORS ( 64 )
#define ARRAY2SH_MAX_GAIN_MIN_VALUE ( 0.0f )
#define ARRAY2SH_MAX_GAIN_MAX_VALUE ( 80.0f )
#define ARRAY2SH_POST_GAIN_MIN_VALUE ( -60.0f )
#define ARRAY2SH_POST_GAIN_MAX_VALUE ( 12.0f )
#define ARRAY2SH_SPEED_OF_SOUND_MIN_VALUE ( 200.0f )
#define ARRAY2SH_SPEED_OF_SOUND_MAX_VALUE ( 2000.0f )
#define ARRAY2SH_ARRAY_RADIUS_MIN_VALUE ( 1.0f )
#define ARRAY2SH_ARRAY_RADIUS_MAX_VALUE ( 200.0f )
#define ARRAY2SH_BAFFLE_RADIUS_MIN_VALUE ( 1.0f )
#define ARRAY2SH_BAFFLE_RADIUS_MAX_VALUE ( 200.0f )
#define ARRAY2SH_PROGRESSBARTEXT_CHAR_LENGTH 256
    
    
/* ========================================================================== */
/*                               Main Functions                               */
/* ========================================================================== */

/*
 * Function: array2sh_create
 * -------------------------
 * Creates an instance of array2sh
 *
 * Input Arguments:
 *     phA2sh - & address of array2sh handle
 */
void array2sh_create(void** const phA2sh);

/*
 * Function: array2sh_destroy
 * --------------------------
 * Destroys an instance of array2sh
 *
 * Input Arguments:
 *     phA2sh - & address of array2sh handle
 */
void array2sh_destroy(void** const phA2sh);

/*
 * Function: array2sh_init
 * -----------------------
 * Initialises an instance of array2sh with default settings
 *
 * Input Arguments:
 *     hA2sh      - array2sh handle
 *     samplerate - host samplerate.
 */
void array2sh_init(void* const hA2sh,
                   int samplerate);
    
/*
 * Function: array2sh_evalEncoder
 * ------------------------------
 * Evaluates the encoder, based on current global/user parameters
 *
 * Input Arguments:
 *     hA2sh - array2sh handle
 */
void array2sh_evalEncoder(void* const hA2sh);

/*
 * Function: array2sh_process
 * --------------------------
 * Spatially encode microphone/hydrophone array signals into spherical harmonic
 * signals
 *
 * Input Arguments:
 *     hA2sh     - array2sh handle
 *     inputs    - input channel buffers; 2-D array: nInputs x nSamples
 *     outputs   - Output channel buffers; 2-D array: nOutputs x nSamples
 *     nInputs   - number of input channels
 *     nOutputs  - number of output channels
 *     nSamples  - number of samples in 'inputs'/'output' matrices
 *     isPlaying - flag to say if there is audio in the input buffers, 0: no
 *                 audio, reduced processing, 1: audio, full processing
 */
void array2sh_process(void* const hA2sh,
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
 * Function: array2sh_refreshSettings
 * ----------------------------------
 * Sets all intialisation flags to 1. i.e. re-initialise all settings/variables
 * as array2sh is currently configured, at next available opportunity.
 *
 * Input Arguments:
 *     hA2sh - array2sh handle
 */
void array2sh_refreshSettings(void* const hA2sh);
 
/*
 * Function: array2sh_setEncodingOrder
 * -----------------------------------
 * Sets the encoding order.
 *
 * Input Arguments:
 *     hA2sh    - array2sh handle
 *     newOrder - new encoding order (see 'ENCODING_ORDERS' enum)
 */
void array2sh_setEncodingOrder(void* const hA2sh, int newOrder);
    
/*
 * Function: array2sh_setRequestEncoderEvalFLAG
 * --------------------------------------------
 * Evaluates the performance of the current encoding filters when applied to a
 * theoretical model of the currently configured array. Two established
 * objective metrics are computed. More information in [1]
 *
 * Input Arguments:
 *     hA2sh - array2sh handle
 *
 * [1] Moreau, S., Daniel, J., Bertet, S., 2006, 3D sound field recording with
 *     higher order ambisonics-objective measurements and validation of
 *     spherical microphone. In Audio Engineering Society Convention 120.
 */
void array2sh_setRequestEncoderEvalFLAG(void* const hA2sh, int newState);
    
/*
 * Function: array2sh_setEvalStatus
 * --------------------------------
 * Sets current eval status.
 *
 * Input Arguments:
 *     hA2sh      - array2sh handle
 *     evalStatus - see 'EVAL_STATUS' enum
 */
void array2sh_setEvalStatus(void* const hA2sh, EVAL_STATUS evalStatus);

/*
 * Function: array2sh_setDiffEQpastAliasing
 * ----------------------------------------
 * Analyses what the theoretical spatial aliasing frequency is, and conducts
 * diffuse-field equalisation above this.
 * Thanks to Dr. Archontis Politis for suggesting and designing this feature.
 *
 * Input Arguments:
 *     hA2sh - array2sh handle
 */
void array2sh_setDiffEQpastAliasing(void* const hA2sh, int newState);
    
/*
 * Function: array2sh_setPreset
 * ----------------------------
 * Sets a pre-defined microphone/hydrophone array preset. See PRESETS enum
 *
 * Input Arguments:
 *     hA2sh  - array2sh handle
 *     preset - see PRESETS enum
 */
void array2sh_setPreset(void* const hA2sh, int preset);
    
/*
 * Function: array2sh_setSensorAzi_rad
 * -----------------------------------
 * Sets a particular sensor's azimuth w.r.t to the origin of the array.
 *
 * Input Arguments:
 *     hA2sh      - array2sh handle
 *     index      - sensor index
 *     newAzi_rad - sensor azimuth in RADIANS
 */
void array2sh_setSensorAzi_rad(void* const hA2sh, int index, float newAzi_rad);
    
/*
 * Function: array2sh_setSensorElev_rad
 * ------------------------------------
 * Sets a particular sensor's elevation w.r.t to the origin of the array.
 *
 * Input Arguments:
 *     hA2sh       - array2sh handle
 *     index       - sensor index
 *     newElev_rad - sensor elevation in RADIANS
 */
void array2sh_setSensorElev_rad(void* const hA2sh, int index, float newElev_rad);
    
/*
 * Function: array2sh_setSensorAzi_deg
 * -----------------------------------
 * Sets a particular sensor's azimuth w.r.t to the origin of the array.
 *
 * Input Arguments:
 *     hA2sh      - array2sh handle
 *     index      - sensor index
 *     newAzi_deg - sensor azimuth in DEGREES
 */
void array2sh_setSensorAzi_deg(void* const hA2sh, int index, float newAzi_deg);
    
/*
 * Function: array2sh_setSensorElev_deg
 * ------------------------------------
 * Sets a particular sensor's elevation w.r.t to the origin of the array.
 *
 * Input Arguments:
 *     hA2sh       - array2sh handle
 *     index       - sensor index
 *     newElev_deg - sensor elevation in DEGREES
 */
void array2sh_setSensorElev_deg(void* const hA2sh, int index, float newElev_deg);
    
/*
 * Function: array2sh_setNumSensors
 * --------------------------------
 * Sets the number of sensors in the array.
 *
 * Input Arguments:
 *     hA2sh - array2sh handle
 *     newQ  - new number of sensors
 */
void array2sh_setNumSensors(void* const hA2sh, int newQ);
    
/*
 * Function: array2sh_setr
 * -----------------------
 * Sets the radius of the array
 *
 * Input Arguments:
 *     hA2sh - array2sh handle
 *     newr  - new array radius
 */
void array2sh_setr(void* const hA2sh, float newr);
    
/*
 * Function: array2sh_setR
 * -----------------------
 * Sets the radius of the scatterer. Only for Rigid arrays.
 * Note: R <= r. i.e. the sensors may protrude of the rigid scattering surface,
 * or be flush with the surface of the array
 *
 * Input Arguments:
 *     hA2sh - array2sh handle
 *     newr  - new scatterer radius
 */
void array2sh_setR(void* const hA2sh, float newR);
    
/*
 * Function: array2sh_setArrayType
 * -------------------------------
 * Sets the type of array. See ARRAY_TYPES enum
 *
 * Input Arguments:
 *     hA2sh   - array2sh handle
 *     newType - new array type. See ARRAY_TYPES enum
 */
void array2sh_setArrayType(void* const hA2sh, int newType);

/*
 * Function: array2sh_setWeightType
 * --------------------------------
 * Sets the type of weights to use. See WEIGHT_TYPES enum
 *
 * Input Arguments:
 *     hA2sh   - array2sh handle
 *     newType - new weight type. See WEIGHT_TYPES enum
 */
void array2sh_setWeightType(void* const hA2sh, int newType);
    
/*
 * Function: array2sh_setFilterType
 * --------------------------------
 * Sets the type filter design to employ for computing the encoding matrices.
 * See FILTER_TYPES enum
 *
 * Input Arguments:
 *     hA2sh   - array2sh handle
 *     newType - new filter type. See FILTER_TYPES enum
 */
void array2sh_setFilterType(void* const hA2sh, int newType);
    
/*
 * Function: array2sh_setRegPar
 * ----------------------------
 * Sets the value of the regurlisation parameter. i.e. the maximum permitted
 * gain provided by the filters
 *
 * Input Arguments:
 *     hA2sh  - array2sh handle
 *     newVal - new filter maximum gain, in DECIBELS
 */
void array2sh_setRegPar(void* const hA2sh, float newVal);
    
    /*
 * Function: array2sh_setChOrder
 * -----------------------------
 * Sets the Ambisonic channel ordering convention to encode with, in order to
 * match the convention employed by the input signals
 *
 * Input Arguments:
 *     hA2sh    - array2sh handle
 *     newOrder - convention to use (see 'CH_ORDER' enum)
 */
void array2sh_setChOrder(void* const hA2sh, int newOrder);
    
/*
 * Function: array2sh_setNormType
 * ------------------------------
 * Sets the Ambisonic normalisation convention to encode with, in order to match
 * with the convention employed by the input signals.
 *
 * Input Arguments:
 *     hA2sh   - array2sh handle
 *     newType - convention to use (see 'NORM_TYPE' enum)
 */
void array2sh_setNormType(void* const hA2sh, int newType);

/*
 * Function: array2sh_setc
 * -----------------------
 * Sets the speed of sound of the medium (~343m/s air, ~1480m/s water).
 *
 * Input Arguments:
 *     hA2sh - array2sh handle
 *     newc  - new speed of sound, in m/s
 */
void array2sh_setc(void* const hA2sh, float newc);
    
/*
 * Function: array2sh_setGain
 * --------------------------
 * Sets the amount of post gain to apply after the encoding
 *
 * Input Arguments:
 *     hA2sh   - array2sh handle
 *     newGain - new post gain, in DECIBELS
 */
void array2sh_setGain(void* const hA2sh, float newGain);


/* ========================================================================== */
/*                                Get Functions                               */
/* ========================================================================== */

/*
 * Function: array2sh_getEvalStatus
 * --------------------------------
 * Returns current eval status.
 *
 * Input Arguments:
 *     hA2sh - array2sh handle
 * Returns:
 *     codec status (see 'EVAL_STATUS' enum)
 */
EVAL_STATUS array2sh_getEvalStatus(void* const hA2sh);

/*
 * Function: array2sh_getProgressBar0_1
 * ------------------------------------
 * (Optional) Returns current intialisation/processing progress, between 0..1
 * 0: intialisation/processing has started
 * 1: intialisation/processing has ended
 *
 * Input Arguments:
 *     hA2sh - array2sh handle
 * Returns:
 *     current progress, 0..1
 */
float array2sh_getProgressBar0_1(void* const hA2sh);

/*
 * Function: array2sh_getProgressBarText
 * -------------------------------------
 * (Optional) Returns current intialisation/processing progress text
 * Note: "text" string should be (at least) of length:
 *     ARRAY2SH_PROGRESSBARTEXT_CHAR_LENGTH
 *
 * Input Arguments:
 *     hA2sh - array2sh handle
 * Output Arguments:
 *     text  - process bar text; ARRAY2SH_PROGRESSBARTEXT_CHAR_LENGTH x 1
 */
void array2sh_getProgressBarText(void* const hA2sh, char* text);

/*
 * Function: array2sh_getDiffEQpastAliasing
 * ----------------------------------------
 *
 * Input Arguments:
 *     hA2sh - array2sh handle
 * Returns:
 *     current flag state, 0: disabled, 1: enabled
 */
int array2sh_getDiffEQpastAliasing(void* const hA2sh);

/*
 * Function: array2sh_getRequestEncoderEvalFLAG
 * --------------------------------------------
 * Returns
 *
 * Input Arguments:
 *     hA2sh - array2sh handle
 * Returns:
 *     0:
 */
int array2sh_getRequestEncoderEvalFLAG(void* const hA2sh);
    
/*
 * Function: array2sh_getEncodingOrder
 * -----------------------------------
 * Returns the encoding order.
 *
 * Input Arguments:
 *     hA2sh - array2sh handle
 * Returns:
 *     current encoding order (see 'ENCODING_ORDERS' enum)
 */
int array2sh_getEncodingOrder(void* const hA2sh);
    
/*
 * Function: array2sh_getSensorAzi_rad
 * -----------------------------------
 * Returns a particular sensor's azimuth w.r.t to the origin of the array.
 *
 * Input Arguments:
 *     hA2sh - array2sh handle
 *     index - sensor index
 * Returns:
 *     sensor azimuth in RADIANS
 */
float array2sh_getSensorAzi_rad(void* const hA2sh, int index);
    
/*
 * Function: array2sh_getSensorElev_rad
 * ------------------------------------
 * Returns a particular sensor's elevation w.r.t to the origin of the array.
 *
 * Input Arguments:
 *     hA2sh - array2sh handle
 *     index - sensor index
 * Returns:
 *     sensor elevation in RADIANS
 */
float array2sh_getSensorElev_rad(void* const hA2sh, int index);
    
/*
 * Function: array2sh_getSensorAzi_deg
 * -----------------------------------
 * Returns a particular sensor's azimuth w.r.t to the origin of the array.
 *
 * Input Arguments:
 *     hA2sh - array2sh handle
 *     index - sensor index
 * Returns:
 *     sensor azimuth in DEGREES
 */
float array2sh_getSensorAzi_deg(void* const hA2sh, int index);
    
/*
 * Function: array2sh_getSensorElev_deg
 * ------------------------------------
 * Returns a particular sensor's elevation w.r.t to the origin of the array.
 *
 * Input Arguments:
 *     hA2sh - array2sh handle
 *     index - sensor index
 * Returns:
 *     sensor elevation in DEGREES
 */
float array2sh_getSensorElev_deg(void* const hA2sh, int index);

/*
 * Function: array2sh_getNumSensors
 * --------------------------------
 * Returns the number of sensors in the array.
 *
 * Input Arguments:
 *     hA2sh - array2sh handle
 * Returns:
 *     current number of sensors
 */
int array2sh_getNumSensors(void* const hA2sh);
    
/*
 * Function: array2sh_getMaxNumSensors
 * -----------------------------------
 * Returns the maximum number of sensors which can be in the array.
 *
 * Returns:
 *     maximum number of sensors (64)
 */
int array2sh_getMaxNumSensors(void);
    
/*
 * Function: array2sh_getMinNumSensors
 * -----------------------------------
 * Returns the minimum number of sensors which can be in the array.
 *
 * Input Arguments:
 *     hA2sh - array2sh handle
 * Returns:
 *     minimum number of sensors [(current_order+1)^2]
 */
int array2sh_getMinNumSensors(void* const hA2sh);
    
/*
 * Function: array2sh_getNSHrequired
 * ---------------------------------
 * Returns the number of spherical harmonic signals required by the current
 * encoding order i.e. (current_order+1)^2
 *
 * Input Arguments:
 *     hA2sh - array2sh handle
 * Returns:
 *     number of required spherical harmonic signals required by current
 *     encoding order
 */
int array2sh_getNSHrequired(void* const hA2sh);
    
/*
 * Function: array2sh_getr
 * -----------------------
 * Returns the radius of the array
 *
 * Input Arguments:
 *     hA2sh - array2sh handle
 * Returns:
 *     current array radius
 */
float array2sh_getr(void* const hA2sh);
    
/*
 * Function: array2sh_getR
 * -----------------------
 * Returns the radius of the scatterer.
 *
 * Input Arguments:
 *     hA2sh - array2sh handle
 * Returns:
 *     current scatterer radius
 */
float array2sh_getR(void* const hA2sh);
    
/*
 * Function: array2sh_getArrayType
 * -------------------------------
 * Returns the type of array. See ARRAY_TYPES enum
 *
 * Input Arguments:
 *     hA2sh   - array2sh handle
 * Returns:
 *     current array type. See ARRAY_TYPES enum
 */
int array2sh_getArrayType(void* const hA2sh);

/*
 * Function: array2sh_getWeightType
 * --------------------------------
 * Returns the type of weights to use. See WEIGHT_TYPES enum
 *
 * Input Arguments:
 *     hA2sh   - array2sh handle
 * Returns:
 *     current weight type. See WEIGHT_TYPES enum
 */
int array2sh_getWeightType(void* const hA2sh);

/*
 * Function: array2sh_getFilterType
 * --------------------------------
 * Returns the type filter design to employ for computing the encoding matrices.
 * See FILTER_TYPES enum
 *
 * Input Arguments:
 *     hA2sh   - array2sh handle
 * Returns:
 *     filter type. See FILTER_TYPES enum
 */
int array2sh_getFilterType(void* const hA2sh);

/*
 * Function: array2sh_getRegPar
 * ----------------------------
 * Returns the value of the regurlisation parameter. i.e. the maximum permitted
 * gain provided by the filters
 *
 * Input Arguments:
 *     hA2sh  - array2sh handle
 * Returns:
 *     current filter maximum gain, in DECIBELS
 */
float array2sh_getRegPar(void* const hA2sh);
    
/*
 * Function: array2sh_getChOrder
 * -----------------------------
 * Returns the Ambisonic channel ordering convention currently being used to
 * decode with, which should match the convention employed by the input signals
 *
 * Input Arguments:
 *     hA2sh - array2sh handle
 * Returns:
 *     convention currently being used (see 'CH_ORDER' enum)
 */
int array2sh_getChOrder(void* const hA2sh);

/*
 * Function: array2sh_getNormType
 * ------------------------------
 * Returns the Ambisonic normalisation convention currently being usedto decode
 * with, which should match the convention employed by the input signals.
 *
 * Input Arguments:
 *     hA2sh - array2sh handle
 * Returns:
 *     convention currently being used (see 'NORM_TYPE' enum)
 */
int array2sh_getNormType(void* const hA2sh);
    
/*
 * Function: array2sh_getc
 * -----------------------
 * Returns the speed of sound of the medium (~343m/s air, ~1480m/s water).
 *
 * Input Arguments:
 *     hA2sh - array2sh handle
 * Returns:
 *      speed of sound, in m/s
 */
float array2sh_getc(void* const hA2sh);

/*
 * Function: array2sh_getGain
 * --------------------------
 * Returns the amount of post gain to apply after the encoding
 *
 * Input Arguments:
 *     hA2sh   - array2sh handle
 * Returns:
 *     post gain, in DECIBELS
 */
float array2sh_getGain(void* const hA2sh);

/*
 * Function: array2sh_getFreqVector
 * --------------------------------
 * Returns a pointer to the frequency vector
 *
 * Input Arguments:
 *     hA2sh       - array2sh handle
 * Output Arguments:
 *     nFreqPoints - & number of frequencies
 * Returns:
 *     vector of centre frequencies; nFreqPoints x 1
 */
float* array2sh_getFreqVector(void* const hA2sh, int* nFreqPoints);
    
/*
 * Function: array2sh_getbN_inv
 * ----------------------------
 * Returns the regularised inversion of the modal coefficients per frequency.
 * May be used for optional plotting purposes.
 *
 * Input Arguments:
 *     hA2sh       - array2sh handle
 * Output Arguments:
 *     nCurves     - & number of equalisation curves (current_order+1)
 *     nFreqPoints - & number of frequencies
 * Returns:
 *     equalisation curves/regularised modal coefficients; nCurves x nFreqPoints
 */
float** array2sh_getbN_inv(void* const hA2sh, int* nCurves, int* nFreqPoints);
    
/*
 * Function: array2sh_getbN_modal
 * ------------------------------
 * Returns the direct inversion of the modal coefficients per frequency.
 * May be used for optional plotting purposes.
 *
 * Input Arguments:
 *     hA2sh       - array2sh handle
 * Output Arguments:
 *     nCurves     - & number of equalisation curves (current_order+1)
 *     nFreqPoints - & number of frequencies
 * Returns:
 *     unregularised modal coefficients; nCurves x nFreqPoints
 */
float** array2sh_getbN_modal(void* const hA2sh, int* nCurves, int* nFreqPoints);

/*
 * Function: array2sh_getSpatialCorrelation_Handle
 * -----------------------------------------------
 * Returns a pointer to the spatial correlation  [1] data. This is given per
 * frequency, and is measure of how similar the encoded spherical harmonics
 * using the current configuration is to ideal spherical harmonics. 1=perfect
 * <1: less good/ aliasing
 * Note: that this objective measure is based on analytical models of the
 * currently configured array, and may differ in practice (i.e. with a real
 * microphone array)
 *
 * Input Arguments:
 *     hA2sh       - array2sh handle
 * Output Arguments:
 *     nCurves     - & number of equalisation curves (current_order+1)
 *     nFreqPoints - & number of frequencies
 * Returns:
 *     spatial correlation per order and frequency; FLAT: nCurves x nFreqPoints
 *
 * [1] Moreau, S., Daniel, J., Bertet, S., 2006, 3D sound field recording with
 *     higher order ambisonics-objective measurements and validation of
 *     spherical microphone. In Audio Engineering Society Convention 120.
 */
float* array2sh_getSpatialCorrelation_Handle(void* const hA2sh, int* nCurves, int* nFreqPoints);

/*
 * Function: array2sh_getLevelDifference_Handle
 * --------------------------------------------
 * Returns a pointer to the level-difference  [1] data. This is given per
 * frequency, and is measure of the mean level difference between the encoded
 * spherical harmonics using the current configuration is to ideal spherical
 * harmonics
 * Note: that this objective measure is based on analytical models of the
 * currently configured array, and may differ in practice (i.e. with a real
 * microphone array)
 *
 * Input Arguments:
 *     hA2sh       - array2sh handle
 * Output Arguments:
 *     nCurves     - & number of equalisation curves (current_order+1)
 *     nFreqPoints - & number of frequencies
 * Returns:
 *     level difference per order and frequency; FLAT: nCurves x nFreqPoints
 *
 * [1] Moreau, S., Daniel, J., Bertet, S., 2006, 3D sound field recording with
 *     higher order ambisonics-objective measurements and validation of
 *     spherical microphone. In Audio Engineering Society Convention 120.
 */
float* array2sh_getLevelDifference_Handle(void* const hA2sh, int* nCurves, int* nFreqPoints);
    
/*
 * Function: array2sh_getSamplingRate
 * ----------------------------------
 * Returns the DAW/Host sample rate
 *
 * Input Arguments:
 *     hA2sh - array2sh handle
 * Returns:
 *     DAW/Host sampling rate
 */
int array2sh_getSamplingRate(void* const hA2sh);
    
/*
 * Function: array2sh_getProcessingDelay
 * -------------------------------------
 * Returns the processing delay in samples. May be used for delay compensation
 * features
 *
 * Input Arguments:
 *     hA2sh - array2sh handle
 * Returns:
 *     processing delay in samples
 */
int array2sh_getProcessingDelay(void);
   
    
#ifdef __cplusplus
} /* extern "C" { */
#endif /* __cplusplus */

#endif /* __ARRAY2SH_H_INCLUDED__ */
