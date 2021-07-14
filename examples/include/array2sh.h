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
 * @example array2sh.h
 * @brief Spatially encodes spherical microphone array signals into spherical
 *        harmonic signals (aka: Ambisonic signals) utilising theoretical
 *        encoding filters.
 *
 * ### Files
 * array2sh.h (include), array2sh_internal.h, array2sh.c, array2sh_internal.c
 * ### Example Usage
 * \code{.c}
 * int main(void) {
 *     void* hA2sh;
 *     int frameSize;
 *
 *     // Create and initialise an instance of array2sh
 *     array2sh_create(&hA2sh);
 *     array2sh_init(hA2sh, hostSamplingRate);
 *
 *     // Call any set functions, e.g.:
 *     array2sh_setPreset(hA2sh, MICROPHONE_ARRAY_PRESET_EIGENMIKE32);
 *     array2sh_setNormType(hA2sh, NORM_N3D);
 *     array2sh_setGain(hA2sh, 6.0f);
 *
 *     // The framesize of this example is fixed, and can be found with
 *     frameSize = array2sh_getFrameSize();
 *
 *     // Processing frame-by-frame
 *     ...
 *     // Load signals into inputSignalBuffer (numberOfInputs x frameSize)
 *     array2sh_process(hA2sh, inputSignalBuffer, outputSignalBuffer,
 *                      numberOfInputs, numberOfOutputs, frameSize);
 *     // Copy signals from outputSignalBuffer (numberOfOutputs x frameSize)
 *     ...
 *
 *     // Destroy this instance of array2sh
 *     array2sh_destroy(&hA2sh);
 * }
 *\endcode
 * ### Include Header
 */

/**
 * @file array2sh.h
 * @brief Spatially encodes spherical microphone array signals into spherical
 *        harmonic signals (aka: Ambisonic signals) utilising theoretical
 *        encoding filters.
 *
 * The algorithms within array2sh were pieced together and developed in
 * collaboration with Symeon Delikaris-Manias and Angelo Farina.
 * A detailed explanation of the algorithms within array2sh can be found in [1].
 * Also included, is a diffuse-field equalisation option for frequencies past
 * aliasing, developed in collaboration with Archontis Politis, 8.02.2019
 *
 * @note Since the algorithms are based on theory, only array designs where
 *       there are analytical solutions available are supported. i.e. only
 *       spherical or cylindrical arrays, which have phase-matched sensors.
 *       For more information, the reader is referred to [2,3].
 * @test test__saf_example_array2sh()
 *
 * @see [1] McCormack, L., Delikaris-Manias, S., Farina, A., Pinardi, D., and
 *          Pulkki, V., "Real-time conversion of sensor array signals into
 *          spherical harmonic signals with applications to spatially localised
 *          sub-band sound-field analysis," in Audio Engineering Society
 *          Convention 144, Audio Engineering Society, 2018.
 * @see [2] Williams EG. Fourier acoustics: sound radiation and nearfield
 *          acoustical holography. Elsevier; 1999 Jun 10.
 * @see [3] Rafaely B. Fundamentals of spherical array processing. Berlin:
 *          Springer; 2015 Feb 18.
 *
 * @author Leo McCormack
 * @date 13.09.2017
 * @license ISC
 */

#ifndef __ARRAY2SH_H_INCLUDED__
#define __ARRAY2SH_H_INCLUDED__

#include "_common.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */
    
/* ========================================================================== */
/*                             Presets + Constants                            */
/* ========================================================================== */

/** Available microphone array presets */
typedef enum {
    MICROPHONE_ARRAY_PRESET_DEFAULT = 1,
    MICROPHONE_ARRAY_PRESET_AALTO_HYDROPHONE,
    MICROPHONE_ARRAY_PRESET_SENNHEISER_AMBEO,
    MICROPHONE_ARRAY_PRESET_CORE_SOUND_TETRAMIC,
    MICROPHONE_ARRAY_PRESET_ZOOM_H3VR_PRESET,
    MICROPHONE_ARRAY_PRESET_SOUND_FIELD_SPS200,
    MICROPHONE_ARRAY_PRESET_ZYLIA_1D,
    MICROPHONE_ARRAY_PRESET_EIGENMIKE32,
    MICROPHONE_ARRAY_PRESET_DTU_MIC

}ARRAY2SH_MICROPHONE_ARRAY_PRESETS;

/**
 * Available encoding filter approaches.
 *
 * @see [1] Bernschutz, B., Porschmann, C., Spors, S., Weinzierl, S.,
 *          Versterkung,  B., 2011. Soft-limiting der modalen
 *          amplitudenverst?rkung bei sph?rischen mikrofonarrays im plane wave
 *          decomposition verfahren. Proceedings of the 37. Deutsche
 *          Jahrestagung fur Akustik (DAGA 2011)
 * @see [2] Moreau, S., Daniel, J., Bertet, S., 2006, 3D sound field recording
 *          with higher order ambisonics-objective measurements and validation
 *          of spherical microphone. In Audio Engineering Society Convention
 *          120.
 * @see [3] Zotter, F. A Linear-Phase Filter-Bank Approach to Process Rigid
 *          Spherical  Microphone Array Recordings.
 */
typedef enum {
    FILTER_SOFT_LIM = 1, /**< Encoding filters based on a 'soft-limiting'
                          *   regularised inversion of the modal responses [1]
                          */
    FILTER_TIKHONOV,     /**< Encoding filters based on a 'Tikhonov' regularised
                          *   inversion of the modal responses [2] */
    FILTER_Z_STYLE,      /**< Encoding filters based on a linear-phase filter-
                          *   bank approach [3] */
    FILTER_Z_STYLE_MAXRE /**< Same as #FILTER_Z_STYLE, only it also has max_rE
                          *   weights baked in */
    
}ARRAY2SH_FILTER_TYPES;
    
/** Number of available filter types */
#define ARRAY2SH_NUM_FILTER_TYPES ( 4 )

/**
 * List of supported array types.
 *
 * @note Although supported, cylindrical arrays have not really been tested as
 *       we don't own one. Just to keep in mind.
 */
typedef enum {
    ARRAY_SPHERICAL = 1, /**< Spherical arrangement of sensors (open/rigid) */
    ARRAY_CYLINDRICAL    /**< Cylindrial arrangement of sensors (open/rigid) */
    
}ARRAY2SH_ARRAY_TYPES;
    
/** Number of supported array types */
#define ARRAY2SH_NUM_ARRAY_TYPES ( 2 )

/** List of supported sensor directivities and array construction types */
typedef enum {
    WEIGHT_RIGID_OMNI = 1, /**< Rigid baffle construction with omni sensors */
    WEIGHT_RIGID_CARD,     /**< Rigid baffle construction with cardioid sensors
                            */
    WEIGHT_RIGID_DIPOLE,   /**< Rigid baffle construction with dipole sensors */
    WEIGHT_OPEN_OMNI,      /**< Open array construction with omni sensors */
    WEIGHT_OPEN_CARD,      /**< Open array construction with cardioid sensors */
    WEIGHT_OPEN_DIPOLE     /**< Open array construction with dipole sensors */
    
}ARRAY2SH_WEIGHT_TYPES;
    
/** Number of supported sensor directivities and array construction types */
#define ARRAY2SH_NUM_WEIGHT_TYPES ( 6 )
    
/**
 * Current status of the encoder evaluation output data
 *
 * These are some objective metrics which you can use to ascertain the
 * performance of the microphone array and the encoding.
 */
typedef enum {
    EVAL_STATUS_EVALUATED = 0,      /**< Encoder has been evaluated */
    EVAL_STATUS_RECENTLY_EVALUATED, /**< Encoder has recently been evaluated */
    EVAL_STATUS_NOT_EVALUATED,      /**< Encoder has not been evaluated */
    EVAL_STATUS_EVALUATING          /**< Encoder is being evaluated */

}ARRAY2SH_EVAL_STATUS;

/** Maximum number of sensors supported */
#define ARRAY2SH_MAX_NUM_SENSORS ( 64 )

/** Minimum gain value used for regularised inverse of modal coeffs, dB */
#define ARRAY2SH_MAX_GAIN_MIN_VALUE ( 0.0f )

/** Maximum gain value used for regularised inverse of modal coeffs, dB */
#define ARRAY2SH_MAX_GAIN_MAX_VALUE ( 80.0f )

/** Minimum post-gain, dB */
#define ARRAY2SH_POST_GAIN_MIN_VALUE ( -60.0f )

/** Maximum post-gain, dB */
#define ARRAY2SH_POST_GAIN_MAX_VALUE ( 12.0f )

/** Minimum speed of sound value, m/s */
#define ARRAY2SH_SPEED_OF_SOUND_MIN_VALUE ( 200.0f )

/** Maximum speed of sound value, m/s */
#define ARRAY2SH_SPEED_OF_SOUND_MAX_VALUE ( 2000.0f )

/** Minimum array radius supported, mm */
#define ARRAY2SH_ARRAY_RADIUS_MIN_VALUE ( 1.0f )

/** Maximum array radius supported, mm */
#define ARRAY2SH_ARRAY_RADIUS_MAX_VALUE ( 400.0f )

/** Minimum baffle radius supported, mm */
#define ARRAY2SH_BAFFLE_RADIUS_MIN_VALUE ( 1.0f )

/** Maximum baffle radius supported, mm */
#define ARRAY2SH_BAFFLE_RADIUS_MAX_VALUE ( 400.0f ) 


/* ========================================================================== */
/*                               Main Functions                               */
/* ========================================================================== */

/**
 * Creates an instance of array2sh
 *
 * @param[in] phA2sh (&) address of array2sh handle
 */
void array2sh_create(void** const phA2sh);

/**
 * Destroys an instance of array2sh
 *
 * @param[in] phA2sh (&) address of array2sh handle
 */
void array2sh_destroy(void** const phA2sh);

/**
 * Initialises an instance of array2sh with default settings
 *
 * @param[in] hA2sh      array2sh handle
 * @param[in] samplerate Host samplerate.
 */
void array2sh_init(void* const hA2sh,
                   int samplerate);

/**
 * Evaluates the encoder, based on current global/user parameters
 *
 * @param[in] hA2sh array2sh handle
 */
void array2sh_evalEncoder(void* const hA2sh);

/**
 * Spatially encode microphone/hydrophone array signals into spherical harmonic
 * signals
 *
 * @param[in] hA2sh     array2sh handle
 * @param[in] inputs    Input channel buffers; 2-D array: nInputs x nSamples
 * @param[in] outputs   Output channel buffers; 2-D array: nOutputs x nSamples
 * @param[in] nInputs   Number of input channels
 * @param[in] nOutputs  Number of output channels
 * @param[in] nSamples  Number of samples in 'inputs'/'output' matrices
 */
void array2sh_process(void* const hA2sh,
                      const float *const * inputs,
                      float** const outputs,
                      int nInputs,
                      int nOutputs,
                      int nSamples);


/* ========================================================================== */
/*                                Set Functions                               */
/* ========================================================================== */

/**
 * Sets all intialisation flags to 1; re-initialising all settings/variables
 * as array2sh is currently configured, at next available opportunity.
 */
void array2sh_refreshSettings(void* const hA2sh);
 
/** Sets the encoding order (see #SH_ORDERS enum) */
void array2sh_setEncodingOrder(void* const hA2sh, int newOrder);

/**
 * Evaluates the performance of the current encoding filters when applied to a
 * theoretical model of the currently configured array; two established
 * objective metrics are then computed; more information in [1]
 *
 * @see [1] Moreau, S., Daniel, J., Bertet, S., 2006, 3D sound field recording
 *          with higher order ambisonics-objective measurements and validation
 *          of spherical microphone. In Audio Engineering Society Convention
 *          120.
 */
void array2sh_setRequestEncoderEvalFLAG(void* const hA2sh, int newState);
    
/** Sets current eval status (see #ARRAY2SH_EVAL_STATUS enum) */
void array2sh_setEvalStatus(void* const hA2sh, ARRAY2SH_EVAL_STATUS evalStatus);

/**
 * Analyses what the theoretical spatial aliasing frequency is, and conducts
 * diffuse-field equalisation above this (enable: 1, disable: 0).
 *
 * Thanks to Dr. Archontis Politis for suggesting and designing this feature.
 */
void array2sh_setDiffEQpastAliasing(void* const hA2sh, int newState);

/**
 * Sets a pre-defined microphone/hydrophone array preset (See
 * #ARRAY2SH_MICROPHONE_ARRAY_PRESETS enum)
 */
void array2sh_setPreset(void* const hA2sh,
                        ARRAY2SH_MICROPHONE_ARRAY_PRESETS preset);
    
/**
 * Sets a particular sensor's azimuth (radians) w.r.t to the origin of the
 * array.
 *
 * @param[in] hA2sh      array2sh handle
 * @param[in] index      Sensor index
 * @param[in] newAzi_rad Sensor azimuth in RADIANS
 */
void array2sh_setSensorAzi_rad(void* const hA2sh, int index, float newAzi_rad);
    
/**
 * Sets a particular sensor's elevation (radians) w.r.t to the origin of the
 * array.
 *
 * @param[in] hA2sh       array2sh handle
 * @param[in] index       Sensor index
 * @param[in] newElev_rad Sensor elevation in RADIANS
 */
void array2sh_setSensorElev_rad(void* const hA2sh, int index, float newElev_rad);
    
/**
 * Sets a particular sensor's azimuth (degrees) w.r.t to the origin of the
 * array.
 *
 * Input Arguments:
 *     hA2sh      - array2sh handle
 *     index      - sensor index
 *     newAzi_deg - sensor azimuth in DEGREES
 */
void array2sh_setSensorAzi_deg(void* const hA2sh, int index, float newAzi_deg);
    
/**
 * Sets a particular sensor's elevation (degrees) w.r.t to the origin of the
 * array.
 *
 * @param[in] hA2sh       array2sh handle
 * @param[in] index       Sensor index
 * @param[in] newElev_deg Sensor elevation in DEGREES
 */
void array2sh_setSensorElev_deg(void* const hA2sh, int index, float newElev_deg);
    
/** Sets the number of sensors in the array. */
void array2sh_setNumSensors(void* const hA2sh, int newQ);
    
/** Sets the radius of the array */
void array2sh_setr(void* const hA2sh, float newr);
    
/**
 * Sets the radius (in meters) of the scatterer (only for Rigid arrays).
 *
 * @note R <= r. i.e. the sensors may protrude from the rigid scattering
 *       surface, or be flush with the surface of the array
 */
void array2sh_setR(void* const hA2sh, float newR);
    
/** Sets the type of array (see #ARRAY2SH_ARRAY_TYPES enum) */
void array2sh_setArrayType(void* const hA2sh, int newType);

/** Sets the type of weights to use (see #ARRAY2SH_WEIGHT_TYPES enum) */
void array2sh_setWeightType(void* const hA2sh, int newType);
    
/**
 * Sets the type filter design to employ for computing the encoding matrices
 * (see #ARRAY2SH_FILTER_TYPES enum)
 */
void array2sh_setFilterType(void* const hA2sh, int newType);
    
/**
 * Sets the value of the regularisation parameter (the maximum permitted gain
 * of the filters), in DECIBELS
 */
void array2sh_setRegPar(void* const hA2sh, float newVal);
    
/**
 * Sets the Ambisonic channel ordering convention to encode with, in order to
 * match the convention employed by the input signals (see #CH_ORDER enum)
 */
void array2sh_setChOrder(void* const hA2sh, int newOrder);
    
/**
 * Sets the Ambisonic normalisation convention to encode with, in order to match
 * with the convention employed by the input signals (see #NORM_TYPES enum)
 */
void array2sh_setNormType(void* const hA2sh, int newType);

/**
 * Sets the speed of sound of the medium (~343m/s air, ~1480m/s water), in m/s
 */
void array2sh_setc(void* const hA2sh, float newc);
    
/** Sets the amount of post gain to apply after the encoding, in DECIBELS */
void array2sh_setGain(void* const hA2sh, float newGain);


/* ========================================================================== */
/*                                Get Functions                               */
/* ========================================================================== */

/**
 * Returns the processing framesize (i.e., number of samples processed with
 * every _process() call )
 */
int array2sh_getFrameSize(void);

/** Returns current eval status (see #ARRAY2SH_EVAL_STATUS enum) */
ARRAY2SH_EVAL_STATUS array2sh_getEvalStatus(void* const hA2sh);

/**
 * Returns 0 if SHT is not be reinitialised, 1: if it is
 */
int array2sh_getReinitSHTmatrixFLAG(void* const hA2sh);

/**
 * (Optional) Returns current intialisation/processing progress, between 0..1
 *  - 0: intialisation/processing has started
 *  - 1: intialisation/processing has ended
 */
float array2sh_getProgressBar0_1(void* const hA2sh);

/**
 * (Optional) Returns current intialisation/processing progress text
 *
 * @note "text" string should be (at least) of length:
 *       #PROGRESSBARTEXT_CHAR_LENGTH
 */
void array2sh_getProgressBarText(void* const hA2sh, char* text);

/**
 * Flag to enabled/disable diffuse equalisation above the spatial aliasing
 * limit of the array (0: disabled, 1: enabled).
 *
 * Developed in collaboration with Archontis Politis.
 *
 * @note In general, theoretical encoding filters have a tendency to boost
 *       the aliased frequencies. Whereas, measurement-based filters (through a
 *       least-squares solution), tend to attenuate them. Neither of these
 *       are correct or incorrect, since, strictly (spatially) speaking, we
 *       should be placing a low-pass filter at the spatial aliasing frequency.
 *       However, since we usually do not want to remove this high frequency
 *       energy from e.g. an Ambisonic reproduction, we would argue that
 *       equalising the aliased components so that they have a flat spectrum,
 *       is probably the way to go; and is exactly what this feature does.
 */
int array2sh_getDiffEQpastAliasing(void* const hA2sh);

/**
 * Returns a flag, which is '1' if there has been a recent request to evaluate
 * the current encoding performance, or '0', if there hasn't.
 */
int array2sh_getRequestEncoderEvalFLAG(void* const hA2sh);
    
/** Returns the current encoding order (see #SH_ORDERS enum) */
int array2sh_getEncodingOrder(void* const hA2sh);
    
/**
 * Returns a particular sensor's azimuth w.r.t to the origin of the array, in
 * RADIANS
 */
float array2sh_getSensorAzi_rad(void* const hA2sh, int index);
    
/**
 * Returns a particular sensor's elevation w.r.t to the origin of the array, in
 * RADIANS
 */
float array2sh_getSensorElev_rad(void* const hA2sh, int index);
    
/**
 * Returns a particular sensor's azimuth w.r.t to the origin of the array, in
 * DEGREES
 */
float array2sh_getSensorAzi_deg(void* const hA2sh, int index);
    
/**
 * Returns a particular sensor's elevation w.r.t to the origin of the array, in
 * DEGREES
 */
float array2sh_getSensorElev_deg(void* const hA2sh, int index);

/** Returns the number of sensors in the array */
int array2sh_getNumSensors(void* const hA2sh);
    
/** Returns the maximum supported number of sensors which can be in the array */
int array2sh_getMaxNumSensors(void);
    
/**
 * Returns the minimum number of sensors which can be in the array
 * [(current_order+1)^2]
 */
int array2sh_getMinNumSensors(void* const hA2sh);
    
/**
 * Returns the number of spherical harmonic signals required by the current
 * encoding order: (current_order+1)^2
 */
int array2sh_getNSHrequired(void* const hA2sh);
    
/** Returns the radius of the array, in meters */
float array2sh_getr(void* const hA2sh);
    
/** Returns the radius of the scatterer, in meters */
float array2sh_getR(void* const hA2sh);
    
/** Returns the type of array. See #ARRAY2SH_ARRAY_TYPES enum */
int array2sh_getArrayType(void* const hA2sh);

/** Returns the type of weights to use see #ARRAY2SH_WEIGHT_TYPES enum */
int array2sh_getWeightType(void* const hA2sh);

/**
 * Returns the type filter design to employ for computing the encoding matrices
 * (see #ARRAY2SH_FILTER_TYPES enum)
 */
int array2sh_getFilterType(void* const hA2sh);

/**
 * Returns the value of the regurlisation parameter; the maximum permitted
 * gain provided by the filters, in DECIBELS
 */
float array2sh_getRegPar(void* const hA2sh);
    
/**
 * Returns the Ambisonic channel ordering convention currently being used to
 * decode with, which should match the convention employed by the input signals
 * (see #CH_ORDER enum)
 */
int array2sh_getChOrder(void* const hA2sh);

/**
 * Returns the Ambisonic normalisation convention currently being usedto decode
 * with, which should match the convention employed by the input signals
 * (see #NORM_TYPES enum)
 */
int array2sh_getNormType(void* const hA2sh);
    
/**
 * Returns the speed of sound of the medium (~343m/s air, ~1480m/s water), in
 * m/s
 */
float array2sh_getc(void* const hA2sh);

/** Returns the amount of post gain to apply after the encoding, in DECIBELS */
float array2sh_getGain(void* const hA2sh);

/**
 * Returns a pointer to the frequency vector
 *
 * @param[in]  hA2sh       array2sh handle
 * @param[out] nFreqPoints (&) number of frequencies
 * @returns                Vector of centre frequencies; nFreqPoints x 1
 */
float* array2sh_getFreqVector(void* const hA2sh, int* nFreqPoints);
    
/**
 * Returns the regularised inversion of the modal coefficients per frequency
 * (may be used for optional plotting purposes)
 *
 * @param[in]  hA2sh       array2sh handle
 * @param[out] nCurves     (&) number of equalisation curves (current_order+1)
 * @param[out] nFreqPoints (&) number of frequencies
 * @returns                Equalisation curves/regularised modal coefficients;
 *                         nCurves x nFreqPoints
 */
float** array2sh_getbN_inv(void* const hA2sh, int* nCurves, int* nFreqPoints);
    
/**
 * Returns the direct inversion of the modal coefficients per frequency
 * (may be used for optional plotting purposes)
 *
 * @param[in]  hA2sh       array2sh handle
 * @param[out] nCurves     (&) number of equalisation curves (current_order+1)
 * @param[out] nFreqPoints (&) number of frequencies
 * @returns                Unregularised modal coefficients;
 *                         nCurves x nFreqPoints
 */
float** array2sh_getbN_modal(void* const hA2sh, int* nCurves, int* nFreqPoints);

/**
 * Returns a pointer to the spatial correlation [1] data. This is given per
 * frequency, and is measure of how similar the encoded spherical harmonics
 * using the current configuration is to ideal spherical harmonics. 1=perfect
 * <1: less good/ aliasing
 *
 * @note This objective measure is based on analytical models of the currently
 *       configured array, and may differ in practice (i.e. with a real
 *       microphone array)
 *
 * @param[in]  hA2sh       array2sh handle
 * @param[out] nCurves     (&) number of equalisation curves (current_order+1)
 * @param[out] nFreqPoints (&) number of frequencies
 * @returns                Spatial correlation per order and frequency;
 *                         FLAT: nCurves x nFreqPoints
 *
 * @see [1] Moreau, S., Daniel, J., Bertet, S., 2006, 3D sound field recording
 *          with higher order ambisonics-objective measurements and validation
 *          of spherical microphone. In Audio Engineering Society Convention
 *          120.
 */
float* array2sh_getSpatialCorrelation_Handle(void* const hA2sh,
                                             int* nCurves,
                                             int* nFreqPoints);

/**
 * Returns a pointer to the level-difference [1] data. This is given per
 * frequency, and is measure of the mean level difference between the encoded
 * spherical harmonics using the current configuration is to ideal spherical
 * harmonics
 *
 * @note This objective measure is based on analytical models of the currently
 *       configured array, and may differ in practice (i.e. with a real
 *       microphone array)
 *
 * @param[in]  hA2sh       array2sh handle
 * @param[out] nCurves     (&) number of equalisation curves (current_order+1)
 * @param[out] nFreqPoints (&) number of frequencies
 * @returns                Level difference per order and frequency;
 *                         FLAT: nCurves x nFreqPoints
 *
 * @see [1] Moreau, S., Daniel, J., Bertet, S., 2006, 3D sound field recording
 *          with higher order ambisonics-objective measurements and validation
 *          of spherical microphone. In Audio Engineering Society Convention
 *          120.
 */
float* array2sh_getLevelDifference_Handle(void* const hA2sh, int* nCurves, int* nFreqPoints);
    
/** Returns the DAW/Host sample rate */
int array2sh_getSamplingRate(void* const hA2sh);
    
/**
 * Returns the processing delay in samples (may be used for delay compensation
 * features) 
 */
int array2sh_getProcessingDelay(void);
   
    
#ifdef __cplusplus
} /* extern "C" { */
#endif /* __cplusplus */

#endif /* __ARRAY2SH_H_INCLUDED__ */
