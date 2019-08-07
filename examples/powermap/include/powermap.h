/*
 * Copyright 2016-2018 Leo McCormack
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
 * Filename: powermap.h (include header)
 * -------------------------------------
 * A sound-field visualiser, which utilises spherical harmonic signals as input.
 *
 * Dependencies:
 *     saf_utilities, afSTFTlib, saf_vbap, saf_sh
 * Author, date created:
 *     Leo McCormack, 26.04.2016
 */

#ifndef __POWERMAP_H_INCLUDED__
#define __POWERMAP_H_INCLUDED__

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */
    
/* ========================================================================== */
/*                             Presets + Constants                            */
/* ========================================================================== */

#define POWERMAP_MAX_NUM_INPUT_CHANNELS ( 64 )

/* Microphone/Hydrophone array options */
#define ENABLE_ZYLIA_MIC_PRESET
#define ENABLE_EIGENMIKE32_MIC_PRESET
#define ENABLE_DTU_MIC_MIC_PRESET
 
/*
 * Enum: MASTER_ORDERS
 * -------------------
 * "Master order" relates to the current maximum order to expect. However, the
 * analysis order can be lower for a given frequency, due to the
 * "analysisOrderPerBand" vector, which can contain lower values than the
 * master order, but not higher.
 *
 * Options:
 *     MASTER_ORDER_FIRST   - First-order analysis (4 channel input)
 *     MASTER_ORDER_SECOND  - Second-order analysis (9 channel input)
 *     MASTER_ORDER_THIRD   - Third-order analysis (16 channel input)
 *     MASTER_ORDER_FOURTH  - Fourth-order analysis (25 channel input)
 *     MASTER_ORDER_FIFTH   - Fifth-order analysis (36 channel input)
 *     MASTER_ORDER_SIXTH   - Sixth-order analysis (49 channel input)
 *     MASTER_ORDER_SEVENTH - Seventh-order analysis (64 channel input)
 */
typedef enum _MASTER_ORDERS{
    MASTER_ORDER_FIRST = 1,
    MASTER_ORDER_SECOND,
    MASTER_ORDER_THIRD,
    MASTER_ORDER_FOURTH,
    MASTER_ORDER_FIFTH,
    MASTER_ORDER_SIXTH,
    MASTER_ORDER_SEVENTH
    
}MASTER_ORDERS;

typedef enum _MIC_PRESETS{
    MIC_PRESET_IDEAL = 1
#ifdef ENABLE_ZYLIA_MIC_PRESET
    , MIC_PRESET_ZYLIA
#endif
#ifdef ENABLE_EIGENMIKE32_MIC_PRESET
    , MIC_PRESET_EIGENMIKE32
#endif
#ifdef ENABLE_DTU_MIC_MIC_PRESET
    , MIC_PRESET_DTU_MIC
#endif
}MIC_PRESETS;

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
typedef enum _NORM_TYPES{
    NORM_N3D = 1,
    NORM_SN3D,
    NORM_FUMA   /* first-order only */
    
}NORM_TYPES;

/*
 * Enum: POWERMAP_MODES
 * --------------------
 * Available power-map/activity-map options
 *
 * Options:
 *     PM_MODE_PWD         - activity-map based on the energy of hyper-cardioid
 *                           [plane-wave decomposition (PWD)] beamformers
 *     PM_MODE_MVDR        - activity-map based on the energy of minimum-
 *                           variance distortionless response (MVDR) beamformers
 *     PM_MODE_CROPAC_LCMV - Experimental! activity-map based on a linearly-
 *                           contrained minimum-variance (LCMV) formulation of
 *                           the Cross-Pattern Coherence (CroPaC) spatial filter
 *     PM_MODE_MUSIC       - activity-map based on the sub-space method:
 *                           multiple signal classification (MUSIC)
 *     PM_MODE_MUSIC_LOG   - same as PM_MODE_MUSIC, but log(out_values)
 *     PM_MODE_MINNORM     - activity-map based on the sub-space method:
 *                           minimum-norm (Min-Norm)
 *     PM_MODE_MINNORM_LOG - same as PM_MODE_MINNORM, but log(out_values)
 */
typedef enum _POWERMAP_MODES {
    PM_MODE_PWD = 1,
    PM_MODE_MVDR,
    PM_MODE_CROPAC_LCMV,
    PM_MODE_MUSIC,
    PM_MODE_MUSIC_LOG, 
    PM_MODE_MINNORM,
    PM_MODE_MINNORM_LOG
    
} POWERMAP_MODES;
    
/*
 * Enum: HFOV_OPTIONS
 * ------------------
 * Available horizontal feild-of-view (FOV) options
 *
 * Options:
 *     HFOV_360 - 360 degrees
 */
typedef enum _HFOV_OPTIONS{
    HFOV_360 = 1
    
}HFOV_OPTIONS;
    
/*
 * Enum: ASPECT_RATIO_OPTIONS
 * --------------------------
 * Available aspect ratios
 *
 * Options:
 *     ASPECT_RATIO_2_1  - 2:1
 */
typedef enum _ASPECT_RATIO_OPTIONS{
    ASPECT_RATIO_2_1 = 1
    
}ASPECT_RATIO_OPTIONS;
    

/* ========================================================================== */
/*                               Main Functions                               */
/* ========================================================================== */

/*
 * Function: powermap_create
 * -------------------------
 * Creates an instance of the powermap
 *
 * Input Arguments:
 *     phPm - & address of powermap handle
 */
void powermap_create(void** const phPm);

/*
 * Function: powermap_destroy
 * --------------------------
 * Destroys an instance of the powermap
 *
 * Input Arguments:
 *     phPm - & address of powermap handle
 */
void powermap_destroy(void** const phPm);

/*
 * Function: powermap_init
 * -----------------------
 * Initialises an instance of powermap with default settings
 *
 * Input Arguments:
 *     hPm        - powermap handle
 *     samplerate - host samplerate.
 */
void powermap_init(void* const hPm,
                   float  samplerate);

/*
 * Function: powermap_process
 * --------------------------
 * Analyses the input spherical harmonic signals to generate an activity-map
 *
 * Input Arguments:
 *     hPm       - powermap handle
 *     inputs    - input channel buffers; 2-D array: nInputs x nSamples
 *     nInputs   - number of input channels
 *     nSamples  - number of samples in 'inputs'/'output' matrices
 *     isPlaying - flag to say if there is audio in the input buffers, 0: no
 *                 audio, reduced processing, 1: audio, full processing
 */
void powermap_analysis(void* const hPm,
                       float** const inputs,
                       int nInputs,
                       int nSamples,
                       int isPlaying);
    
   
/* ========================================================================== */
/*                                Set Functions                               */
/* ========================================================================== */

/*
 * Function: powermap_refreshSettings
 * ----------------------------------
 * Sets all intialisation flags to 1. i.e. re-initialise all settings/variables
 * as powermap is currently configured, at next available opportunity.
 *
 * Input Arguments:
 *     hPm - powermap handle
 */
void powermap_refreshSettings(void* const hPm);

/*
 * Function: powermap_checkReInit
 * ------------------------------
 * Check if any reInit Flags are active, and reinitialise if they are.
 * Note: Only call when playback has stopped.
 *
 * Input Arguments:
 *     hPm - powermap handle
 */
void powermap_checkReInit(void* const hPm);

/*
 * Function: powermap_setPowermapMode
 * ----------------------------------
 * Sets the powermap/activity-map approach, See "POWERMAP_MODES" enum
 *
 * Input Arguments:
 *     hPm     - powermap handle
 *     newMode - new mode. See "POWERMAP_MODES" enum.
 */
void powermap_setPowermapMode(void* const hPm, int newMode);

/*
 * Function: powermap_setMasterOrder
 * ---------------------------------
 * Sets the maximum input/analysis order
 *
 * Input Arguments:
 *     hPm      - powermap handle
 *     newValue - new maximum analysis order (see 'MASTER_ORDERS' enum)
 */
void powermap_setMasterOrder(void* const hPm,  int newValue);

/*
 * Function: powermap_setAnaOrder
 * ------------------------------
 * Sets the input/analysis order for one specific frequency band.
 *
 * Input Arguments:
 *     hPm      - powermap handle
 *     newValue - new analysis order
 *     bandIdx  - frequency band index
 */
void powermap_setAnaOrder(void* const hPm,  int newValue, int bandIdx);

/*
 * Function: powermap_setAnaOrderAllBands
 * --------------------------------------
 * Sets the input/analysis order for all frequency bands.
 *
 * Input Arguments:
 *     hPm      - powermap handle
 *     newValue - new analysis order
 */
void powermap_setAnaOrderAllBands(void* const hPm, int newValue);
    
/*
 * Function: powermap_setPowermapEQ
 * --------------------------------
 * Sets the weighting coefficient for a particular frequency band, allowing
 * one to "equalise" the activity-map.
 *
 * Input Arguments:
 *     hPm      - powermap handle
 *     newValue - new weight value
 *     bandIdx  - frequency band index
 */
void powermap_setPowermapEQ(void* const hPm,  float newValue, int bandIdx);

/*
 * Function: powermap_setPowermapEQAllBands
 * ----------------------------------------
 * Sets the weighting coefficient for all frequency bands.
 *
 * Input Arguments:
 *     hPm      - powermap handle
 *     newValue - new weight value
 */
void powermap_setPowermapEQAllBands(void* const hPm,  float newValue);
    
/*
 * Function: powermap_setCovAvgCoeff
 * ---------------------------------
 * Sets the covariance matrix averaging coefficient
 *
 * Input Arguments:
 *     hPm    - powermap handle
 *     newAvg - new averaging coefficient
 */
void powermap_setCovAvgCoeff(void* const hPm, float newAvg);

/*
 * Function: powermap_setChOrder
 * -----------------------------
 * Sets the Ambisonic channel ordering convention to decode with, in order to
 * match the convention employed by the input signals
 *
 * Input Arguments:
 *     hPm      - powermap handle
 *     newOrder - convention to use (see 'CH_ORDER' enum)
 */
void powermap_setChOrder(void* const hPm, int newOrder);

/*
 * Function: powermap_setNormType
 * ------------------------------
 * Sets the Ambisonic normalisation convention to decode with, in order to match
 * with the convention employed by the input signals.
 *
 * Input Arguments:
 *     hPm     - powermap handle
 *     newType - convention to use (see 'NORM_TYPE' enum)
 */
void powermap_setNormType(void* const hPm, int newType);

/*
 * Function: powermap_setSourcePreset
 * ----------------------------------
 * Sets an input preset, i.e. the microphone/hyrophone array used to capture
 * the input signals.
 *
 * Input Arguments:
 *     hPm         - powermap handle
 *     newPresetID - preset to use (see 'MIC_PRESETS' enum)
 */
void powermap_setSourcePreset(void* const hPm, int newPresetID);

/*
 * Function: powermap_setNumSources
 * --------------------------------
 * Sets the number of sources present in the input sound scene.
 *
 * Input Arguments:
 *     hPm      - powermap handle
 *     newValue - new number of sources in the sound scene
 */
void powermap_setNumSources(void* const hPm, int newValue);

/*
 * Function: powermap_setDispFOV
 * -----------------------------
 * Sets the visualisation display window horizontal field-of-view (FOV).
 *
 * Input Arguments:
 *     hPm       - powermap handle
 *     newOption - horizontal FOV to use (see 'HFOV_OPTIONS' enum)
 */
void powermap_setDispFOV(void* const hPm, int newOption);

/*
 * Function: powermap_setAspectRatio
 * ---------------------------------
 * Sets the visualisation display window aspect-ratio
 *
 * Input Arguments:
 *     hPm       - powermap handle
 *     newOption - aspect ratio to use (see 'ASPECT_RATIO_OPTIONS' enum)
 */
void powermap_setAspectRatio(void* const hPm, int newOption);

/*
 * Function: powermap_setPowermapAvgCoeff
 * --------------------------------------
 * Sets the activity-map averaging coefficient.
 *
 * Input Arguments:
 *     hPm     - powermap handle
 *     newValue - new averaging coefficient, 0..1
 */
void powermap_setPowermapAvgCoeff(void* const hPm, float newValue);

/*
 * Function: powermap_requestPmapUpdate
 * ------------------------------------
 * Informs powermap that it should compute a new activity-map at its own
 * convenience, if it would be so kind. Thank you, god bless.
 *
 * Input Arguments:
 *     hPm - powermap handle
 */
void powermap_requestPmapUpdate(void* const hPm);


/* ========================================================================== */
/*                                Get Functions                               */
/* ========================================================================== */

/*
 * Function: powermap_getMasterOrder
 * ---------------------------------
 * Returns the current maximum analysis/input order.
 *
 * Input Arguments:
 *     hPm - powermap handle
 * Returns:
 *     current maximum analysis order (see 'MASTER_ORDERS' enum)
 */
int powermap_getMasterOrder(void* const hPm);

/*
 * Function: powermap_getPowermapMode
 * ----------------------------------
 * Returns the powermap/activity-map mode to employed for the analysis.
 *
 * Input Arguments:
 *     hPm - powermap handle
 * Returns:
 *     current type. See "BEAM_TYPE" enum.
 */
int powermap_getPowermapMode(void* const hPm);

/*
 * Function: powermap_getSamplingRate
 * ----------------------------------
 * Returns the current sampling rate
 *
 * Input Arguments:
 *     hPm - powermap handle
 *  Returns:
 *     current sampling rate, in Hz
 */
int powermap_getSamplingRate(void* const hPm);

/*
 * Function: powermap_getCovAvgCoeff
 * ---------------------------------
 * Returns the current covariance averaging coefficient value
 *
 * Input Arguments:
 *     hPm - powermap handle
 *  Returns:
 *     current sampling rate, in Hz
 */
float powermap_getCovAvgCoeff(void* const hPm);

/*
 * Function: powermap_getNumberOfBands
 * -----------------------------------
 * Returns the number of frequency bands used for the analysis
 *
 *  Returns:
 *     Returns the number of frequency bands
 */
int powermap_getNumberOfBands(void);
    
/*
 * Function: powermap_getNSHrequired
 * -------------------------------
 * Returns the number of spherical harmonic signals required by the current
 * analysis order i.e. (current_order + 1)^2
 *
 * Input Arguments:
 *     hPm - powermap handle
 * Returns:
 *     number of required spherical harmonic signals required by current
 *     analysis order
 */
int powermap_getNSHrequired(void* const hPm);

/*
 * Function: powermap_getPowermapEQ
 * --------------------------------
 * Returns the weighting coefficient for a particular frequency band, allowing
 * one to "equalise" the activity-map.
 *
 * Input Arguments:
 *     hPm     - powermap handle
 *     bandIdx - frequency band index
 * Returns:
 *     current weight value
 */
float powermap_getPowermapEQ(void* const hPm, int bandIdx);

/*
 * Function: powermap_getPowermapEQAllBands
 * ----------------------------------------
 * Returns the weighting coefficient for the first frequency band
 *
 * Input Arguments:
 *     hPm - powermap handle
 * Returns:
 *     current weight value
 */
float powermap_getPowermapEQAllBands(void* const hPm);

/*
 * Function: powermap_getPowermapEQHandle
 * --------------------------------------
 * Returns the weighting coefficient for all frequency bands
 *
 * Input Arguments:
 *     hPm       - powermap handle
 * Output Arguments:
 *     pX_vector - & frequency vector; pNpoints x 1
 *     pY_values - & weighting coefficients; pNpoints x 1
 *     pNpoints  - & number of frequency bands
 */
void powermap_getPowermapEQHandle(void* const hPm,
                                  float** pX_vector,
                                  float** pY_values,
                                  int* pNpoints);

/*
 * Function: powermap_getAnaOrder
 * ------------------------------
 * Returns the input/analysis order for one specific frequency band.
 *
 * Input Arguments:
 *     hPm     - powermap handle
 *     bandIdx - frequency band index
 * Returns:
 *     current analysis order at this band index
 */
int powermap_getAnaOrder(void* const hPm, int bandIdx);

/*
 * Function: powermap_getAnaOrderAllBands
 * --------------------------------------
 * Returns the input/analysis order for the first frequency band
 *
 * Input Arguments:
 *     hPm - powermap handle
 * Returns:
 *     current analysis order
 */
int powermap_getAnaOrderAllBands(void* const hPm);

/*
 * Function: powermap_getAnaOrderHandle
 * ------------------------------------
 * Returns the input/analysis order for all frequency bands
 *
 * Input Arguments:
 *     hPm       - powermap handle
 * Output Arguments:
 *     pX_vector - & frequency vector; pNpoints x 1
 *     pY_values - & input/analysis orders; pNpoints x 1
 *     pNpoints  - & number of frequency bands
 */
void powermap_getAnaOrderHandle(void* const hPm,
                                float** pX_vector,
                                int** pY_values,
                                int* pNpoints);

/*
 * Function: powermap_getChOrder
 * -----------------------------
 * Returns the Ambisonic channel ordering convention currently being used to
 * decode with, which should match the convention employed by the input signals
 *
 * Input Arguments:
 *     hPm - powermap handle
 * Returns:
 *     convention currently being used (see 'CH_ORDER' enum)
 */
int powermap_getChOrder(void* const hPm);

/*
 * Function: powermap_getNormType
 * ------------------------------
 * Returns the Ambisonic normalisation convention currently being usedto decode
 * with, which should match the convention employed by the input signals.
 *
 * Input Arguments:
 *     hPm - powermap handle
 * Returns:
 *     convention currently being used (see 'NORM_TYPE' enum)
 */
int powermap_getNormType(void* const hPm);

/*
 * Function: powermap_getNumSources
 * --------------------------------
 * Returns the number of sources present in the input sound scene.
 *
 * Input Arguments:
 *     hPm      - powermap handle
 * Returns:
 *     number of sources in the sound scene
 */
int powermap_getNumSources(void* const hPm);

/*
 * Function: powermap_getDispFOV
 * -----------------------------
 * Returns the current visualisation display window horizontal field-of-view
 * (FOV).
 *
 * Input Arguments:
 *     hPm - powermap handle
 * Returns:
 *     horizontal FOV currently in use (see 'HFOV_OPTIONS' enum)
 */
int powermap_getDispFOV(void* const hPm);

/*
 * Function: powermap_getAspectRatio
 * ---------------------------------
 * Returns the current visualisation display window aspect-ratio
 *
 * Input Arguments:
 *     hPm - powermap handle
 * Returns:
 *     aspect ratio currently in use (see 'ASPECT_RATIO_OPTIONS' enum)
 */
int powermap_getAspectRatio(void* const hPm);

/*
 * Function: powermap_getPowermapAvgCoeff
 * --------------------------------------
 * Returns the current activity-map averaging coefficient.
 *
 * Input Arguments:
 *     hPm - powermap handle
 * Returns:
 *     current averaging coefficient, 0..1
 */
float powermap_getPowermapAvgCoeff(void* const hPm);
    
/*
 * Function: powermap_getPmap
 * --------------------------
 * Returns the latest computed activity-map if it is ready. Otherwise it returns
 * 0, and you'll just have to wait a bit
 *
 * Input Arguments:
 *     hPm         - powermap handle
 * Output Arguments:
 *     grid_dirs   - & scanning grid directions, in DEGREES; nDirs x 1
 *     pmap        - & activity-map values; nDirs x 1
 *     nDirs       - & number of directions
 *     pmapWidth   - & activity-map width in pixels
 *     hfov        - & horizontal FOV used to generate activity-map
 *     aspectRatio - & aspect ratio used to generate activity-map
 * Returns:
 *     flag, if activity-map is ready, 1: it is, 0: it is NOT
 */
int powermap_getPmap(void* const hPm,
                     float** grid_dirs,
                     float** pmap,
                     int* nDirs,
                     int* pmapWidth,
                     int* hfov,
                     int* aspectRatio);


#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */

#endif /* __POWERMAP_H_INCLUDED__ */
