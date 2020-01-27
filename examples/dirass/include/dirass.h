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
 * Filename: dirass.h (include header)
 * -----------------------------------
 * A sound-field visualiser based on the directional re-assignment of beamformer
 * energy, utilising the DoA estimates extracted from spatially-localised
 * active-intensity (SLAI) vectors; which are centred around each of the
 * corresponding scanning grid directions [1].
 *
 * Dependencies:
 *     saf_utilities, saf_vbap, saf_sh
 * Author, date created:
 *     Leo McCormack, 21.02.2019
 *
 * [1] McCormack, L., Politis, A., and Pulkki, V. (2019). "Sharpening of angular
 *     spectra based on a directional re-assignment approach for ambisonic
 *     sound-field visualisation". IEEE International Conference on Acoustics,
 *     Speech and Signal Processing (ICASSP).
 */

#ifndef __DIRASS_H_INCLUDED__
#define __DIRASS_H_INCLUDED__

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */
    
/* ========================================================================== */
/*                             Presets + Constants                            */
/* ========================================================================== */
 
#define DIRASS_MAX_NUM_INPUT_CHANNELS ( 64 )
    
/*
 * Enum: INPUT_ORDERS
 * ------------------
 * Available analysis orders.
 *
 * Options:
 *     INPUT_ORDER_FIRST   - First-order analysis (4 channel input)
 *     INPUT_ORDER_SECOND  - Second-order analysis (9 channel input)
 *     INPUT_ORDER_THIRD   - Third-order analysis (16 channel input)
 *     INPUT_ORDER_FOURTH  - Fourth-order analysis (25 channel input)
 *     INPUT_ORDER_FIFTH   - Fifth-order analysis (36 channel input)
 *     INPUT_ORDER_SIXTH   - Sixth-order analysis (49 channel input)
 *     INPUT_ORDER_SEVENTH - Seventh-order analysis (64 channel input)
 */
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
 * Enum: UPSCALE_ORDERS
 * --------------------
 * Available upscaling orders.
 *
 * Options:
 *     UPSCALE_ORDER_FIRST       - First-order upscaling
 *     UPSCALE_ORDER_SECOND      - Second-order upscaling
 *     UPSCALE_ORDER_THIRD       - Third-order upscaling
 *     UPSCALE_ORDER_FOURTH      - Fourth-order upscaling
 *     UPSCALE_ORDER_FIFTH       - Fifth-order upscaling
 *     UPSCALE_ORDER_SIXTH       - Sixth-order upscaling
 *     UPSCALE_ORDER_SEVENTH     - Seventh-order upscaling
 *     UPSCALE_ORDER_EIGHTH      - Eighth-order upscaling
 *     UPSCALE_ORDER_NINTH       - Ninth-order upscaling
 *     UPSCALE_ORDER_TENTH       - Tenth-order upscaling
 *     UPSCALE_ORDER_ELEVENTH    - Eleventh-order upscaling
 *     UPSCALE_ORDER_TWELFTH     - Twelfth-order upscaling
 *     UPSCALE_ORDER_THIRTEENTH  - Thirteenth-order upscaling
 *     UPSCALE_ORDER_FOURTEENTH  - Fourteenth-order upscaling
 *     UPSCALE_ORDER_FIFTEENTH   - Fifteenth-order upscaling
 *     UPSCALE_ORDER_SIXTHTEENTH - Sixthteenth-order upscaling
 *     UPSCALE_ORDER_SEVENTEENTH - Seventeenth-order upscaling
 *     UPSCALE_ORDER_EIGHTEENTH  - Eighteenth-order upscaling
 *     UPSCALE_ORDER_NINETEENTH  - Ninteenth-order upscaling
 *     UPSCALE_ORDER_TWENTIETH   - Twentieth-order upscaling
 */
typedef enum _UPSCALE_ORDERS{
    UPSCALE_ORDER_FIRST = 1,
    UPSCALE_ORDER_SECOND,
    UPSCALE_ORDER_THIRD,
    UPSCALE_ORDER_FOURTH,
    UPSCALE_ORDER_FIFTH,
    UPSCALE_ORDER_SIXTH,
    UPSCALE_ORDER_SEVENTH,
    UPSCALE_ORDER_EIGHTH,
    UPSCALE_ORDER_NINTH,
    UPSCALE_ORDER_TENTH,
    UPSCALE_ORDER_ELEVENTH,
    UPSCALE_ORDER_TWELFTH,
    UPSCALE_ORDER_THIRTEENTH,
    UPSCALE_ORDER_FOURTEENTH,
    UPSCALE_ORDER_FIFTEENTH,
    UPSCALE_ORDER_SIXTHTEENTH,
    UPSCALE_ORDER_SEVENTEENTH,
    UPSCALE_ORDER_EIGHTEENTH,
    UPSCALE_ORDER_NINETEENTH,
    UPSCALE_ORDER_TWENTIETH
    
}UPSCALE_ORDERS;
   
/*
 * Enum: GRID_OPTIONS
 * ------------------
 * Available scanning grid options
 *
 * Options:
 *     T_DESIGN_3        - 6 points
 *     T_DESIGN_4        - 12 points
 *     T_DESIGN_6        - 24 points
 *     T_DESIGN_9        - 48 points
 *     T_DESIGN_13       - 94 points
 *     T_DESIGN_18       - 180 points
 *     GRID_GEOSPHERE_6  - 362 points
 *     T_DESIGN_30       - 480 points
 *     GRID_GEOSPHERE_8  - 642 points
 *     GRID_GEOSPHERE_9  - 812 points
 *     GRID_GEOSPHERE_10 - 1002 points
 *     GRID_GEOSPHERE_12 - 1442 points
 */
typedef enum _GRID_OPTIONS {
    T_DESIGN_3 = 1,
    T_DESIGN_4,
    T_DESIGN_6,
    T_DESIGN_9,
    T_DESIGN_13,
    T_DESIGN_18,
    GRID_GEOSPHERE_6,
    T_DESIGN_30,
    GRID_GEOSPHERE_8,
    GRID_GEOSPHERE_9,
    GRID_GEOSPHERE_10,
    GRID_GEOSPHERE_12
    
}GRID_OPTIONS;

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
 * Enum: BEAM_TYPES
 * ----------------
 * Available sector beamforming patterns
 *
 * Options:
 *     BEAM_TYPE_CARD      - Cardioid
 *     BEAM_TYPE_HYPERCARD - Hyper-cardioid
 *     BEAM_TYPE_MAX_EV    - Hyper-cardioid with max_rE weighting
 */
typedef enum _BEAM_TYPES {
    BEAM_TYPE_CARD = 1,
    BEAM_TYPE_HYPERCARD,
    BEAM_TYPE_MAX_EV
    
} BEAM_TYPES;

/*
 * Enum: REASS_MODES
 * -----------------
 * Available processing modes. More information can be found in [1]
 *
 * Options:
 *     REASS_MODE_OFF - Re-assignment is disabled. i.e. dirass generates a
 *                      standard (beamformer)energy-based map
 *     REASS_NEAREST  - Each sector beamformer energy is re-assigned to the
 *                      nearest interpolation grid point, based on the analysed
 *                      DoA
 *     REASS_UPSCALE  - Each sector beamformer is re-encoded into spherical
 *                      harmonics of a higher order. The map is then derived
 *                      from the upscaled SHs as normal.
 *
 * [1] McCormack, L., Politis, A., and Pulkki, V. (2019). "Sharpening of angular
 *     spectra based on a directional re-assignment approach for ambisonic
 *     sound-field visualisation". IEEE International Conference on Acoustics,
 *     Speech and Signal Processing (ICASSP).
 */
typedef enum _REASS_MODE {
    REASS_MODE_OFF = 1,
    REASS_NEAREST,
    REASS_UPSCALE
    
} REASS_MODES;
    
/*
 * Enum: HFOV_OPTIONS
 * ------------------
 * Available horizontal feild-of-view (FOV) options
 *
 * Options:
 *     HFOV_360 - 360 degrees
 *     HFOV_180 - 180 degrees
 *     HFOV_90  - 90 degrees
 *     HFOV_60  - 60 degrees
 */
typedef enum _HFOV_OPTIONS{
    HFOV_360 = 1,
    HFOV_180,
    HFOV_90,
    HFOV_60
    
}HFOV_OPTIONS;
    
/*
 * Enum: ASPECT_RATIO_OPTIONS
 * --------------------------
 * Available aspect ratios
 *
 * Options:
 *     ASPECT_RATIO_2_1  - 2:1
 *     ASPECT_RATIO_16_9 - 16:9
 *     ASPECT_RATIO_4_3  - 4:3
 */
typedef enum _ASPECT_RATIO_OPTIONS{
    ASPECT_RATIO_2_1 = 1,
    ASPECT_RATIO_16_9,
    ASPECT_RATIO_4_3

}ASPECT_RATIO_OPTIONS;
    
/*
 * Enum: CODEC_STATUS
 * ------------------
 * Current status of the codec.
 *
 * Options:
 *     CODEC_STATUS_INITIALISED     - Codec is initialised and ready to process
 *                                    input audio.
 *     CODEC_STATUS_NOT_INITIALISED - Codec has not yet been initialised, or
 *                                    the codec configuration has changed. Input
 *                                    audio should not be processed.
 *     CODEC_STATUS_INITIALISING    - Codec is currently being initialised,
 *                                    input audio should not be processed.
 */
typedef enum _CODEC_STATUS{
    CODEC_STATUS_INITIALISED = 0,
    CODEC_STATUS_NOT_INITIALISED,
    CODEC_STATUS_INITIALISING
}CODEC_STATUS;

#define DIRASS_PROGRESSBARTEXT_CHAR_LENGTH 256
    

/* ========================================================================== */
/*                               Main Functions                               */
/* ========================================================================== */

/*
 * Function: dirass_create
 * -----------------------
 * Creates an instance of the dirass
 *
 * Input Arguments:
 *     phDir - & address of dirass handle
 */
void dirass_create(void** const phDir);

/*
 * Function: dirass_destroy
 * ------------------------
 * Destroys an instance of the dirass
 *
 * Input Arguments:
 *     phDir - & address of dirass handle
 */
void dirass_destroy(void** const phDir);

/*
 * Function: dirass_init
 * ---------------------
 * Initialises an instance of dirass with default settings
 *
 * Input Arguments:
 *     hDir       - dirass handle
 *     samplerate - host samplerate.
 */
void dirass_init(void* const hDir,
                 float  samplerate);
    
/*
 * Function: dirass_initCodec
 * --------------------------
 * Intialises the codec variables, based on current global/user parameters
 *
 * Input Arguments:
 *     hDir - dirass handle
 */
void dirass_initCodec(void* const hDir);

/*
 * Function: dirass_process
 * ------------------------
 * Analyses the input spherical harmonic signals to generate an activity-map as
 * in [1]
 *
 * Input Arguments:
 *     hDir      - dirass handle
 *     inputs    - input channel buffers; 2-D array: nInputs x nSamples
 *     nInputs   - number of input channels
 *     nSamples  - number of samples in 'inputs'/'output' matrices
 *     isPlaying - flag to say if there is audio in the input buffers, 0: no
 *                 audio, reduced processing, 1: audio, full processing
 *
 * [1] McCormack, L., Politis, A., and Pulkki, V. (2019). "Sharpening of angular
 *     spectra based on a directional re-assignment approach for ambisonic
 *     sound-field visualisation". IEEE International Conference on Acoustics,
 *     Speech and Signal Processing (ICASSP).
 */
void dirass_analysis(void* const hDir,
                     float** const inputs,
                     int nInputs,
                     int nSamples,
                     int isPlaying);
    
   
/* ========================================================================== */
/*                                Set Functions                               */
/* ========================================================================== */

/*
 * Function: dirass_refreshSettings
 * --------------------------------
 * Sets all intialisation flags to 1. i.e. re-initialise all settings/variables
 * as dirass is currently configured, at next available opportunity.
 *
 * Input Arguments:
 *     hDir - dirass handle
 */
void dirass_refreshSettings(void* const hDir);
    
/*
 * Function: dirass_setBeamType
 * ----------------------------
 * Sets the sector beamforming pattern to employ for the analysis
 *
 * Input Arguments:
 *     hDir    - dirass handle
 *     newType - new type. See "BEAM_TYPE" enum.
 */
void dirass_setBeamType(void* const hDir, int newType);
    
/*
 * Function: dirass_setInputOrder
 * ------------------------------
 * Sets the input/analysis order
 *
 * Input Arguments:
 *     hDir     - dirass handle
 *     newValue - new analysis order (see 'INPUT_ORDERS' enum)
 */
void dirass_setInputOrder(void* const hDir,  int newValue);
    
/*
 * Function: dirass_setDisplayGridOption
 * -------------------------------------
 * Sets a new display grid option. see "GRID_OPTIONS" enum
 * Note: not safe to call while simultaneously calling "dirass_analysis"!
 *
 * Input Arguments:
 *     hDir      - dirass handle
 *     newOption - new grid option (see 'GRID_OPTIONS' enum)
 */
void dirass_setDisplayGridOption(void* const hDir,  int newOption);
   
/*
 * Function: dirass_setDispWidth
 * -----------------------------
 * Sets the output display width in pixels
 * Note: not safe to call while simultaneously calling "dirass_analysis"!
 *
 * Input Arguments:
 *     hDir     - dirass handle
 *     newValue - new number of pixels
 */
void dirass_setDispWidth(void* const hDir,  int newValue);
    
/*
 * Function: dirass_setUpscaleOrder
 * --------------------------------
 * Sets the upscale order, if REASS_MODE is set to REASS_UPSCALE
 *
 * Input Arguments:
 *     hDir     - dirass handle
 *     newOrder - new order. See "UPSCALE_ORDER" enum.
 */
void dirass_setUpscaleOrder(void* const hDir,  int newOrder);
    
/*
 * Function: dirass_setDiRAssMode
 * ------------------------------
 * Sets the analysis directional re-assignment mode. See "REASS_MODE" enum.
 *
 * Input Arguments:
 *     hDir    - dirass handle
 *     newMode - new mode. See "REASS_MODE" enum.
 */
void dirass_setDiRAssMode(void* const hDir,  int newMode);
    
/*
 * Function: dirass_setMinFreq
 * ---------------------------
 * Sets the minimum analysis frequency
 *
 * Input Arguments:
 *     hDir     - dirass handle
 *     newValue - new minimum analysis frequency, in Hz
 */
void dirass_setMinFreq(void* const hDir,  float newValue);

/*
 * Function: dirass_setMaxFreq
 * ---------------------------
 * Sets the maximum analysis frequency
 *
 * Input Arguments:
 *     hDir     - dirass handle
 *     newValue - new maximum analysis frequency, in Hz
 */
void dirass_setMaxFreq(void* const hDir,  float newValue);

/*
 * Function: dirass_setChOrder
 * ---------------------------
 * Sets the Ambisonic channel ordering convention to decode with, in order to
 * match the convention employed by the input signals
 *
 * Input Arguments:
 *     hDir     - dirass handle
 *     newOrder - convention to use (see 'CH_ORDER' enum)
 */
void dirass_setChOrder(void* const hDir, int newOrder);

/*
 * Function: dirass_setNormType
 * ----------------------------
 * Sets the Ambisonic normalisation convention to decode with, in order to match
 * with the convention employed by the input signals.
 *
 * Input Arguments:
 *     hDir    - dirass handle
 *     newType - convention to use (see 'NORM_TYPE' enum)
 */
void dirass_setNormType(void* const hDir, int newType);

/*
 * Function: dirass_setDispFOV
 * ---------------------------
 * Sets the visualisation display window horizontal field-of-view (FOV).
 *
 * Input Arguments:
 *     hDir      - dirass handle
 *     newOption - horizontal FOV to use (see 'HFOV_OPTIONS' enum)
 */
void dirass_setDispFOV(void* const hDir, int newOption);
    
/*
 * Function: dirass_setAspectRatio
 * -------------------------------
 * Sets the visualisation display window aspect-ratio
 *
 * Input Arguments:
 *     hDir      - dirass handle
 *     newOption - aspect ratio to use (see 'ASPECT_RATIO_OPTIONS' enum)
 */
void dirass_setAspectRatio(void* const hDir, int newOption);
    
/*
 * Function: dirass_setMapAvgCoeff
 * -------------------------------
 * Sets the activity-map averaging coefficient.
 *
 * Input Arguments:
 *     hDir     - dirass handle
 *     newValue - new averaging coefficient, 0..1
 */
void dirass_setMapAvgCoeff(void* const hDir, float newValue);
    
/*
 * Function: dirass_requestPmapUpdate
 * ----------------------------------
 * Informs dirass that it should compute a new activity-map
 *
 * Input Arguments:
 *     hDir - dirass handle
 */
void dirass_requestPmapUpdate(void* const hDir);
    

/* ========================================================================== */
/*                                Get Functions                               */
/* ========================================================================== */

/*
 * Function: dirass_getCodecStatus
 * -------------------------------
 * Returns current codec status.
 *
 * Input Arguments:
 *     hDir - dirass handle
 * Returns:
 *     codec status (see 'CODEC_STATUS' enum)
 */
CODEC_STATUS dirass_getCodecStatus(void* const hDir);

/*
 * Function: dirass_getProgressBar0_1
 * ----------------------------------
 * (Optional) Returns current intialisation/processing progress, between 0..1
 * 0: intialisation/processing has started
 * 1: intialisation/processing has ended
 *
 * Input Arguments:
 *     hDir - dirass handle
 * Returns:
 *     current progress, 0..1
 */
float dirass_getProgressBar0_1(void* const hDir);

/*
 * Function: dirass_getProgressBarText
 * -----------------------------------
 * (Optional) Returns current intialisation/processing progress text
 * Note: "text" string should be (at least) of length:
 *     DIRASS_PROGRESSBARTEXT_CHAR_LENGTH
 *
 * Input Arguments:
 *     hDir - dirass handle
 * Output Arguments:
 *     text - process bar text; DIRASS_PROGRESSBARTEXT_CHAR_LENGTH x 1
 */
void dirass_getProgressBarText(void* const hDir, char* text);

/*
 * Function: dirass_getInputOrder
 * ------------------------------
 * Returns the current analysis/input order.
 *
 * Input Arguments:
 *     hDir - dirass handle
 * Returns:
 *     current maximum analysis order (see 'INPUT_ORDERS' enum)
 */
int dirass_getInputOrder(void* const hDir);
    
/*
 * Function: dirass_getBeamType
 * ----------------------------
 * Returns the sector beamforming pattern to employed for the analysis
 *
 * Input Arguments:
 *     hDir - dirass handle
 * Returns:
 *     current type. See "BEAM_TYPE" enum.
 */
int dirass_getBeamType(void* const hDir);

/*
 * Function: dirass_getDisplayGridOption
 * -------------------------------------
 * Returns the current display grid option. see "GRID_OPTIONS" enum
 *
 * Input Arguments:
 *     hDir - dirass handle
 * Returns:
 *     current grid option (see 'GRID_OPTIONS' enum)
 */
int dirass_getDisplayGridOption(void* const hDir);

/*
 * Function: dirass_getDispWidth
 * -----------------------------
 * Returns the current output display width in pixels
 *
 * Input Arguments:
 *     hDir - dirass handle
 * Returns:
 *     current number of horizontal display pixels
 */
int dirass_getDispWidth(void* const hDir);

/*
 * Function: dirass_getUpscaleOrder
 * --------------------------------
 * Returns the current upscale order
 *
 * Input Arguments:
 *     hDir - dirass handle
 * Returns:
 *     current order. See "UPSCALE_ORDER" enum.
 */
int dirass_getUpscaleOrder(void* const hDir);

/*
 * Function: dirass_getDiRAssMode
 * ------------------------------
 * Returns the current analysis directional re-assignment mode. See "REASS_MODE"
 * enum.
 *
 * Input Arguments:
 *     hDir - dirass handle
 * Returns:
 *     current mode. See "REASS_MODE" enum.
 */
int dirass_getDiRAssMode(void* const hDir); 

/*
 * Function: dirass_getMinFreq
 * ---------------------------
 * Returns the current minimum analysis frequency
 *
 * Input Arguments:
 *     hDir - dirass handle
 * Returns:
 *     current minimum analysis frequency, in Hz
 */
float dirass_getMinFreq(void* const hDir);

/*
 * Function: dirass_getMaxFreq
 * ---------------------------
 * Returns the current maximum analysis frequency
 *
 * Input Arguments:
 *     hDir - dirass handle
 *  Returns:
 *     current maximum analysis frequency, in Hz
 */
float dirass_getMaxFreq(void* const hDir);

/*
 * Function: dirass_getSamplingRate
 * --------------------------------
 * Returns the current sampling rate
 *
 * Input Arguments:
 *     hDir - dirass handle
 *  Returns:
 *     current sampling rate, in Hz
 */
int dirass_getSamplingRate(void* const hDir); 
    
/*
 * Function: dirass_getNSHrequired
 * -------------------------------
 * Returns the number of spherical harmonic signals required by the current
 * analysis order i.e. (current_order + 1)^2
 *
 * Input Arguments:
 *     hDir - dirass handle
 * Returns:
 *     number of required spherical harmonic signals required by current
 *     analysis order
 */
int dirass_getNSHrequired(void* const hDir);

/*
 * Function: dirass_getChOrder
 * ---------------------------
 * Returns the Ambisonic channel ordering convention currently being used to
 * decode with, which should match the convention employed by the input signals
 *
 * Input Arguments:
 *     hDir - dirass handle
 * Returns:
 *     convention currently being used (see 'CH_ORDER' enum)
 */
int dirass_getChOrder(void* const hDir);

/*
 * Function: dirass_getNormType
 * ----------------------------
 * Returns the Ambisonic normalisation convention currently being usedto decode
 * with, which should match the convention employed by the input signals.
 *
 * Input Arguments:
 *     hDir - dirass handle
 * Returns:
 *     convention currently being used (see 'NORM_TYPE' enum)
 */
int dirass_getNormType(void* const hDir);

/*
 * Function: dirass_getDispFOV
 * ---------------------------
 * Returns the current visualisation display window horizontal field-of-view
 * (FOV).
 *
 * Input Arguments:
 *     hDir - dirass handle
 * Returns:
 *     horizontal FOV currently in use (see 'HFOV_OPTIONS' enum)
 */
int dirass_getDispFOV(void* const hDir);

/*
 * Function: dirass_getAspectRatio
 * -------------------------------
 * Returns the current visualisation display window aspect-ratio
 *
 * Input Arguments:
 *     hDir - dirass handle
 * Returns:
 *     aspect ratio currently in use (see 'ASPECT_RATIO_OPTIONS' enum)
 */
int dirass_getAspectRatio(void* const hDir);

/*
 * Function: dirass_getMapAvgCoeff
 * -------------------------------
 * Returns the current activity-map averaging coefficient.
 *
 * Input Arguments:
 *     hDir - dirass handle
 * Returns:
 *     current averaging coefficient, 0..1
 */
float dirass_getMapAvgCoeff(void* const hDir);
    
/*
 * Function: dirass_getPmap
 * ------------------------
 * Returns the latest computed activity-map if it is ready. Otherwise it returns
 * 0, and you'll just have to wait a bit  
 *
 * Input Arguments:
 *     hDir        - dirass handle
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
int dirass_getPmap(void* const hDir,
                   float** grid_dirs,
                   float** pmap,
                   int* nDirs,
                   int* pmapWidth,
                   int* hfov,
                   float* aspectRatio);


#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */

#endif /* __DIRASS_H_INCLUDED__ */
