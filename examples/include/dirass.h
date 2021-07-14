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

/**
 * @example dirass.h
 * @brief A sound-field visualiser based on the directional re-assignment of
 *        beamformer energy based on local DoA estimates
 *
 * ### Files
 * dirass.h (include), dirass_internal.h, dirass.c, dirass_internal.c
 * ### Include Header
 */

/**
 * @file dirass.h
 * @brief A sound-field visualiser based on the directional re-assignment of
 *        beamformer energy based on local DoA estimates [1,2]
 *
 * @see [1] McCormack, L., Delikaris-Manias, S., Politis, A., Pavlidi, D.,
 *          Farina, A., Pinardi, D. and Pulkki, V., 2019. Applications of
 *          Spatially Localized Active-Intensity Vectors for Sound-Field
 *          Visualization. Journal of the Audio Engineering Society, 67(11),
 *          pp.840-854.
 * @see [2] McCormack, L., Politis, A., and Pulkki, V. (2019). "Sharpening of
 *          angular spectra based on a directional re-assignment approach for
 *          ambisonic sound-field visualisation". IEEE International Conference
 *          on Acoustics, Speech and Signal Processing (ICASSP).
 *
 * @author Leo McCormack
 * @date 21.02.2019
 * @license ISC
 */

#ifndef __DIRASS_H_INCLUDED__
#define __DIRASS_H_INCLUDED__

#include "_common.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/* ========================================================================== */
/*                             Presets + Constants                            */
/* ========================================================================== */

/** Available upscaling orders */
typedef enum {
    UPSCALE_ORDER_FIRST = 1,   /**< First-order upscaling */
    UPSCALE_ORDER_SECOND,      /**< Second-order upscaling */
    UPSCALE_ORDER_THIRD,       /**< Third-order upscaling */
    UPSCALE_ORDER_FOURTH,      /**< Fourth-order upscaling */
    UPSCALE_ORDER_FIFTH,       /**< Fifth-order upscaling */
    UPSCALE_ORDER_SIXTH,       /**< Sixth-order upscaling */
    UPSCALE_ORDER_SEVENTH,     /**< Seventh-order upscaling */
    UPSCALE_ORDER_EIGHTH,      /**< Eighth-order upscaling */
    UPSCALE_ORDER_NINTH,       /**< Ninth-order upscaling */
    UPSCALE_ORDER_TENTH,       /**< Tenth-order upscaling */
    UPSCALE_ORDER_ELEVENTH,    /**< Eleventh-order upscaling */
    UPSCALE_ORDER_TWELFTH,     /**< Twelfth-order upscaling */
    UPSCALE_ORDER_THIRTEENTH,  /**< Thirteenth-order upscaling */
    UPSCALE_ORDER_FOURTEENTH,  /**< Fourteenth-order upscaling */
    UPSCALE_ORDER_FIFTEENTH,   /**< Fifteenth-order upscaling */
    UPSCALE_ORDER_SIXTHTEENTH, /**< Sixthteenth-order upscaling */
    UPSCALE_ORDER_SEVENTEENTH, /**< Seventeenth-order upscaling */
    UPSCALE_ORDER_EIGHTEENTH,  /**< Eighteenth-order upscaling */
    UPSCALE_ORDER_NINETEENTH,  /**< Ninteenth-order upscaling */
    UPSCALE_ORDER_TWENTIETH    /**< Twentieth-order upscaling */
    
}DIRASS_UPSCALE_ORDERS;
   
/** Available scanning grid options */
typedef enum {
    T_DESIGN_3 = 1,    /**< T_DESIGN_3        - 6 points */
    T_DESIGN_4,        /**< T_DESIGN_4        - 12 points */
    T_DESIGN_6,        /**< T_DESIGN_6        - 24 points */
    T_DESIGN_9,        /**< T_DESIGN_9        - 48 points */
    T_DESIGN_13,       /**< T_DESIGN_13       - 94 points */
    T_DESIGN_18,       /**< T_DESIGN_18       - 180 points */
    GRID_GEOSPHERE_6,  /**< GRID_GEOSPHERE_6  - 362 points */
    T_DESIGN_30,       /**< T_DESIGN_30       - 480 points */
    GRID_GEOSPHERE_8,  /**< GRID_GEOSPHERE_8  - 642 points */
    GRID_GEOSPHERE_9,  /**< GRID_GEOSPHERE_9  - 812 points */
    GRID_GEOSPHERE_10, /**< GRID_GEOSPHERE_10 - 1002 points */
    GRID_GEOSPHERE_12  /**< GRID_GEOSPHERE_12 - 1442 points */
    
}DIRASS_GRID_OPTIONS;

/**
 * Available processing modes. More information can be found in [1]
 *
 * @see [1] McCormack, L., Politis, A., and Pulkki, V. (2019). "Sharpening of
 *          angular spectra based on a directional re-assignment approach for
 *          ambisonic sound-field visualisation". IEEE International Conference
 *          on Acoustics, Speech and Signal Processing (ICASSP).
 */
typedef enum {
    REASS_MODE_OFF = 1, /**< Re-assignment is disabled. i.e. dirass generates a
                         *   standard (beamformer)energy-based map */
    REASS_NEAREST,      /**< Each sector beamformer energy is re-assigned to the
                         *   nearest interpolation grid point, based on the
                         *   analysed DoA */
    REASS_UPSCALE       /**< Each sector beamformer is re-encoded into spherical
                         *   harmonics of a higher order. The map is then
                         *   derived from the upscaled SHs as normal. */
    
} DIRASS_REASS_MODES;


/* ========================================================================== */
/*                               Main Functions                               */
/* ========================================================================== */

/**
 * Creates an instance of the dirass
 *
 * @param[in] phDir (&) address of dirass handle
 */
void dirass_create(void** const phDir);

/**
 * Destroys an instance of the dirass
 *
 * @param[in] phDir (&) address of dirass handle
 */
void dirass_destroy(void** const phDir);

/**
 * Initialises an instance of dirass with default settings
 *
 * @warning This should not be called while _process() is on-going!
 *
 * @param[in] hDir       dirass handle
 * @param[in] samplerate Host samplerate.
 */
void dirass_init(void* const hDir,
                 float  samplerate);
    
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
 * @param[in] hDir dirass handle
 */
void dirass_initCodec(void* const hDir);

/**
 * Analyses the input spherical harmonic signals to generate an activity-map as
 * in [1,2]
 *
 * @param[in] hDir      dirass handle
 * @param[in] inputs    Input channel buffers; 2-D array: nInputs x nSamples
 * @param[in] nInputs   Number of input channels
 * @param[in] nSamples  Number of samples in 'inputs'/'output' matrices
 * @param[in] isPlaying Flag to indicate if there is audio in the input buffers,
 *                      0: no audio, reduced processing, 1: audio, full
 *                      processing
 *
 * @see [1] McCormack, L., Politis, A., and Pulkki, V. (2019). "Sharpening of
 *          angular spectra based on a directional re-assignment approach for
 *          ambisonic sound-field visualisation". IEEE International Conference
 *          on Acoustics, Speech and Signal Processing (ICASSP).
 * @see [2] McCormack, L., Delikaris-Manias, S., Politis, A., Pavlidi, D.,
 *          Farina, A., Pinardi, D. and Pulkki, V., 2019. Applications of
 *          Spatially Localized Active-Intensity Vectors for Sound-Field
 *          Visualization. Journal of the Audio Engineering Society, 67(11),
 *          pp.840-854.
 */
void dirass_analysis(void* const hDir,
                     const float *const * inputs,
                     int nInputs,
                     int nSamples,
                     int isPlaying);
    
   
/* ========================================================================== */
/*                                Set Functions                               */
/* ========================================================================== */

/**
 * Sets all intialisation flags to 1; re-initialising all settings/variables
 * as dirass is currently configured, at next available opportunity.
 */
void dirass_refreshSettings(void* const hDir);
    
/**
 * Sets the sector beamforming pattern to employ for the analysis (see
 * #STATIC_BEAM_TYPES enum).
 */
void dirass_setBeamType(void* const hDir, int newType);
    
/** Sets the input/analysis order (see #SH_ORDERS enum) */
void dirass_setInputOrder(void* const hDir,  int newValue);
    
/**
 * Sets a new display grid option (see #DIRASS_GRID_OPTIONS enum)
 *
 * @warning Not safe to call while simultaneously calling dirass_analysis()!
 */
void dirass_setDisplayGridOption(void* const hDir,  int newOption);
   
/**
 * Sets the output display width in pixels
 *
 * @warning Not safe to call while simultaneously calling dirass_analysis()!
 */
void dirass_setDispWidth(void* const hDir,  int newValue);
    
/**
 * Sets the upscale order, only if #DIRASS_REASS_MODES is set to #REASS_UPSCALE,
 * (see #DIRASS_UPSCALE_ORDERS enum).
 */
void dirass_setUpscaleOrder(void* const hDir,  int newOrder);
    
/**
 * Sets the analysis directional re-assignment mode (see #DIRASS_REASS_MODES
 * enum)
 */
void dirass_setDiRAssMode(void* const hDir,  int newMode);
    
/** Sets the minimum analysis frequency, in Hz */
void dirass_setMinFreq(void* const hDir,  float newValue);

/** Sets the maximum analysis frequency, in Hz */
void dirass_setMaxFreq(void* const hDir,  float newValue);

/**
 * Sets the Ambisonic channel ordering convention to decode with, in order to
 * match the convention employed by the input signals (see #CH_ORDER enum)
 */
void dirass_setChOrder(void* const hDir, int newOrder);

/**
 * Sets the Ambisonic normalisation convention to decode with, in order to match
 * with the convention employed by the input signals (see #NORM_TYPES enum)
 */
void dirass_setNormType(void* const hDir, int newType);

/**
 * Sets the visualisation display window horizontal field-of-view (FOV)
 * (see #HFOV_OPTIONS enum)
 */
void dirass_setDispFOV(void* const hDir, int newOption);
    
/**
 * Sets the visualisation display window aspect-ratio (see
 * #ASPECT_RATIO_OPTIONS enum)
 */
void dirass_setAspectRatio(void* const hDir, int newOption);
    
/** Sets the activity-map averaging coefficient, 0..1 */
void dirass_setMapAvgCoeff(void* const hDir, float newValue);
    
/** Informs dirass that it should compute a new activity-map */
void dirass_requestPmapUpdate(void* const hDir);
    

/* ========================================================================== */
/*                                Get Functions                               */
/* ========================================================================== */

/**
 * Returns the processing framesize (i.e., number of samples processed with
 * every _process() call )
 */
int dirass_getFrameSize(void);

/** Returns current codec status (see #CODEC_STATUS enum) */
CODEC_STATUS dirass_getCodecStatus(void* const hDir);

/**
 * (Optional) Returns current intialisation/processing progress, between 0..1
 * - 0: intialisation/processing has started
 * - 1: intialisation/processing has ended
 */
float dirass_getProgressBar0_1(void* const hDir);

/**
 * (Optional) Returns current intialisation/processing progress text
 *
 * @note "text" string should be (at least) of length:
 *       #PROGRESSBARTEXT_CHAR_LENGTH
 *
 * @param[in]  hDir dirass handle
 * @param[out] text Process bar text; #PROGRESSBARTEXT_CHAR_LENGTH x 1
 */
void dirass_getProgressBarText(void* const hDir, char* text);

/**
 * Returns the current analysis/input order (see #SH_ORDERS enum)
 */
int dirass_getInputOrder(void* const hDir);
    
/**
 * Returns the sector beamforming pattern to employed for the analysis (see
 * #STATIC_BEAM_TYPES enum)
 */
int dirass_getBeamType(void* const hDir);

/** Returns the current display grid option (see #DIRASS_GRID_OPTIONS enum) */
int dirass_getDisplayGridOption(void* const hDir);

/** Returns the current output display width in pixels */
int dirass_getDispWidth(void* const hDir);

/** Returns the current upscale order (see #DIRASS_UPSCALE_ORDERS enum) */
int dirass_getUpscaleOrder(void* const hDir);

/**
 * Returns the current analysis directional re-assignment mode (see
 * #DIRASS_REASS_MODES enum)
 */
int dirass_getDiRAssMode(void* const hDir); 

/** Returns the current minimum analysis frequency, in Hz */
float dirass_getMinFreq(void* const hDir);

/** Returns the current maximum analysis frequency, in Hz */
float dirass_getMaxFreq(void* const hDir);

/** Returns the current sampling rate, in Hz */
int dirass_getSamplingRate(void* const hDir); 
    
/**
 * Returns the number of spherical harmonic signals required by the current
 * analysis order: (current_order + 1)^2
 */
int dirass_getNSHrequired(void* const hDir);

/**
 * Returns the Ambisonic channel ordering convention currently being used to
 * decode with, which should match the convention employed by the input signals
 * (see #CH_ORDER enum)
 */
int dirass_getChOrder(void* const hDir);

/**
 * Returns the Ambisonic normalisation convention currently being usedto decode
 * with, which should match the convention employed by the input signals (see
 * #NORM_TYPES enum)
 */
int dirass_getNormType(void* const hDir);

/**
 * Returns the current visualisation display window horizontal field-of-view
 * (FOV) (see #HFOV_OPTIONS enum)
 */
int dirass_getDispFOV(void* const hDir);

/**
 * Returns the current visualisation display window aspect-ratio (see
 * #ASPECT_RATIO_OPTIONS enum)
 */
int dirass_getAspectRatio(void* const hDir);

/** Returns the current activity-map averaging coefficient, 0..1 */
float dirass_getMapAvgCoeff(void* const hDir);
    
/**
 * Returns the latest computed activity-map if it is ready; otherwise it returns
 * 0, and you'll just have to wait a bit  
 *
 * @param[in]  hDir        (&) dirass handle
 * @param[out] grid_dirs   (&) scanning grid directions, in DEGREES; nDirs x 1
 * @param[out] pmap        (&) activity-map values; nDirs x 1
 * @param[out] nDirs       (&) number of directions
 * @param[out] pmapWidth   (&) activity-map width in pixels
 * @param[out] hfov        (&) horizontal FOV used to generate activity-map
 * @param[out] aspectRatio (&) aspect ratio used to generate activity-map
 * @returns                flag, if activity-map is ready, 1: it is, 0: it is
 *                         NOT
 */
int dirass_getPmap(void* const hDir,
                   float** grid_dirs,
                   float** pmap,
                   int* nDirs,
                   int* pmapWidth,
                   int* hfov,
                   float* aspectRatio);

/**
 * Returns the processing delay in samples (may be used for delay compensation
 * features)
 */
int dirass_getProcessingDelay(void);


#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */

#endif /* __DIRASS_H_INCLUDED__ */
