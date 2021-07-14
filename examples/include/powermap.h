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

/**
 * @example powermap.h
 * @brief A sound-field visualiser, which utilises spherical harmonic signals as
 *        input
 *
 * ### Files
 * powermap.h (include), powermap_internal.h, powermap.c, powermap_internal.c
 * ### Include Header
 */

/**
 * @file powermap.h
 * @brief A sound-field visualiser, which utilises spherical harmonic signals as
 *        input; note this code is a remnant from the work conducted in [1]
 *
 * @see [1] McCormack, L., Delikaris-Manias, S. and Pulkki, V., 2017. Parametric
 *          acoustic camera for real-time sound capture, analysis and tracking.
 *          In Proceedings of the 20th International Conference on Digital Audio
 *          Effects (DAFx-17) (pp. 412-419)
 *
 * @author Leo McCormack
 * @date 26.04.2016
 * @license ISC
 */

#ifndef __POWERMAP_H_INCLUDED__
#define __POWERMAP_H_INCLUDED__

#include "_common.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */
    
/* ========================================================================== */
/*                             Presets + Constants                            */
/* ========================================================================== */

/**
 * Available power-map/activity-map options
 */
typedef enum {
    PM_MODE_PWD = 1,     /**< Activity-map based on the energy of hyper-cardioid
                          *   [plane-wave decomposition (PWD)] beamformers*/
    PM_MODE_MVDR,        /**< Activity-map based on the energy of minimum-
                          *   variance distortionless response (MVDR)
                          *   beamformers */
    PM_MODE_CROPAC_LCMV, /**< Experimental! activity-map based on a linearly-
                          *   contrained minimum-variance (LCMV) formulation of
                          *   the Cross-Pattern Coherence (CroPaC) spatial
                          *   filter */
    PM_MODE_MUSIC,       /**< Activity-map based on the sub-space method:
                          *   multiple signal classification (MUSIC) */
    PM_MODE_MUSIC_LOG,   /**< Same as #PM_MODE_MUSIC, but log(out_values) */
    PM_MODE_MINNORM,     /**< Activity-map based on the sub-space method:
                          *   minimum-norm (Min-Norm) */
    PM_MODE_MINNORM_LOG  /**< Same as #PM_MODE_MINNORM, but log(out_values) */
    
} POWERMAP_MODES;


/* ========================================================================== */
/*                               Main Functions                               */
/* ========================================================================== */

/**
 * Creates an instance of the powermap
 *
 * @param[in] phPm (&) address of powermap handle
 */
void powermap_create(void** const phPm);

/**
 * Destroys an instance of the powermap
 *
 * @param[in] phPm (&) address of powermap handle
 */
void powermap_destroy(void** const phPm);

/**
 * Initialises an instance of powermap with default settings
 *
 * @warning This should not be called while _process() is on-going!
 *
 * @param[in] hPm        powermap handle
 * @param[in] samplerate Host samplerate.
 */
void powermap_init(void* const hPm,
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
 * @param[in] hPm powermap handle
 */
void powermap_initCodec(void* const hPm);

/**
 * Analyses the input spherical harmonic signals to generate an activity-map
 *
 * @param[in] hPm       powermap handle
 * @param[in] inputs    Input channel buffers; 2-D array: nInputs x nSamples
 * @param[in] nInputs   Number of input channels
 * @param[in] nSamples  Number of samples in 'inputs'/'output' matrices
 * @param[in] isPlaying flag to say if there is audio in the input buffers, 0:
 *                      no audio, reduced processing, 1: audio, full processing
 */
void powermap_analysis(void* const hPm,
                       const float *const * inputs,
                       int nInputs,
                       int nSamples,
                       int isPlaying);


/* ========================================================================== */
/*                                Set Functions                               */
/* ========================================================================== */

/**
 * Sets all intialisation flags to 1; re-initialising all settings/variables
 * as powermap is currently configured, at next available opportunity.
 */
void powermap_refreshSettings(void* const hPm);

/** Sets the powermap/activity-map approach, (see #POWERMAP_MODES enum) */
void powermap_setPowermapMode(void* const hPm, int newMode);

/** Sets the maximum input/analysis order (see #SH_ORDERS enum) */
void powermap_setMasterOrder(void* const hPm,  int newValue);

/** Sets the input/analysis order for one specific frequency band index */
void powermap_setAnaOrder(void* const hPm,  int newValue, int bandIdx);

/** Sets the input/analysis order for all frequency bands */
void powermap_setAnaOrderAllBands(void* const hPm, int newValue);
    
/**
 * Sets the weighting coefficient for a particular frequency band, allowing
 * one to "equalise" the activity-map.
 */
void powermap_setPowermapEQ(void* const hPm,  float newValue, int bandIdx);

/** Sets the weighting coefficient for all frequency bands */
void powermap_setPowermapEQAllBands(void* const hPm,  float newValue);
    
/** Sets the covariance matrix averaging coefficient, 0..1 */
void powermap_setCovAvgCoeff(void* const hPm, float newAvg);

/**
 * Sets the Ambisonic channel ordering convention to decode with, in order to
 * match the convention employed by the input signals (see #CH_ORDER enum)
 */
void powermap_setChOrder(void* const hPm, int newOrder);

/**
 * Sets the Ambisonic normalisation convention to decode with, in order to match
 * with the convention employed by the input signals (see #NORM_TYPES enum)
 */
void powermap_setNormType(void* const hPm, int newType);

/**
 * Sets an input preset; the microphone/hyrophone array used to capture
 * the input signals, (see #MIC_PRESETS enum)
 */
void powermap_setSourcePreset(void* const hPm, int newPresetID);

/** Sets the number of sources present in the input sound scene */
void powermap_setNumSources(void* const hPm, int newValue);

/**
 * Sets the visualisation display window horizontal field-of-view (FOV)
 * (see #HFOV_OPTIONS enum)
 */
void powermap_setDispFOV(void* const hPm, int newOption);

/**
 * Sets the visualisation display window aspect-ratio (see
 * #ASPECT_RATIO_OPTIONS enum)
 */
void powermap_setAspectRatio(void* const hPm, int newOption);

/** Sets the activity-map averaging coefficient, 0..1 */
void powermap_setPowermapAvgCoeff(void* const hPm, float newValue);

/**
 * Informs powermap that it should compute a new activity-map at its own
 * convenience, if it would be so kind; thank you, God bless.
 */
void powermap_requestPmapUpdate(void* const hPm);


/* ========================================================================== */
/*                                Get Functions                               */
/* ========================================================================== */

/**
 * Returns the processing framesize (i.e., number of samples processed with
 * every _process() call )
 */
int powermap_getFrameSize(void);

/** Returns current codec status (see #CODEC_STATUS enum) */
CODEC_STATUS powermap_getCodecStatus(void* const hPm);

/**
 * (Optional) Returns current intialisation/processing progress, between 0..1
 *  - 0: intialisation/processing has started
 *  - 1: intialisation/processing has ended
 */
float powermap_getProgressBar0_1(void* const hPm);

/**
 * (Optional) Returns current intialisation/processing progress text
 * @note "text" string should be (at least) of length:
 *        #PROGRESSBARTEXT_CHAR_LENGTH
 */
void powermap_getProgressBarText(void* const hPm, char* text);

/**
 * Returns the current maximum analysis/input order (see #SH_ORDERS enum)
 */
int powermap_getMasterOrder(void* const hPm);

/**
 * Returns the powermap/activity-map mode to employed for the analysis
 * see #POWERMAP_MODES enum.
 */
int powermap_getPowermapMode(void* const hPm);

/** Returns the current sampling rate, in Hz */
int powermap_getSamplingRate(void* const hPm);

/** Returns the current covariance averaging coefficient value, in Hz */
float powermap_getCovAvgCoeff(void* const hPm);

/** Returns the number of frequency bands used for the analysis */
int powermap_getNumberOfBands(void);
    
/**
 * Returns the number of spherical harmonic signals required by the current
 * analysis order: (current_order + 1)^2
 */
int powermap_getNSHrequired(void* const hPm);

/**
 * Returns the weighting coefficient for a particular frequency band index,
 * allowing one to "equalise" the activity-map.
 */
float powermap_getPowermapEQ(void* const hPm, int bandIdx);

/**
 * Returns the weighting coefficient for the first frequency band
 */
float powermap_getPowermapEQAllBands(void* const hPm);

/**
 * Returns the weighting coefficient for all frequency bands
 *
 * @param[in]  hPm       powermap handle
 * @param[out] pX_vector (&) frequency vector; pNpoints x 1
 * @param[out] pY_values (&) weighting coefficients; pNpoints x 1
 * @param[out] pNpoints  (&) number of frequency bands
 */
void powermap_getPowermapEQHandle(void* const hPm,
                                  float** pX_vector,
                                  float** pY_values,
                                  int* pNpoints);

/** Returns the input/analysis order for one specific frequency band */
int powermap_getAnaOrder(void* const hPm, int bandIdx);

/** Returns the input/analysis order for the first frequency band */
int powermap_getAnaOrderAllBands(void* const hPm);

/**
 * Returns the input/analysis order for all frequency bands
 *
 * @param[in]  hPm       powermap handle
 * @param[out] pX_vector (&) frequency vector; pNpoints x 1
 * @param[out] pY_values (&) input/analysis orders; pNpoints x 1
 * @param[out] pNpoints  (&) number of frequency bands
 */
void powermap_getAnaOrderHandle(void* const hPm,
                                float** pX_vector,
                                int** pY_values,
                                int* pNpoints);

/**
 * Returns the Ambisonic channel ordering convention currently being used to
 * decode with, which should match the convention employed by the input signals
 * (see #CH_ORDER enum)
 */
int powermap_getChOrder(void* const hPm);

/**
 * Returns the Ambisonic normalisation convention currently being usedto decode
 * with, which should match the convention employed by the input signals (see
 * #NORM_TYPES enum)
 */
int powermap_getNormType(void* const hPm);

/** Returns the number of sources present in the input sound scene */
int powermap_getNumSources(void* const hPm);

/**
 * Returns the current visualisation display window horizontal field-of-view
 * (FOV) (see #HFOV_OPTIONS enum)
 */
int powermap_getDispFOV(void* const hPm);

/**
 * Returns the current visualisation display window aspect-ratio (see
 * #ASPECT_RATIO_OPTIONS enum)
 */
int powermap_getAspectRatio(void* const hPm);

/** Returns the current activity-map averaging coefficient, 0..1 */
float powermap_getPowermapAvgCoeff(void* const hPm);

/**
 * Returns the latest computed activity-map if it is ready. Otherwise it returns
 * 0, and you'll just have to wait a bit
 *
 * @param[in]  hPm         powermap handle
 * @param[out] grid_dirs   (&) scanning grid directions, in DEGREES; nDirs x 1
 * @param[out] pmap        (&) activity-map values; nDirs x 1
 * @param[out] nDirs       (&) number of directions
 * @param[out] pmapWidth   (&) activity-map width in pixels
 * @param[out] hfov        (&) horizontal FOV used to generate activity-map
 * @param[out] aspectRatio (&) aspect ratio used to generate activity-map
 * @returns flag, if activity-map is ready, 1: it is, 0: it is NOT
 */
int powermap_getPmap(void* const hPm,
                     float** grid_dirs,
                     float** pmap,
                     int* nDirs,
                     int* pmapWidth,
                     int* hfov,
                     int* aspectRatio);

/**
 * Returns the processing delay in samples (may be used for delay compensation
 * features)
 */
int powermap_getProcessingDelay(void);


#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */

#endif /* __POWERMAP_H_INCLUDED__ */
