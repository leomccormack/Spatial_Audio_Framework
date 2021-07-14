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
 * @example sldoa.h
 * @brief A spatially-localised active-intensity (SLAI) based direction-of-
 *        arrival estimator (SLDoA)
 * 
 * ### Files
 * sldoa.h (include), sldoa_internal.h, sldoa_database.h, sldoa.c,
 * sldoa_internal.c, sldoa_database.c
 * ### Include Header
 */

/**
 * @file sldoa.h
 * @brief A spatially-localised active-intensity (SLAI) based direction-of-
 *        arrival estimator (SLDoA)
 *
 * VBAP gain patterns are imposed on the spherical harmonic signals, such that
 * the DoA can be estimated in a spatially-constrained region; thus mitigating
 * the effect of interferes and reflections arriving from other directions.
 * The DoA is estimated per sector for each frequency band.
 *
 * The algorithms within sldoa were developed in collaboration with Symeon
 * Delikaris-Manias and Angelo Farina, and are explained in more detail in [1,2]
 *
 * @see [1] McCormack, L., Delikaris-Manias, S., Politis, A., Pavlidi, D.,
 *          Farina, A., Pinardi, D. and Pulkki, V., 2019. Applications of
 *          Spatially Localized Active-Intensity Vectors for Sound-Field
 *          Visualization. Journal of the Audio Engineering Society, 67(11),
 *          pp.840-854.
 * @see [2] McCormack, L., Delikaris-Manias, S., Farina, A., Pinardi, D., and
 *          Pulkki, V., "Real-time conversion of sensor array signals into
 *          spherical harmonic signals with applications to spatially localised
 *          sub-band sound-field analysis," in Audio Engineering Society
 *          Convention 144, Audio Engineering Society, 2018.
 *
 * @author Leo McCormack
 * @date 18.10.2017
 * @license ISC
 */

#ifndef __SLDOA_H_INCLUDED__
#define __SLDOA_H_INCLUDED__

#include "_common.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/* ========================================================================== */
/*                               Main Functions                               */
/* ========================================================================== */

/**
 * Creates an instance of the sldoa
 *
 * @param[in] phSld (&) address of sldoa handle
 */
void sldoa_create(void** const phSld);

/**
 * Destroys an instance of the sldoa
 *
 * @param[in] phSld (&) address of sldoa handle
 */
void sldoa_destroy(void** const phSld);

/**
 * Initialises an instance of sldoa with default settings
 *
 * @warning This should not be called while _process() is on-going!
 *
 * @param[in] hSld       sldoa handle
 * @param[in] samplerate Host samplerate.
 */
void sldoa_init(void* const hSld,
                float samplerate);
    
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
 * @param[in] hSld - sldoa handle
 */
void sldoa_initCodec(void* const hSld);
    
/**
 * Applies the spatially-localised active-intensity based direction-of-arrival
 * estimator (SLDoA) onto the input signals [1,2].
 *
 * @param[in] hSld      sldoa handle
 * @param[in] inputs    Input channel buffers; 2-D array: nInputs x nSamples
 * @param[in] nInputs   Number of input channels
 * @param[in] nSamples  Number of samples in 'inputs'/'output' matrices
 * @param[in] isPlaying Flag to say if there is audio in the damn input buffers,
 *                      0: no audio, reduced processing, 1: audio, full processing
 *
 * @see [1] McCormack, L., Delikaris-Manias, S., Farina, A., Pinardi, D., and
 *          Pulkki, V., “Real-time conversion of sensor array signals into
 *          spherical harmonic signals with applications to spatially localised
 *          sub-band sound-field analysis,” in Audio Engineering Society
 *          Convention 144, Audio Engineering Society, 2018.
 * @see [2] McCormack, L., Delikaris-Manias, S., Politis, A., Pavlidi, D.,
 *          Farina, A., Pinardi, D. and Pulkki, V., 2019. Applications of
 *          Spatially Localized Active-Intensity Vectors for Sound-Field
 *          Visualization. Journal of the Audio Engineering Society, 67(11),
 *          pp.840-854.
 */
void sldoa_analysis(void* const hSld,
                    const float *const * inputs,
                    int nInputs,
                    int nSamples,
                    int isPlaying);


/* ========================================================================== */
/*                                Set Functions                               */
/* ========================================================================== */

/** Sets the maximum input/analysis order (see #SH_ORDERS enum) */
void sldoa_setMasterOrder(void* const hSld,  int newValue);

/**
 * Sets all intialisation flags to 1; re-initialising all settings/variables
 * as sldoa is currently configured, at next available opportunity.
 */
void sldoa_refreshSettings(void* const hSld);
   
/** Sets the maximum analysis frequency, in Hz */
void sldoa_setMaxFreq(void* const hSld, float newFreq);

/** Sets the minimum analysis frequency, in Hz */
void sldoa_setMinFreq(void* const hSld, float newFreq);
    
/** Sets the DoA averaging coefficient, 0..1 */
void sldoa_setAvg(void* const hSld, float newAvg);

/** Sets the input/analysis order for one specific frequency band */
void sldoa_setAnaOrder(void* const hSld, int newValue, int bandIdx);

/** Sets the input/analysis order for all frequency bands */
void sldoa_setAnaOrderAllBands(void* const hSld,  int newValue);

/**
 * Sets the Ambisonic channel ordering convention to decode with, in order to
 * match the convention employed by the input signals (see #CH_ORDER enum)
 */
void sldoa_setChOrder(void* const hSld, int newOrder);

/**
 * Sets the Ambisonic normalisation convention to decode with, in order to match
 * with the convention employed by the input signals (see #NORM_TYPES enum)
 */
void sldoa_setNormType(void* const hSld, int newType);

/**
 * Sets an input preset, the microphone/hyrophone array used to capture
 * the input signals (see #MIC_PRESETS enum)
 */
void sldoa_setSourcePreset(void* const hSld, int newPresetID);


/* ========================================================================== */
/*                                Get Functions                               */
/* ========================================================================== */

/**
 * Returns the processing framesize (i.e., number of samples processed with
 * every _process() call )
 */
int sldoa_getFrameSize(void);

/** Returns current codec status (see #CODEC_STATUS enum) */
CODEC_STATUS sldoa_getCodecStatus(void* const hSld);

/**
 * (Optional) Returns current intialisation/processing progress, between 0..1
 *  - 0: intialisation/processing has started
 *  - 1: intialisation/processing has ended
 */
float sldoa_getProgressBar0_1(void* const hSld);

/**
 * (Optional) Returns current intialisation/processing progress text
 *
 * @note "text" string should be (at least) of length:
 *       #PROGRESSBARTEXT_CHAR_LENGTH
 */
void sldoa_getProgressBarText(void* const hSld, char* text);

/** Returns the current maximum analysis/input order (see #SH_ORDERS enum) */
int sldoa_getMasterOrder(void* const hSld);

/** Returns the current sampling rate, in Hz */
int sldoa_getSamplingRate(void* const hSld);

/** Returns the maximum analysis frequency, in Hz */
float sldoa_getMaxFreq(void* const hSld);

/** Returns the minimum analysis frequency, in Hz */
float sldoa_getMinFreq(void* const hSld);

/** Returns the current DoA averaging coefficient value, 0..1 */
float sldoa_getAvg(void* const hSld);

/** Returns the number frequency bands employed by sldoa */
int sldoa_getNumberOfBands(void);

/**
 * Returns the number of spherical harmonic signals required by the current
 * analysis order: (current_order + 1)^2
 */
int sldoa_getNSHrequired(void* const hSld);
    
/**
 * Returns the analysis output data. Including the DoAs per frequency, and per
 * sector, accompanied by colour coefficients (red: high frequencies, blue:
 * low frequencies), and alpha coefficients (more opaque: higher energy, more
 * transpararent: less energy).
 *
 * @note nBands can be found by using sldoa_getNumberOfBands()
 *
 * @param[in]  hSld             sldoa handle
 * @param[out] pAzi_deg         (&) azimuth of estimated DoAs;
 *                              FLAT: pNsectorsPerBand x nBands
 * @param[out] pElev_deg        (&) elevation of estimated DoAs;
 *                              FLAT: pNsectorsPerBand x nBands
 * @param[out] pColourScale     (&) colour scale, 0..1, 1:red, 0: blue
 *                              FLAT: pNsectorsPerBand x nBands
 * @param[out] pAlphaScale      (&) alpha scale, 0..1, 1: opaque, 0: transparent
 *                              FLAT: pNsectorsPerBand x nBands
 * @param[out] pNsectorsPerBand (&) number of sectors per frequency;
 *                              pNsectorsPerBand x 1
 * @param[out] pMaxNumSectors   (&) maximum number of sectors
 * @param[out] pStartBand       (&) band index corresponding to lowest frequency
 * @param[out] pEndBand         (&) band index corresponding to highest
 *                              frequency
 */
void sldoa_getDisplayData(void *  const hSld,
                          float** pAzi_deg,
                          float** pElev_deg,
                          float** pColourScale,
                          float** pAlphaScale,
                          int** pNsectorsPerBand,
                          int* pMaxNumSectors,
                          int* pStartBand,
                          int* pEndBand);
    
/** Returns the input/analysis order for one specific frequency band */
int sldoa_getAnaOrder(void* const hSld, int bandIdx);

/** Returns the input/analysis order for the first frequency band */
int sldoa_getAnaOrderAllBands(void* const hSld);

/**
 * Returns the input/analysis order for all frequency bands
 *
 * @param[in]  hSld      sldoa handle
 * @param[out] pX_vector (&) frequency vector; pNpoints x 1
 * @param[out] pY_values (&) input/analysis orders; pNpoints x 1
 * @param[out] pNpoints  (&) number of frequency bands
 */
void sldoa_getAnaOrderHandle(void* const hSld,
                             float** pX_vector,
                             int** pY_values,
                             int* pNpoints);

/**
 * Returns the Ambisonic channel ordering convention currently being used to
 * decode with, which should match the convention employed by the input signals
 * (see #CH_ORDER enum)
 */
int sldoa_getChOrder(void* const hSld);

/**
 * Returns the Ambisonic normalisation convention currently being usedto decode
 * with, which should match the convention employed by the input signals.
 * (see #NORM_TYPES enum)
 */
int sldoa_getNormType(void* const hSld);

/**
 * Returns the processing delay in samples (may be used for delay compensation
 * features)
 */
int sldoa_getProcessingDelay(void);

    
#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */

#endif /* __SLDOA_H_INCLUDED__ */
