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
 * @file sldoa.h
 * @brief A spatially-localised active-intensity based direction-of-arrival
 *        estimator (SLDoA).
 *
 * VBAP gain patterns are imposed on the spherical harmonic signals, such that
 * the DoA can be estimated in a spatially-constrained region; thus mitigating
 * the effect of interferes and reflections arriving from other directions.
 * The DoA is estimated per sector for each frequency band.
 *
 * The algorithms within sldoa were developed in collaboration with Symeon
 * Delikaris-Manias and Angelo Farina, and are explained in more detail in [1,2]
 *
 * @see [1] McCormack, L., Delikaris-Manias, S., Farina, A., Pinardi, D., and
 *          Pulkki, V., "Real-time conversion of sensor array signals into
 *          spherical harmonic signals with applications to spatially localised
 *          sub-band sound-field analysis," in Audio Engineering Society
 *          Convention 144, Audio Engineering Society, 2018.
 * @see [2] McCormack, L., Delikaris-Manias, S., Politis, A., Pavlidi, D.,
 *          Farina, A., Pinardi, D. and Pulkki, V., 2019. Applications of
 *          Spatially Localized Active-Intensity Vectors for Sound-Field
 *          Visualization. Journal of the Audio Engineering Society, 67(11),
 *          pp.840-854.
 *
 * @author Leo McCormack
 * @date 18.10.2017
 */

#ifndef __SLDOA_H_INCLUDED__
#define __SLDOA_H_INCLUDED__

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */
    
/* ========================================================================== */
/*                             Presets + Constants                            */
/* ========================================================================== */

#define SLDOA_MAX_NUM_INPUT_CHANNELS ( 64 )
    
/* Microphone/Hydrophone array options */
#define ENABLE_ZYLIA_MIC_PRESET
#define ENABLE_EIGENMIKE32_MIC_PRESET
#define ENABLE_DTU_MIC_MIC_PRESET

/**
 * "Master order" relates to the current maximum order to expect. However, the
 * analysis order can be lower for a given frequency, due to the
 * "analysisOrderPerBand" vector, which can contain lower values than the
 * master order, but not higher.
 */
typedef enum _SLDOA_MASTER_ORDERS{
    MASTER_ORDER_FIRST = 1, /**< First-order analysis (4 channel input) */
    MASTER_ORDER_SECOND,    /**< Second-order analysis (9 channel input) */
    MASTER_ORDER_THIRD,     /**< Third-order analysis (16 channel input) */
    MASTER_ORDER_FOURTH,    /**< Fourth-order analysis (25 channel input) */
    MASTER_ORDER_FIFTH,     /**< Fifth-order analysis (36 channel input) */
    MASTER_ORDER_SIXTH,     /**< Sixth-order analysis (49 channel input) */
    MASTER_ORDER_SEVENTH    /**< Seventh-order analysis (64 channel input) */
    
}SLDOA_MASTER_ORDERS;
   
/**
 * Available microphone array presets. These determine the frequency ranges
 * where the microphone array provides usable spherical harmonic components
 * at each order.
 */
typedef enum _SLDOA_MIC_PRESETS{
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
}SLDOA_MIC_PRESETS;
    
/**
 * Available Ambisonic channel ordering conventions
 *
 * @note CH_FUMA only supported for 1st order input.
 */
typedef enum _SLDOA_H_ORDER {
    CH_ACN = 1, /**< Ambisonic Channel Numbering (ACN) */
    CH_FUMA     /**< (Obsolete) Furse-Malham/B-format (WXYZ) */
    
} SLDOA_CH_ORDER;

/** Number of channel ordering options */
#define SLDOA_NUM_CH_ORDERINGS ( 2 )

/**
 * Available Ambisonic normalisation conventions
 *
 * @note NORM_FUMA only supported for 1st order input and does NOT have the
 *       1/sqrt(2) scaling on the omni.
 */
typedef enum _SLDOA_NORM_TYPES {
    NORM_N3D = 1, /**< orthonormalised (N3D) */
    NORM_SN3D,    /**< Schmidt semi-normalisation (SN3D) */
    NORM_FUMA     /**< (Obsolete) Same as NORM_SN3D for 1st order */
    
} SLDOA_NORM_TYPES;

/** Number of normalisation options */
#define SLDOA_NUM_NORM_TYPES ( 3 )

/**
 * Current status of the codec.
 */
typedef enum _SLDOA_CODEC_STATUS {
    CODEC_STATUS_INITIALISED = 0, /**< Codec is initialised and ready to process
                                   *   input audio. */
    CODEC_STATUS_NOT_INITIALISED, /**< Codec has not yet been initialised, or
                                   *   the codec configuration has changed.
                                   *   Input audio should not be processed. */
    CODEC_STATUS_INITIALISING     /**< Codec is currently being initialised,
                                   *   input audio should not be processed. */
} SLDOA_CODEC_STATUS;

#define SLDOA_PROGRESSBARTEXT_CHAR_LENGTH ( 256 )
    

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
 * @param[in] hSld       sldoa handle
 * @param[in] samplerate Host samplerate.
 */
void sldoa_init(void* const hSld,
                float samplerate);
    
/**
 * Intialises the codec variables, based on current global/user parameters
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
                    float** const inputs,
                    int nInputs,
                    int nSamples,
                    int isPlaying);


/* ========================================================================== */
/*                                Set Functions                               */
/* ========================================================================== */

/**
 * Sets the maximum input/analysis order (see 'SLDOA_MASTER_ORDERS' enum)
 */
void sldoa_setMasterOrder(void* const hSld,  int newValue);

/**
 * Sets all intialisation flags to 1; re-initialising all settings/variables
 * as sldoa is currently configured, at next available opportunity.
 */
void sldoa_refreshSettings(void* const hSld);
   
/**
 * Sets the maximum analysis frequency, in Hz
 */
void sldoa_setMaxFreq(void* const hSld, float newFreq);

/**
 * Sets the minimum analysis frequency, in Hz
 */
void sldoa_setMinFreq(void* const hSld, float newFreq);
    
/**
 * Sets the DoA averaging coefficient, 0..1
 */
void sldoa_setAvg(void* const hSld, float newAvg);

/**
 * Sets the input/analysis order for one specific frequency band.
 */
void sldoa_setAnaOrder(void* const hSld, int newValue, int bandIdx);

/**
 * Sets the input/analysis order for all frequency bands.
 */
void sldoa_setAnaOrderAllBands(void* const hSld,  int newValue);

/**
 * Sets the Ambisonic channel ordering convention to decode with, in order to
 * match the convention employed by the input signals (see 'SLDOA_CH_ORDER'
 * enum)
 */
void sldoa_setChOrder(void* const hSld, int newOrder);

/**
 * Sets the Ambisonic normalisation convention to decode with, in order to match
 * with the convention employed by the input signals (see 'SLDOA_NORM_TYPE'
 * enum)
 */
void sldoa_setNormType(void* const hSld, int newType);

/**
 * Sets an input preset, the microphone/hyrophone array used to capture
 * the input signals (see 'SLDOA_MIC_PRESETS' enum)
 */
void sldoa_setSourcePreset(void* const hSld, int newPresetID);


/* ========================================================================== */
/*                                Get Functions                               */
/* ========================================================================== */

/**
 * Returns current codec status (see 'SLDOA_CODEC_STATUS' enum)
 */
SLDOA_CODEC_STATUS sldoa_getCodecStatus(void* const hSld);

/**
 * (Optional) Returns current intialisation/processing progress, between 0..1
 *  - 0: intialisation/processing has started
 *  - 1: intialisation/processing has ended
 */
float sldoa_getProgressBar0_1(void* const hSld);

/**
 * (Optional) Returns current intialisation/processing progress text
 * @note "text" string should be (at least) of length:
 *       SLDOA_PROGRESSBARTEXT_CHAR_LENGTH 
 */
void sldoa_getProgressBarText(void* const hSld, char* text);

/**
 * Returns the current maximum analysis/input order (see 'SLDOA_MASTER_ORDERS'
 * enum)
 */
int sldoa_getMasterOrder(void* const hSld);

/**
 * Returns the current sampling rate, in Hz
 */
int sldoa_getSamplingRate(void* const hSld);

/**
 * Returns the maximum analysis frequency, in Hz
 */
float sldoa_getMaxFreq(void* const hSld);

/**
 * Returns the minimum analysis frequency, in Hz
 */
float sldoa_getMinFreq(void* const hSld);

/**
 * Returns the current DoA averaging coefficient value, 0..1
 */
float sldoa_getAvg(void* const hSld);

/**
 * Returns the number frequency bands employed by sldoa
 */
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
 * @param[in]  hSld             sldoa handle
 * @param[out] pAzi_deg         (&) azimuth of estimated DoAs;
 *                              FLAT: pNsectorsPerBand*sldoa_getNumberOfBands()
 * @param[out] pElev_deg        (&) elevation of estimated DoAs;
 *                              FLAT: pNsectorsPerBand*sldoa_getNumberOfBands()
 * @param[out] pColourScale     (&) colour scale, 0..1, 1:red, 0: blue
 *                              FLAT: pNsectorsPerBand*sldoa_getNumberOfBands()
 * @param[out] pAlphaScale      (&) alpha scale, 0..1, 1: opaque, 0: transparent;
 *                              FLAT: pNsectorsPerBand*sldoa_getNumberOfBands()
 * @param[out] pNsectorsPerBand (&) number of sectors per frequency;
 *                              pNsectorsPerBand x 1
 * @param[out] pMaxNumSectors   (&) maximum number of sectors
 * @param[out] pStartBand       (&) band index corresponding to lowest frequency
 * @param[out] pEndBand         (&) band index corresponding to highest frequency
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
    
/**
 * Returns the input/analysis order for one specific frequency band.
 */
int sldoa_getAnaOrder(void* const hSld, int bandIdx);

/**
 * Returns the input/analysis order for the first frequency band
 */
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
 * (see 'SLDOA_CH_ORDER' enum)
 */
int sldoa_getChOrder(void* const hSld);

/**
 * Returns the Ambisonic normalisation convention currently being usedto decode
 * with, which should match the convention employed by the input signals.
 * (see 'SLDOA_NORM_TYPE' enum)
 */
int sldoa_getNormType(void* const hSld);

    
#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */

#endif /* __SLDOA_H_INCLUDED__ */
