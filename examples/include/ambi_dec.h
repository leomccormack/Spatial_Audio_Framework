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
 * @example ambi_dec.h
 * @brief A frequency-dependent Ambisonic decoder for reproducing Ambisonic
 *        sound scenes over loudspeakers
 *
 * ### Files
 * ambi_dec.h (include), ambi_dec_internal.h, ambi_dec.c, ambi_dec_internal.c
 * ### Include Header
 */

/**
 * @file ambi_dec.h
 * @brief A frequency-dependent Ambisonic decoder for reproducing Ambisonic
 *        sound scenes over loudspeakers
 *
 * Different decoder settings can be specified for the low and high frequencies.
 * A number of decoding options are also offered, including [1,2]. When
 * utilising spherical harmonic signals derived from real microphone arrays,
 * this implementation also allows the decoding order to be specified per
 * frequency band; of course, this may also be used creatively. An optional,
 * loudspeaker channel binauraliser is included, along with with SOFA file
 * loading, for headphone listening.
 *
 * The algorithms utilised in this Ambisonic decoder were pieced together and
 * developed in collaboration with Archontis Politis.
 *
 * @test test__saf_example_ambi_dec()
 *
 * @see [1] Zotter F, Pomberger H, Noisternig M. Energy--preserving ambisonic
 *          decoding. Acta Acustica united with Acustica. 2012 Jan 1;
 *          98(1):37-47.
 * @see [2] Zotter F, Frank M. All-round ambisonic panning and decoding. Journal
 *          of the audio engineering society. 2012 Nov 26; 60(10):807-20.
 *
 * @author Leo McCormack
 * @date 07.12.2017
 * @license ISC
 */

#ifndef __AMBI_DEC_H_INCLUDED__
#define __AMBI_DEC_H_INCLUDED__

#include "_common.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/* ========================================================================== */
/*                             Presets + Constants                            */
/* ========================================================================== */

/**
 * Available decoding methods. See saf_hoa_internal.h for a more indepth
 * description of each decoding approach.
 */
typedef enum {
    DECODING_METHOD_SAD = 1, /**< Sampling Ambisonic Decoder (SAD) */
    DECODING_METHOD_MMD,     /**< Mode-Matching Decoder (MMD) */
    DECODING_METHOD_EPAD,    /**< Energy-Preserving Ambisonic Decoder (EPAD) */
    DECODING_METHOD_ALLRAD   /**< All-Round Ambisonic Decoder (AllRAD) */
    
} AMBI_DEC_DECODING_METHODS;

/** Number of decoding method options */
#define AMBI_DEC_NUM_DECODING_METHODS ( 4 )

/**
 * When using mixed order decoding (i.e. different decoding orders for
 * different frequencies), this equalisation helps maintain equal perceived
 * "loudness"
 *
 * At low frequencies, preserving amplitude is more favourable, whereas for high
 * frequencies, preserving energy is better.
 */
typedef enum {
    AMPLITUDE_PRESERVING=1, /**< preserve omni amplitude */
    ENERGY_PRESERVING       /**< preserve omni energy */
    
} AMBI_DEC_DIFFUSE_FIELD_EQ_APPROACH;

/** Minimum transition value between low/high frequency decoders, in Hz */
#define AMBI_DEC_TRANSITION_MIN_VALUE ( 500.0f )

/** Maximum transition value between low/high frequency decoders, in Hz */
#define AMBI_DEC_TRANSITION_MAX_VALUE ( 2000.0f )


/* ========================================================================== */
/*                               Main Functions                               */
/* ========================================================================== */

/**
 * Creates an instance of the ambi_dec
 *
 * @param[in] phAmbi (&) address of ambi_dec handle
 */
void ambi_dec_create(void** const phAmbi);

/**
 * Destroys an instance of the ambi_dec
 *
 * @param[in] phAmbi (&) address of ambi_dec handle
 */
void ambi_dec_destroy(void** const phAmbi);

/**
 * Initialises an instance of ambi_dec with default settings
 *
 * @warning This should not be called while _process() is on-going!
 *
 * @param[in] hAmbi      ambi_dec handle
 * @param[in] samplerate Host samplerate.
 */
void ambi_dec_init(void* const hAmbi,
                   int samplerate);

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
 * @param[in] hAmbi      ambi_dec handle
 */
void ambi_dec_initCodec(void* const hAmbi);

/**
 * Decodes input spherical harmonic signals to the loudspeaker channels
 *
 * @param[in] hAmbi    ambi_dec handle
 * @param[in] inputs   Input channel buffers; 2-D array: nInputs x nSamples
 * @param[in] outputs  Output channel buffers; 2-D array: nOutputs x nSamples
 * @param[in] nInputs  Number of input channels
 * @param[in] nOutputs Number of output channels
 * @param[in] nSamples Number of samples in 'inputs'/'output' matrices
 */
void ambi_dec_process(void* const hAmbi,
                      const float *const * inputs,
                      float** const outputs,
                      int nInputs,
                      int nOutputs,
                      int nSamples);


/* ========================================================================== */
/*                                Set Functions                               */
/* ========================================================================== */

/**
 * Sets all intialisation flags to 1. Re-initialising all settings/variables,
 * as ambi_dec is currently configured, at next available opportunity.
 */
void ambi_dec_refreshSettings(void* const hAmbi);
    
/**
 * Sets the master decoding order. However, the decoding order may be lower than
 * this for any given frequency, this is just the maximum.
 *
 * If decoding order is higher than the input signal order, the extra required
 * channels are filled with zeros. If the decoding order is lower than the input
 * signal order, the number input signals is truncated accordingly. (see
 * #SH_ORDERS enum)
 */
void ambi_dec_setMasterDecOrder(void* const hAmbi,  int newValue);
    
/**
 * Sets the decoding order for a given frequency band
 *
 * @note The maximum order is dictated by ambi_dec_setMasterDecOrder()
 *
 * @param[in] hAmbi    ambi_dec handle
 * @param[in] newValue New decoding order
 * @param[in] bandIdx  Frequency band index
 */
void ambi_dec_setDecOrder(void* const hAmbi,  int newValue, int bandIdx);

/**
 * Sets the decoding order for all frequency bands
 *
 * @note the maximum order is dictated by ambi_dec_setMasterDecOrder()
 *
 * @param[in] hAmbi    ambi_dec handle
 * @param[in] newValue New decoding order
 */
void ambi_dec_setDecOrderAllBands(void* const hAmbi,  int newValue);

/**
 * Sets the azimuth of a specific loudspeaker
 *
 * @param[in] hAmbi      ambi_dec handle
 * @param[in] index      Loudspeaker index
 * @param[in] newAzi_deg New azimuth in DEGREES
 */
void ambi_dec_setLoudspeakerAzi_deg(void* const hAmbi,
                                    int index,
                                    float newAzi_deg);

/**
 * Sets the elevation of a specific loudspeaker
 *
 * @param[in] hAmbi       ambi_dec handle
 * @param[in] index       Loudspeaker index
 * @param[in] newElev_deg New elevation in DEGREES
 */
void ambi_dec_setLoudspeakerElev_deg(void* const hAmbi,
                                     int index,
                                     float newElev_deg);

/**
 * Sets the number of loudspeakers to decode to
 */
void ambi_dec_setNumLoudspeakers(void* const hAmbi, int new_nLoudspeakers);

/**
 * Sets flag to dictate whether the output loudspeaker signals should be
 * binauralised
 *
 * @param[in] hAmbi    ambi_dec handle
 * @param[in] newState '0' output loudspeaker signals, '1' output binaural
 *                     signals
 */
void ambi_dec_setBinauraliseLSflag(void* const hAmbi, int newState);

/**
 * Sets flag to dictate whether the default HRIRs in the Spatial_Audio_Framework
 * should be used (1), or a custom HRIR set loaded via a SOFA file (0).
 *
 * @note If the custom set fails to load correctly, ambi_dec will revert to the
 *       default set. Use ambi_dec_getUseDefaultHRIRsflag() to check if loading
 *       was successful.
 *
 * @param[in] hAmbi     ambi_dec handle
 * @param[in] newState  '0' use custom HRIR set, '1' use default HRIR set
 */
void ambi_dec_setUseDefaultHRIRsflag(void* const hAmbi, int newState);

/**
 * Sets the file path for a .sofa file, in order to employ a custom HRIR set for
 * the decoding.
 *
 * @note If the custom set failes to load correctly, ambi_dec will revert to the
 *       default set. Use ambi_dec_getUseDefaultHRIRsflag() to check if loading
 *       was successful.
 *
 * @param[in] hAmbi ambi_dec handle
 * @param[in] path  File path to .sofa file (WITH file extension)
 */
void ambi_dec_setSofaFilePath(void* const hAmbi, const char* path);

/** Enable (1) or disable (0) the pre-processing applied to the HRTFs. */
void ambi_dec_setEnableHRIRsPreProc(void* const hAmbi, int newState);

/**
 * Sets the source preset (ideal SH or SH signals derived from mic arrays)
 *
 * By default the decoder will decode at the maximum order for all frequencies.
 * However, in the case of spherical harmonic input derived from microphone
 * arrays, the available order is frequency dependent, therefore, different
 * bands require different decoding orders.
 * For conveinience, presets for a handful of comerically available microphone
 * array are included (see #MIC_PRESETS enum).
 */
void ambi_dec_setSourcePreset(void* const hAmbi, int newPresetID);

/**
 * Sets the output loudspeaker preset.
 *
 * For conveinience, presets for several popular arrangements are included (see
 * #LOUDSPEAKER_ARRAY_PRESETS enum).
 */
void ambi_dec_setOutputConfigPreset(void* const hAmbi, int newPresetID);

/**
 * Sets the Ambisonic channel ordering convention to decode with, in order to
 * match the convention employed by the input signals (see #CH_ORDER enum)
 */
void ambi_dec_setChOrder(void* const hAmbi, int newOrder);

/**
 * Sets the Ambisonic normalisation convention to decode with, in order to match
 * with the convention employed by the input signals (see #NORM_TYPES enum)
 */
void ambi_dec_setNormType(void* const hAmbi, int newType);

/**
 * Sets the decoding method for a specific decoder. ambi_dec employs two
 * decoders, one for low frequencies and one for high frequencies.
 * (use ambi_dec_setTransitionFreq() to dictate the transition frequency)
 *
 * @param[in] hAmbi ambi_dec handle
 * @param[in] index Index for low (0) or high (1) frequency decoder
 * @param[in] newID New decoding method (see #AMBI_DEC_DECODING_METHODS enum)
 */
void ambi_dec_setDecMethod(void* const hAmbi, int index, int newID);

/**
 * Sets a flag to enable/disable (1 or 0) the max_rE weighting for one of the
 * decoders.
 *
 * @param[in] hAmbi ambi_dec handle
 * @param[in] index For low (0) or high (1) frequency decoder
 * @param[in] newID '0' disabled, '1' enabled
 */
void ambi_dec_setDecEnableMaxrE(void* const hAmbi, int index, int newID);

/**
 * Sets the equalisation approach for one of the decoders. This is used to help
 * keep the perceived loudness consistent, when using mixed decoding orders
 * (i.e. different decoding orders for different frequency bands)
 * ambi_dec either to preserves amplitude or energy for each order.
 *
 * @note It is suggested to preserve amplitude at low-frequencies and energy
 *       at high-frequencies.
 *
 * @param[in] hAmbi ambi_dec handle
 * @param[in] index For low (0) or high (1) frequency decoder
 * @param[in] newID see #AMBI_DEC_DIFFUSE_FIELD_EQ_APPROACH enum
 */
void ambi_dec_setDecNormType(void* const hAmbi, int index, int newID);

/**
 * Sets the frequeny at which to transition from the low frequency decoder to
 * the high frequency decoder.
 *
 * @param[in] hAmbi    ambi_dec handle
 * @param[in] newValue New transition frequency, in Hz
 */
void ambi_dec_setTransitionFreq(void* const hAmbi, float newValue);


/* ========================================================================== */
/*                                Get Functions                               */
/* ========================================================================== */

/**
 * Returns the processing framesize (i.e., number of samples processed with
 * every _process() call )
 */
int ambi_dec_getFrameSize(void);

/** Returns current codec status (see #CODEC_STATUS enum) */
CODEC_STATUS ambi_dec_getCodecStatus(void* const hAmbi);

/**
 * (Optional) Returns current intialisation/processing progress, between 0..1
 * - 0: intialisation/processing has started
 * - 1: intialisation/processing has ended
 */
float ambi_dec_getProgressBar0_1(void* const hAmbi);

/**
 * (Optional) Returns current intialisation/processing progress text
 *
 * @note "text" string should be (at least) of length:
 *       #PROGRESSBARTEXT_CHAR_LENGTH
 */
void ambi_dec_getProgressBarText(void* const hAmbi, char* text);

/** Returns the master/maximum decoding order (see #SH_ORDERS enum) */
int ambi_dec_getMasterDecOrder(void* const hAmbi);
    
/**
 * Returns the decoding order for a given frequency band index (see #SH_ORDERS
 * enum)
 */
int ambi_dec_getDecOrder(void* const hAmbi, int bandIdx);
  
/** Returns the decoding order for the first band (see #SH_ORDERS enum) */
int ambi_dec_getDecOrderAllBands(void* const hAmbi);

/**
 * Returns handles for the decoding orders and frequency vector.
 *
 * @param[in]  hAmbi     ambi_dec handle
 * @param[out] pX_vector (&) frequency vector; pNpoints x 1
 * @param[out] pY_values (&) decoding order per frequency; pNpoints x 1
 * @param[out] pNpoints  (&) number of grid points.
 */
void ambi_dec_getDecOrderHandle(void* const hAmbi,
                                float** pX_vector,
                                int** pY_values,
                                int* pNpoints);

/** Returns the number of frequency bands employed by ambi_dec */
int ambi_dec_getNumberOfBands(void);

/** Returns the loudspeaker azimuth in degrees for a given index */
float ambi_dec_getLoudspeakerAzi_deg(void* const hAmbi, int index);

/** Returns the loudspeaker elevation in degrees for a given index */
float ambi_dec_getLoudspeakerElev_deg(void* const hAmbi, int index);

/** Returns the number of loudspeakers in the current layout */
int ambi_dec_getNumLoudspeakers(void* const hAmbi);

/** Returns the maximum number of loudspeakers supported by ambi_dec */
int ambi_dec_getMaxNumLoudspeakers(void);

/**
 * Returns the number of spherical harmonic signals required by the current
 * decoding order: (current_order + 1)^2
 */
int  ambi_dec_getNSHrequired(void* const hAmbi); 

/**
 * Returns the value of a flag used to dictate whether the loudspeaker signals
 * should be binauralised (0: output loudspeaker signals, 1: output binaural
 * signals).
 */
int ambi_dec_getBinauraliseLSflag(void* const hAmbi);

/**
 * Returns the value of a flag used to dictate whether the default HRIRs in the
 * Spatial_Audio_Framework should be used (1), or a custom HRIR set loaded via a
 * SOFA file (0).
 *
 * @note If the custom set failes to load correctly, ambi_dec will revert to the
 *       default set.
 */
int ambi_dec_getUseDefaultHRIRsflag(void* const hAmbi);

/**
 * Returns the file path for a .sofa file (WITH file extension)
 *
 * @note If the custom set failes to load correctly, ambi_dec will revert to the
 *       default set. Use ambi_dec_getUseDefaultHRIRsflag() to check if loading
 *       was successful.
 */
char* ambi_dec_getSofaFilePath(void* const hAmbi);

/**
 * Returns the flag indicating whether the pre-processing applied to the HRTFs
 * is enabled (1) or disabled (0)
 */
int ambi_dec_getEnableHRIRsPreProc(void* const hAmbi);

/**
 * Returns the Ambisonic channel ordering convention currently being used to
 * decode with, which should match the convention employed by the input signals
 * (see #CH_ORDER enum)
 */
int ambi_dec_getChOrder(void* const hAmbi);

/**
 * Returns the Ambisonic normalisation convention currently being usedto decode
 * with, which should match the convention employed by the input signals (see
 * #NORM_TYPES enum).
 */
int ambi_dec_getNormType(void* const hAmbi);

/**
 * Returns the currently selected decoding method (see
 * #AMBI_DEC_DECODING_METHODS enum)
 */
int ambi_dec_getDecMethod(void* const hAmbi, int index);
    
/**
 * Returns the value of a flag used to dictate whether the max_rE weighting is
 * being applied by a given decoder
 *
 * @param[in] hAmbi ambi_dec handle
 * @param[in] index Index for low (0) or high (1) frequency decoder
 * @returns '0' if enabled, '1' if disabled
 */
int ambi_dec_getDecEnableMaxrE(void* const hAmbi, int index);

/**
 * Returns the current equalisation approach for one of the decoders (see
 * #AMBI_DEC_DIFFUSE_FIELD_EQ_APPROACH enum)
 *
 * @note It is suggested to preserve amplitude at low-frequencies and energy
 *       at high-frequencies
 *
 * @param[in] hAmbi ambi_dec handle
 * @param[in] index Index for low (0) or high (1) frequency decoder
 * @returns (see #AMBI_DEC_DIFFUSE_FIELD_EQ_APPROACH enum)
 */
int ambi_dec_getDecNormType(void* const hAmbi, int index);

/**
 * Returns the frequency (in Hz) at which to transition from the low frequency
 * decoder to the high frequency decoder.
 */
float ambi_dec_getTransitionFreq(void* const hAmbi);
    
/** Returns the HRIR sample rate */
int ambi_dec_getHRIRsamplerate(void* const hAmbi);

/** Returns the DAW/Host sample rate */
int ambi_dec_getDAWsamplerate(void* const hAmbi);
    
/**
 * Returns the processing delay in samples; may be used for delay compensation
 * features
 */
int ambi_dec_getProcessingDelay(void);


#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */

#endif /* __AMBI_DEC_H_INCLUDED__ */
