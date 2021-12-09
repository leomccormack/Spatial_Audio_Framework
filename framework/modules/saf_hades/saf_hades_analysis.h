/*
 * This file is part of the saf_hades module.
 * Copyright (c) 2021 - Leo McCormack & Janani Fernandez
 *
 * The saf_hades module is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at your
 * option) any later version.
 *
 * The saf_hades module is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 * or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
 * more details.
 *
 * See <http://www.gnu.org/licenses/> for a copy of the GNU General Public
 * License.
 */

/**
 * @file saf_hades_analysis.h
 * @ingroup HADES
 * @brief Header for the HADES analysis (#SAF_HADES_MODULE)
 *
 * The framework for binaural rendering of Hearing-Assistive/Augmented-reality
 * Devices (HADES) is described further in [1].
 *
 * @see [1] paper submitted for review.
 *
 * @author Leo McCormack and Janani Fernandez
 * @date 01.02.2021
 * @license GNU GPLv2
 */

#ifndef __SAF_HADES_ANALYSIS_H_INCLUDED__
#define __SAF_HADES_ANALYSIS_H_INCLUDED__

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#ifdef  SAF_ENABLE_HADES_MODULE

/** Maximum number of microphones */
#define HADES_MAX_NMICS ( 64 )

/** Handle for the hades analysis data */
typedef struct _hades_analysis_data* hades_analysis_handle;

/** Handle for the hades parameter container data */
typedef struct _hades_param_container_data* hades_param_container_handle;

/** Handle for the hades signal container data */
typedef struct _hades_signal_container_data* hades_signal_container_handle;

/* ========================================================================== */
/*                   HADES Analysis Configurations Options                    */
/* ========================================================================== */

/**
 * Options for diffuseness estimation for hades_analysis
 *
 * @see [1] Epain, N. and Jin, C.T., 2016. Spherical harmonic signal covariance
 *        and sound field diffuseness. IEEE/ACM Transactions on Audio, Speech,
 *        and Language Processing, 24(10), pp.1796-1807.
 */
typedef enum {
    HADES_USE_COMEDIE /**< As in [1], after spatially whitening the array SCM */
}HADES_DIFFUSENESS_ESTIMATORS;

/** Options for DoA estimation for hades_analysis */
typedef enum {
    HADES_USE_MUSIC      /**< Use MUSIC */
}HADES_DOA_ESTIMATORS;

/** Filterbank options */
typedef enum {
    HADES_USE_AFSTFT_LD, /**< Alias-free STFT filterbank (low delay) */
    HADES_USE_AFSTFT     /**< Alias-free STFT filterbank */
}HADES_FILTERBANKS;


/* ========================================================================== */
/*                              HADES Analysis                                */
/* ========================================================================== */

/**
 * Creates and returns a handle to an instance of a hades analysis object
 *
 * @param[in] phAna         (&) address of hades analysis handle
 * @param[in] fs            Samplerate, Hz
 * @param[in] fbOption      Filterbank to use (see #HADES_FILTERBANKS)
 * @param[in] hopsize       Filterbank hopsize
 * @param[in] blocksize     Number of time-domain samples to process at a time
 * @param[in] hybridmode    1: enable filterbank hybrid-filtering, 0: disable
 * @param[in] h_array       Array impulse responses; FLAT: nGrid x nMics x h_len
 * @param[in] grid_dirs_deg Measurement dirs [azi elev] degrees; FLAT: nGrid x 2
 * @param[in] nGrid         Number of mesurement directions
 * @param[in] nMics         Number of microphones
 * @param[in] h_len         Length of impulse responses, in samples
 * @param[in] diffOption    Diffusness parameter estimator to use (see
 *                          #HADES_DIFFUSENESS_ESTIMATORS)
 * @param[in] doaOption     DoA estimator to use (see #HADES_DOA_ESTIMATORS)
 */
void hades_analysis_create(/* Input Arguments */
                           hades_analysis_handle* const phAna,
                           float fs,
                           HADES_FILTERBANKS fbOption,
                           int hopsize,
                           int blocksize,
                           int hybridmode,
                           float* h_array,
                           float* grid_dirs_deg,
                           int nGrid,
                           int nMics,
                           int h_len,
                           HADES_DIFFUSENESS_ESTIMATORS diffOption,
                           HADES_DOA_ESTIMATORS doaOption);

/**
 * Destroys an instance of a hades analysis object
 *
 * @param[in] phAna (&) address of hades analysis handle
 */
void hades_analysis_destroy(/* Input Arguments */
                            hades_analysis_handle* const phAna);

/**
 * Flushes run-time buffers with zeros
 *
 * Call this ONCE before calling hades_analysis_apply()
 *
 * @param[in] hAna hades analysis handle
 */
void hades_analysis_reset(hades_analysis_handle const hAna);

/**
 * Performs hades encoding: forward time-frequency transform, diffuseness
 * and DoA estimation per band
 *
 * @note See hades_param_container_create() and
 *       hades_signal_container_create() for creating the parameter and signal
 *       containers, respectively. The former contains the estimated spatial
 *       parameters (a diffuseness measure, DoA for each source), while
 *       the latter contains the input signals in the time-frequency domain,
 *       and their spatial covariance matrices per band. These containers can
 *       then be passed to hades_synthesis_apply() to reproduce the encoded
 *       scene over the target setup.
 *
 * @param[in]  hAna      hades analysis handle
 * @param[in]  input     Input buffer; nChannels x blocksize
 * @param[in]  nChannels Number of channels in input buffer
 * @param[in]  blocksize Number of samples in input buffer
 * @param[out] hPCon     hades parameter container handle
 * @param[out] hSCon     hades signal container handle
 */
void hades_analysis_apply(/* Input Arguments */
                          hades_analysis_handle const hAna,
                          float** input,
                          int nChannels,
                          int blocksize,
                          /* Output Arguments */
                          void* const hPCon,
                          void* const hSCon);

/**
 * Returns a pointer to the frequency vector (read-only)
 *
 * @param[in]  hAna   hades analysis handle
 * @param[out] nBands (&) Number of bands (set to NULL if not needed)
 * @returns pointer to freqVector (or NULL if hAna is not initialised);
 *          nBands x 1
 */
const float* hades_analysis_getFrequencyVectorPtr(/* Input Arguments */
                                                  hades_analysis_handle const hAna,
                                                  /* Output Arguments */
                                                  int* nBands);

/** Returns number of frequency bands (0 if hAna is not initialised) */
int hades_analysis_getNbands(hades_analysis_handle const hAna);

/**
 * Returns a pointer to the covariance matrix averaging scalar [0..1], which can
 * be changed at run-time
 *
 * @param[in] hAna hades analysis handle
 * @returns pointer to the covariance matrix averaging scalar (or NULL if hAna
 *          is not initialised); 1 x 1
 */
float* hades_analysis_getCovarianceAvagingCoeffPtr(hades_analysis_handle const hAna);

/**
 * Returns the analyser processing delay, in samples
 *
 * @note The total delay for an analyser -> synthesiser configuration is
 *       computed as:
 *           hades_analysis_getProcDelay() + hades_synthesis_getProcDelay().
 */
int hades_analysis_getProcDelay(hades_analysis_handle const hAna);


/* ========================================================================== */
/*                      Parameter and Signal Containers                       */
/* ========================================================================== */

/**
 * Creates an instance of a container used for storing the parameters
 * estimated by an analyser for one 'blocksize'
 *
 * @note There should be one container per analyser, but this container can be
 *       passed to multiple different synthesisers. You may also create multiple
 *       containers, fill them using an analyser, store them, and pass them to
 *       the synthesiser(s) later.
 *
 * @param[in] phPCon (&) address of hades parameter container handle
 * @param[in] hAna   hades analysis handle
 */
void hades_param_container_create(/* Input Arguments */
                                  hades_param_container_handle* const phPCon,
                                  hades_analysis_handle const hAna);

/**
 * Destroys an instance of a hades parameter container
 *
 * @param[in] phPCon (&) address of hades parameter container handle
 */
void hades_param_container_destroy(/* Input Arguments */
                                   hades_param_container_handle* const phPCon);

/**
 * Creates an instance of a container used for storing the TF-domain audio
 * returned by an analyser for one 'blocksize'
 *
 * @note There should be one container per analyser, but this container can be
 *       passed to multiple different synthesisers. You may also create multiple
 *       containers, fill them using an analyser, store them, and pass them to
 *       the synthesiser(s) later.
 *
 * @param[in] phSCon (&) address of hades signal container handle
 * @param[in] hAna   hades analysis handle
 */
void hades_signal_container_create(hades_signal_container_handle* const phSCon,
                                   hades_analysis_handle const hAna);

/**
 * Destroys an instance of a hades signal container
 *
 * @param[in] phSCon (&) address of hades signal container handle
 */
void hades_signal_container_destroy(/* Input Arguments */
                                    hades_signal_container_handle* const phSCon);

#endif /* SAF_ENABLE_HADES_MODULE */


#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */

#endif /* __SAF_HADES_ANALYSIS_H_INCLUDED__ */
