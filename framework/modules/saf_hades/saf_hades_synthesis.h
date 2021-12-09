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
 * @file saf_hades_synthesis.h
 * @ingroup HADES
 * @brief Header for the HADES synthesis (#SAF_HADES_MODULE)
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

#ifndef __SAF_HADES_SYNTHESIS_H_INCLUDED__
#define __SAF_HADES_SYNTHESIS_H_INCLUDED__

#include "saf_hades_analysis.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#ifdef  SAF_ENABLE_HADES_MODULE

/** Handle for the hades synthesis data */
typedef struct _hades_synthesis_data* hades_synthesis_handle;

/** Handle for the hades radial editor data */
typedef struct _hades_radial_editor_data* hades_radial_editor_handle;

/* ========================================================================== */
/*                   HADES Synthesis Configurations Options                   */
/* ========================================================================== */

/** Binaural configuration struct */
typedef struct _hades_binaural_config{
    int lHRIR;                  /**< Length of HRIRs in samples */
    int nHRIR;                  /**< Number of HRIRs */
    int hrir_fs;                /**< HRIR sample rate */
    float* hrirs;               /**< Matrix of HRIR data;
                                 *   FLAT: nHRIR x #NUM_EARS x lHRIR */
    float* hrir_dirs_deg;       /**< HRTF directions in [azimuth elevation]
                                 *   format, in degrees; FLAT: nHRIR x 2 */
}hades_binaural_config;

/** Beamforming options for hades_synthesis */
typedef enum {
    HADES_BEAMFORMER_NONE,            /**< No beamforming (ref sensors only) */
    HADES_BEAMFORMER_FILTER_AND_SUM,  /**< Filter-and-sum beamforming */
    HADES_BEAMFORMER_BMVDR            /**< Binaural minimum-variance distortion-
                                       *   less response (MVDR) beamforming */
}HADES_BEAMFORMER_TYPE;

/** HRTF interpolation options for hades_synthesis */
typedef enum {
    HADES_HRTF_INTERP_NEAREST,    /**< Quantise to nearest measurement */
    HADES_HRTF_INTERP_TRIANGULAR  /**< Triangular interpolation */
}HADES_HRTF_INTERP_OPTIONS;


/* ========================================================================== */
/*                            HADES Radial Editor                             */
/* ========================================================================== */

/**
 * Creates and returns a handle to an instance of a hades radial editor object,
 * which allows for direction-dependent (360degree) manipulation of gains
 *
 * @param[in] phREd (&) address of hades radial editor handle
 * @param[in] hAna  hades analysis handle
 */
void hades_radial_editor_create(/* Input Arguments */
                                hades_radial_editor_handle* const phREd,
                                hades_analysis_handle const hAna);

/**
 * Destroys an instance of a hades radial editor object
 *
 * @param[in] phREd (&) address of hades radial editor handle
 */
void hades_radial_editor_destroy(/* Input Arguments */
                                 hades_radial_editor_handle* const phREd);

/**
 * Applies the radial (360 degree) parameter editing
 *
 * @param[in] hREd         hades radial editor handle
 * @param[in] hPCon        hades parameter container handle
 * @param[in] dirGain_dB   Extra directional gains for the direct stream, in dB
 */
void hades_radial_editor_apply(/* Input Arguments */
                               hades_radial_editor_handle const hREd,
                               hades_param_container_handle  const hPCon,
                               float dirGain_dB[360]);


/* ========================================================================== */
/*                              HADES Synthesis                               */
/* ========================================================================== */

/**
 * Creates and returns a handle to an instance of a hades synthesis object
 *
 * @param[in] phSyn        (&) address of hades synthesis handle
 * @param[in] hAna         hades analysis handle
 * @param[in] beamOption   see #HADES_BEAMFORMER_TYPE
 * @param[in] enableCM     0: disabled, 1: enable covariance matching
 * @param[in] binConfig    Binaural configuration
 * @param[in] interpOption see #HADES_HRTF_INTERP_OPTIONS
 */
void hades_synthesis_create(/* Input Arguments */
                            hades_synthesis_handle* const phSyn,
                            hades_analysis_handle const hAna,
                            HADES_BEAMFORMER_TYPE beamOption,
                            int enableCM,
                            int refIndices[2],
                            hades_binaural_config* binConfig,
                            HADES_HRTF_INTERP_OPTIONS interpOption);

/**
 * Destroys an instance of hades synthesis
 *
 * @param[in] phSyn (&) address of hades synthesis handle
 */
void hades_synthesis_destroy(/* Input Arguments */
                             hades_synthesis_handle* const phSyn);

/**
 * Flushes run-time buffers with zeros
 *
 * Call this ONCE before calling hades_synthesis_apply()
 *
 * @param[in] hSyn hades synthesis handle
 */
void hades_synthesis_reset(hades_synthesis_handle const hSyn);

/**
 * Performs hades synthesis
 *
 * @note If nChannels is higher than the number required by the configuration,
 *       then these extra channels are zero'd. If there are too few, then
 *       the channels are truncated.
 *
 * @param[in]  hSyn      hades synthesis handle
 * @param[in]  hPCon     hades parameter container handle
 * @param[in]  hSCon     hades signal container handle
 * @param[in]  nChannels Number of channels in output buffer
 * @param[in]  blocksize Number of samples in output buffer
 * @param[out] output    Output buffer; nChannels x blocksize
 */
void hades_synthesis_apply(/* Input Arguments */
                           hades_synthesis_handle const hSyn,
                           hades_param_container_handle  const hPCon,
                           hades_signal_container_handle const hSCon,
                           int nChannels,
                           int blocksize,
                           /* Output Arguments */
                           float** output);

/**
 * Returns a pointer to the eq vector, which can be changed at run-time
 *
 * @param[in]  hSyn   hades synthesis handle
 * @param[out] nBands (&) Number of bands (set to NULL if not needed)
 * @returns pointer to the eq vector (or NULL if hSyn is not initialised);
 *          nBands x 1
 */
float* hades_synthesis_getEqPtr(hades_synthesis_handle const hSyn,
                                int* nBands);

/**
 * Returns a pointer to the stream balance vector [0..2], which can be changed
 * at run-time
 *
 * @param[in]  hSyn   hades synthesis handle
 * @param[out] nBands (&) Number of bands (set to NULL if not needed)
 * @returns pointer to the stream balance vector (or NULL if hSyn is not
 *          initialised); nBands x 1
 */
float* hades_synthesis_getStreamBalancePtr(hades_synthesis_handle const hSyn,
                                           int* nBands);

/**
 * Returns a pointer to the synthesis averaging coefficient scalar [0..1], which
 * can be changed at run-time
 *
 * @param[in]  hSyn   hades synthesis handle
 * @returns pointer to the mixing matrix averaging coeff scalar (or NULL if hSyn
 *          is not initialised); 1 x 1
 */
float* hades_synthesis_getSynthesisAveragingCoeffPtr(hades_synthesis_handle const hSyn);

/**
 * Returns the synthesiser processing delay, in samples
 *
 * @note This is not inclusive of the time-frequency transform delay, as you
 *       may get this using hades_analysis_getProcDelay(). The total delay is:
 *       hades_analysis_getProcDelay() + hades_synthesis_getProcDelay().
 */
int hades_synthesis_getProcDelay(hades_synthesis_handle const hSyn);

#endif /* SAF_ENABLE_HADES_MODULE */


#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */

#endif /* __SAF_HADES_SYNTHESIS_H_INCLUDED__ */
