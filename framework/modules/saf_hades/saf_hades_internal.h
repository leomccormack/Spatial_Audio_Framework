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
 * @file saf_hades_internal.h
 * @ingroup HADES
 * @brief Internal header for the HADES module (#SAF_HADES_MODULE)
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

#ifndef __SAF_HADES_INTERNAL_H_INCLUDED__
#define __SAF_HADES_INTERNAL_H_INCLUDED__

#include "saf_hades_analysis.h"
#include "saf_hades_synthesis.h"
#include "saf.h"
#include "saf_externals.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#ifdef  SAF_ENABLE_HADES_MODULE

/** Maximum supported blocksize */
#define HADES_MAX_BLOCKSIZE ( 4096 )

/** Helper struct for averaging covariance matrices (block-wise) */
typedef struct _CxMic{
    float_complex Cx[HADES_MAX_NMICS*HADES_MAX_NMICS];
}CxMic;

/* ========================================================================== */
/*                           Main Internal Structs                            */
/* ========================================================================== */

/** Main structure for hades analysis */
typedef struct _hades_analysis_data
{
    /* User parameters (defined at intialisation stage) */
    float fs;                             /**< Host samplerate, Hz */
    HADES_FILTERBANKS fbOpt;              /**< see #HADES_FILTERBANKS */
    int hopsize;                          /**< Filterbank hop size (blocksize must be divisable by this */
    int blocksize;                        /**< Number of samples to process at a time (note that 1 doa and diffuseness estimate is made per block) */
    int hybridmode;                       /**< Optionally, the lowest TF bands may be subdivided to improve low-freq resolution. */
    float* h_array;                       /**< Array impulse responses; FLAT: nGrid x nMics x h_len */
    float* grid_dirs_deg;                 /**< Array grid dirs in degrees; FLAT: nGrid x 2 */
    int nGrid;                            /**< Number of grid/scanning directions */
    int nMics;                            /**< Number of microphones */
    int h_len;                            /**< Length of impulse responses, in samples */
    HADES_DIFFUSENESS_ESTIMATORS diffOpt; /**< see #HADES_DIFFUSENESS_ESTIMATORS */
    HADES_DOA_ESTIMATORS doaOpt;          /**< see #HADES_DOA_ESTIMATORS */

    /* Optional user parameters (that can also be manipulated at run-time) */
    float covAvgCoeff;                    /**< Temporal averaging coefficient [0 1] */

    /* Time-frequency transform and array data */
    void* hFB_enc;                        /**< Time-frequency transform handle */
    int nBands;                           /**< Number of frequency bands */
    int timeSlots;                        /**< Number of time slots */
    int filterbankDelay;                  /**< Filterbank delay, in time-domain samples */
    float* freqVector;                    /**< Centre frequencies; nBands x 1 */
    float_complex* DCM_array;             /**< Diffuse covariance matrix (computed over all grid directions and weighted); FLAT: nBands x nMics x nMics */
    float_complex* H_array;               /**< Array IRs in the frequency domain; FLAT: nBands x nMics x nDirs */
    float_complex* H_array_w;             /**< Array IRs in the frequency domain spatially weightend; FLAT: nBands x nMics x nDirs */

    /* DoA and diffuseness estimator data */
    void* hEig;                           /**< handle for the eigen solver */
    float_complex** T;                    /**< for covariance whitening; nBands x (nMics x nMics) */
    void* hDoA;                           /**< DoA estimator handle */
    float* grid_dirs_xyz;                 /**< Scanning grid coordinates (unit vectors and only used by grid-based estimators); FLAT: nGrid x 3 */
    float_complex* W;                     /**< Diffuse integration weighting matrix; FLAT: nGrid x nGrid */

    /* Run-time variables */
    float** inputBlock;                   /**< Input frame; nMics x blocksize */
    CxMic* Cx;                            /**< Current (time-averaged) covariance matrix per band; nBands x 1 */
    float_complex* V;                     /**< Eigen vectors; FLAT: nMics x nMics */
    float_complex* Vn;                    /**< Noise subspace; FLAT: nMics x (nMics-1) */
    float* lambda;                        /**< Eigenvalues; nMics x 1 */

}hades_analysis_data;

/** Main structure for hades synthesis */
typedef struct _hades_synthesis_data
{
    /* User parameters */
    HADES_BEAMFORMER_TYPE beamOption; /**< see #HADES_BEAMFORMER_TYPE */
    int enableCM;                     /**< Flag: whether the spatial covariance matching is enabled (1) or disabled (0) */
    hades_binaural_config* binConfig; /**< Internal copy of user configuration */
    int refIndices[2];                /**< Indices into [0 nMics-1], defining the reference sensors */
    HADES_HRTF_INTERP_OPTIONS interpOption; /**< HRIR interpolation option, see #HADES_HRTF_INTERP_OPTIONS */

    /* Optional user parameters (that can also be manipulated at run-time) */
    float* eq;                       /**< Gain factor per band; nBands x 1 */
    float* streamBalance;            /**< Stream balance per band (0:fully diffuse, 1:balanced, 2:fully direct); nBands x 1 */
    float synAvgCoeff;               /**< Mixing matrix averaging coefficent [0..1] */

    /* Things relevant to the synthesiser, which are copied from the hades_analysis_create() to keep everything aligned */
    HADES_FILTERBANKS fbOpt;         /**< Filterbank option, see #HADES_FILTERBANKS */
    int nBands;                      /**< Number of bands in the time-frequency transform domain */
    int hopsize;                     /**< hopsize in samples */
    int blocksize;                   /**< blocksize in samples */
    int nGrid;                       /**< Number of grid/scanning directions */
    int nMics;                       /**< Number of microphones */
    float_complex* H_array;          /**< Array IRs in the frequency domain; FLAT: nBands x nMics x nGrid */
    float* grid_dirs_deg;            /**< Array grid dirs in degrees; FLAT: nGrid x 2 */
    float** grid_dirs_xyz;           /**< Grid dirs as Cartesian coordinates of unit length; nGrid x 3 */
    int timeSlots;                   /**< Number of time frames in the time-frequency transform domain */ 
    float* freqVector;               /**< Frequency vector (band centre frequencies); nBands x 1 */
    float_complex* DCM_array;        /**< Diffuse coherence matrix for the array; FLAT: nBands x nMics x nMics */
    float_complex* W;                /**< Diffuse integration weighting matrix; FLAT: nGrid x nGrid */
 
    /* Time-frequency transform */
    void* hFB_dec;                   /**< Filterbank handle */

    /* HRTF and diffuse rendering variables */
    float_complex* H_bin;            /**< To spatialise the source beamformers; FLAT: nBands x #NUM_EARS x nGrid */
    float_complex* DCM_bin_norm;     /**< Diffuse coherence matrix for the HRTF set, normalised with 1/trace(DCM_bin); FLAT: nBands x nMics x nMics */
    float* diffEQ;                   /**< EQ curve to bring the overall diffuse-field magnitude response of the array to that of the HRTFs instead; nBands x 1 */

    /* Run-time variables */
    void* hPinv;                     /**< Handle for computing the Moore-Penrose pseudo inverse */
    void* hLinSolve;                 /**< Handle for solving linear equations (Ax=b) */
    void* hCDF;                      /**< Handle for solving the covariance matching problem */
    float_complex* As;               /**< Array steering vector for DoA; FLAT: nMics x 1 */
    float_complex* As_l;             /**< Array steering vector relative to left reference sensor; FLAT: nMics x 1 */
    float_complex* As_r;             /**< Array steering vector relative to right reference sensor; FLAT: nMics x 1 */
    float_complex* Q_diff;           /**< Mixing matrix for the diffuse stream; FLAT: #NUM_EARS x nMics */
    float_complex* Q_dir;            /**< Mixing matrix for the direct stream; FLAT: #NUM_EARS x nMics */
    float_complex* Q;                /**< Mixing matrix for the direct and diffuse streams combined (based on the diffuseness value); FLAT: #NUM_EARS x nMics */
    float_complex* Cy;               /**< Target binaural spatial covariance matrix; FLAT: #NUM_EARS x #NUM_EARS */
    float_complex* new_M;            /**< New mixing matrix (not yet temporally averaged); FLAT: #NUM_EARS x nMics */
    float_complex** M;               /**< Mixing matrix per band; nBands x FLAT: (#NUM_EARS x nMics) */

    /* Run-time audio buffers */
    float_complex*** outTF;          /**< nBands x #NUM_EARS x timeSlots */
    float** outTD;                   /**< output time-domain buffer; #NUM_EARS x blocksize */
 
}hades_synthesis_data;

/** Parameter container to store the data from an analyser for one blocksize of audio */
typedef struct _hades_param_container_data {
    int nBands;                      /**< Number of bands */

    /* Estimated Parameters */
    float* diffuseness;              /**< Diffuseness value per band; nBands x 1 */
    int* doa_idx;                    /**< Beamforming direction index per band; nBands x 1 */
    int* gains_idx;                  /**< Reproduction direction index per band; nBands x 1 */

    /* Optional parameters */
    float* gains_dir;                /**< Extra direct reproduction gain per band (default=1.0f); nBands x 1  */
    float* gains_diff;               /**< Extra diffuse reproduction gain per band (default=1.0f); nBands x 1  */

} hades_param_container_data;

/** Main structure for hades radial (360degree) gain and direct-to-diffuse ratio editor */
typedef struct _hades_radial_editor_data {
    int nBands;                      /**< Number of bands */
    int nGrid;                       /**< Number of grid/scanning directions */
    float* pGrid_dirs_deg;           /**< Pointer to grid dirs in degrees; FLAT: nGrid x 2 */
    float* pGrid_dirs_xyz;           /**< Pointer to grid dirs as Cartesian coordinates of unit length; FLAT: nGrid x 3 */

} hades_radial_editor_data;

/** Signal container to store one block of TF-domain audio data */
typedef struct _hades_signal_container_data {
    int nMics;                       /**< Number of spherical harmonic components */
    int nBands;                      /**< Number of bands in the time-frequency transform */
    int timeSlots;                   /**< Number of time frames in time-frequency transform */

    /* Covariance matrices and signal statistics computed during the analysis */
    CxMic* Cx;                       /**< NON-time-averaged covariance matrix per band; nBands x .Cx(nMics x nMics) */

    /* TF frame to carry over to a decoder */
    float_complex*** inTF;           /**< Input frame in TF-domain; nBands x nMics x timeSlots */

} hades_signal_container_data;


/* ========================================================================== */
/*                             Internal Functions                             */
/* ========================================================================== */

/**
 * Binaural filter interpolator
 *
 * @param[in]  hAna            hades analysis handle
 * @param[in]  interpOption    see #HADES_HRTF_INTERP_OPTIONS
 * @param[in]  binConfig       Binaural configuration
 * @param[in]  target_dirs_deg Target/interpolation dirs, in degrees;
 *                             FLAT: nTargetDirs x 2
 * @param[in]  nTargetDirs     Number of target/interpolation directions
 * @param[out] hrtf_interp     The interpolated HRTFs;
 *                             nBands x #NUM_EARS x nTargetDirs
 */
void hades_getInterpolatedHRTFs(/* Input Arguments */
                                hades_analysis_handle const hAna,
                                HADES_HRTF_INTERP_OPTIONS interpOption,
                                hades_binaural_config* binConfig,
                                float* target_dirs_deg,
                                int nTargetDirs,
                                /* Output Arguments */
                                float_complex* hrtf_interp);

/**
 * Creates an instance of the space-domain MUSIC implementation
 *
 * @param[in] phMUSIC       (&) address of the sdMUSIC handle
 * @param[in] nMics         Number of microphones in the array
 * @param[in] grid_dirs_deg Scanning grid directions; FLAT: nDirs x 2
 * @param[in] nDirs         Number of scanning directions
 */
void hades_sdMUSIC_create(void ** const phMUSIC,
                          int nMics,
                          float* grid_dirs_deg,
                          int nDirs);

/**
 * Destroys an instance of the spherical harmonic domain MUSIC implementation,
 * which may be used for computing pseudo-spectrums for visualisation/DoA
 * estimation purposes
 *
 * @param[in] phMUSIC    (&) address of the sdMUSIC handle
 */
void hades_sdMUSIC_destroy(void ** const phMUSIC);

/**
 * Computes a pseudo-spectrum based on the MUSIC algorithm optionally returning
 * the grid indices corresponding to the N highest peaks (N=nSrcs)
 *
 * @warning The number of sources should not exceed: floor(nMics/2)!
 *
 * @param[in] hMUSIC    sdMUSIC handle
 * @param[in] A_grid    Scanning steering vectors; nMics x nGrid
 * @param[in] Vn        Noise subspace; FLAT: nSH x (nSH - nSrcs)
 * @param[in] nSrcs     Number of sources
 * @param[in] P_music   Pseudo-spectrum (set to NULL if not wanted); nDirs x 1
 * @param[in] peak_inds Indices corresponding to the "nSrcs" highest peaks in
 *                      the pseudo-spectrum (set to NULL if not wanted);
 *                      nSrcs x 1
 */
void hades_sdMUSIC_compute(/* Input arguments */
                           void* const hMUSIC,
                           float_complex* A_grid,
                           float_complex* Vn,
                           int nSrcs,
                           /* Output arguments */
                           float* P_music,
                           int* peak_inds);

/**
 * Returns an estimate of the diffuseness, based on [1]
 *
 * @param[in] lambda Eigenvalues; N x 1
 * @param[in] N      Number of eigenvalues
 * @returns an estimate of the diffuseness
 *
 * @see [1] Epain, N. and Jin, C.T., 2016. Spherical harmonic signal covariance
 *          and sound field diffuseness. IEEE/ACM Transactions on Audio, Speech,
 *          and Language Processing, 24(10), pp.1796-1807.
 */
float hades_comedie(float* lambda,
                    int N);

#endif /* SAF_ENABLE_HADES_MODULE */


#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */

#endif /* __SAF_HADES_INTERNAL_H_INCLUDED__ */
