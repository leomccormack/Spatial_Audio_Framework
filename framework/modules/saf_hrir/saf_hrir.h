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
 *@addtogroup HRIR
 *@{
 * @file saf_hrir.h
 * @brief Main header for the HRIR/HRTF processing module (#SAF_HRIR_MODULE)
 *
 * A collection functions for processing head-related impulse-response (HRIR).
 * Including: estimation of the interaural time differences (ITDs), conversion
 * of HRIRs to HRTFs or filterbank coefficients; diffuse-field equalisation and
 * phase simplication; and HRTF interpolation.
 *
 * @author Leo McCormack
 * @date 12.12.2016
 * @license ISC
 */

#ifndef __SAF_HRIR_H_INCLUDED__
#define __SAF_HRIR_H_INCLUDED__

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */
    
#include "../saf_utilities/saf_utility_complex.h"
    
/* ========================================================================== */
/*                               Default HRIRs                                */
/* ========================================================================== */
/* The default HRIR dataset is a Genelec Aural ID of a KEMAR Dummy Head (@48kHz)
 * Kindly provided by Aki Mäkivirta and Jaan Johansson */

/** The default HRIR data for SAF */
extern const float __default_hrirs[836][2][256];

/** The measurement directions used for the default HRIR dataset */
extern const float __default_hrir_dirs_deg[836][2];

/** The number of directions/measurements in the default HRIR dataset */
extern const int __default_N_hrir_dirs;

/** The length of the filters, in samples, for the default HRIR dataset */
extern const int __default_hrir_len;

/** The samplerate used to measure the default HRIR filters  */
extern const int __default_hrir_fs;


/* ========================================================================== */
/*                               Main Functions                               */
/* ========================================================================== */

/**
 * Estimates the interaural time-differences (ITDs) for each HRIR based on the
 * cross-correlation between the left and right channels, which are first
 * low-pass filtered at 750Hz
 *
 * @param[in]  hrirs    HRIRs; FLAT: N_dirs x #NUM_EARS x hrir_len
 * @param[in]  N_dirs   Number of HRIRs
 * @param[in]  hrir_len Length of the HRIRs in samples
 * @param[in]  fs       Sampling rate of the HRIRs
 * @param[out] itds_s   ITDs in seconds; N_dirs x 1
 */
void estimateITDs(/* Input Arguments */
                  float* hrirs,
                  int N_dirs,
                  int hrir_len,
                  int fs,
                  /* Output Arguments */
                  float* itds_s);

/**
 * Passes zero padded HRIRs through the afSTFT filterbank
 *
 * The filterbank coefficients are then normalised with the energy of an
 * impulse, which is centered at approximately the beginning of the median HRIR
 * peak.
 *
 * @warning This function is NOT suitable for binaural room impulse responses
 *          (BRIRs)!
 *
 * @param[in]  hrirs      HRIRs; FLAT: N_dirs x #NUM_EARS x hrir_len
 * @param[in]  N_dirs     Number of HRIRs
 * @param[in]  hrir_len   Length of the HRIRs in samples
 * @param[in]  hopsize    Hop size in samples
 * @param[in]  LDmode     Low-Delay mode, 0:disabled, 1:enabled
 * @param[in]  hybridmode Hybrid-filtering, 0:disabled, 1:enabled
 * @param[out] hrtf_fb    HRTFs as filterbank coeffs; FLAT:
 *                        (hybrid ? hopsize+5 : hopsize+1) x #NUM_EARS x N_dirs
 */
void HRIRs2HRTFs_afSTFT(/* Input Arguments */
                        float* hrirs,
                        int N_dirs,
                        int hrir_len,
                        int hopsize,
                        int LDmode,
                        int hybridmode,
                        /* Output Arguments */
                        float_complex* hrtf_fb);

/**
 * Passes zero padded HRIRs through the qmf filterbank
 *
 * The filterbank coefficients are then normalised with the energy of an
 * impulse, which is centered at approximately the beginning of the median HRIR
 * peak.
 *
 * @warning This function is NOT suitable for binaural room impulse responses
 *          (BRIRs)!
 *
 * @param[in]  hrirs      HRIRs; FLAT: N_dirs x #NUM_EARS x hrir_len
 * @param[in]  N_dirs     Number of HRIRs
 * @param[in]  hrir_len   Length of the HRIRs in samples
 * @param[in]  hopsize    Hop size in samples
 * @param[in]  hybridmode 0:disabled, 1:enabled
 * @param[out] hrtf_fb    HRTFs as filterbank coeffs;
 *                        FLAT:
 *                        (hybrid ? hopsize+7 : hopsize+1) x #NUM_EARS x N_dirs
 */
void HRIRs2HRTFs_qmf(/* Input Arguments */
                     float* hrirs,
                     int N_dirs,
                     int hrir_len,
                     int hopsize,
                     int hybridmode,
                     /* Output Arguments */
                     float_complex* hrtf_fb);

/**
 * Converts HRIRs to HRTFs for a given FFT size
 *
 * @note If the HRIRs are shorter than the FFT size (hrir_len<fftSize), then the
 *       HRIRs are zero-padded. If they are longer, then they are truncated.
 *
 * @param[in]  hrirs    HRIRs; FLAT: N_dirs x #NUM_EARS x hrir_len
 * @param[in]  N_dirs   Number of HRIRs
 * @param[in]  hrir_len Length of the HRIRs in samples
 * @param[in]  fftSize  FFT size
 * @param[out] hrtfs    HRTFs; FLAT: (fftSize/2+1) x #NUM_EARS x N_dirs
 */
void HRIRs2HRTFs(/* Input Arguments */
                 float* hrirs,
                 int N_dirs,
                 int hrir_len,
                 int fftSize,
                 /* Output Arguments */
                 float_complex* hrtfs);

/**
 * Applies pre-processing to a set of HRTFs, which can either be diffuse-field
 * EQ of an (optionally weighted) average of all HRTFs (CTF), phase
 * simplification based on ITDs, or both.
 *
 * @note 'weights' (if used) should sum to 4pi. 'itds_s' and 'centreFreq' are
 *       only required if applyPhase==1, so can be set to NULL otherwise.
 * @warning This function is NOT suitable for binaural room impulse responses
 *          (BRIRs)!
 *
 * @param[in]     N_dirs     Number of HRTFs
 * @param[in]     itds_s     HRIR ITDs (set to NULL if not needed); N_dirs x 1
 * @param[in]     centreFreq Frequency vector (set to NULL if not needed);
 *                           N_bands x 1
 * @param[in]     N_bands    Number of frequency bands/bins
 * @param[in]     weights    Grid weights (set to NULL if not available);
 *                           N_dirs x 1
 * @param[in]     applyEQ    Diffuse-field EQ / CTF; 0:disabled, 1:enabled
 * @param[in]     applyPhase Phase simplification; 0:disabled, 1:enabled
 * @param[in,out] hrtfs      The HRTFs; FLAT: N_bands x #NUM_EARS x N_dirs
 */
void diffuseFieldEqualiseHRTFs(/* Input Arguments */
                               int N_dirs,
                               float* itds_s,
                               float* centreFreq,
                               int N_bands,
                               float* weights,
                               int applyEQ,
                               int applyPhase,
                               /* Input/Output Arguments */
                               float_complex* hrtfs);

/**
 * Interpolates a set of HRTFs based on a specified interpolation table
 *
 * @note For 'interp_table' you can use e.g. VBAPgainTable2InterpTable() to take
 *       a conventional energy-normalised VBAP gain table, and convert it to an
 *       amplitude-normalised interpolation table. Note: amplitude-normalised
 *       VBAP gains are the same as triangular interpolation weights.
 * @note If itds!=NULL && freqVector!=NULL, then the interpolation is performed
 *       by applying interpolation gains to the HRTF magnitudes and HRIR inter-
 *       aural time differences (ITDs) separately. The inter-aural phase
 *       differences (IPDs) are then reintroduced for each frequency band.
 *       If itds==NULL || freqVector==NULL then the interpolatation is applied
 *       directly on the complex spectra.
 * @warning This function is NOT suitable for binaural room impulse responses
 *          (BRIRs)!
 *
 * @param[in]  hrtfs         HRTFs as filterbank coeffs;
 *                           FLAT: N_bands x #NUM_EARS x N_hrtf_dirs
 * @param[in]  itds          The inter-aural time difference (ITD) for each
 *                           HRIR (set to NULL if you do not want phase
 *                           simplication to be applied); N_hrtf_dirs x 1
 * @param[in]  freqVector    Frequency vector (set to NULL if you do not want
 *                           phase simplication to be applied); N_bands x 1
 * @param[in]  interp_table  Amplitude-Normalised VBAP gain table;
 *                           FLAT: N_interp_dirs x N_hrtf_dirs
 * @param[in]  N_hrtf_dirs   Number of HRTF directions
 * @param[in]  N_bands       Number of frequency bands
 * @param[in]  N_interp_dirs Number of interpolated hrtf positions
 * @param[out] hrtf_interp   interpolated HRTFs;
 *                           FLAT: N_bands x #NUM_EARS x N_interp_dirs
 */
void interpHRTFs(/* Input Arguments */
                 float_complex* hrtfs,
                 float* itds,
                 float* freqVector,
                 float* interp_table,
                 int N_hrtf_dirs,
                 int N_bands,
                 int N_interp_dirs,
                 /* Output Arguments */
                 float_complex* hrtf_interp);

/**
 * Computes the binaural diffuse coherence per frequency for a given HRTF set,
 * as described in [1]
 *
 * @param[in]  hrtfs       HRTFs as filterbank coeffs
 *                         FLAT: N_bands x #NUM_EARS x N_hrtf_dirs
 * @param[in]  itds        The inter-aural time difference (ITD) for each HRIR;
 *                         N_hrtf_dirs x 1
 * @param[in]  freqVector  Frequency vector; N_bands x 1
 * @param[in]  N_hrtf_dirs Number of HRTF directions
 * @param[in]  N_bands     Number of frequency bands
 * @param[out] HRTFcoh     Binaural coherence per frequency; N_bands x 1
 *
 * @see [1] A. Politis, "Diffuse-field coherence of sensors with arbitrary
 *          directional responses," arXiv preprint arXiv:1608.07713,2016.
 */
void binauralDiffuseCoherence(/* Input Arguments */
                              float_complex* hrtfs,
                              float* itds,
                              float* freqVector,
                              int N_hrtf_dirs,
                              int N_bands,
                              /* Output Arguments */
                              float* HRTFcoh);

/**
 * Resamples a set of HRIRs from its original samplerate to a new samplerate
 *
 * @note The included speex_resampler.h is used by default. If SAF_USE_INTEL_IPP
 *       is defined, then the IPP resampler is employed instead.
 *
 * @param[in]  hrirs_in      Input HRIRs;
 *                           FLAT: hrirs_N_dirs x #NUM_EARS x hrirs_in_len
 * @param[in]  hrirs_N_dirs  Number of HRIRs
 * @param[in]  hrirs_in_len  Length of input HRIRs, in samples
 * @param[in]  hrirs_in_fs   Original sampling rate, in Hz
 * @param[in]  hrirs_out_fs  New sampling rate, in Hz
 * @param[in]  padToNextPow2 1: length of output HRIRs padded to next 2^x, 0: no
 * @param[out] hrirs_out     Resampled HRIRs;
 *                           FLAT: hrirs_N_dirs x #NUM_EARS x hrirs_out_len
 * @param[out] hrirs_out_len (&) New HRIR length
 */
void resampleHRIRs(/* Input Arguments */
                   float* hrirs_in,
                   int hrirs_N_dirs,
                   int hrirs_in_len,
                   int hrirs_in_fs,
                   int hrirs_out_fs,
                   int padToNextPow2,
                   /* Output Arguments */
                   float** hrirs_out,
                   int* hrirs_out_len);


#ifdef __cplusplus
} /* extern "C" */
#endif  /* __cplusplus */

#endif /* __SAF_HRIR_H_INCLUDED__ */

/**@} */ /* doxygen addtogroup HRIR */
