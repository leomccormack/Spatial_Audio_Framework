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

/*
 * Filename: saf_hrir.h (include header)
 * -------------------------------------
 * A collection of head-related impulse-response (HRIR) functions. Including
 * estimation of the interaural time differences (ITDs), conversion of HRIRs to
 * HRTF filterbank coefficients, and HRTF interpolation utilising amplitude-
 * normalised VBAP gains.
 *
 * Dependencies:
 *     saf_utilities, afSTFTlib
 * Author, date created:
 *     Leo McCormack, 12.12.2016
 */

#ifndef __SAF_HRIR_H_INCLUDED__
#define __SAF_HRIR_H_INCLUDED__

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */
    
#include "../saf_utilities/saf_utilities.h"
    
/* misc. */
float matlab_fmodf(float x,  float y);
    
/* ========================================================================== */
/*                               Default HRIRs                                */
/* ========================================================================== */
    
/* Default HRIRs, measured at one of Aalto University's anechoic chambers, using
 * a KEMAR Dummy Head. 217 measurements in total @48kHz. */
extern const double __default_hrirs[217][2][1024];
extern const double __default_hrir_dirs_deg[217][2];
extern const int __default_N_hrir_dirs;
extern const int __default_hrir_len;
extern const int __default_hrir_fs;
    

/* ========================================================================== */
/*                               Main Functions                               */
/* ========================================================================== */
    
/*
 * Function: estimateITDs
 * ----------------------
 * Estimates the interaural time-differences (ITDs) for each HRIR in a set via
 * the cross-correlation between the left and right IRs
 *
 * Input Arguments:
 *     hrirs    - HRIRs; FLAT: N_dirs x 2 x hrir_len
 *     N_dirs   - number of HRIRs
 *     hrir_len - length of the HRIRs in samples
 *     fs       - sampling rate of the HRIRs
 * Output Arguments:
 *     itds_s   - ITDs in seconds; N_dirs x 1
 */
void estimateITDs(/* Input Arguments */
                  float* hrirs,
                  int N_dirs,
                  int hrir_len,
                  int fs,
                  /* Output Arguments */
                  float* itds_s);

/*
 * Function: HRIRs2FilterbankHRTFs
 * -------------------------------
 * Passes zero padded HRIRs through the afSTFT filterbank. The filterbank
 * coefficients are then normalised with the energy of an impulse, which is
 * centered at approximately the beginning of the HRIR peak.
 * Note: this function is NOT suitable for binaural room impulse responses
 * (BRIRs).
 * Further note: Currently, this is hard-coded for 128 hop size with hybrid mode
 * enabled. (133 bands in total)
 *
 * Input Arguments:
 *     hrirs    - HRIRs; FLAT: N_dirs x 2 x hrir_len
 *     N_dirs   - number of HRIRs
 *     hrir_len - length of the HRIRs in samples
 * Output Arguments:
 *     hrtf_fb  - HRTFs as filterbank coeffs; FLAT: 133 x 2 x N_dirs
 */
void HRIRs2FilterbankHRTFs(/* Input Arguments */
                           float* hrirs,
                           int N_dirs,
                           int hrir_len,
                           /* Output Arguments */
                           float_complex* hrtf_fb);

/*
 * Function: HRIRs2HRTFs
 * ---------------------
 * Converts a HRIR set to HRTFs, with a given FFT size.
 * Note: if the HRIRs are shorter than the FFT size (hrir_len<fftSize), then the
 * HRIRs are zero-padded. If they are longer, then they are truncated.
 *
 * Input Arguments:
 *     hrirs    - HRIRs; FLAT: N_dirs x 2 x hrir_len
 *     N_dirs   - number of HRIRs
 *     hrir_len - length of the HRIRs in samples
 *     fftSize  - FFT size
 * Output Arguments:
 *     hrtf     - HRTFs; FLAT: (fftSize/2+1) x 2 x N_dirs
 */
void HRIRs2HRTFs(/* Input Arguments */
                 float* hrirs,
                 int N_dirs,
                 int hrir_len,
                 int fftSize,
                 /* Output Arguments */
                 float_complex* hrtfs);

/*
 * Function: diffuseFieldEqualiseHRTFs
 * -----------------------------------
 * Applies diffuse-field equalisation to a set of HRTFs
 * Note: this function is NOT suitable for binaural room impulse responses
 * (BRIRs).
 *
 * Input Arguments:
 *     N_dirs      - number of HRTFs
 *     itds_s      - HRIR ITDs; N_dirs x 1
 *     centreFreq  - frequency vector; N_bands x 1
 *     N_bands     - number of frequency bands/bins
 * Input/Output Arguments:
 *     hrtfs       - the HRTFs; FLAT: N_bands x 2 x N_dirs
 */
void diffuseFieldEqualiseHRTFs(/* Input Arguments */
                               int N_dirs,
                               float* itds_s,
                               float* centreFreq,
                               int N_bands,
                               /* Input/Output Arguments */
                               float_complex* hrtfs);

/*
 * Function: interpFilterbankHRTFs
 * -------------------------------
 * Interpolates a set of HRTFs for specified directions; defined by an amplitude
 * normalised VBAP interpolation table (see saf_vbap.h). The interpolation is
 * performed by applying interpolation gains to the HRTF magnitudes and HRIR
 * inter-aural time differences separately. The inter-aural phase differences
 * are then reintroduced for each frequency band
 * Note: use "VBAPgainTable2InterpTable" to take a conventional energy-
 * normalised VBAP gain table, and convert it to an amplitude-normalised
 * interpolation table
 *
 * Input Arguments:
 *     hrtfs                - HRTFs as filterbank coeffs;
 *                            FLAT: N_bands x 2 x N_hrtf_dirs
 *     itds                 - the inter-aural time difference (ITD) for each
 *                            HRIR; N_hrtf_dirs x 1
 *     freqVector           - frequency vector; N_bands x 1
 *     vbap_gtable          - Amplitude-Normalised VBAP gain table;
 *                            FLAT: N_interp_dirs x N_hrtf_dirs
 *     N_hrtf_dirs          - number of HRTF directions
 *     N_bands              - number of frequency bands
 *     N_interp_dirs        - number of interpolated hrtf positions
 * Output Arguments:
 *     hrtf_interp          - pre-alloc, interpolated HRTFs;
 *                            FLAT: N_bands x 2 x N_interp_dirs
 */
void interpFilterbankHRTFs(/* Input Arguments */
                           float_complex* hrtfs,
                           float* itds,
                           float* freqVector,
                           float* vbap_gtable,
                           int N_hrtf_dirs,
                           int N_bands,
                           int N_interp_dirs,
                           /* Output Arguments */
                           float_complex* hrtf_interp);
   
/*
 * Function: binauralDiffuseCoherence
 * ----------------------------------
 * Computes the binaural diffuse coherence per frequency for a given HRTF set,
 * as in [1].
 *
 * Input Arguments:
 *     hrtfs       - HRTFs as filterbank coeffs; FLAT: N_bands x 2 x N_hrtf_dirs
 *     itds        - the inter-aural time difference (ITD) for each HRIR;
 *                   N_hrtf_dirs x 1
 *     freqVector  - frequency vector; N_bands x 1
 *     N_hrtf_dirs - number of HRTF directions
 *     N_bands     - number of frequency bands
 * Output Arguments:
 *     HRTFcoh     -  binaural coeherence per frequency; N_bands x 1
 *
 * [1] A. Politis, "Diffuse-field coherence of sensors with arbitrary
 *     directional responses," arXiv preprint arXiv:1608.07713,2016.
 */
void binauralDiffuseCoherence(/* Input Arguments */
                              float_complex* hrtfs,
                              float* itds,
                              float* freqVector,
                              int N_hrtf_dirs,
                              int N_bands,
                              /* Output Arguments */
                              float* HRTFcoh);
    

#ifdef __cplusplus
} /* extern "C" */
#endif  /* __cplusplus */

#endif /* __SAF_HRIR_H_INCLUDED__ */
