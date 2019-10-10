/*
 * Copyright 2018 Leo McCormack
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
 * Filename: saf_hoa_internal.h
 * ----------------------------
 * A collection of higher-order Ambisonics related functions. Some of which are
 * derived from the Matlab library by Archontis Politis, found here:
 *     https://github.com/polarch/Higher-Order-Ambisonics
 *
 * Dependencies:
 *     saf_utilities, saf_sh, saf_vbap
 * Author, date created:
 *     Leo McCormack, 19.03.2018
 */

#if defined(SAF_ENABLE_HOA) && !defined(__SAF_HOA_INTERNAL_H_INCLUDED__)
#define __SAF_HOA_INTERNAL_H_INCLUDED__

#include <stdio.h>
#include <math.h> 
#include <string.h>
#include "saf_hoa.h"
#include "../saf_sh/saf_sh.h" /* for computing spherical harmonics */
#include "../saf_vbap/saf_vbap.h" /* for VBAP gains utilised by AllRAD */
#include "../saf_utilities/saf_utilities.h" /* for linear algebra speed-ups */

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */
    
/* ========================================================================== */
/*                       Loudspeaker Ambisonic Decoders                       */
/* ========================================================================== */
    
/*
 * Function: getEPAD
 * -----------------
 * Computes the "Energy preserving Ambisonic decoder", as detailed in [1]
 * Note: the function has been written to also work when the number of spherical
 * harmonic components exceeds the number of loudspeakers. In which case, the
 * 'U' matrix from the SVD is truncated instead. However, ideally, nLS > nSH,
 * like in the paper
 *
 * Input Arguments:
 *     order       - decoding order
 *     ls_dirs_deg - loudspeaker directions in DEGREES [azi elev]; FLAT: nLS x 2
 *     nLS         - number of loudspeakers
 * Output Arguments:
 *     decMtx      - & decoding matrix; FLAT: nLS x (order+1)^2
 *
 * [1] Zotter, F., Pomberger, H., Noisternig, M. (2012). Energy-Preserving
 *     Ambisonic Decoding. Acta Acustica United with Acustica, 98(1), 37:47.
 */
void getEPAD(/* Input Arguments */
             int order,
             float* ls_dirs_deg,
             int nLS,
             /* Output Arguments */
             float** decMtx);

/*
 * Function: getAllRAD
 * -------------------
 * Computes the "All-round Ambisonics decoder", as detailed in [1].
 *
 * Input Arguments:
 *     order       - decoding order
 *     ls_dirs_deg - loudspeaker directions in DEGREES [azi elev]; FLAT: nLS x 2
 *     nLS         - number of loudspeakers
 * Output Arguments:
 *     decMtx      - & decoding matrix; FLAT: nLS x (order+1)^2
 *
 * [1] Zotter, F., Frank, M. (2012). All-Round Ambisonic Panning and Decoding.
 *     Journal of the Audio Engineering Society, 60(10), 807:820
 */
void getAllRAD(/* Input Arguments */
               int order,
               float* ls_dirs_deg,
               int nLS,
               /* Output Arguments */
               float **decMtx);

    
/* ========================================================================== */
/*                         Binaural Ambisonic Decoders                        */
/* ========================================================================== */

/*
 * Function: getBinDecoder_LS
 * --------------------------
 * Computes a standard least-squares (LS) binaural ambisonic decoder for each
 * frequency, ready to be applied to input SH signals in the time-frequency
 * domain, or, take the inverse-FFT and apply it via matrix convolution.
 * Note: This standard LS decoder typically exhibits strong timbral colourations
 * in the output when using lower-order input. This is due to input order
 * truncation, as the HRTF grid is typically of much higher modal order than
 * that of the input. This colouration especially affects high-frequencies,
 * since high-frequency energy is predominantly concentrated in the higher-order
 * components.
 *
 * Input Arguments:
 *     hrtfs         - the HRTFs; FLAT: N_bands x 2 x N_dirs
 *     hrtf_dirs_deg - HRTF directions; FLAT: N_dirs x 2
 *     N_dirs        - number of HRTF directions in set
 *     N_bands       - number of frequency bands/bins
 *     order         - decoding order
 *     weights       - integration weights (set to NULL if not available);
 *                     N_dirs x 1
 * Output Arguments:
 *     decMtx        - decoding matrix; FLAT: N_bands x 2 x (order+1)^2
 */
void getBinDecoder_LS(/* Input Arguments */
                      float_complex* hrtfs,
                      float* hrtf_dirs_deg,
                      int N_dirs,
                      int N_bands,
                      int order,
                      float* weights,
                      /* Output Arguments */
                      float_complex* decMtx);

/*
 * Function: getBinDecoder_LSDIFFEQ
 * --------------------------------
 * Computes a least-squares (LS) binaural ambisonic decoder for each
 * frequency with diffuse-field equalisation [1], ready to be applied to input
 * SH signals in the time-frequency domain, or, take the inverse-FFT and apply
 * it via matrix convolution.
 * Note: this equalisation mitagates some of the timbral colourations exhibited
 * by standard LS decoding at lower input orders.
 *
 * Input Arguments:
 *     hrtfs         - the HRTFs; FLAT: N_bands x 2 x N_dirs
 *     hrtf_dirs_deg - HRTF directions; FLAT: N_dirs x 2
 *     N_dirs        - number of HRTF directions in set
 *     N_bands       - number of frequency bands/bins
 *     order         - decoding order
 *     weights       - integration weights (set to NULL if not available);
 *                     N_dirs x 1
 * Output Arguments:
 *     decMtx        - decoding matrix; FLAT: N_bands x 2 x (order+1)^2
 *
 * [1] Z. Ben-Hur, F. Brinkmann, J. Sheaffer, S. Weinzierl, and B. Rafaely,
 *     "Spectral equalization in binaural signals represented by order-
 *     truncated spherical harmonics," The Jour- nal of the Acoustical Society
 *     of America, vol. 141, no. 6, pp. 4087–4096, 2017.
 */
void getBinDecoder_LSDIFFEQ(/* Input Arguments */
                            float_complex* hrtfs,
                            float* hrtf_dirs_deg,
                            int N_dirs,
                            int N_bands,
                            int order,
                            float* weights,
                            /* Output Arguments */
                            float_complex* decMtx);

/*
 * Function: getBinDecoder_SPR
 * ---------------------------
 * Computes a binaural ambisonic decoder based on spatial resampling (aka:
 * virtual loudspeaker decoding) [1]
 * Note: like "getBinDecoder_LSDIFFEQ" this method mitagates some of the timbral
 * colourations exhibited by standard LS decoding at lower input orders.
 * However, it operates without equalisation. Instead, the modal order of the
 * HRTF grid is brought closer to the decoding order, by simply reducing the
 * number of HRTF points used, and calculating the LS decoder with this reduced
 * number of HRTFs. Therefore, rather than assigning high-frequency energy to
 * higher-order components and subsequently discarding it, due to order
 * truncation, the energy is instead aliased back into the lower-order
 * components and preserved.
 *
 * Input Arguments:
 *     hrtfs         - the HRTFs; FLAT: N_bands x 2 x N_dirs
 *     hrtf_dirs_deg - HRTF directions; FLAT: N_dirs x 2
 *     N_dirs        - number of HRTF directions in set
 *     N_bands       - number of frequency bands/bins
 *     order         - decoding order
 *     weights       - integration weights (set to NULL if not available);
 *                     N_dirs x 1
 * Output Arguments:
 *     decMtx        - decoding matrix; FLAT: N_bands x 2 x (order+1)^2
 *
 * [1] B. Bernschu ̈tz, A. V. Giner, C. Pörschmann, and J. Arend, “Binaural
 *     reproduction of plane waves with reduced modal order,” Acta Acustica
 *     united with Acustica, vol. 100, no. 5, pp. 972–983, 2014.
 */
void getBinDecoder_SPR(/* Input Arguments */
                       float_complex* hrtfs,
                       float* hrtf_dirs_deg,
                       int N_dirs,
                       int N_bands,
                       int order,
                       float* weights,
                       /* Output Arguments */
                       float_complex* decMtx);

/*
 * Function: getBinDecoder_TA
 * --------------------------
 * Computes a binaural ambisonic decoder based on the time-alignment (TA) method
 * described in [1].
 * Note: Since the standard LS decoder is unable to sufficiently fit lower-order
 * spherical harmonics to the highly directive HRTF patterns, this approach
 * addresses this by conducting preliminary time-alignment of the Head-related
 * impulse responses (HRIRs), which aids the LS fitting. This method essentially
 * exploits prior knowledge of the bandwidth in which the inter-aural level
 * differences (ILDs) are the dominant localisation cues; which is above
 * approximately 1.5 kHz. By discarding the phase information of the HRTFs at
 * frequencies above 1.5 kHz, the LS fitting instead prioritises the delivery of
 * the correct magnitude responses; rather than the phase. Thus it ultimately
 * yields improved ILD cues and diminished inter-aural time difference (ITD)
 * cues; but in a frequency range where ILD cues are more important for
 * localisation. This method, therefore, mitagates some of the localisation
 * deficiencies compared with the standard LS decoding at lower input orders.
 * Further note: The paper also detailed a diffuse-field covariance contraint,
 * and the original name was TAC (C=contrained), however, in this framework,
 * this constraint is a seperate independent operation. One may impose it on any
 * binaural decoder using the "applyDiffCovMatching" function.
 *
 * Input Arguments:
 *     hrtfs         - the HRTFs; FLAT: N_bands x 2 x N_dirs
 *     hrtf_dirs_deg - HRTF directions; FLAT: N_dirs x 2
 *     N_dirs        - number of HRTF directions in set
 *     N_bands       - number of frequency bands/bins
 *     order         - decoding order
 *     freqVector    - frequency vector; N_bands x 1
 *     itd_s         - interaural time differences (ITDs), seconds; N_dirs x 1
 *     weights       - integration weights (set to NULL if not available);
 *                     N_dirs x 1
 * Output Arguments:
 *     decMtx        - decoding matrix; FLAT: N_bands x 2 x (order+1)^2
 *
 * [1] Zaunschirm M, Schörkhuber C, Höldrich R. Binaural rendering of
 *     Ambisonic signals by head-related impulse response time alignment and
 *     a diffuseness constraint. The Journal of the Acoustical Society of
 *     America. 2018 Jun 19;143(6):3616-27
 */
void getBinDecoder_TA(/* Input Arguments */
                      float_complex* hrtfs,
                      float* hrtf_dirs_deg,
                      int N_dirs,
                      int N_bands,
                      int order,
                      float* freqVector,
                      float* itd_s,
                      float* weights,
                      /* Output Arguments */
                      float_complex* decMtx);

/*
 * Function: getBinDecoder_MAGLS
 * -----------------------------
 * Computes a binaural ambisonic decoder based on the magnitude least-squares
 * (MagLS) method described in [1].
 * Note: Mag-LS operates under the same principles held by the TA/TAC decoder,
 * differing in the manner in which the phase is neglected at frequencies above
 * 1.5kHz.
 *
 * Input Arguments:
 *     hrtfs         - the HRTFs; FLAT: N_bands x 2 x N_dirs
 *     hrtf_dirs_deg - HRTF directions; FLAT: N_dirs x 2
 *     N_dirs        - number of HRTF directions in set
 *     N_bands       - number of frequency bands/bins
 *     order         - decoding order
 *     freqVector    - frequency vector; N_bands x 1
 *     weights       - integration weights (set to NULL if not available);
 *                     N_dirs x 1
 * Output Arguments:
 *     decMtx        - decoding matrix; FLAT: N_bands x 2 x (order+1)^2
 *
 * [1] Schörkhuber C, Zaunschirm M, Höldrich R. Binaural Rendering of
 *     Ambisonic Signals via Magnitude Least Squares. InProceedings of the
 *     DAGA 2018 (Vol. 44, pp. 339-342).
 */
void getBinDecoder_MAGLS(/* Input Arguments */
                         float_complex* hrtfs,
                         float* hrtf_dirs_deg,
                         int N_dirs,
                         int N_bands,
                         int order,
                         float* freqVector,
                         float* weights,
                         /* Output Arguments */
                         float_complex* decMtx);


#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */

#endif /* defined(SAF_ENABLE_HOA) && !defined(__SAF_HOA_INTERNAL_H_INCLUDED__) */
