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

/**
 * @file saf_hoa_internal.h
 * @ingroup HOA
 * @brief Internal header for the higher-order Ambisonics module
 *        (#SAF_HOA_MODULE)
 *
 * A collection of Ambisonics related functions. Many of which are derived from
 * the MATLAB library found in [1].
 *
 * @see [1] https://github.com/polarch/Higher-Order-Ambisonics
 *          Copyright (c) 2015, Archontis Politis, BSD-3-Clause License
 *
 * @author Leo McCormack
 * @date 19.03.2018
 * @license ISC
 */

#ifndef __SAF_HOA_INTERNAL_H_INCLUDED__
#define __SAF_HOA_INTERNAL_H_INCLUDED__

#include "saf_hoa.h"
#include "../saf_sh/saf_sh.h"       /* for computing the legendre polynamials */
#include "../saf_vbap/saf_vbap.h"        /* for VBAP gains utilised by AllRAD */
#include "../saf_utilities/saf_utilities.h"   /* for linear algebra speed-ups */
#include "saf_externals.h" 

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/* ========================================================================== */
/*                       Loudspeaker Ambisonic Decoders                       */
/* ========================================================================== */

/**
 * Computes the Energy preserving Ambisonic decoder (EPAD), as detailed in [1]
 *
 * @note The function has been written to also work when the number of spherical
 *       harmonic components exceeds the number of loudspeakers. In which case,
 *       the 'U' matrix from the SVD is truncated instead. However, ideally,
 *       nLS > nSH.
 * @note Additional scaling is applied so that when the loudspeakers are
 *       uniformly arranged, the decoding matrix gains are equivalent to those
 *       produced by SAD/MMD.
 *
 * @param[in]  order       Decoding order
 * @param[in]  ls_dirs_deg Loudspeaker directions in DEGREES [azi elev];
 *                         FLAT: nLS x 2
 * @param[in]  nLS         Number of loudspeakers
 * @param[out] decMtx      Decoding matrix; FLAT: nLS x (order+1)^2
 *
 * @see [1] Zotter, F., Pomberger, H., Noisternig, M. (2012). Energy-Preserving
 *          Ambisonic Decoding. Acta Acustica United with Acustica, 98(1), 37:47
 */
void getEPAD(/* Input Arguments */
             int order,
             float* ls_dirs_deg,
             int nLS,
             /* Output Arguments */
             float* decMtx);

/**
 * Computes the All-round Ambisonics decoder (AllRAD), as detailed in [1], which
 * is essentially a spherical harmonic approximation of VBAP patterns for the
 * target loudspeaker setup.
 *
 * @param[in]  order       Decoding order
 * @param[in]  ls_dirs_deg Loudspeaker directions in DEGREES [azi elev];
 *                         FLAT: nLS x 2
 * @param[in]  nLS         Number of loudspeakers
 * @param[out] decMtx      Decoding matrix; FLAT: nLS x (order+1)^2
 *
 * @see [1] Zotter, F., Frank, M. (2012). All-Round Ambisonic Panning and
 *          Decoding. Journal of the Audio Engineering Society, 60(10), 807:820
 */
void getAllRAD(/* Input Arguments */
               int order,
               float* ls_dirs_deg,
               int nLS,
               /* Output Arguments */
               float* decMtx);


/* ========================================================================== */
/*                         Binaural Ambisonic Decoders                        */
/* ========================================================================== */

/**
 * Computes a standard least-squares (LS) binaural ambisonic decoder
 *
 * The binaural ambisonic decoder is computed for each frequency bin/band,
 * ready to be applied to input SH signals in the time-frequency domain, or,
 * take the inverse-FFT and apply it via matrix convolution.
 *
 * @note This standard LS decoder typically produces strong timbral colourations
 *       in the output when using lower-order input. This is due to input order
 *       truncation, since the HRTF grid is typically of much higher modal order
 *       than that of the input order. This colouration especially affects high-
 *       frequencies, since high-frequency energy is predominantly concentrated
 *       in the higher-order components so is then lost by truncating the input
 *       order. This phenomenon therefore gets worse when increasing the number
 *       of HRTFs in the set.
 *
 * @param[in]  hrtfs         The HRTFs; FLAT: N_bands x #NUM_EARS x N_dirs
 * @param[in]  hrtf_dirs_deg HRTF directions; FLAT: N_dirs x 2
 * @param[in]  N_dirs        Number of HRTF directions in set
 * @param[in]  N_bands       Number of frequency bands/bins
 * @param[in]  order         Decoding order
 * @param[in]  weights       Integration weights (set to NULL if not available);
 *                           N_dirs x 1
 * @param[out] decMtx        Decoding matrix;
 *                           FLAT: N_bands x #NUM_EARS x (order+1)^2
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

/**
 * Computes a least-squares (LS) binaural ambisonic decoder with diffuse-field
 * equalisation [1]
 *
 * The binaural ambisonic decoder is computed for each frequency bin/band,
 * ready to be applied to input SH signals in the time-frequency domain, or,
 * take the inverse-FFT and apply it via matrix convolution.
 *
 * @note This equalisation mitagates some of the timbral colourations exhibited
 *       by standard LS decoding; especially at lower input orders.
 *
 * @param[in]  hrtfs         The HRTFs; FLAT: N_bands x #NUM_EARS x N_dirs
 * @param[in]  hrtf_dirs_deg HRTF directions; FLAT: N_dirs x 2
 * @param[in]  N_dirs        Number of HRTF directions in set
 * @param[in]  N_bands       Number of frequency bands/bins
 * @param[in]  order         Decoding order
 * @param[in]  weights       Integration weights (set to NULL if not available);
 *                           N_dirs x 1
 * @param[out] decMtx        Decoding matrix;
 *                           FLAT: N_bands x #NUM_EARS x (order+1)^2
 *
 * @see [1] Z. Ben-Hur, F. Brinkmann, J. Sheaffer, S. Weinzierl, and B. Rafaely,
 *          "Spectral equalization in binaural signals represented by order-
 *          truncated spherical harmonics," The Journal of the Acoustical
 *          Society of America, vol. 141, no. 6, pp. 4087--4096, 2017.
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

/**
 * Computes a binaural ambisonic decoder based on spatial resampling (i.e,
 * virtual loudspeaker decoding) [1]
 *
 * The binaural ambisonic decoder is computed for each frequency bin/band,
 * ready to be applied to input SH signals in the time-frequency domain, or,
 * take the inverse-FFT and apply it via matrix convolution.
 *
 * @note Like getBinDecoder_LSDIFFEQ() this method mitagates some of the timbral
 *       colourations exhibited by standard LS decoding at lower input orders.
 *       However, it operates without equalisation. Instead, the modal order of
 *       the HRTF grid is brought closer to the decoding order by simply
 *       reducing the number of HRTF points. The LS decoder is then computed
 *       using this reduced HRTF set. Therefore, rather than assigning high-
 *       frequency energy to higher-order components and subsequently discarding
 *       it due to order truncation, the energy is instead aliased back into the
 *       lower-order components and preserved.
 *
 * @param[in]  hrtfs         The HRTFs; FLAT: N_bands x #NUM_EARS x N_dirs
 * @param[in]  hrtf_dirs_deg HRTF directions; FLAT: N_dirs x 2
 * @param[in]  N_dirs        Number of HRTF directions in set
 * @param[in]  N_bands       Number of frequency bands/bins
 * @param[in]  order         Decoding order
 * @param[in]  weights       Integration weights (set to NULL if not available);
 *                           N_dirs x 1
 * @param[out] decMtx        Decoding matrix;
 *                           FLAT: N_bands x #NUM_EARS x (order+1)^2
 *
 * @see [1] B. Bernschu"tz, A. V. Giner, C. Po"rschmann, and J. Arend, "Binaural
 *          reproduction of plane waves with reduced modal order" Acta Acustica
 *          united with Acustica, vol. 100, no. 5, pp. 972--983, 2014.
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

/**
 * Computes a binaural ambisonic decoder based on the time-alignment (TA) method
 * described in [1]
 *
 * The binaural ambisonic decoder is computed for each frequency bin/band,
 * ready to be applied to input SH signals in the time-frequency domain, or,
 * take the inverse-FFT and apply it via matrix convolution.
 *
 * @note Since the standard LS decoder is unable to sufficiently fit lower-order
 *       spherical harmonics to the highly directive HRTF patterns, this
 *       approach addresses this by conducting preliminary time-alignment of
 *       the Head-related impulse responses (HRIRs), which aids the LS fitting.
 *       This method essentially exploits prior knowledge of the bandwidth in
 *       which the inter-aural level differences (ILDs) are the dominant
 *       localisation cues; which is above approximately 1.5 kHz. By discarding
 *       the phase information of the HRTFs at frequencies above 1.5 kHz, the LS
 *       fitting instead prioritises the delivery of the correct magnitude
 *       responses; rather than the phase. Thus it ultimately yields improved
 *       ILD cues and diminished inter-aural time difference (ITD) cues, but in
 *       a frequency range where ILD cues are more important for localisation.
 *       This method, therefore, mitagates many of the localisation deficiencies
 *       compared with the standard LS decoding at lower input orders.
 * @note The paper [1] also detailed a diffuse-field covariance constraint,
 *       and the original accronym was TAC (C=contrained), however, in SAF, this
 *       constraint is implemented as an independent operation. One may impose
 *       this constraint on any binaural decoder using the
 *       applyDiffCovMatching() function.
 *
 * @param[in]  hrtfs         The HRTFs; FLAT: N_bands x #NUM_EARS x N_dirs
 * @param[in]  hrtf_dirs_deg HRTF directions; FLAT: N_dirs x 2
 * @param[in]  N_dirs        Number of HRTF directions in set
 * @param[in]  N_bands       Number of frequency bands/bins
 * @param[in]  order         Decoding order
 * @param[in]  freqVector    Frequency vector; N_bands x 1
 * @param[in]  itd_s         Interaural time differences (ITDs), seconds;
 *                           N_dirs x 1
 * @param[in]  weights       Integration weights (set to NULL if not available);
 *                           N_dirs x 1
 * @param[out] decMtx        Decoding matrix;
 *                           FLAT: N_bands x #NUM_EARS x (order+1)^2
 *
 * @see [1] Zaunschirm M, Scho"rkhuber C, Ho"ldrich R. Binaural rendering of
 *          Ambisonic signals by head-related impulse response time alignment
 *          and a diffuseness constraint. The Journal of the Acoustical Society
 *          of America. 2018 Jun 19;143(6):3616--27
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

/**
 * Computes a binaural ambisonic decoder based on the magnitude least-squares
 * (MagLS) method, first described in [1], with the algorithm given in [2]
 *
 * The binaural ambisonic decoder is computed for each frequency bin/band,
 * ready to be applied to input SH signals in the time-frequency domain, or,
 * take the inverse-FFT and apply it via matrix convolution.
 *
 * @note Mag-LS operates under similar principles held by the TA/TAC decoder,
 *       differing in the manner in which the phase is neglected at frequencies
 *       above 1.5kHz.
 *
 * @param[in]  hrtfs         The HRTFs; FLAT: N_bands x #NUM_EARS x N_dirs
 * @param[in]  hrtf_dirs_deg HRTF directions; FLAT: N_dirs x 2
 * @param[in]  N_dirs        Number of HRTF directions in set
 * @param[in]  N_bands       Number of frequency bands/bins
 * @param[in]  order         Decoding order
 * @param[in]  freqVector    Frequency vector; N_bands x 1
 * @param[in]  weights       Integration weights (set to NULL if not available);
 *                           N_dirs x 1
 * @param[out] decMtx        Decoding matrix;
 *                           FLAT: N_bands x #NUM_EARS x (order+1)^2
 *
 * @see [1] Scho"rkhuber C, Zaunschirm M, Ho"ldrich R. Binaural Rendering of
 *          Ambisonic Signals via Magnitude Least Squares. InProceedings of the
 *          DAGA 2018 (Vol. 44, pp. 339--342).
 * @see [2] Zotter, F., & Frank, M. (2019). Ambisonics. Springer Open.
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

#endif /* __SAF_HOA_INTERNAL_H_INCLUDED__ */
