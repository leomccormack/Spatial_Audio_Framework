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
 *@addtogroup HOA
 *@{
 * @file saf_hoa.h
 * @brief Main header for the higher-order Ambisonics module (#SAF_HOA_MODULE)
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

#ifndef __SAF_HOA_H_INCLUDED__
#define __SAF_HOA_H_INCLUDED__

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include "../saf_utilities/saf_utility_complex.h"

/* ========================================================================== */
/*                                    Enums                                   */
/* ========================================================================== */

/**
 * Ambisonic decoding options for loudspeaker playback
 *
 * Note that the MMD and EPAD decoding options revert back to "SAD" if the
 * loudspeakers are uniformly distributed on the sphere. The benefits afforded
 * by MMD, EPAD [1], and AllRAD [2] relate to their improved performance when
 * using irregular loudspeaker arrangements.
 *
 * @see [1] Zotter F, Pomberger H, Noisternig M. Energy--preserving ambisonic
 *          decoding. Acta Acustica united with Acustica. 2012 Jan 1;
 *          98(1):37-47.
 * @see [2] Zotter F, Frank M. All-round ambisonic panning and decoding. Journal
 *          of the audio engineering society. 2012 Nov 26; 60(10):807-20.
 */
typedef enum {
    /**
     * The default decoder is #LOUDSPEAKER_DECODER_SAD
     */
    LOUDSPEAKER_DECODER_DEFAULT,
    /**
     * Sampling Ambisonic Decoder (SAD): transpose of the loudspeaker spherical
     * harmonic matrix, scaled by the number of loudspeakers. This is the
     * simplest decoding approach, as it essentially just generates hyper-
     * cardioid beamformers (aka virtual microphones) towards each loudspeaker
     * direction. This approach is numerically robust to irregular loudspeaker
     * arrangements. However, it does not preserve the energy of a source (or
     * localisation cues) as it is panned around in different directions over
     * irregular setups.
     */
    LOUDSPEAKER_DECODER_SAD,
    /**
     * Mode-Matching Decoder (MMD): pseudo-inverse of the loudspeaker spherical
     * harmonic matrix. Due to the pseudo-inverse, more signal energy is lent to
     * regions on the surface of the sphere that are more sparsely populated
     * with loudspeakers; (this is essentially a least-squares solution).
     * Therefore, this approach can help balance out directional loudness
     * differences when using slightly irregular setups. However, one must also
     * be careful since loudspeakers that are very far way from all the other
     * loudspeakers (e.g. voice-of-god) may be given significantly more signal
     * energy than expected. Therefore, this approach is not recommended for
     * highly irregular loudspeaker arrangements!
     */
    LOUDSPEAKER_DECODER_MMD,
    /**
     * Energy-Preserving Ambisonic Decoder (EPAD) [1]. This decoder aims to
     * preserve the energy of a source, as it panned around to directions of the
     * sphere; essentially, addressing the energy-preserving issues of the
     * SAD and MMD decoding approaches for irregular layouts.
     */
    LOUDSPEAKER_DECODER_EPAD,
    /**
     * All-Round Ambisonic Decoder (AllRAD): SAD decoding to a t-design, panned
     * for the target loudspeaker directions using VBAP [2]. Perhaps the
     * Ambisonic decoder we would most recommend for irregular loudspeaker
     * layouts. Note, given a high (well... technically infinite) order, AllRAD
     * will converge to VBAP. However, since lower-orders are employed in
     * practice, AllRAD is not as spatially "sharp" as VBAP, but it will yield
     * more consistent source spread when panning a source inbetween the
     * loudspeakers.
     * The approach is highly robust to irregular loudspeaker setups, and
     * exhibits low directional error and good energy-preserving properties.
     */
    LOUDSPEAKER_DECODER_ALLRAD
    
} LOUDSPEAKER_AMBI_DECODER_METHODS;

/**
 * Ambisonic decoding options for binaural/headphone playback
 *
 * @note A more detailed description of each method may be found in
 *       saf_hoa_internal.h.
 *
 * @see [1] Z. Ben-Hur, F. Brinkmann, J. Sheaffer, S. Weinzierl, and B. Rafaely,
 *          "Spectral equalization in binaural signals represented by order-
 *          truncated spherical harmonics" The Journal of the Acoustical
 *          Society of America, vol. 141, no. 6, pp. 4087--4096, 2017.
 * @see [2] Zaunschirm M, Scho"rkhuber C, Ho"ldrich R. Binaural rendering of
 *          Ambisonic signals by head-related impulse response time alignment
 *          and a diffuseness constraint. The Journal of the Acoustical Society
 *          of America. 2018 Jun 19 143(6) 3616-27
 * @see [3] Scho"rkhuber C, Zaunschirm M, Ho"ldrich R. Binaural Rendering of
 *          Ambisonic Signals via Magnitude Least Squares. InProceedings of the
 *          DAGA 2018 (Vol. 44, pp. 339-342).
 * @see [4] B. Bernschutz, A. V. Giner, C. Po"rschmann, and J. Arend, "Binaural
 *          reproduction of plane waves with reduced modal order" Acta Acustica
 *          united with Acustica, vol. 100, no. 5, pp. 972--983, 2014
 */
typedef enum {
    /**
     * The default decoder is #BINAURAL_DECODER_LS
     */
    BINAURAL_DECODER_DEFAULT,
    /**
     * Least-squares (LS) decoder. The simplest binaural decoder, which is based
     * on a least-squares fit of the spherical harmonic patterns onto the HRTF
     * directivity patterns.
     */
    BINAURAL_DECODER_LS,
    /**
     * Least-squares (LS) decoder with diffuse-field spectral equalisation [1].
     * Note that the diffuse-field EQ is applied in the spherical harmonic
     * domain (to account for the truncation error/loss of high frequencies), so
     * this is not the same as applying diffuseFieldEqualiseHRTFs() on the HRTFs
     * followed by BINAURAL_DECODER_LS.
     */
    BINAURAL_DECODER_LSDIFFEQ,
    /**
     * Spatial resampling decoder (on the same lines as the virtual loudspeaker
     * approach) [4].
     */
    BINAURAL_DECODER_SPR,
    /**
     * Time-alignment decoder [2]. Relies on discarding the phase information of
     * the HRTFs, past the frequency at which humans are less sensitive to
     * inter-aural time difference cues. Therefore, the least-squares fitting
     * priorites matching the interaural level differences (ILDs), rather than
     * the interaural time differences (ITDs).
     */
    BINAURAL_DECODER_TA,
    /**
     * Magnitude least-squares decoder [3]. On similar lines to the time-
     * alignment decoder, but differing slightly in its execution.
     */
    BINAURAL_DECODER_MAGLS
    
} BINAURAL_AMBI_DECODER_METHODS;

/**
 * Available Ambisonic channel ordering conventions
 *
 * @note ACN channel ordering with SN3D normalisation is often collectively
 *       referred to as the 'AmbiX' format.
 * @warning FuMa is a deprecated legacy format and is only supported for first-
 *          order! The recommended Ambisonic conventions are ACN with N3D or
 *          SN3D normalisation.
 */
typedef enum {
    HOA_CH_ORDER_ACN,  /**< Ambisonic Channel numbering (ACN) convention, which
                        *   is employed by all spherical harmonic related
                        *   functions in SAF */
    HOA_CH_ORDER_FUMA  /**< Furse-Malham (FuMa) convention, often used by older
                        *   recordings. The convention follows the WXYZ ordering
                        *   of the omni and dipoles, and is suitable only for
                        *   1st order. */

}HOA_CH_ORDER;

/**
 * Available Ambisonic normalisation conventions
 *
 * @note ACN channel ordering with SN3D normalisation is often collectively
 *       referred to as the 'AmbiX' format.
 * @warning FuMa is a deprecated legacy format and is only supported for first-
 *          order! The recommended Ambisonic conventions are ACN with N3D/SN3D
 *          normalisation.
 */
typedef enum {
    HOA_NORM_N3D,  /**< Orthonormalised (N3D) convention, which is the default
                    *   convention used by SAF */
    HOA_NORM_SN3D, /**< Schmidt semi-normalisation (SN3D) convention, as used
                    *   by the AmbiX standard */
    HOA_NORM_FUMA  /**< Furse-Malham (FuMa) convention. This is similar to SN3D
                    *   (at first-order), except theres an additional 1/sqrt(2)
                    *   scaling applied to the omni. This is also known as maxN
                    *   normalisation. */

}HOA_NORM;


/* ========================================================================== */
/*                               Misc. Functions                              */
/* ========================================================================== */

/**
 * Converts an Ambisonic signal from one channel ordering convention to another
 *
 * @warning If one of the in/out conventions is FuMa, then only the first 4
 *          channels are converted, and any remaining channels of 'insig' are
 *          set to zeros (i.e. FuMa is strictly first-order only in SAF).
 * @note insig is converted "in-place". Also, if the in/out conventions are the
 *       same, then the function is bypassed.
 *
 * @param[in,out] insig         Input signal with the channel ordering
 *                              convention of: inConvention;
 *                              FLAT: (order+1)^2 x signalLength
 * @param[in]     order         Ambisonic order
 * @param[in]     signalLength  Signal length in samples
 * @param[in]     inConvention  Channel order convention of input signals
 * @param[in]     outConvention Channel order convention of output signals
 */
void convertHOAChannelConvention(/* Input Arguments */
                                 float* insig,
                                 int order,
                                 int signalLength,
                                 HOA_CH_ORDER inConvention,
                                 HOA_CH_ORDER outConvention);

/**
 * Converts an Ambisonic signal from one normalisation convention to another
 *
 * @warning If one of the in/out conventions is FuMa, then only the first 4
 *          channels are converted, and any remaining channels of 'insig' are
 *          set to zeros (FuMa is strictly first-order only in SAF).
 * @note insig is converted "in-place". Also, if the in/out conventions are the
 *       same, then the function is bypassed.
 *
 * @param[in,out] insig         Input signal with the channel ordering
 *                              convention of: inConvention, which should be
 *                              converted to: outConvention, "in-place";
 *                              FLAT: (order+1)^2 x signalLength
 * @param[in]     order         Ambisonic order
 * @param[in]     signalLength  Signal length in samples
 * @param[in]     inConvention  Normalisation convention of the input signals
 * @param[in]     outConvention Normalisation convention of the output signals
 */
void convertHOANormConvention(/* Input Arguments */
                              float* insig,
                              int order,
                              int signalLength,
                              HOA_NORM inConvention,
                              HOA_NORM outConvention);

/**
 * Computes real-valued spherical harmonics [1] for each given direction on the
 * unit sphere
 *
 * The spherical harmonic values are computed WITHOUT the 1/sqrt(4*pi) term.
 * Compared to getRSH_recur(), this function uses unnorm_legendreP() and double
 * precision, so is more suitable for being computed in an initialisation stage.
 * This version is indeed slower, but more precise (especially for high orders).
 *
 * @note This function is mainly intended for Ambisonics, due to the omission of
 *       the 1/sqrt(4*pi) scaling, and the directions are given in
 *       [azimuth elevation] (degrees). In Ambisonics literature, the format
 *       convention of 'Y' is referred to as ACN/N3D
 *
 * @param[in]  order    Order of spherical harmonic expansion
 * @param[in]  dirs_deg Directions on the sphere [azi, ELEVATION] convention, in
 *                      DEGREES; FLAT: nDirs x 2
 * @param[in]  nDirs    Number of directions
 * @param[out] Y        The SH weights [WITHOUT the 1/sqrt(4*pi)];
 *                      FLAT: (order+1)^2 x nDirs
 *
 * @see [1] Rafaely, B. (2015). Fundamentals of spherical array processing
 *          (Vol. 8). Berlin: Springer.
 */
void getRSH(/* Input Arguments */
            int order,
            float* dirs_deg,
            int nDirs,
            /* Output Arguments */
            float* Y);

/**
 * Computes real-valued spherical harmonics [1] for each given direction on the
 * unit sphere
 *
 * The real spherical harmonics are computed WITHOUT the 1/sqrt(4*pi) term.
 * Compared to getRSH(), this function uses unnorm_legendreP_recur() and single
 * precision, so is more suitable for being computed in a real-time loop. It
 * sacrifices some precision, and numerical error propogates through the
 * recursion, but it is much faster.
 *
 * The function also uses static memory buffers for single direction and up to
 * 7th order, which speeds things up considerably for such use cases.
 *
 * @note This function is mainly intended for Ambisonics, due to the omission of
 *       the 1/sqrt(4*pi) scaling, and the directions are given in
 *       [azimuth elevation] (degrees). In Ambisonics literature, the format
 *       convention of 'Y' is referred to as ACN/N3D
 *
 * @param[in]  order    Order of spherical harmonic expansion
 * @param[in]  dirs_deg Directions on the sphere [azi, ELEVATION] convention, in
 *                      DEGREES; FLAT: nDirs x 2
 * @param[in]  nDirs    Number of directions
 * @param[out] Y        The SH weights [WITHOUT the 1/sqrt(4*pi)];
 *                      FLAT: (order+1)^2 x nDirs
 *
 * @see [1] Rafaely, B. (2015). Fundamentals of spherical array processing
 *          (Vol. 8). Berlin: Springer.
 */
void getRSH_recur(/* Input Arguments */
                  int order,
                  float* dirs_deg,
                  int nDirs,
                  /* Output Arguments */
                  float* Y);


/* ========================================================================== */
/*                               Main Functions                               */
/* ========================================================================== */

/**
 * Computes the weights required to manipulate a hyper-cardioid beam-pattern,
 * such that it has maximum energy in the given look-direction
 *
 * Due to the side and back lobes of the beamformers employed by the Ambisonic
 * decoder, when panning a source there can be unwanted energy given to
 * loudspeakers directly opposite the true source direction. This max_rE
 * weighting [1] essentially spatially tapers the spherical harmonic components
 * used to generate the beamformers, thus reducing the contribution of the
 * higher order components. This results in worse spatial selectivity, as
 * the width of the beam pattern main lobe is widened, however, the back lobes
 * are also reduced, thus mitigating perceptual issues that may arise due to the
 * aforementioned problem.
 *
 * @param[in]  order       Order of spherical harmonic expansion
 * @param[in]  diagMtxFlag Set to '0' if you want the weights to be returned as
 *                         a vector, or to '1' as a diagonal matrix instead.
 * @param[out] a_n         The max_rE weights, as a vector/diagonal matrix;
 *                         (order+1)^2 x 1 OR FLAT: (order+1)^2 x (order+1)^2
 *
 * @see [1] Zotter F, Frank M. All-round ambisonic panning and decoding. Journal
 *          of the audio engineering society. 2012 Nov 26; 60(10):807-20.
 */
void getMaxREweights(/* Input Arguments */
                     int order,
                     int diagMtxFlag,
                     /* Output Arguments */
                     float* a_n);

/**
 * Filter that equalises the high frequency roll-off due to SH truncation and
 * tapering; as described in [1].
 *
 * @param[in]  w_n             Tapering weights; (order_truncated + 1) x 1
 *                             E.g. maxRE, or all ones for truncation only
 * @param[in]  order_truncated Input SH order
 * @param[in]  order_target    Target SH order, (should be higher, e.g. 38)
 * @param[in]  kr              kr vector, r e.g. 0.085 m; nBands x 1
 * @param[in]  nBands          Number of frequency bins
 * @param[in]  softThreshold   Threshold in dB, soft limiting above to +6dB
 * @param[out] gain            Gain factor for compensation filter; nBands x 1
 *
 * @see [1] Hold, C., Gamper, H., Pulkki, V., Raghuvanshi, N., & Tashev, I. J. 
 *          (2019). Improving Binaural Ambisonics Decoding by Spherical 
 *          Harmonics Domain Tapering and Coloration Compensation. ICASSP, 
 *          IEEE International Conference on Acoustics, Speech and Signal 
 *          Processing - Proceedings.
 */
void truncationEQ(/* Input arguments */
                   float* w_n,
                   int order_truncated,
                   int order_target,
                   double* kr,
                   int nBands,
                   float softThreshold,
                   /* Output arguments */
                   float* gain);

/**
 * Computes an ambisonic decoding matrix of a specific order, for a given
 * loudspeaker layout
 *
 * @test test__getLoudspeakerDecoderMtx()
 *
 * @param[in]  ls_dirs_deg Loudspeaker directions in DEGREES [azi elev];
 *                         FLAT: nLS x 2
 * @param[in]  nLS         Number of loudspeakers
 * @param[in]  method      Decoding method (see
 *                         #LOUDSPEAKER_AMBI_DECODER_METHODS enum)
 * @param[in]  order       Decoding order
 * @param[in]  enableMaxrE Set to '0' to disable, '1' to enable
 * @param[out] decMtx      Decoding matrix; FLAT: nLS x (order+1)^2
 */
void getLoudspeakerDecoderMtx(/* Input Arguments */
                              float* ls_dirs_deg,
                              int nLS,
                              LOUDSPEAKER_AMBI_DECODER_METHODS method,
                              int order,
                              int enableMaxrE,
                              /* Output Arguments */
                              float* decMtx);

/**
 * Computes binaural ambisonic decoding matrices (one per frequency) at a
 * specific order, for a given HRTF set
 *
 * @param[in]  hrtfs         The HRTFs; FLAT: N_bands x #NUM_EARS x N_dirs
 * @param[in]  hrtf_dirs_deg HRTF directions; FLAT: N_dirs x 2
 * @param[in]  N_dirs        Number of HRTF directions
 * @param[in]  N_bands       Number of frequency bands/bins
 * @param[in]  method        Decoder method (see #BINAURAL_AMBI_DECODER_METHODS
 *                           enum)
 * @param[in]  order         Decoding order
 * @param[in]  freqVector    Only needed for #BINAURAL_DECODER_TA or
 *                           #BINAURAL_DECODER_MAGLS decoders (set to NULL if
 *                           using a different method); N_bands x 1
 * @param[in]  itd_s         Only needed for #BINAURAL_DECODER_TA decoder (set
 *                           to NULL if using different method); N_dirs x 1
 * @param[in]  weights       Integration weights (set to NULL if not available);
 *                           N_dirs x 1
 * @param[in]  enableDiffCM  Set to '0' to disable diffuse correction, '1' to
 *                           enable
 * @param[in]  enableMaxrE   Set to '0' to disable maxRE weighting, '1' to
 *                           enable
 * @param[out] decMtx        Decoding matrices (one per frequency);
 *                           FLAT: N_bands x #NUM_EARS x (order+1)^2
 */
void getBinauralAmbiDecoderMtx(/* Input Arguments */
                               float_complex* hrtfs,
                               float* hrtf_dirs_deg,
                               int N_dirs,
                               int N_bands,
                               BINAURAL_AMBI_DECODER_METHODS method,
                               int order,
                               float* freqVector,
                               float* itd_s,
                               float* weights,
                               int enableDiffCM,
                               int enableMaxrE,
                               /* Output Arguments */
                               float_complex* decMtx);

/**
 * Computes binaural ambisonic decoding filters for a given HRTF set
 *
 * @param[in]  hrtfs         The HRTFs; FLAT: (fftSize/2+1) x #NUM_EARS x N_dirs
 * @param[in]  hrtf_dirs_deg HRTF directions; FLAT: N_dirs x 2
 * @param[in]  N_dirs        Number of HRTF directions
 * @param[in]  fftSize       FFT size
 * @param[in]  fs            Sampling rate
 * @param[in]  method        Decoder method (see #BINAURAL_AMBI_DECODER_METHODS
 *                           enum)
 * @param[in]  order         Decoding order
 * @param[in]  itd_s         Only needed for #BINAURAL_DECODER_TA decoder (can
 *                           set to NULL if using different method); N_dirs x 1
 * @param[in]  weights       Integration weights (set to NULL if not available);
 *                           N_dirs x 1
 * @param[in]  enableDiffCM  Set to '0' to disable diffuse correction, '1' to
 *                           enable
 * @param[in]  enableMaxrE   Set to '0' to disable maxRE weighting, '1' to
 *                           enable
 * @param[out] decFilters    Decoding filters;
 *                           FLAT: #NUM_EARS x (order+1)^2 x fftSize
 */
void getBinauralAmbiDecoderFilters(/* Input Arguments */
                                   float_complex* hrtfs,
                                   float* hrtf_dirs_deg,
                                   int N_dirs,
                                   int fftSize,
                                   float fs,
                                   BINAURAL_AMBI_DECODER_METHODS method,
                                   int order,
                                   float* itd_s,
                                   float* weights,
                                   int enableDiffCM,
                                   int enableMaxrE,
                                   /* Output Arguments */
                                   float* decFilters);

/**
 * Imposes a diffuse-field covariance constraint on a given binaural decoding
 * matrix, as described in [1]
 *
 * @note decMtx is altered in-place.
 *
 * @param[in]     hrtfs         The HRTFs; FLAT: N_bands x #NUM_EARS x N_dirs
 * @param[in]     hrtf_dirs_deg HRTF directions; FLAT: N_dirs x 2
 * @param[in]     N_dirs        Number of HRTF directions
 * @param[in]     N_bands       Number of frequency bands/bins
 * @param[in]     order         Decoding order
 * @param[in]     weights       Integration weights (set to NULL if not
 *                              available); N_dirs x 1
 * @param[in,out] decMtx        Decoding matrix;
 *                              FLAT: N_bands x #NUM_EARS x (order+1)^2
 *
 * @see [1] Zaunschirm M, Scho"rkhuber C, Ho"ldrich R. Binaural rendering of
 *          Ambisonic signals by head-related impulse response time alignment
 *          and a diffuseness constraint. The Journal of the Acoustical Society
 *          of America. 2018 Jun 19;143(6):3616-27
 */
void applyDiffCovMatching(/* Input Arguments */
                          float_complex* hrtfs,
                          float* hrtf_dirs_deg,
                          int N_dirs,
                          int N_bands,
                          int order,
                          float* weights,
                          /* Input/Output Arguments */
                          float_complex* decMtx);


#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */

#endif /* __SAF_HOA_H_INCLUDED__ */

/**@} */ /* doxygen addtogroup HOA */
