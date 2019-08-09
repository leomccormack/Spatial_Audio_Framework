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
 * Filename: saf_hoa.h (include header)
 * ------------------------------------
 * A collection of higher-order Ambisonics related functions. Some of which are
 * derived from the Matlab library by Archontis Politis, found here:
 *     https://github.com/polarch/Higher-Order-Ambisonics
 *
 * Dependencies:
 *     saf_utilities, saf_sh, saf_vbap
 * Author, date created:
 *     Leo McCormack, 19.03.2018
 */

#ifndef __SAF_HOA_H_INCLUDED__
#define __SAF_HOA_H_INCLUDED__

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */
    
#include "../saf_utilities/saf_complex.h"
#include "../saf_utilities/saf_error.h"
    
/* ========================================================================== */
/*                                    Enums                                   */
/* ========================================================================== */
    
/*
 * Enum: AMBI_DECODER_METHODS
 * --------------------------
 * Ambisonic decoding options for loudspeaker playback
 * Note: all of these decoding options revert to "SAD" if the loudspeakers are
 * uniformly distributed on the sphere. The benefits afforded by MMD and AllRAD,
 * relate to their improved performance when using irregular loudspeaker
 * arrangements.
 *
 * Options:
 *     DECODER_DEFAULT - DECODER_DEFAULT is "DECODER_SAD"
 *     DECODER_SAD     - Sampling Ambisonic Decoder (SAD): transpose of the
 *                       loudspeaker spherical harmonic matrix, scaled by the
 *                       number of loudspeakers. This is the simplest decoding
 *                       approach, as it simply relies on generating hyper-
 *                       cardioid beamformers for each loudspeaker direction.
 *     DECODER_MMD     - Mode-Matching Decoder (MMD): pseudo-inverse of the
 *                       loudspeaker spherical harmonic matrix. Due to the
 *                       pseudo-inverse, more signal energy is lent to regions
 *                       on the surface of the sphere that are more sparsely
 *                       populated with loudspeakers. Therefore, one must also
 *                       be careful, as some loudspeakers may be given a huge
 *                       amount of signal energy and wake the dead.
 *     DECODER_EPAD    - Energy-Preserving Ambisonic Decoder (EPAD) [1]
 *     DECODER_ALLRAD  - All-Round Ambisonic Decoder (AllRAD): SAD decoding to
 *                       t-design, panned for the target loudspeaker directions
 *                       using VBAP [2]. Perhaps the only Ambisonic decoder we
 *                       would actually recommend for irregular loudspeaker
 *                       layouts. It just simply makes the most sense and sounds
 *                       the most reasonable. Interestingly, prior to its
 *                       proposal, we had almost 40 years! of people insisting
 *                       Ambisonics could only work for uniformly distributed
 *                       loudspeaker arrays... So thank you Franz & Frank!
 *
 * [1] Zotter F, Pomberger H, Noisternig M. Energy- preserving ambisonic
 *     decoding. Acta Acustica united with Acustica. 2012 Jan 1; 98(1):37-47.
 * [2] Zotter F, Frank M. All-round ambisonic panning and decoding. Journal of
 *     the audio engineering society. 2012 Nov 26; 60(10):807-20.
 */
typedef enum _AMBI_DECODER_METHODS {
    DECODER_DEFAULT,
    DECODER_SAD,
    DECODER_MMD,
    DECODER_EPAD,
    DECODER_ALLRAD
    
} AMBI_DECODER_METHODS;
    
/*
 * Enum: BINAURAL_AMBI_DECODER_METHODS
 * -----------------------------------
 * Ambisonic decoding options for headphone playback.
 * Note: more detailed descriptions of each method may be found in:
 * "saf_hoa_internal.h"
 *
 * Options:
 *     BINAURAL_DECODER_DEFAULT  - DECODER_DEFAULT is "DECODER_SAD"
 *     BINAURAL_DECODER_LS       - Least-squares (LS) decoder. The simplest
 '                                 binaural decoder.
 *     BINAURAL_DECODER_LSDIFFEQ - Least-squares (LS) decoder with diffuse-field
 *                                 spectral equalisation  [1]
 *     BINAURAL_DECODER_SPR      - Spatial resampling decoder (on the same lines
 *                                 as the virtual loudspeaker approach) [2]
 *     BINAURAL_DECODER_TA       - Time-alignment decoder [3]. Relies on
 *                                 discarding the phase information in HRTFs,
 *                                 past the frequency at which humans are less
 *                                 sensitive to inter-aural time differences.
 *     BINAURAL_DECODER_MAGLS    - Magnitude least-squares decoder [4]. On
 *                                 similar lines to the time-alignment decoder,
 *                                 but differing in its execution.
 *
 * [1] Z. Ben-Hur, F. Brinkmann, J. Sheaffer, S. Weinzierl, and B. Rafaely,
 *     "Spectral equalization in binaural signals represented by order-
 *     truncated spherical harmonics," The Journal of the Acoustical Society
 *     of America, vol. 141, no. 6, pp. 4087–4096, 2017.
 * [2] B. Bernschutz, A. V. Giner, C. Pörschmann, and J. Arend, “Binaural
 *     reproduction of plane waves with reduced modal order,” Acta Acustica
 *     united with Acustica, vol. 100, no. 5, pp. 972–983, 2014.
 * [3] Zaunschirm M, Schörkhuber C, Höldrich R. Binaural rendering of
 *     Ambisonic signals by head-related impulse response time alignment and
 *     a diffuseness constraint. The Journal of the Acoustical Society of
 *     America. 2018 Jun 19;143(6):3616-27
 * [4] Schörkhuber C, Zaunschirm M, Höldrich R. Binaural Rendering of
 *     Ambisonic Signals via Magnitude Least Squares. InProceedings of the
 *     DAGA 2018 (Vol. 44, pp. 339-342).
 */
typedef enum _BINAURAL_AMBI_DECODER_METHODS {
    BINAURAL_DECODER_DEFAULT,
    BINAURAL_DECODER_LS,
    BINAURAL_DECODER_LSDIFFEQ,
    BINAURAL_DECODER_SPR,
    BINAURAL_DECODER_TA,
    BINAURAL_DECODER_MAGLS
    
} BINAURAL_AMBI_DECODER_METHODS;
    
    
/* ========================================================================== */
/*                               Main Functions                               */
/* ========================================================================== */
    
/*
 * Function: getMaxREweights
 * -------------------------
 * Returns the weights required to manipulate a hyper-cardioid beam-pattern,
 * such that it has maximum energy in the given look-direction [1].
 * Traditionally, due to the back lobes of beamformers when panning a source
 * via Ambisonics encoding/decoding, there is unwanted energy given to
 * loudspeakers directly opposite the true source direction. This max_rE
 * weighting essentially spatially "tapers" the spherical harmonic patterns
 * used to generate said beams, by reducing the contribution of the higher
 * orders to the beam patterns. This results in worse spatial selectivity, as
 * the width of the beam pattern main lobe is widened. However, the back lobes
 * are also reduced, thus mitigating the aforementioned problem.
 *
 * Input Arguments:
 *     order - decoding order
 * Output Arguments:
 *     a_n   - the max_rE weights, as a diagonal matrix;
 *             FLAT: (order+1)^2 x (order+1)^2
 *
 * [1] Zotter, F., Frank, M. (2012). All-Round Ambisonic Panning and Decoding.
 *     Journal of the Audio Engineering Society, 60(10), 807--820.
 */
void getMaxREweights(/* Input Arguments */
                     int order,
                     /* Output Arguments */
                     float* a_n);

/*
 * Function: getAmbiDecoder
 * ------------------------
 * Returns an ambisonic decoding matrix of a specific order, for a specific
 * loudspeaker layout
 *
 * Input Arguments:
 *     ls_dirs_deg - loudspeaker directions in DEGREES [azi elev]; FLAT: nLS x 2
 *     nLS         - number of loudspeakers
 *     method      - decoding method (see "AMBI_DECODER_METHODS" enum)
 *     order       - decoding order
 * Output Arguments:
 *     decMtx      - & decoding matrix; FLAT: nLS x (order+1)^2
 */
void getAmbiDecoder(/* Input Arguments */
                    float* ls_dirs_deg,
                    int nLS,
                    AMBI_DECODER_METHODS method,
                    int order,
                    /* Output Arguments */
                    float** decMtx);
    
/*
 * Function: getBinauralAmbiDecoder
 * --------------------------------
 * Returns an ambisonic decoding matrix of a specific order, for a specific set
 * of HRTFs
 *
 * Input Arguments:
 *     hrtfs         - the HRTFs; FLAT: N_bands x 2 x N_dirs
 *     hrtf_dirs_deg - HRTF directions; FLAT: N_dirs x 2
 *     N_dirs        - number of HRTF directions
 *     N_bands       - number of frequency bands/bins
 *     method        - decoding method (see BINAURAL_AMBI_DECODER_METHODS enum)
 *     order         - decoding order
 *     freqVector    - only needed for TAC/Mag-LS decoders (can set to NULL if
 *                     not needed); N_bands x 1
 *     itd_s         - only needed for TAC  decoder (can set to NULL if
 *                     not needed); N_dirs x 1
 *     weights       - integration weights (set to NULL if not available);
 *                     N_dirs x 1
 * Output Arguments:
 *     decMtx        - decoding matrix; FLAT: N_bands x 2 x (order+1)^2
 */
void getBinauralAmbiDecoder(/* Input Arguments */
                            float_complex* hrtfs,
                            float* hrtf_dirs_deg,
                            int N_dirs,
                            int N_bands,
                            BINAURAL_AMBI_DECODER_METHODS method,
                            int order,
                            float* freqVector,
                            float* itd_s,
                            float* weights,
                            /* Output Arguments */
                            float_complex* decMtx);

/*
 * Function: applyDiffCovMatching
 * ------------------------------
 * Imposes a diffuse-field covariance constraint on a given decoding matrix [1]
 * Note: decMtx is altered in-place.
 *
 * Input Arguments:
 *     hrtfs         - the HRTFs; FLAT: N_bands x 2 x N_dirs
 *     hrtf_dirs_deg - HRTF directions; FLAT: N_dirs x 2
 *     N_dirs        - number of HRTF directions
 *     N_bands       - number of frequency bands/bins
 *     order         - decoding order
 *     weights       - integration weights (set to NULL if not available);
 *                     N_dirs x 1
 * Input/Output Arguments:
 *     decMtx        - decoding matrix; FLAT: N_bands x 2 x (order+1)^2
 *
 * [1] Zaunschirm M, Schörkhuber C, Höldrich R. Binaural rendering of
 *     Ambisonic signals by head-related impulse response time alignment and
 *     a diffuseness constraint. The Journal of the Acoustical Society of
 *     America. 2018 Jun 19;143(6):3616-27
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
