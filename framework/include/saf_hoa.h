/*
 Copyright 2018 Leo McCormack
 
 Permission to use, copy, modify, and/or distribute this software for any purpose with or
 without fee is hereby granted, provided that the above copyright notice and this permission
 notice appear in all copies.
 
 THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES WITH REGARD TO
 THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS. IN NO EVENT
 SHALL THE AUTHOR BE LIABLE FOR ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR
 ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF
 CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE
 OR PERFORMANCE OF THIS SOFTWARE.
*/
/*
 * Filename:
 *     saf_hoa.h (include header)
 * Description:
 *     A Higher-order Ambisonics C library; largely derived from the MatLab library by
 *     Archontis Politis: https://github.com/polarch/Higher-Order-Ambisonics
 * Dependencies:
 *     saf_utilities, saf_sh, saf_vbap
 * Author, date created:
 *     Leo McCormack, 19.03.2018
 */

#ifndef __SAF_HOA_H_INCLUDED__
#define __SAF_HOA_H_INCLUDED__

#ifdef __cplusplus
extern "C" {
#endif
    
/*************************/
/* Processing modes tags */
/*************************/

typedef enum _AMBI_DECODER_METHODS {
    DECODER_DEFAULT,        /* default is "DECODER_SAD" */
    DECODER_SAD,            /* Sampling Ambisonic Decoder (SAD); transpose of loudspeaker spherical harmonic matrix, scaled by number of loudspeakers */
    DECODER_MMD,            /* Mode-Matching Decoder (MMD); pseudo-inverse of the loudspeaker spherical harmonic matrix */
    DECODER_EPAD,           /* Energy-Preserving Ambisonic Decoder (EPAD) */
    DECODER_ALLRAD          /* All-Round Ambisonic Decoder (AllRAD); the most VBAP-like decoder ;) */
    
} AMBI_DECODER_METHODS;

/******************/
/* Main Functions */
/******************/
    
/* returns the weights required to manipulate the beam-patterns, such that they have maximum energy towards a given look-direction */
/* Zotter, F., Frank, M. (2012). All-Round Ambisonic Panning and Decoding. Journal of the Audio Engineering Society, 60(10), 807?820. */
void getMaxREweights(/* Input arguments */
                     int order,                  /* decoding order */
                     /* Output arguments */
                     float* a_n);                /* the max_rE weights, as a diagonal matrix; (order+1)^2 x (order+1)^2 */

/* returns an ambisonic decoding matrix of a specific order, for a specific loudspeaker set-up */
void getAmbiDecoder(/* Input arguments */
                    float* ls_dirs_deg,          /* loudspeaker directions in degrees [azi elev]; FLAT: nLS x 2 */
                    int nLS,                     /* number of loudspeakers */
                    AMBI_DECODER_METHODS method, /* decoding method to use (see AMBI_DECODER_METHODS enum) */
                    int order,                   /* decoding order */
                    /* Output arguments */
                    float** decMtx);             /* & decoding matrix; FLAT: nLS x (order+1)^2 */
 
#ifdef __cplusplus
}
#endif

#endif /* __SAF_HOA_H_INCLUDED__ */




