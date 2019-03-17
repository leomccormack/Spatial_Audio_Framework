/*
 Copyright 2017-2018 Leo McCormack
 
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
 *     saf_hrir.h (include header)
 * Description:
 *     A collection of head-related impulse-response (HRIR)- related functions.
 * Dependencies:
 *     saf_utilities, afSTFTlib
 * Author, date created:
 *     Leo McCormack, 12.12.2016
 */

#ifndef __SAF_HRIR_H_INCLUDED__
#define __SAF_HRIR_H_INCLUDED__

#ifdef __cplusplus
extern "C" {
#endif
    
#include "../saf_utilities/saf_complex.h"
    
#ifndef M_PI
  #define M_PI ( 3.14159265359f )
#endif

/* estimates the interaural time-differences (ITDs) for each HRIR via the cross-correlation between the left and right IRs */
void estimateITDs(/* Input arguments */
                  float* hrirs,                           /* HRIRs; FLAT: N_dirs x 2 x hrir_len */
                  int N_dirs,                             /* number of HRIRs */
                  int hrir_len,                           /* length of the HRIRs in samples */
                  int fs,                                 /* sampling rate of the HRIRs */
                  /* Output arguments */
                  float** itds_s);                        /* & ITDs in seconds; N_dirs x 1 */

/* passes zero padded HRIRs through the afSTFT filterbank. The filterbank coefficients are then normalised with the energy
 * of an impulse, which is centered at approximately the beginning of the HRIR peak. The HRTF FB coefficients are then
 * diffuse-field equalised before reintroducing the interaural phase differences (IPDs) per frequency band.
 * Please note that this function is NOT suitable for binaural room impulse responses (BRIRs). */
void HRIRs2FilterbankHRTFs(/* Input arguments */
                           float* hrirs,                  /* HRIRs; FLAT: N_dirs x 2 x hrir_len */
                           int N_dirs,                    /* number of HRIRs */
                           int hrir_len,                  /* length of the HRIRs in samples */
                           float* itds_s,                 /* HRIR ITDs; N_dirs x 1 */
                           float* centreFreq,             /* filterbank centre frequencies; N_bands x 1 */
                           int N_bands,                   /* number of frequency bands */
                           int enablePhaseManipFLAG,      /* 0: off, 1: on */
                           /* Output arguments */
                           float_complex** hrtf_fb);      /* & HRTFs as filterbank coeffs; FLAT: N_bands x 2 x N_dirs */

/* Interpolates a set of HRTFs for specified directions; defined by a amplitude normalised vbap interpolation table (see saf_vbap).
 * The interpolation is performed by applying interpolation gains to the HRTF magnitudes and HRIR inter-aural time differences separately.
 * The inter-aural phase differences are then reintroduced for each frequency band */
void interpFilterbankHRTFs(/* Input arguments */
                           float_complex* hrtfs,          /* HRTFs as filterbank coeffs; FLAT: N_bands x 2 x N_hrtf_dirs */
                           float* itds,                   /* the inter-aural time difference for each HRIR; N_hrtf_dirs x 1 */
                           float* freqVector,             /* frequency vector; N_bands x 1 */
                           float* vbap_gtable,            /* vbap gain table; FLAT: N_interp_dirs x N_hrtf_dirs */
                           int N_hrtf_dirs,               /* number of HRTF directions */
                           int N_bands,                   /* number of frequency bands */
                           int N_interp_dirs,             /* number of interpolated hrtf positions  */
                           int enablePhaseManipFLAG,      /* 0: off, 1: on */
                           /* Output arguments */
                           float_complex* hrtf_interp);   /* pre-alloc, interpolated HRTFs; FLAT: N_bands x 2 x N_interp_dirs */
    
/* computes the binaural diffuse coherence per frequency */
void binauralDiffuseCoherence(/* Input arguments */
                              float_complex* hrtfs,       /* HRTFs as filterbank coeffs; FLAT: N_bands x 2 x N_hrtf_dirs */
                              float* itds,                /* the inter-aural time difference for each HRIR; N_hrtf_dirs x 1 */
                              float* freqVector,          /* frequency vector; N_bands x 1 */
                              int N_hrtf_dirs,            /* number of HRTF directions */
                              int N_bands,                /* number of frequency bands */
                              /* Output arguments */
                              float* HRTFcoh);            /* binaural coeherence per frequency; N_bands x 1 */
    

#ifdef __cplusplus
}
#endif


#endif /* __SAF_HRIR_H_INCLUDED__ */




