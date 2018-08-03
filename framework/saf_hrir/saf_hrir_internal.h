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
 *     saf_hrir_internal.h
 * Description:
 *     A collection of head-related impulse-response (HRIR)- related functions.
 * Dependencies:
 *     saf_utilities, afSTFTlib
 * Author, date created:
 *     Leo McCormack, 12.12.2016
 */

#ifndef __HRIR_INTERNAL_H_INCLUDED__
#define __HRIR_INTERNAL_H_INCLUDED__

#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <string.h>
#include <float.h>
#include "saf_hrir.h"
#include "afSTFTlib.h"
#include "saf_utilities.h"

#ifdef __cplusplus
extern "C" {
#endif
 
#ifndef NUM_EARS
  #define NUM_EARS 2
#endif
    
/* Calculates the cross correlation between two vectors */
void cxcorr(float* a,                               /* vector a */
            float* b,                               /* vector b */
            float* x_ab,                            /* cross-correlation result between a and b */
            size_t la,                              /* length of vector a */
            size_t lb);                             /* length of vector b */
 
/* Converts and FIR filter into Filterbank Coefficients
 * It is currently hard coded for a 128 hop size with hybrid mode enabled (see afSTFTlib) */
void FIRtoFilterbankCoeffs(float* hIR               /* time-domain FIR; N_dirs x nCH x ir_len */,
                           int N_dirs,              /* number of FIR sets */
                           int nCH,                 /* number of channels per FIR set */
                           int ir_len,              /* length of the FIR */
                           int N_bands,             /* number of time-frequency domain bands */
                           float_complex** hFB);    /* & the FIRs as Filterbank coefficients; N_bands x nCH x N_dirs */
    
    
#ifdef __cplusplus
}
#endif


#endif /* __HRIR_INTERNAL_H_INCLUDED__ */




















