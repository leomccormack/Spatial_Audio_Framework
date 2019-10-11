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
 * Filename: saf_hrir_internal.h
 * -----------------------------
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

#ifndef __HRIR_INTERNAL_H_INCLUDED__
#define __HRIR_INTERNAL_H_INCLUDED__

#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <string.h>
#include "saf_hrir.h"
#include "../../resources/afSTFT/afSTFTlib.h"
#include "../saf_utilities/saf_utilities.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */
 
#ifndef NUM_EARS
# define NUM_EARS 2
#endif
    
/* ========================================================================== */
/*                             Internal Functions                             */
/* ========================================================================== */
    
/*
 * Function: cxcorr
 * ----------------
 * Calculates the cross correlation between two vectors
 *
 * Input Arguments:
 *     a    - vector a; la x 1
 *     b    - vector b; lb x 1
 *     la   - length of vector a
 *     lb   - length of vector b
 * Output Arguments:
 *     x_ab - cross-correlation between a and b; (la + lb - 1) x 1
 */
void cxcorr(float* a,
            float* b,
            float* x_ab,
            size_t la,
            size_t lb);
 
/*
 * Function: FIRtoFilterbankCoeffs
 * -------------------------------
 * Converts and FIR filter into Filterbank Coefficients
 * Note: This is currently hard coded for a 128 hop size with hybrid mode
 * enabled (see afSTFTlib)
 *
 * Input Arguments:
 *     hIR     - time-domain FIR; FLAT: N_dirs x nCH x ir_len
 *     N_dirs  - number of FIR sets
 *     nCH     - number of channels per FIR set
 *     ir_len  - length of the FIR
 *     N_bands - number of time-frequency domain bands
 * Output Arguments:
 *     hFB     - & the FIRs as Filterbank coefficients;
 *               FLAT: N_bands x nCH x N_dirs
 */
void FIRtoFilterbankCoeffs(/* Input Arguments */
                           float* hIR,
                           int N_dirs,
                           int nCH,
                           int ir_len,
                           int N_bands,
                           /* Output Arguments */
                           float_complex** hFB);
    
    
#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */

#endif /* __HRIR_INTERNAL_H_INCLUDED__ */
