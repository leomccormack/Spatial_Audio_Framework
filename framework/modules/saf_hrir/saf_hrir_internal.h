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
 * @file saf_hrir_internal.h
 * @ingroup HRIR
 * @brief Internal header for the HRIR/HRTF processing module (#SAF_HRIR_MODULE)
 *
 * A collection of head-related impulse-response (HRIR) functions. Including
 * estimation of the interaural time differences (ITDs), conversion of HRIRs to
 * HRTF filterbank coefficients, and HRTF interpolation utilising amplitude-
 * normalised VBAP gains.
 *
 * @author Leo McCormack
 * @date 12.12.2016
 */

#ifndef __HRIR_INTERNAL_H_INCLUDED__
#define __HRIR_INTERNAL_H_INCLUDED__

#include <stdio.h>
#include <math.h> 
#include <string.h>

#include "saf_hrir.h"
#include "../saf_utilities/saf_utilities.h"
#include "saf_externals.h" 

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/* ========================================================================== */
/*                             Internal Functions                             */
/* ========================================================================== */

/**
 * Converts FIR filters into Filterbank Coefficients by passing them through
 * afSTFT
 *
 * @param[in]  hIR        Time-domain FIR; FLAT: N_dirs x nCH x ir_len
 * @param[in]  N_dirs     Number of FIR sets
 * @param[in]  nCH        Number of channels per FIR set
 * @param[in]  ir_len     Length of the FIR
 * @param[in]  hopSize    Hop size
 * @param[in]  hybridmode 0: disabled, 1:enabled
 * @param[out] hFB        The FIRs as Filterbank coefficients;
 *                        FLAT: N_bands x nCH x N_dirs
 */
void FIRtoFilterbankCoeffs(/* Input Arguments */
                           float* hIR,
                           int N_dirs,
                           int nCH,
                           int ir_len,
                           int hopSize,
                           int hybridmode,
                           /* Output Arguments */
                           float_complex* hFB);


#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */

#endif /* __HRIR_INTERNAL_H_INCLUDED__ */
