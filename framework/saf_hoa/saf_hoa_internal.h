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
 *     saf_hoa_internal.h
 * Description:
 *     A Higher-order Ambisonics C library; largely derived from the MatLab library by
 *     Archontis Politis: https://github.com/polarch/Higher-Order-Ambisonics
 * Dependencies:
 *     saf_utilities, saf_sh, saf_vbap
 * Author, date created:
 *     Leo McCormack, 19.03.2018
 */

#ifndef __SAF_HOA_INTERNAL_H_INCLUDED__
#define __SAF_HOA_INTERNAL_H_INCLUDED__

#include <stdio.h>
#include <math.h> 
#include <string.h>
#include "saf_hoa.h"
#include "saf_sh.h"        /* for computing the spherical harmonics */
#include "saf_vbap.h"      /* for vbap gains utilised by allrad */
#include "saf_utilities.h" /* for linear algebra speed-ups */

#ifdef __cplusplus
extern "C" {
#endif
 
/* returns the "Energy preserving Ambisonic decoder" detailed in:
 * Zotter, F., Pomberger, H., Noisternig, M. (2012). Energy-Preserving Ambisonic Decoding. Acta Acustica United with Acustica, 98(1), 37:47.
 * The function has been written to also work when the number of spherical harmonic components exceeds the number of loudspeakers.
 * In which case, the 'U' matrix from the SVD is truncated instead. However, ideally... nLS > nSH, like in the paper */
void getEPAD(int order,                   /* decoding order */
             float* ls_dirs_deg,          /* loudspeaker directions in degrees [azi elev]; FLAT: nLS x 2 */
             int nLS,                     /* number of loudspeakers */
             float **decMtx);             /* & decoding matrix; FLAT: nLS x (order+1)^2 */
    
    
/* returns the "All-round Ambisonics decoder" detailed in:
 * Zotter, F., Frank, M. (2012). All-Round Ambisonic Panning and Decoding. Journal of the Audio Engineering Society, 60(10), 807:820 */
void getAllRAD(int order,                 /* decoding order */
               float* ls_dirs_deg,        /* loudspeaker directions in degrees [azi elev]; FLAT: nLS x 2 */
               int nLS,                   /* number of loudspeakers */
               float **decMtx);           /* & decoding matrix; FLAT: nLS x (order+1)^2 */
    
#ifdef __cplusplus
}
#endif

#endif /* __SAF_HOA_INTERNAL_H_INCLUDED__ */




















