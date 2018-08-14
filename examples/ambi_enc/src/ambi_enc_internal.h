/*
 Copyright 2016-2018 Leo McCormack
 
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
 *     ambi_enc_internal.h
 * Description:
 *     A simple, but flexible, Ambisonic encoder (aka: Ambisonic Panner).
 * Dependencies:
 *     saf_utilities, saf_sh
 * Author, date created:
 *     Leo McCormack, 07.10.2016
 */

#ifndef __AMBI_ENC_INTERNAL_H_INCLUDED__
#define __AMBI_ENC_INTERNAL_H_INCLUDED__

#include <stdio.h>
#include <math.h>
#include <string.h>
#include "ambi_enc.h"
#define SAF_ENABLE_SH /* for spherical harmonic weights */
#include "saf.h"

#ifdef __cplusplus
extern "C" {
#endif
  
#define MAX_ORDER ( 7 )
#define MAX_NUM_INPUTS ( 64 )
#define MAX_NUM_SH_SIGNALS ( (MAX_ORDER + 1)*(MAX_ORDER + 1) )    /* (L+1)^2 */
    
typedef struct _ambi_enc
{
    float inputFrameTD[MAX_NUM_INPUTS][FRAME_SIZE];
    float outputFrameTD[MAX_NUM_SH_SIGNALS][FRAME_SIZE];
    float fs;
    int recalc_SH_FLAG[MAX_NUM_INPUTS];
    float Y[MAX_NUM_SH_SIGNALS][MAX_NUM_INPUTS];
    int order;
    
    /* user parameters */
    int nSources;
    int new_nSources;
    float src_dirs_deg[MAX_NUM_INPUTS][2];
    CH_ORDER chOrdering;
    NORM_TYPES norm;
    OUTPUT_ORDERS outputOrderPreset;
    
} ambi_enc_data;
    
/* Loads directions from preset */
void ambi_enc_loadPreset(PRESETS preset,                         /* PRESET enum */
                         float dirs_deg[MAX_NUM_INPUTS][2],      /* source/loudspeaker directions */
                         int* newNCH);                           /* & new number of channels */


#ifdef __cplusplus
}
#endif


#endif /* __AMBI_ENC_INTERNAL_H_INCLUDED__ */




















