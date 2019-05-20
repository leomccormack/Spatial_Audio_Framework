/*
 Copyright 2019 Leo McCormack
 
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
 *     beamformer_internal.h  
 * Description:
 *     Generates beamformers/virtual microphones in arbitrary directions. Several
 *     different beam pattern types are included.
 * Dependencies:
 *     saf_utilities, afSTFTlib, saf_sh
 * Author, date created:
 *     Leo McCormack, 17.05.2019
 */

#ifndef __BEAMFORMER_INTERNAL_H_INCLUDED__
#define __BEAMFORMER_INTERNAL_H_INCLUDED__

#include <stdio.h>
#include <math.h>
#include <string.h>
#include "beamformer.h"
#define SAF_ENABLE_AFSTFT   /* for time-frequency transform */
#define SAF_ENABLE_SH       /* for spherical harmonic weights */
#include "saf.h"

#ifdef __cplusplus
extern "C" {
#endif
    
/***************/
/* Definitions */
/***************/
    
#define HOP_SIZE ( 128 )                        /* STFT hop size = nBands */
#define HYBRID_BANDS ( HOP_SIZE + 5 )           /* hybrid mode incurs an additional 5 bands  */
#define TIME_SLOTS ( FRAME_SIZE / HOP_SIZE )    /* 4/8/16 */
#define MAX_SH_ORDER ( 7 )
#define MAX_NUM_SH_SIGNALS ( (MAX_SH_ORDER+1)*(MAX_SH_ORDER+1) ) /* Maximum number of spherical harmonic components */
#define MAX_NUM_BEAMS ( 64 )                    /* Maximum permitted channels for the VST standard */
    
    
/***********/
/* Structs */
/***********/
    
typedef struct _beamformer
{
    /* audio buffers + afSTFT time-frequency transform handle */
    float SHFrameTD[MAX_NUM_SH_SIGNALS][FRAME_SIZE];
    //float_complex SHframeTF[HYBRID_BANDS][MAX_NUM_SH_SIGNALS][TIME_SLOTS];
    //float_complex outputFrameTF[HYBRID_BANDS][MAX_NUM_BEAMS][TIME_SLOTS];
    float outputFrameTD[MAX_NUM_BEAMS][FRAME_SIZE];
    //complexVector* STFTInputFrameTF;
    //complexVector* STFTOutputFrameTF;
    //void* hSTFT;                             /* afSTFT handle */
    //int afSTFTdelay;                         /* for host delay compensation */
    //float** tempHopFrameTD;                  /* temporary multi-channel time-domain buffer of size "HOP_SIZE". */
    int fs;                                  /* host sampling rate */
    //float freqVector[HYBRID_BANDS];          /* frequency vector for time-frequency transform, in Hz */
    
    /* internal variables */
    int nSH;
    int new_nBeams;                          /* if new_nLoudpkrs != nLoudpkrs, afSTFT is reinitialised */
    int new_nSH; 
    float beamWeights[MAX_NUM_BEAMS][MAX_NUM_SH_SIGNALS];
    float_complex beamWeights_cmplx[MAX_NUM_BEAMS][MAX_NUM_SH_SIGNALS];
    
    /* flags */
    int recalc_beamWeights[MAX_NUM_BEAMS];   /* 0: no init required, 1: init required */ 
    int reInitTFT;                           /* 0: no init required, 1: init required, 2: init in progress */
    
    /* user parameters */
    int beamOrder;                           /* beam order */
    int nBeams;                              /* number of loudspeakers/virtual loudspeakers */
    float beam_dirs_deg[MAX_NUM_BEAMS][2];   /* beam directions in degrees [azi, elev] */
    BEAM_TYPES beamType;                     /* see 'BEAM_TYPES' enum */
    CH_ORDER chOrdering;                     /* only ACN is supported */
    NORM_TYPES norm;                         /* N3D or SN3D */
    
} beamformer_data;


/**********************/
/* Internal functions */
/**********************/
    
/* Initialise the filterbank */
void beamformer_initTFT(void* const hBeam);  /* beamformer handle */
    

#ifdef __cplusplus
}
#endif


#endif /* __BEAMFORMER_INTERNAL_H_INCLUDED__ */




















