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
 *     panner_internal.h
 * Description:
 *     A frequency-dependent 3D panner, based on the Vector-base Amplitude Panning (VBAP)
 *     method. Depending on the room, it may be beneficial to utilise amplitude-normalised
 *     gains for low frequencies, and energy-normalised gains for high frequencies; which
 *     this implemenation takes into account with one parameter "DTT". Set "DTT" to 0 for a
 *     normal room, 0.5 for listening room, and 1 for anechoic.
 * Dependencies:
 *     saf_utilities, saf_vbap, afSTFTlib
 * Author, date created:
 *     Leo McCormack, 25.09.2017
 */

#ifndef __PANNER_INTERNAL_H_INCLUDED__
#define __PANNER_INTERNAL_H_INCLUDED__

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "panner.h"
#include "panner_database.h"
#define SAF_ENABLE_AFSTFT /* for time-frequency transform */
#define SAF_ENABLE_VBAP   /* for VBAP gains */
#include "saf.h"

#ifdef __cplusplus
extern "C" {
#endif

/***************/
/* Definitions */
/***************/ 
    
#define HOP_SIZE ( 128 )                                    /* STFT hop size = nBands */
#define HYBRID_BANDS ( HOP_SIZE + 5 )                       /* hybrid mode incurs an additional 5 bands  */
#define TIME_SLOTS ( FRAME_SIZE / HOP_SIZE )                /* 4/8/16 */
#define MAX_NUM_INPUTS ( 64 )                               /* Maximum permited channels for the VST standard */
#define MAX_NUM_OUTPUTS ( 64 )                              /* Maximum permited channels for the VST standard */
#define NUM_EARS ( 2 )                                      /* true for most humans */
 
    
/***********/
/* Structs */
/***********/

typedef struct _panner
{
    /* audio buffers */
    float inputFrameTD[MAX_NUM_INPUTS][FRAME_SIZE];
    float outframeTD[NUM_EARS][FRAME_SIZE];
    float_complex inputframeTF[HYBRID_BANDS][MAX_NUM_INPUTS][TIME_SLOTS];
    float_complex outputframeTF[HYBRID_BANDS][MAX_NUM_OUTPUTS][TIME_SLOTS];
    complexVector** STFTInputFrameTF;
    complexVector** STFTOutputFrameTF;
    float** tempHopFrameTD;
    int fs;
    
    /* time-frequency transform */
    float freqVector[HYBRID_BANDS];
    void* hSTFT;
    
    /* Loudspeaker version */
    int vbapTableRes[2];
    float* vbap_gtable; /* N_hrtf_vbap_gtable x nLoudpkrs */
    int N_vbap_gtable;
    int reInitGainTables;
    int reInitTFT;
    
    /* misc. */
    int nTriangles;
    int input_nDims; /* both 2D and 3D setups are supported, however, triangulation can fail if LS directions are shady */
    int output_nDims;
    
    /* pValue */
    float pValue[HYBRID_BANDS];
    
    /* user parameters */
    int nSources;
    int new_nSources;
    float src_dirs_deg[MAX_NUM_INPUTS][2];
    float DTT;
    int nLoudpkrs;
    int new_nLoudpkrs;
    float loudpkrs_dirs_deg[MAX_NUM_INPUTS][2];
    
} panner_data;
     

/**********************/
/* Internal functions */
/**********************/
    
/* Generate a VBAP gain table for current loudspeaker configuration. */
void panner_initGainTables(void* const hPan);                /* panner handle */
    
/* Initialise the filterbank used by panner */
void panner_initTFT(void* const hPan);                       /* panner handle */
    
/* Loads directions from preset */
void panner_loadPreset(PRESETS preset,                       /* PRESET enum */
                       float dirs_deg[MAX_NUM_INPUTS][2],    /* source/loudspeaker directions */
                       int* newNCH,                          /* & new number of channels */
                       int* nDims);                          /* & estimate of the number of dimensions (2 or 3) */
    
/* Calculates pValue per frequency, DTT = 1 for anechoic conditions, ~0.5 for listening rooms, 0 for standard power normalisation */
void panner_getPvalue(float DTT,                             /* pValue coefficient 0..1 */
                      float f[HYBRID_BANDS],                 /* frequency vector */
                      float p[HYBRID_BANDS] );               /* pValues per frequency */

#ifdef __cplusplus
}
#endif

#endif /* __PANNER_INTERNAL_H_INCLUDED__ */






