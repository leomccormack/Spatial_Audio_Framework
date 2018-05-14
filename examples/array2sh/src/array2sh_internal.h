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
 *     array2sh_internal.h
 * Description:
 *     Spatially encodes spherical or cylindrical sensor array signals into spherical harmonic
 *     signals utilising theoretical encoding filters.
 *     The algorithms within array2sh were pieced together and developed in collaboration
 *     with Symeon Delikaris-Manias.
 *     A detailed explanation of the algorithms in array2sh can be found in:
 *     McCormack, L., Delikaris-Manias, S., Farina, A., Pinardi, D., and Pulkki, V.,
 *     “Real-time conversion of sensor array signals into spherical harmonic signals with
 *     applications to spatially localised sub-band sound-field analysis,” in Audio
 *     Engineering Society Convention 144, Audio Engineering Society, 2018.
 * Dependencies:
 *     saf_utilities, afSTFTlib, saf_sh
 * Author, date created:
 *     Leo McCormack, 13.09.2017
 */

#ifndef __ARRAY2SH_INTERNAL_H_INCLUDED__
#define __ARRAY2SH_INTERNAL_H_INCLUDED__

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include "array2sh.h"
#include "array2sh_database.h"
#define SAF_ENABLE_AFSTFT /* for time-frequency transform */
#define SAF_ENABLE_SH     /* for spherical harmonic weights */
#include "saf.h"

#ifdef __cplusplus
extern "C" {
#endif

/***************/
/* Definitions */
/***************/
    
#define NUM_SH_SIGNALS ( (SH_ORDER + 1)*(SH_ORDER + 1)  )   /* (L+1)^2 */
#define HOP_SIZE ( 128 )                                    /* STFT hop size = nBands */
#define HYBRID_BANDS ( HOP_SIZE + 5 )                       /* hybrid mode incurs an additional 5 bands  */
#define TIME_SLOTS ( FRAME_SIZE / HOP_SIZE )                /* 4/8/16 */
#define MAX_NUM_SENSORS ( 64 )                              /* Maximum permited channels for the VST standard */
#define MAX_EVAL_FREQ_HZ ( 20e3f )                          /* Up to which frequency should the evaluation be accurate */
    
/***********/
/* Structs */
/***********/
    
typedef struct _arrayPars {
    int Q, newQ;
    float r;
    float R;
    float admittance;
    ARRAY_TYPES arrayType;
    WEIGHT_TYPES weightType;
    float sensorCoords_rad[MAX_NUM_SENSORS][2];
    float sensorCoords_deg[MAX_NUM_SENSORS][2];
        
}arrayPars;

typedef struct _array2sh
{
    /* audio buffers */
    float inputFrameTD[MAX_NUM_SENSORS][FRAME_SIZE];
    float SHframeTD[NUM_SH_SIGNALS][FRAME_SIZE];
    float_complex inputframeTF[HYBRID_BANDS][MAX_NUM_SENSORS][TIME_SLOTS];
    float_complex SHframeTF[HYBRID_BANDS][NUM_SH_SIGNALS][TIME_SLOTS];
    complexVector** STFTInputFrameTF;
    complexVector** STFTOutputFrameTF;
    float** tempHopFrameTD;
    
    /* intermediates */
    double_complex bN_modal[HYBRID_BANDS][SH_ORDER + 1];
    double_complex bN[HYBRID_BANDS][SH_ORDER + 1];
    double_complex bN_inv[HYBRID_BANDS][SH_ORDER + 1];
    double_complex bN_inv_R[HYBRID_BANDS][NUM_SH_SIGNALS];
    double Y[NUM_SH_SIGNALS][MAX_NUM_SENSORS];
    double YYT[NUM_SH_SIGNALS][NUM_SH_SIGNALS];
    float_complex Y_cmplx[NUM_SH_SIGNALS][MAX_NUM_SENSORS];
    float_complex W[HYBRID_BANDS][NUM_SH_SIGNALS][MAX_NUM_SENSORS];
    
    /* for displaying the bNs */
    float** bN_modal_dB;
    float** bN_inv_dB;
    float* cSH;
    float* lSH; 
    float* disp_freqVector;
    
    /* time-frequency transform and array details */
    float freqVector[HYBRID_BANDS];
    void* hSTFT;
    void* arraySpecs;
    
    /* additional user parameters that are not included in the array presets */
    PRESETS preset;
    REG_TYPES regType;
    float regPar;
    CH_ORDER chOrdering;
    NORM_TYPES norm;
    float c;
    float gain_dB;
    float maxFreq; 
    
    int reinitSHTmatrixFLAG;
    int reinitTFTFLAG;
    int recalcEvalFLAG;

} array2sh_data;
     

/**********************/
/* Internal functions */
/**********************/
    
void array2sh_initTFT(void* const hA2sh);
    
void array2sh_calculate_sht_matrix(void* const hA2sh);
    
void array2sh_calculate_mag_curves(void* const hA2sh);
    
void array2sh_evaluateSHTfilters(void* hA2sh);
    
void array2sh_createArray(void ** const hPars);

void array2sh_destroyArray(void ** const hPars);
    
void array2sh_initArray(void * const hPars,
                        PRESETS preset,
                        int firstInitFLAG);

#ifdef __cplusplus
}
#endif

#endif /* __ARRAY2SH_INTERNAL_H_INCLUDED__ */

