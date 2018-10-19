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
 *     array2sh.h (include header)
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

#ifndef __ARRAY2SH_H_INCLUDED__
#define __ARRAY2SH_H_INCLUDED__

#ifdef __cplusplus
extern "C" {
#endif
    
#define ENABLE_AALTO_HYDROPHONE_PRESET
#define ENABLE_SENNHEISER_AMBEO_PRESET
#define ENABLE_CORE_SOUND_TETRAMIC_PRESET
#define ENABLE_SOUND_FIELD_SPS200_PRESET
#define ENABLE_ZYLIA_1D_PRESET
#define ENABLE_EIGENMIKE32_PRESET
#define ENABLE_DTU_MIC_PRESET
    
typedef enum _ENCODING_ORDERS{
    ENCODING_ORDER_FIRST = 1,
    ENCODING_ORDER_SECOND,
    ENCODING_ORDER_THIRD,
    ENCODING_ORDER_FOURTH,
    ENCODING_ORDER_FIFTH,
    ENCODING_ORDER_SIXTH,
    ENCODING_ORDER_SEVENTH
    
}ENCODING_ORDERS;
    
typedef enum _PRESETS{
    PRESET_DEFAULT = 1
#ifdef ENABLE_AALTO_HYDROPHONE_PRESET
    ,PRESET_AALTO_HYDROPHONE
#endif
#ifdef ENABLE_SENNHEISER_AMBEO_PRESET
    ,PRESET_SENNHEISER_AMBEO
#endif
#ifdef ENABLE_CORE_SOUND_TETRAMIC_PRESET
    ,PRESET_CORE_SOUND_TETRAMIC
#endif
#ifdef ENABLE_SOUND_FIELD_SPS200_PRESET
    ,PRESET_SOUND_FIELD_SPS200
#endif
#ifdef ENABLE_ZYLIA_1D_PRESET
    ,PRESET_ZYLIA_1D
#endif
#ifdef ENABLE_EIGENMIKE32_PRESET
    ,PRESET_EIGENMIKE32
#endif
#ifdef ENABLE_DTU_MIC_PRESET
    ,PRESET_DTU_MIC
#endif
}PRESETS;

typedef enum _REG_TYPES{
    REG_DAS = 1,
    REG_SOFT_LIM,
    REG_TIKHONOV
}REG_TYPES;

typedef enum _CH_ORDER{
    CH_ACN = 1,
    CH_SID
}CH_ORDER;

typedef enum _NORM_TYPES{
    NORM_N3D = 1,
    NORM_SN3D
}NORM_TYPES;

typedef enum _ARRAY_TYPES{
    ARRAY_SPHERICAL = 1,
    ARRAY_CYLINDRICAL
}ARRAY_TYPES;

typedef enum _WEIGHT_TYPES{
    WEIGHT_RIGID = 1,
    WEIGHT_OPEN_OMNI,
    WEIGHT_OPEN_CARD,
    WEIGHT_OPEN_DIPOLE
}WEIGHT_TYPES;
    
/******************/
/* Main Functions */
/******************/  

/* creates an instance of the array2sh */
void array2sh_create(void** const phA2sh);              /* address of hA2sh handle */

/* destroys an instance of array2sh */
void array2sh_destroy(void** const phA2sh);             /* address of hA2sh handle */

/* initialises an instance of array2sh */
void array2sh_init(void* const hA2sh,                   /* hA2sh handle */
                   int samplerate);                     /* host sample rate */
    
/* spatially encode microphone/hydrophone array signals into spherical harmonic signals */
void array2sh_process(void* const hA2sh,                /* hA2sh handle */
                      float** const inputs,             /* input channels [NUM_INPUTS][FRAME_SIZE] */
                      float** const outputs,            /* spherical harmonic signals */
                      int nInputs,                      /* number of channels in 'inputs' matrix */
                      int nOutputs,                     /* number of channels in 'outputs' matrix */
                      int nSamples,                     /* number of samples in 'inputs' and 'outputs' matrices */
                      int isPlaying);                   /* flag, 1: if inputs actually has audio in it */
    
    
/*****************/
/* Set Functions */
/*****************/

/* Set reInit Flags to 1 */
void array2sh_refreshSettings(void* const hA2sh);

/* Check if any reInit Flags are active, and reinitialise if they are. Only call when playback has stopped. */
void array2sh_checkReInit(void* const hA2sh);
    
void array2sh_setEncodingOrder(void* const hA2sh, int newOrder);
    
void array2sh_evaluateFilters(void* const hA2sh);
    
void array2sh_setPreset(void* const hA2sh, int preset);
    
void array2sh_setSensorAzi_rad(void* const hA2sh, int index, float newAzi_rad);
    
void array2sh_setSensorElev_rad(void* const hA2sh, int index, float newElev_rad);
    
void array2sh_setSensorAzi_deg(void* const hA2sh, int index, float newAzi_deg);
    
void array2sh_setSensorElev_deg(void* const hA2sh, int index, float newElev_deg);
    
void array2sh_setNumSensors(void* const hA2sh, int newQ);
    
void array2sh_setr(void* const hA2sh, float newr);
    
void array2sh_setR(void* const hA2sh, float newR);
    
void array2sh_setArrayType(void* const hA2sh, int newType);

void array2sh_setWeightType(void* const hA2sh, int newType);
    
void array2sh_setRegType(void* const hA2sh, int newType);
    
void array2sh_setRegPar(void* const hA2sh, float newVal);
    
void array2sh_setChOrder(void* const hA2sh, int newOrder);
    
void array2sh_setNormType(void* const hA2sh, int newType);

void array2sh_setc(void* const hA2sh, float newc);
    
void array2sh_setGain(void* const hA2sh, float newGain);
    
void array2sh_setMaxFreq(void* const hA2sh, float newF);

    
/*****************/
/* Get Functions */
/*****************/

/* returns 1 if there are new evaluation curves */
int array2sh_getEvalReady(void* const hA2sh);
    
int array2sh_getEncodingOrder(void* const hA2sh);
    
float array2sh_getSensorAzi_rad(void* const hA2sh, int index);
    
float array2sh_getSensorElev_rad(void* const hA2sh, int index);
    
float array2sh_getSensorAzi_deg(void* const hA2sh, int index);
    
float array2sh_getSensorElev_deg(void* const hA2sh, int index);

int array2sh_getNumSensors(void* const hA2sh);
    
int array2sh_getMaxNumSensors(void);
    
int array2sh_getMinNumSensors(void* const hA2sh);
    
int array2sh_getNSHrequired(void* const hA2sh);
    
float array2sh_getr(void* const hA2sh);
    
float array2sh_getR(void* const hA2sh);
    
int array2sh_getArrayType(void* const hA2sh);
    
int array2sh_getWeightType(void* const hA2sh);
    
int array2sh_getRegType(void* const hA2sh);
    
float array2sh_getRegPar(void* const hA2sh);
    
int array2sh_getChOrder(void* const hA2sh);
    
int array2sh_getNormType(void* const hA2sh);
    
float array2sh_getc(void* const hA2sh);
    
float array2sh_getGain(void* const hA2sh);
    
float array2sh_getMaxFreq(void* const hA2sh);
    
float* array2sh_getFreqVector(void* const hA2sh, int* nFreqPoints);
    
float** array2sh_getbN_inv(void* const hA2sh, int* nCurves, int* nFreqPoints);
    
float** array2sh_getbN_modal(void* const hA2sh, int* nCurves, int* nFreqPoints);
    
float* array2sh_getSpatialCorrelation_Handle(void* const hA2sh, int* nCurves, int* nFreqPoints);

float* array2sh_getLevelDifference_Handle(void* const hA2sh, int* nCurves, int* nFreqPoints);
    
#ifdef __cplusplus
}
#endif


#endif /* __ARRAY2SH_H_INCLUDED__ */





