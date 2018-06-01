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
 *     powermap.h (include header)
 * Description:
 *     A powermap-based sound-field visualiser, which utilises spherical harmonic
 *     signals as input.
 * Dependencies:
 *     saf_utilities, afSTFTlib, saf_vbap, saf_sh
 * Author, date created:
 *     Leo McCormack, 26.04.2016
 */

#ifndef __POWERMAP_H_INCLUDED__
#define __POWERMAP_H_INCLUDED__

#ifdef __cplusplus
extern "C" {
#endif
 
/***********/
/* Presets */
/***********/

/* Microphone/Hydrophone array options */
#define ENABLE_ZYLIA_MIC_PRESET
#define ENABLE_EIGENMIKE32_MIC_PRESET
#define ENABLE_DTU_MIC_MIC_PRESET

typedef enum _MIC_PRESETS{
    MIC_PRESET_IDEAL = 1
#ifdef ENABLE_ZYLIA_MIC_PRESET
    , MIC_PRESET_ZYLIA
#endif
#ifdef ENABLE_EIGENMIKE32_MIC_PRESET
    , MIC_PRESET_EIGENMIKE32
#endif
#ifdef ENABLE_DTU_MIC_MIC_PRESET
    , MIC_PRESET_DTU_MIC
#endif
}MIC_PRESETS;

typedef enum _CH_ORDER{
    CH_ACN = 1
}CH_ORDER;

typedef enum _NORM_TYPES{
    NORM_N3D = 1,
    NORM_SN3D
}NORM_TYPES;

typedef enum _POWERMAP_MODES {
    PM_MODE_PWD = 1,
    PM_MODE_MVDR,
    PM_MODE_CROPAC_LCMV,
    PM_MODE_MUSIC,
    PM_MODE_MUSIC_LOG, 
    PM_MODE_MINNORM,
    PM_MODE_MINNORM_LOG
    
} POWERMAP_MODES;
    
typedef enum _HFOV_OPTIONS{
    HFOV_360 = 1
    
}HFOV_OPTIONS;
    
typedef enum _ASPECT_RATIO_OPTIONS{
    ASPECT_RATIO_2_1 = 1
    
}ASPECT_RATIO_OPTIONS;
    

/******************/
/* Main Functions */
/******************/  

/* creates an instance of powermap */
void powermap_create(void** const phPm);                     /* address of powermap handle */

/* destroys an instance of accropac */
void powermap_destroy(void** const phPm);                    /* address of powermap handle */

/* initialises an instance of powermap */
void powermap_init(void* const hPm,                          /* powermap handle */
                   float  samplerate);                       /* host sample rate */
    
/* applies powermap analysis to input frame  */
void powermap_analysis(void* const hPm,                      /* powermap handle */
                       float** const inputs,                 /* input channels [NUM_SH_SIGNALS][FRAME_SIZE] */
                       int nInputs,                          /* number of channels in 'inputs' matrix */
                       int nSamples,                         /* number of samples in 'inputs' and 'outputs' matrices */
                       int isPlaying);                       /* flag, 0: no audio in buffer, 1: buffers have been filled */
    
   
/*****************/
/* Set Functions */
/*****************/
    
void powermap_setPowermapMode(void* const hPm, int newMode);

void powermap_setAnaOrder(void* const hPm,  int newValue, int bandIdx);

void powermap_setAnaOrderAllBands(void* const hPm, int newValue);
    
void powermap_setPowermapEQ(void* const hPm,  float newValue, int bandIdx);
    
void powermap_setPowermapEQAllBands(void* const hPm,  float newValue);
    
void powermap_setCovAvgCoeff(void* const hPm, float newAvg);

void powermap_setChOrder(void* const hPm, int newOrder);

void powermap_setNormType(void* const hPm, int newType);

void powermap_setSourcePreset(void* const hPm, int newPresetID);
    
void powermap_setNumSources(void* const hPm, int newValue);
    
void powermap_setDispFOV(void* const hPm, int newOption);
    
void powermap_setAspectRatio(void* const hPm, int newOption);
    
void powermap_setPowermapAvgCoeff(void* const hPm, float newValue);
    
void powermap_requestPmapUpdate(void* const hPm);

void powermap_refreshSettings(void* const hPm);
    
    
/*****************/
/* Get Functions */
/*****************/
    
int powermap_getPowermapMode(void* const hPm);

float powermap_getSamplingRate(void* const hPm);

float powermap_getCovAvgCoeff(void* const hPm);

int powermap_getNumberOfBands(void);

float powermap_getPowermapEQ(void* const hPm, int bandIdx);

float powermap_getPowermapEQAllBands(void* const hPm);

void powermap_getPowermapEQHandle(void* const hPm,
                                  float** pX_vector,
                                  float** pY_values,
                                  int* pNpoints);
    
int powermap_getAnaOrder(void* const hPm, int bandIdx);

int powermap_getAnaOrderAllBands(void* const hPm);

void powermap_getAnaOrderHandle(void* const hPm,
                                float** pX_vector,
                                int** pY_values,
                                int* pNpoints);

int powermap_getChOrder(void* const hPm);

int powermap_getNormType(void* const hPm);

int powermap_getNumSources(void* const hPm);

int powermap_getDispFOV(void* const hPm);

int powermap_getAspectRatio(void* const hPm);
    
float powermap_getPowermapAvgCoeff(void* const hPm);
    
//TODO: hfov and aspectRatio should be float, if 16:9 etc options are added
int powermap_getPmap(void* const hPm,
                     float** grid_dirs,
                     float** pmap,
                     int* nDirs,
                     int* pmapWidth,
                     int* hfov,
                     int* aspectRatio);


#ifdef __cplusplus
}
#endif


#endif /* __POWERMAP_H_INCLUDED__ */





