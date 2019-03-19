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
 *     sldoa.h (include header)
 * Description:
 *     A spatially-localised pressure-intensity based direction-of-arrival estimator (SLDoA).
 *     VBAP gain patterns are imposed on the spherical harmonic signals, such that the DoA
 *     can be estimated in a spatially-constrained region; thus mitigating interferes and
 *     reflections arriving from other directions. The DoA is estimated per sector for
 *     each frequency band.
 *     The algorithms within sldoa were developed in collaboration with Symeon Delikaris-
 *     Manias, and are explained in more detail in:
 *         McCormack, L., Delikaris-Manias, S., Farina, A., Pinardi, D., and Pulkki, V.,
 *         “Real-time conversion of sensor array signals into spherical harmonic signals with
 *         applications to spatially localised sub-band sound-field analysis,” in Audio
 *         Engineering Society Convention 144, Audio Engineering Society, 2018.
 * Dependencies:
 *     saf_utilities, afSTFTlib, saf_vbap, saf_sh
 * Author, date created:
 *     Leo McCormack, 18.10.2017
 */

#ifndef __SLDOA_H_INCLUDED__
#define __SLDOA_H_INCLUDED__

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
 
/* "Master order" relates to the current maximum order to expect. However,
 * the analysis order can be lower for a given frequency, due to the
 * "analysisOrderPerBand" vector, which can contain lower values than the
 * master order, but not higher. */
typedef enum _MASTER_ORDERS{
    MASTER_ORDER_FIRST = 1,
    MASTER_ORDER_SECOND,
    MASTER_ORDER_THIRD,
    MASTER_ORDER_FOURTH,
    MASTER_ORDER_FIFTH,
    MASTER_ORDER_SIXTH,
    MASTER_ORDER_SEVENTH
    
}MASTER_ORDERS;
    
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

    
/******************/
/* Main Functions */
/******************/

/* creates an instance of the mighty sldoa */
void sldoa_create(void** const phSld);               /* address of sldoa handle */

/* destroys an instance of sldoa */
void sldoa_destroy(void** const phSld);              /* address of sldoa handle */

/* initialises an instance of sldoa */
void sldoa_init(void* const hSld,                    /* sldoa handle */
                float samplerate);                   /* host sample rate */
    
/* applies sldoa analysis on input frame  */
void sldoa_analysis(void* const hSld,                /* sldoa handle */
                    float** const inputs,            /* input channels [NUM_SH_SIGNALS][FRAME_SIZE] */
                    int nInputs,                     /* number of channels in 'inputs' matrix */
                    int nSamples,                    /* number of samples in 'inputs' matrix */
                    int isPlaying);                  /* 0or1, flag to say whether playback is ongoing */

   
/*****************/
/* Set Functions */
/*****************/
    
void sldoa_setMasterOrder(void* const hSld,  int newValue);

void sldoa_setAnalysisOrder(void* const hSld, int newOrder);
    
void sldoa_setNumSectors(void* const hSld, int newSecs);
    
void sldoa_setMaxFreq(void* const hSld, float newFreq);
    
void sldoa_setMinFreq(void* const hSld, float newFreq);
    
void sldoa_setAvg(void* const hSld, float newAvg);
    
void sldoa_setAnaOrder(void* const hSld,  int newValue, int bandIdx);

void sldoa_setAnaOrderAllBands(void* const hSld,  int newValue);
    
void sldoa_setChOrder(void* const hSld, int newOrder);

void sldoa_setNormType(void* const hSld, int newType);

void sldoa_setSourcePreset(void* const hSld, int newPresetID);


/*****************/
/* Get Functions */
/*****************/
    
int sldoa_getMasterOrder(void* const hSld);
    
void sldoa_refreshSettings(void* const hSld);
    
int sldoa_getSamplingRate(void* const hSld);

float sldoa_getMaxFreq(void* const hSld);
    
float sldoa_getMinFreq(void* const hSld);

float sldoa_getAvg(void* const hSld);
    
int sldoa_getNumberOfBands(void);
    
int sldoa_getNSHrequired(void* const hSld);
    
void sldoa_getDisplayData(void *  const hSld,
                          float** pAzi_deg,
                          float** pElev_deg,
                          float** pColourScale,
                          float** pAlphaScale,
                          int** pNsectorsPerBand,
                          int* pMaxNumSectors,
                          int* pStartBand,
                          int* pEndBand);
    
int sldoa_getAnaOrder(void* const hSld, int bandIdx);

int sldoa_getAnaOrderAllBands(void* const hSld);

void sldoa_getAnaOrderHandle(void* const hSld,
                             float** pX_vector,
                             int** pY_values,
                             int* pNpoints);

int sldoa_getChOrder(void* const hSld);

int sldoa_getNormType(void* const hSld);

#ifdef __cplusplus
}
#endif

#endif /* __SLDOA_H_INCLUDED__ */





