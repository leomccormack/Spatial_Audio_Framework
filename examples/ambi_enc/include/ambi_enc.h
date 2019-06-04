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
 *     ambi_enc.h (include header)
 * Description:
 *     A simple, but flexible, Ambisonic encoder (aka: Ambisonic Panner).
 * Dependencies:
 *     saf_utilities, saf_sh
 * Author, date created:
 *     Leo McCormack, 07.10.2016
 */

#ifndef __AMBI_ENC_H_INCLUDED__
#define __AMBI_ENC_H_INCLUDED__

#ifdef __cplusplus
extern "C" {
#endif
    
/***********************/
/* Presets + Constants */
/***********************/
    
#define ENABLE_MONO_PRESET
#define ENABLE_STEREO_PRESET
#define ENABLE_5PX_PRESET
#define ENABLE_7PX_PRESET
#define ENABLE_8PX_PRESET
#define ENABLE_9PX_PRESET
#define ENABLE_10PX_PRESET
#define ENABLE_11PX_PRESET
#define ENABLE_11PX_7_4_PRESET
#define ENABLE_13PX_PRESET
#define ENABLE_22PX_PRESET
#define ENABLE_AALTO_MCC_PRESET
#define ENABLE_AALTO_APAJA_PRESET
#define ENABLE_AALTO_APAJA2_PRESET
#define ENABLE_AALTO_LR_PRESET
#define ENABLE_DTU_AVIL_PRESET
#define ENABLE_T_DESIGN_4_PRESET
#define ENABLE_T_DESIGN_12_PRESET
#define ENABLE_T_DESIGN_24_PRESET
#define ENABLE_T_DESIGN_36_PRESET
#define ENABLE_T_DESIGN_48_PRESET
#define ENABLE_T_DESIGN_60_PRESET

typedef enum _PRESETS{
    PRESET_DEFAULT = 1
#ifdef ENABLE_MONO_PRESET
    , PRESET_MONO
#endif
#ifdef ENABLE_STEREO_PRESET
    , PRESET_STEREO
#endif
#ifdef ENABLE_5PX_PRESET
    , PRESET_5PX
#endif
#ifdef ENABLE_7PX_PRESET
    , PRESET_7PX
#endif
#ifdef ENABLE_8PX_PRESET
    , PRESET_8PX
#endif
#ifdef ENABLE_9PX_PRESET
    , PRESET_9PX
#endif
#ifdef ENABLE_10PX_PRESET
    , PRESET_10PX
#endif
#ifdef ENABLE_11PX_PRESET
    , PRESET_11PX
#endif
#ifdef ENABLE_11PX_7_4_PRESET
    , PRESET_11PX_7_4
#endif
#ifdef ENABLE_13PX_PRESET
    , PRESET_13PX
#endif
#ifdef ENABLE_22PX_PRESET
    , PRESET_22PX
#endif
#ifdef ENABLE_AALTO_MCC_PRESET
    , PRESET_AALTO_MCC
#endif
#ifdef ENABLE_AALTO_APAJA_PRESET
    , PRESET_AALTO_APAJA
#endif
#ifdef ENABLE_AALTO_APAJA2_PRESET
    , PRESET_AALTO_APAJA2
#endif
#ifdef ENABLE_AALTO_LR_PRESET
    , PRESET_AALTO_LR
#endif
#ifdef ENABLE_DTU_AVIL_PRESET
    , PRESET_DTU_AVIL
#endif
#ifdef ENABLE_T_DESIGN_4_PRESET
    , PRESET_T_DESIGN_4
#endif
#ifdef ENABLE_T_DESIGN_12_PRESET
    , PRESET_T_DESIGN_12
#endif
#ifdef ENABLE_T_DESIGN_24_PRESET
    , PRESET_T_DESIGN_24
#endif
#ifdef ENABLE_T_DESIGN_36_PRESET
    , PRESET_T_DESIGN_36
#endif
#ifdef ENABLE_T_DESIGN_48_PRESET
    , PRESET_T_DESIGN_48
#endif
#ifdef ENABLE_T_DESIGN_60_PRESET
    , PRESET_T_DESIGN_60
#endif
    
}PRESETS;
    
#define AMBI_ENC_MAX_SH_ORDER ( 7 )
typedef enum _OUTPUT_ORDERS{ 
    OUTPUT_ORDER_FIRST=1,
    OUTPUT_ORDER_SECOND,
    OUTPUT_ORDER_THIRD,
    OUTPUT_ORDER_FOURTH,
    OUTPUT_ORDER_FIFTH,
    OUTPUT_ORDER_SIXTH,
    OUTPUT_ORDER_SEVENTH
    
}OUTPUT_ORDERS;
    
#define AMBI_ENC_NUM_CH_ORDERINGS ( 2 )
typedef enum _CH_ORDER{
    CH_ACN = 1,
    CH_FUMA     /* first-order only */
}CH_ORDER;

#define AMBI_ENC_NUM_NORM_TYPES ( 3 )
typedef enum _NORM_TYPES{
    NORM_N3D = 1,
    NORM_SN3D,
    NORM_FUMA   /* first-order only */
}NORM_TYPES;
    
#define AMBI_ENC_MAX_NUM_INPUTS ( 64 )
    
    
/******************/
/* Main Functions */
/******************/

/* creates an instance of ambi_enc */
void ambi_enc_create(void** const phAmbi);              /* address of ambi_enc handle */


/* destroys an instance of ambi_enc */
void ambi_enc_destroy(void** const phAmbi);             /* address of ambi_enc handle */


/* initialises an instance of ambi_enc */
void ambi_enc_init(void* const hAmbi,                   /* ambi_enc handle */
                     int samplerate);                   /* host sample rate */
    
    
/* performs ambisonic encoding of the input signals */
void ambi_enc_process(void* const hAmbi,                /* ambi_enc handle */
                        float** const inputs,           /* input channels, [nInputs][nSampes] */
                        float** const outputs,          /* output channels, [nOutputs][nSampes] */
                        int nInputs,                    /* number of channels in 'inputs' matrix */
                        int nOutputs,                   /* number of channels in 'outputs' matrix */
                        int nSamples,                   /* number of samples in 'inputs' matrix */
                        int isPlaying);                 /* flag, 1: inputs is not empty */

    
/*****************/
/* Set Functions */
/*****************/
    
void ambi_enc_refreshParams(void* const hAmbi);
    
void ambi_enc_setOutputOrder(void* const hAmbi, int newType);

void ambi_enc_setSourceAzi_deg(void* const hAmbi, int index, float newAzi_deg);

void ambi_enc_setSourceElev_deg(void* const hAmbi, int index, float newElev_deg);

void ambi_enc_setNumSources(void* const hAmbi, int new_nSources);
 
void ambi_enc_setInputConfigPreset(void* const hAmbi, int newPresetID);
 
void ambi_enc_setChOrder(void* const hAmbi, int newOrder);
    
void ambi_enc_setNormType(void* const hAmbi, int newType);

    
/*****************/
/* Get Functions */
/*****************/
    
int ambi_enc_getOutputOrder(void* const hAmbi);

float ambi_enc_getSourceAzi_deg(void* const hAmbi, int index);

float ambi_enc_getSourceElev_deg(void* const hAmbi, int index);

int ambi_enc_getNumSources(void* const hAmbi);

int ambi_enc_getMaxNumSources(void);
    
int ambi_enc_getNSHrequired(void* const hAmbi);
    
int ambi_enc_getChOrder(void* const hAmbi);
    
int ambi_enc_getNormType(void* const hAmbi);

    
#ifdef __cplusplus
}
#endif


#endif /* __AMBI_ENC_H_INCLUDED__ */


