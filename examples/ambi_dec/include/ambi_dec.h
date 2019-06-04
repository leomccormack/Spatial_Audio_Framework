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
 *     ambi_dec.h (include header)
 * Description:
 *     A frequency-dependent Ambisonic decoder for loudspeakers or headphones. Different
 *     decoder settings can be specified for the low and high frequencies. When utilising
 *     spherical harmonic signals derived from real microphone arrays, this implementation
 *     also allows the decoding order per frequency band to be specified. Optionally, a SOFA
 *     file may be loaded for personalised headphone listening.
 *     The algorithms utilised in this Ambisonic decoder were pieced together and developed
 *     in collaboration with Archontis Politis.
 * Dependencies:
 *     saf_utilities, afSTFTlib, saf_hoa, saf_vbap, saf_hrir, saf_sh
 * Author, date created:
 *     Leo McCormack, 07.12.2017
 */

#ifndef __AMBI_DEC_H_INCLUDED__
#define __AMBI_DEC_H_INCLUDED__

#ifdef __cplusplus
extern "C" {
#endif
    
/***********/
/* Presets */
/***********/
    
/* Microphone array options */
#define ENABLE_ZYLIA_MIC_PRESET
#define ENABLE_EIGENMIKE32_MIC_PRESET
#define ENABLE_DTU_MIC_MIC_PRESET
    
/* Loudspeaker configurations */
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
#define ENABLE_ZYLIA_LAB_PRESET
#define ENABLE_T_DESIGN_4_PRESET
#define ENABLE_T_DESIGN_12_PRESET
#define ENABLE_T_DESIGN_24_PRESET
#define ENABLE_T_DESIGN_36_PRESET
#define ENABLE_T_DESIGN_48_PRESET
#define ENABLE_T_DESIGN_60_PRESET
    
#define AMBI_DEC_MAX_SH_ORDER ( 7 )
typedef enum _MASTER_ORDERS{
    MASTER_ORDER_FIRST = 1,
    MASTER_ORDER_SECOND,
    MASTER_ORDER_THIRD,
    MASTER_ORDER_FOURTH,
    MASTER_ORDER_FIFTH,
    MASTER_ORDER_SIXTH,
    MASTER_ORDER_SEVENTH
    
}MASTER_ORDERS;
    
#define AMBI_DEC_NUM_DECODING_METHODS ( 4 )
typedef enum _DECODING_METHODS {
    DECODING_METHOD_SAD = 1,   /* Sampling Ambisonic Decoder (SAD)*/
    DECODING_METHOD_MMD,       /* Mode-Matching Decoder (MMD) */
    DECODING_METHOD_EPAD,      /* Energy-Preserving Ambisonic Decoder (EPAD) */
    DECODING_METHOD_ALLRAD     /* All-Round Ambisonic Decoder (AllRAD) */
    
} DECODING_METHODS;
    
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

typedef enum _PRESETS{
    PRESET_DEFAULT = 1 
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
#ifdef ENABLE_ZYLIA_LAB_PRESET
    , PRESET_ZYLIA_LAB
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
    
typedef enum _DIFFUSE_FIELD_EQ_APPROACH {
    AMPLITUDE_PRESERVING=1,      /* preserve omni amplitude */
    ENERGY_PRESERVING            /* preserve omni energy */
    
} DIFFUSE_FIELD_EQ_APPROACH;

#define AMBI_DEC_NUM_CH_ORDERINGS ( 2 )
typedef enum _CH_ORDER{
    CH_ACN = 1,
    CH_FUMA     /* first-order only */
}CH_ORDER;

#define AMBI_DEC_NUM_NORM_TYPES ( 3 )
typedef enum _NORM_TYPES{
    NORM_N3D = 1,
    NORM_SN3D,
    NORM_FUMA   /* first-order only */
}NORM_TYPES;
    

/******************/
/* Main Functions */
/******************/
    
/* creates an instance of ambi_dec */
void ambi_dec_create(void** const phAmbi);              /* address of ambi_dec handle */

/* destroys an instance of ambi_dec */
void ambi_dec_destroy(void** const phAmbi);             /* address of ambi_dec handle */

/* initialises an instance of ambi_dec */
void ambi_dec_init(void* const hAmbi,                   /* ambi_dec handle */
                   int samplerate);                     /* host sample rate */
    
/* decodes the input spherical harmonic signals to loudspeaker or headphone signals */
void ambi_dec_process(void* const hAmbi,                /* ambi_dec handle */
                      float** const inputs,             /* input channels, [nInputs][nSampes] */
                      float** const outputs,            /* output channels, [nOutputs][nSampes] */
                      int nInputs,                      /* number of channels in 'inputs' matrix */
                      int nOutputs,                     /* number of channels in 'outputs' matrix */
                      int nSamples,                     /* number of samples in 'inputs' matrix */
                      int isPlaying);                   /* flag, 1: if there is signal in the buffers */

    
/*****************/
/* Set Functions */
/*****************/

/* Set reInit Flags to 1 */
void ambi_dec_refreshSettings(void* const hAmbi);

/* Check if any reInit Flags are active, and reinitialise if they are. Only call when playback has stopped. */
void ambi_dec_checkReInit(void* const hAmbi);
    
void ambi_dec_setMasterDecOrder(void* const hAmbi,  int newValue);
    
void ambi_dec_setDecOrder(void* const hAmbi,  int newValue, int bandIdx);

void ambi_dec_setDecOrderAllBands(void* const hAmbi,  int newValue);

void ambi_dec_setLoudspeakerAzi_deg(void* const hAmbi, int index, float newAzi_deg);

void ambi_dec_setLoudspeakerElev_deg(void* const hAmbi, int index, float newElev_deg);

void ambi_dec_setNumLoudspeakers(void* const hAmbi, int new_nLoudspeakers);

void ambi_dec_setUseDefaultHRIRsflag(void* const hAmbi, int newState);

void ambi_dec_setBinauraliseLSflag(void* const hAmbi, int newState);
    
void ambi_dec_setSofaFilePath(void* const hAmbi, const char* path);

void ambi_dec_setSourcePreset(void* const hAmbi, int newPresetID);

void ambi_dec_setOutputConfigPreset(void* const hAmbi, int newPresetID);

void ambi_dec_setChOrder(void* const hAmbi, int newOrder);

void ambi_dec_setNormType(void* const hAmbi, int newType);
    
void ambi_dec_setDecMethod(void* const hAmbi, int index, int newID);

void ambi_dec_setDecEnableMaxrE(void* const hAmbi, int index, int newID);

void ambi_dec_setDecNormType(void* const hAmbi, int index, int newID);

void ambi_dec_setTransitionFreq(void* const hAmbi, float newValue);

    
/*****************/
/* Get Functions */
/*****************/
    
int ambi_dec_getMasterDecOrder(void* const hAmbi);
    
int ambi_dec_getDecOrder(void* const hAmbi, int bandIdx);
    
int ambi_dec_getDecOrderAllBands(void* const hAmbi);

void ambi_dec_getDecOrderHandle(void* const hAmbi,
                                float** pX_vector,
                                int** pY_values,
                                int* pNpoints);

int ambi_dec_getNumberOfBands(void);

float ambi_dec_getLoudspeakerAzi_deg(void* const hAmbi, int index);

float ambi_dec_getLoudspeakerElev_deg(void* const hAmbi, int index);

int ambi_dec_getNumLoudspeakers(void* const hAmbi);

int ambi_dec_getMaxNumLoudspeakers(void);
    
int  ambi_dec_getNSHrequired(void* const hAmbi);

float ambi_dec_getSourceAzi_deg(void* const hAmbi, int index);

float ambi_dec_getSourceElev_deg(void* const hAmbi, int index);
    
int ambi_dec_getUseDefaultHRIRsflag(void* const hAmbi);
    
int ambi_dec_getBinauraliseLSflag(void* const hAmbi);

char* ambi_dec_getSofaFilePath(void* const hAmbi);

int ambi_dec_getChOrder(void* const hAmbi);

int ambi_dec_getNormType(void* const hAmbi);
    
int ambi_dec_getDecMethod(void* const hAmbi, int index);
    
int ambi_dec_getDecEnableMaxrE(void* const hAmbi, int index);
    
int ambi_dec_getDecNormType(void* const hAmbi, int index);
    
float ambi_dec_getTransitionFreq(void* const hAmbi);
    
int ambi_dec_getHRIRsamplerate(void* const hAmbi);

int ambi_dec_getDAWsamplerate(void* const hAmbi);
    
int ambi_dec_getProcessingDelay(void);
    
#ifdef __cplusplus
}
#endif


#endif /* __SAF_AMBIDEC_H_INCLUDED__ */


