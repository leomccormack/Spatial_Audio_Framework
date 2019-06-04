/*
 Copyright 2017-2018 Leo McCormack
 
 Permission to use, copy, modify, and/or distribute this software for any purpose with or
 without fee is hereby granted, provided that the above copyright notice and this permission
 notice appear in all copies.
 
 THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES WITH REGARD TO
 THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS. IN NO EVENT
 SHALL THE AUTHOR BE LIABLE FOR ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR
 ANY DAMAGES WHATSOEVER RESULTING FROM LOSS O    F USE, DATA OR PROFITS, WHETHER IN AN ACTION OF
 CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE
 OR PERFORMANCE OF THIS SOFTWARE.
 */
/*
 * Filename:
 *     ambi_drc.h (include header)
 * Description:
 *     A frequency-dependent spherical harmonic domain dynamic range compressor (DRC). The
 *     implementation can also keep track of the frequency-dependent gain factors for
 *     the omnidirectional component over time, for optional plotting. The design utilises
 *     a similar approach as in:
 *         McCormack, L., & Välimäki, V. (2017). "FFT-Based Dynamic Range Compression". in
 *         Proceedings of the 14th Sound and Music Computing Conference, July 5-8, Espoo,
 *         Finland.
 *     The DRC gain factors are determined based on analysing the omnidirectional component.
 *     These gain factors are then applied to the higher-order components, in a such a manner
 *     as to retain the spatial information within them.
 * Dependencies:
 *     saf_utilities, afSTFTlib
 * Author, date created:
 *     Leo McCormack, 07.01.2017
 */

#ifndef __AMBI_DRC_H_INCLUDED__
#define __AMBI_DRC_H_INCLUDED__

#ifdef __cplusplus
extern "C" {
#endif
    
#define ENABLE_TF_DISPLAY

/***********************/
/* Presets + Constants */
/***********************/

#define HOP_SIZE ( 128 )                   /* STFT hop size, can be flexible, but only 'hybrid' mode afSTFT is supported (i.e. non uniform) */
#define TIME_SLOTS ( FRAME_SIZE/HOP_SIZE ) /* time-frequency domain frame size */
#define HYBRID_BANDS ( HOP_SIZE + 5 )      /* hybrid mode incurs an additional 5 bands  */
#define SPECTRAL_FLOOR (0.1585)            /* -16dB, maximum gain reduction for a given frequency band */
#define AMBI_DRC_MAX_SH_ORDER ( 7 )
#define MAX_ORDER ( AMBI_DRC_MAX_SH_ORDER )
#define MAX_NUM_SH_SIGNALS ( (MAX_ORDER+1)*(MAX_ORDER+1) )
#ifdef ENABLE_TF_DISPLAY
  #define NUM_DISPLAY_SECONDS   ( 8 )      /* How many seconds the display will show historic TF data */
  #define NUM_DISPLAY_TIME_SLOTS ( (int)(NUM_DISPLAY_SECONDS*48000.0f/(float)HOP_SIZE) )
  #define READ_OFFSET ( 200 )
#endif
    
#define AMBI_DRC_NUM_CH_ORDERINGS ( 2 )
typedef enum _CH_ORDER{
    CH_ACN = 1,
    CH_FUMA     /* first-order only */
}CH_ORDER;

#define AMBI_DRC_NUM_NORM_TYPES ( 3 )
typedef enum _NORM_TYPES{
    NORM_N3D = 1,
    NORM_SN3D,
    NORM_FUMA   /* first-order only */
}NORM_TYPES;
    
typedef enum _INPUT_ORDER{
    INPUT_ORDER_1 = 1,
    INPUT_ORDER_2,
    INPUT_ORDER_3,
    INPUT_ORDER_4,
    INPUT_ORDER_5,
    INPUT_ORDER_6,
    INPUT_ORDER_7
    
}INPUT_ORDER;
    
#define AMBI_DRC_IN_GAIN_MIN_VAL ( -40.0f )
#define AMBI_DRC_IN_GAIN_MAX_VAL ( 20.0f )
#define AMBI_DRC_THRESHOLD_MIN_VAL ( -60.0f )
#define AMBI_DRC_THRESHOLD_MAX_VAL ( 0.0f )
#define AMBI_DRC_RATIO_MIN_VAL ( 1.0f )
#define AMBI_DRC_RATIO_MAX_VAL ( 30.0f )
#define AMBI_DRC_KNEE_MIN_VAL ( 0.0f )
#define AMBI_DRC_KNEE_MAX_VAL ( 10.0f )
#define AMBI_DRC_ATTACK_MIN_VAL ( 10.0f )
#define AMBI_DRC_ATTACK_MAX_VAL ( 200.0f )
#define AMBI_DRC_RELEASE_MIN_VAL ( 50.0f )
#define AMBI_DRC_RELEASE_MAX_VAL ( 1000.0f )
#define AMBI_DRC_OUT_GAIN_MIN_VAL ( -20.0f )
#define AMBI_DRC_OUT_GAIN_MAX_VAL ( 40.0f )
    
    
/******************/
/* Main Functions */
/******************/

/* creates an instance of ambi_drc */
void ambi_drc_create(void** const phAmbi);    /* address of ambi_drc handle */

/* destroys an instance of ambi_drc */
void ambi_drc_destroy(void** const phAmbi);   /* address of ambi_drc handle */

/* initialises an instance of ambi_drc */
void ambi_drc_init(void* const hAmbi,         /* ambi_drc handle */
                   int samplerate);           /* host sample rate */
    
/* performs the frequency-dependent multi-channel compression on input signals */
void ambi_drc_process(void* const hAmbi,      /* ambi_drc handle */
                      float** const inputs,   /* input channels, [nCh][nSampes] */
                      float** const outputs,  /* output channels, [nOutputs][nSampes] */
                      int nCH,                /* number of channels in 'inputs'/'outputs' matrix */
                      int nSamples,           /* number of samples in 'inputs'/'outputs' matrix */
                      int isPlaying);         /* Flag, 1: if there is audio in buffers */


/*****************/
/* Set Functions */
/*****************/
    
void ambi_drc_refreshSettings(void* const hAmbi);

void ambi_drc_setThreshold(void* const hAmbi, float newValue);

void ambi_drc_setRatio(void* const hAmbi, float newValue);

void ambi_drc_setKnee(void* const hAmbi, float newValue);

void ambi_drc_setInGain(void* const hAmbi, float newValue);
    
void ambi_drc_setOutGain(void* const hAmbi, float newValue);
    
void ambi_drc_setAttack(void* const hAmbi, float newValue);

void ambi_drc_setRelease(void* const hAmbi, float newValue);
    
void ambi_drc_setChOrder(void* const hAmbi, int newOrder);

void ambi_drc_setNormType(void* const hAmbi, int newType);
    
void ambi_drc_setInputPreset(void* const hAmbi, INPUT_ORDER newPreset);


/*****************/
/* Get Functions */
/*****************/

#ifdef ENABLE_TF_DISPLAY
float** ambi_drc_getGainTF(void* const hAmbi);

int ambi_drc_getGainTFwIdx(void* const hAmbi);
    
int ambi_drc_getGainTFrIdx(void* const hAmbi);
    
float* ambi_drc_getFreqVector(void* const hAmbi, int* nFreqPoints);
#endif
    
float ambi_drc_getThreshold(void* const hAmbi);

float ambi_drc_getRatio(void* const hAmbi);

float ambi_drc_getKnee(void* const hAmbi);

float ambi_drc_getInGain(void* const hAmbi);
    
float ambi_drc_getOutGain(void* const hAmbi);

float ambi_drc_getAttack(void* const hAmbi);

float ambi_drc_getRelease(void* const hAmbi);

int ambi_drc_getChOrder(void* const hAmbi);

int ambi_drc_getNormType(void* const hAmbi);
    
INPUT_ORDER ambi_drc_getInputPreset(void* const hAmbi);
    
int ambi_drc_getNSHrequired(void* const hAmbi);
    
int ambi_drc_getSamplerate(void* const hAmbi);
    
int ambi_drc_getProcessingDelay(void);
    
    
#ifdef __cplusplus
}
#endif


#endif /* __AMBI_DRC_H_INCLUDED__ */


