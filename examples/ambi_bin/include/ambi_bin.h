/*
 Copyright 2018 Leo McCormack
 
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
 *     ambi_bin.h (include header)
 * Description:
 *     A binaural Ambisonic decoder for reproducing ambisonic signals over headphones. 
 * Dependencies:
 *     saf_utilities, afSTFTlib, saf_hrir, saf_vbap, saf_sh
 * Author, date created:
 *     Leo McCormack, 14.04.2018
 */

#ifndef __AMBI_BIN_H_INCLUDED__
#define __AMBI_BIN_H_INCLUDED__

#ifdef __cplusplus
extern "C" {
#endif
 
/***********************/
/* Presets + Constants */
/***********************/

#define AMBI_BIN_MAX_SH_ORDER ( 7 )
typedef enum _INPUT_ORDERS{
    INPUT_ORDER_FIRST = 1,
    INPUT_ORDER_SECOND,
    INPUT_ORDER_THIRD,
    INPUT_ORDER_FOURTH,
    INPUT_ORDER_FIFTH,
    INPUT_ORDER_SIXTH,
    INPUT_ORDER_SEVENTH
    
}INPUT_ORDERS;
    
#define AMBI_BIN_NUM_DECODING_METHODS ( 5 )
typedef enum _DECODING_METHODS{
    DECODING_METHOD_LS = 1,         /* Least-squares (LS) decoder */
    DECODING_METHOD_LSDIFFEQ,       /* Least-squares (LS) decoder with diffuse-field spectral equalisation */
    DECODING_METHOD_SPR,            /* Spatial resampling decoder (on the same lines as the virtual loudspeaker approach) */
    DECODING_METHOD_TA,             /* Time-alignment */
    DECODING_METHOD_MAGLS,          /* Magnitude least-squares decoder */
    
}DECODING_METHODS;

#define AMBI_BIN_NUM_CH_ORDERINGS ( 2 )
typedef enum _CH_ORDER{
    CH_ACN = 1,
    CH_FUMA     /* first-order only */
}CH_ORDER;

#define AMBI_BIN_NUM_NORM_TYPES ( 3 )
typedef enum _NORM_TYPES{
    NORM_N3D = 1,
    NORM_SN3D,
    NORM_FUMA   /* first-order only */
}NORM_TYPES;
    
    
/******************/
/* Main Functions */
/******************/
    
/* creates an instance of ambi_bin */
void ambi_bin_create(void** const phAmbi);              /* address of ambi_bin handle */

/* destroys an instance of ambi_bin */
void ambi_bin_destroy(void** const phAmbi);             /* address of ambi_bin handle */

/* initialises an instance of ambi_bin */
void ambi_bin_init(void* const hAmbi,                   /* ambi_bin handle */
                   int samplerate);                     /* host sample rate */
    
/* decodes the input spherical harmonic signals to loudspeaker or headphone signals */
void ambi_bin_process(void* const hAmbi,                /* ambi_bin handle */
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
void ambi_bin_refreshParams(void* const hAmbi);

/* Check if any reInit Flags are active, and reinitialise if they are. Only call when playback has stopped. */
void ambi_bin_checkReInit(void* const hAmbi);

void ambi_bin_setUseDefaultHRIRsflag(void* const hAmbi, int newState);
    
void ambi_bin_setSofaFilePath(void* const hAmbi, const char* path);

void ambi_bin_setInputOrderPreset(void* const hAmbi, INPUT_ORDERS newPreset);
    
void ambi_bin_setDecodingMethod(void* const hAmbi, DECODING_METHODS newMethod);

void ambi_bin_setChOrder(void* const hAmbi, int newOrder);

void ambi_bin_setNormType(void* const hAmbi, int newType);

void ambi_bin_setEnableMaxRE(void* const hAmbi, int newState);
    
void ambi_bin_setEnableDiffuseMatching(void* const hAmbi, int newState);
    
void ambi_bin_setEnablePhaseWarping(void* const hAmbi, int newState);
    
void ambi_bin_setEnableRotation(void* const hAmbi, int newState);
    
void ambi_bin_setYaw(void* const hAmbi, float newYaw);

void ambi_bin_setPitch(void* const hAmbi, float newPitch);

void ambi_bin_setRoll(void* const hAmbi, float newRoll);

void ambi_bin_setFlipYaw(void* const hAmbi, int newState);

void ambi_bin_setFlipPitch(void* const hAmbi, int newState);

void ambi_bin_setFlipRoll(void* const hAmbi, int newState);

void ambi_bin_setRPYflag(void* const hAmbi, int newState);
    
    
/*****************/
/* Get Functions */
/*****************/

int ambi_bin_getUseDefaultHRIRsflag(void* const hAmbi);
    
int ambi_bin_getInputOrderPreset(void* const hAmbi);
    
int ambi_bin_getDecodingMethod(void* const hAmbi);

char* ambi_bin_getSofaFilePath(void* const hAmbi);

int ambi_bin_getChOrder(void* const hAmbi);

int ambi_bin_getNormType(void* const hAmbi); 
    
int ambi_bin_getNumEars(void);
    
int ambi_bin_getNSHrequired(void* const hAmbi);
    
int ambi_bin_getEnableMaxRE(void* const hAmbi);

int ambi_bin_getEnableDiffuseMatching(void* const hAmbi);

int ambi_bin_getEnablePhaseWarping(void* const hAmbi);
    
int ambi_bin_getEnableRotation(void* const hAmbi);
    
float ambi_bin_getYaw(void* const hAmbi);

float ambi_bin_getPitch(void* const hAmbi);

float ambi_bin_getRoll(void* const hAmbi);

int ambi_bin_getFlipYaw(void* const hAmbi);

int ambi_bin_getFlipPitch(void* const hAmbi);

int ambi_bin_getFlipRoll(void* const hAmbi);
    
int ambi_bin_getRPYflag(void* const hAmbi);
    
int ambi_bin_getNDirs(void* const hAmbi);

int ambi_bin_getHRIRlength(void* const hAmbi);

int ambi_bin_getHRIRsamplerate(void* const hAmbi);
 
int ambi_bin_getDAWsamplerate(void* const hAmbi);
    
int ambi_bin_getProcessingDelay(void);
    
    
#ifdef __cplusplus
}
#endif

#endif /* __SAF_AMBI_BIN_H_INCLUDED__ */


