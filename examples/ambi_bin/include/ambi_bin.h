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
 *     Optionally, a SOFA file may be loaded for personalised headphone listening.
 * Dependencies:
 *     saf_utilities, afSTFTlib, saf_hrir, saf_sh
 * Author, date created:
 *     Leo McCormack, 14.04.2018
 */

#ifndef __AMBI_BIN_H_INCLUDED__
#define __AMBI_BIN_H_INCLUDED__

#ifdef __cplusplus
extern "C" {
#endif

/***********/
/* Presets */
/***********/
    
typedef enum _INPUT_ORDERS{
    INPUT_OMNI = 1,
    INPUT_ORDER_FIRST,
    INPUT_ORDER_SECOND,
    INPUT_ORDER_THIRD,
    INPUT_ORDER_FOURTH,
    INPUT_ORDER_FIFTH,
    INPUT_ORDER_SIXTH,
    INPUT_ORDER_SEVENTH
    
}INPUT_ORDERS;

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
    
void ambi_bin_refreshParams(void* const hAmbi);

void ambi_bin_setUseDefaultHRIRsflag(void* const hAmbi, int newState);
    
void ambi_bin_setSofaFilePath(void* const hAmbi, const char* path);

void ambi_bin_setInputOrderPreset(void* const hAmbi, INPUT_ORDERS newPreset);

void ambi_bin_setChOrder(void* const hAmbi, int newOrder);

void ambi_bin_setNormType(void* const hAmbi, int newType);

void ambi_bin_setDecEnableMaxrE(void* const hAmbi, int newState);
    
void ambi_bin_setEnableEQ(void* const hAmbi, int newState);
    
void ambi_bin_setYaw(void* const hAmbi, float newYaw);

void ambi_bin_setPitch(void* const hAmbi, float newPitch);

void ambi_bin_setRoll(void* const hAmbi, float newRoll);

void ambi_bin_setFlipYaw(void* const hAmbi, int newState);

void ambi_bin_setFlipPitch(void* const hAmbi, int newState);

void ambi_bin_setFlipRoll(void* const hAmbi, int newState);

    
/*****************/
/* Get Functions */
/*****************/
 
int ambi_bin_getUseDefaultHRIRsflag(void* const hAmbi);
    
int ambi_bin_getInputOrderPreset(void* const hAmbi);

char* ambi_bin_getSofaFilePath(void* const hAmbi);

int ambi_bin_getChOrder(void* const hAmbi);

int ambi_bin_getNormType(void* const hAmbi);
    
int ambi_bin_getDecEnableMaxrE(void* const hAmbi);
    
int ambi_bin_getEnableEQ(void* const hAmbi);
    
float ambi_bin_getYaw(void* const hAmbi);

float ambi_bin_getPitch(void* const hAmbi);

float ambi_bin_getRoll(void* const hAmbi);

int ambi_bin_getFlipYaw(void* const hAmbi);

int ambi_bin_getFlipPitch(void* const hAmbi);

int ambi_bin_getFlipRoll(void* const hAmbi);
    
int ambi_bin_getNDirs(void* const hAmbi);

int ambi_bin_getHRIRlength(void* const hAmbi);

int ambi_bin_getHRIRsamplerate(void* const hAmbi);
 
int ambi_bin_getDAWsamplerate(void* const hAmbi);
    
    
#ifdef __cplusplus
}
#endif


#endif /* __SAF_AMBI_BIN_H_INCLUDED__ */


