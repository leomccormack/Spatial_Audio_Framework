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
 *     binauraliser.h (include header)
 * Description:
 *     Convolves input audio (up to 64 channels) with interpolated HRTFs in the time-frequency
 *     domain. The HRTFs are interpolated by applying amplitude-preserving VBAP gains to the
 *     HRTF magnitude responses and inter-aural time differences (ITDs) individually, before
 *     being re-combined. The example allows the user to specify an external SOFA file for the
 *     convolution.
 * Dependencies:
 *     saf_utilities, saf_hrir, saf_vbap, afSTFTlib
 * Author, date created:
 *     Leo McCormack, 25.09.2017
 */

#ifndef __BINAURALISER_H_INCLUDED__
#define __BINAURALISER_H_INCLUDED__

#ifdef __cplusplus
extern "C" {
#endif
    
/***********/
/* Presets */
/***********/
     
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
#define ENABLE_ZYLIA_LAB_PRESET
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

    
/******************/
/* Main Functions */
/******************/

/* creates an instance of the binauraliser */
void binauraliser_create(void** const phBin);            /* address of binauraliser handle */


/* destroys an instance of binauraliser */
void binauraliser_destroy(void** const phBin);           /* address of binauraliser handle */


/* initialises an instance of binauraliser */
void binauraliser_init(void* const hBin,                 /* binauraliser handle */
                       int samplerate);                  /* host sample rate */
    
    
/* pan input sources to HRIR directions using amplitude-normalised VBAP gains  */
void binauraliser_process(void* const hBin,              /* binauraliser handle */
                          float** const inputs,          /* input channels [NUM_INPUTS][FRAME_SIZE] */
                          float** const outputs,         /* spherical harmonic signals */
                          int nInputs,                   /* number of channels in 'inputs' matrix */
                          int nOutputs,                  /* number of channels in 'outputs' matrix */
                          int nSamples,                  /* number of samples in 'inputs' and 'outputs' matrices */
                          int isPlaying);                /* flag; set to 1 if there really is audio */
    
    
/*****************/
/* Set Functions */
/*****************/
    
void binauraliser_refreshSettings(void* const hBin);
    
void binauraliser_setSourceAzi_deg(void* const hBin, int index, float newAzi_deg);

void binauraliser_setSourceElev_deg(void* const hBin, int index, float newElev_deg);

void binauraliser_setNumSources(void* const hBin, int new_nSources);
 
void binauraliser_setUseDefaultHRIRsflag(void* const hBin, int newState);
    
void binauraliser_setSofaFilePath(void* const hBin, const char* path);

void binauraliser_setInputConfigPreset(void* const hBin, int newPresetID);
    
void binauraliser_setEnableRotation(void* const hBin, int newState);
    
void binauraliser_setYaw(void* const hBin, float newYaw);

void binauraliser_setPitch(void* const hBin, float newPitch);

void binauraliser_setRoll(void* const hBin, float newRoll);

void binauraliser_setFlipYaw(void* const hBin, int newState);

void binauraliser_setFlipPitch(void* const hBin, int newState);

void binauraliser_setFlipRoll(void* const hBin, int newState);

void binauraliser_setRPYflag(void* const hBin, int newState);
    

/*****************/
/* Get Functions */
/*****************/
    
float binauraliser_getSourceAzi_deg(void* const hBin, int index);

float binauraliser_getSourceElev_deg(void* const hBin, int index);

int binauraliser_getNumSources(void* const hBin);

int binauraliser_getMaxNumSources(void);
    
int binauraliser_getNDirs(void* const hBin);
    
int binauraliser_getNTriangles(void* const hBin);
    
float binauraliser_getHRIRAzi_deg(void* const hBin, int index);
    
float binauraliser_getHRIRElev_deg(void* const hBin, int index);
    
int binauraliser_getHRIRlength(void* const hBin);
    
int binauraliser_getHRIRsamplerate(void* const hBin);
    
int binauraliser_getUseDefaultHRIRsflag(void* const hBin);
    
char* binauraliser_getSofaFilePath(void* const hBin);
 
int binauraliser_getDAWsamplerate(void* const hBin);
    
int binauraliser_getEnableRotation(void* const hBin);
    
float binauraliser_getYaw(void* const hBin);

float binauraliser_getPitch(void* const hBin);

float binauraliser_getRoll(void* const hBin);

int binauraliser_getFlipYaw(void* const hBin);

int binauraliser_getFlipPitch(void* const hBin);

int binauraliser_getFlipRoll(void* const hBin);

int binauraliser_getRPYflag(void* const hBin);


#ifdef __cplusplus
}
#endif


#endif /* __BINAURALISER_H_INCLUDED__ */





