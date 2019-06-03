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
 *     panner.h (include header)
 * Description:
 *     A frequency-dependent 3D panner, based on the Vector-base Amplitude Panning (VBAP)
 *     method. Depending on the room, it may be beneficial to utilise amplitude-normalised
 *     gains for low frequencies, and energy-normalised gains for high frequencies; which
 *     this implemenation takes into account with one parameter "DTT". Set "DTT" to 0 for a
 *     normal room, 0.5 for listening room, and 1 for anechoic.
 * Dependencies:
 *     saf_utilities, saf_vbap, afSTFTlib
 * Author, date created:
 *     Leo McCormack, 25.09.2017
 */

#ifndef __PANNER_H_INCLUDED__
#define __PANNER_H_INCLUDED__

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


/******************/
/* Main Functions */
/******************/

/* creates an instance of the panner */
void panner_create(void** const phPan);                  /* address of panner handle */


/* destroys an instance of panner */
void panner_destroy(void** const phBp);                  /* address of panner handle */


/* initialises an instance of panner */
void panner_init(void* const hPan,                       /* panner handle */
                 int samplerate);                        /* host sample rate */
    
    
/* pan input sources using VBAP */
void panner_process(void* const hPan,                    /* panner handle */
                    float** const inputs,                /* input channels [NUM_INPUTS][FRAME_SIZE] */
                    float** const outputs,               /* spherical harmonic signals */
                    int nInputs,                         /* number of channels in 'inputs' matrix */
                    int nOutputs,                        /* number of channels in 'outputs' matrix */
                    int nSamples,                        /* number of samples in 'inputs' and 'outputs' matrices */
                    int isPlaying);                      /* flag; set to 1 if there really is audio */
    
    
/*****************/
/* Set Functions */
/*****************/

/* Set reInit Flags to 1 */
void panner_refreshSettings(void* const hPan);

/* Check if any reInit Flags are active, and reinitialise if they are. Only call when playback has stopped. */
void panner_checkReInit(void* const hPan);
    
void panner_setSourceAzi_deg(void* const hPan, int index, float newAzi_deg);

void panner_setSourceElev_deg(void* const hPan, int index, float newElev_deg);

void panner_setNumSources(void* const hPan, int new_nSources);
    
void panner_setLoudspeakerAzi_deg(void* const hPan, int index, float newAzi_deg);

void panner_setLoudspeakerElev_deg(void* const hPan, int index, float newElev_deg);

void panner_setNumLoudspeakers(void* const hPan, int new_nLoudspeakers);
    
void panner_setOutputConfigPreset(void* const hPan, int newPresetID);

void panner_setInputConfigPreset(void* const hPan, int newPresetID);
    
void panner_setDTT(void* const hPan, float newValue);
    
void panner_setSpread(void* const hPan, float newValue);
    
void panner_setYaw(void* const hPan, float newYaw);

void panner_setPitch(void* const hPan, float newPitch);

void panner_setRoll(void* const hPan, float newRoll);

void panner_setFlipYaw(void* const hPan, int newState);

void panner_setFlipPitch(void* const hPan, int newState);

void panner_setFlipRoll(void* const hPan, int newState);
    
    
/*****************/
/* Get Functions */
/*****************/
    
float panner_getSourceAzi_deg(void* const hPan, int index);

float panner_getSourceElev_deg(void* const hPan, int index);

int panner_getNumSources(void* const hPan);

int panner_getMaxNumSources(void);
    
float panner_getLoudspeakerAzi_deg(void* const hPan, int index);

float panner_getLoudspeakerElev_deg(void* const hPan, int index);

int panner_getNumLoudspeakers(void* const hPan);

int panner_getMaxNumLoudspeakers(void);

int panner_getDAWsamplerate(void* const hPan);
    
float panner_getDTT(void* const hPan);
    
float panner_getSpread(void* const hPan);
    
float panner_getYaw(void* const hPan);

float panner_getPitch(void* const hPan);

float panner_getRoll(void* const hPan);

int panner_getFlipYaw(void* const hPan);

int panner_getFlipPitch(void* const hPan);

int panner_getFlipRoll(void* const hPan);
    
int panner_getProcessingDelay(void);
    

#ifdef __cplusplus
}
#endif


#endif /* __PANNER_H_INCLUDED__ */





