/*
 Copyright 2019 Leo McCormack
 
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
 *     beamformer.h (include header)
 * Description:
 *     Generates beamformers/virtual microphones in arbitrary directions. Several
 *     different beam pattern types are included.
 * Dependencies:
 *     saf_utilities, afSTFTlib, saf_sh
 * Author, date created:
 *     Leo McCormack, 17.05.2019
 */

#ifndef __BEAMFORMER_H_INCLUDED__
#define __BEAMFORMER_H_INCLUDED__

#ifdef __cplusplus
extern "C" {
#endif
    
/***********************/
/* Presets + Constants */
/***********************/
    
#define BEAMFORMER_MAX_SH_ORDER ( 7 )
typedef enum _BEAM_ORDERS{
    BEAM_ORDER_FIRST = 1,
    BEAM_ORDER_SECOND,
    BEAM_ORDER_THIRD,
    BEAM_ORDER_FOURTH,
    BEAM_ORDER_FIFTH,
    BEAM_ORDER_SIXTH,
    BEAM_ORDER_SEVENTH
    
}BEAM_ORDERS;
    
#define BEAMFORMER_NUM_BEAM_TYPES ( 3 )
typedef enum _BEAM_TYPES {
    BEAM_TYPE_CARDIOID = 1,
    BEAM_TYPE_HYPERCARDIOID,
    BEAM_TYPE_MAX_EV
    
} BEAM_TYPES;
 
#define BEAMFORMER_NUM_CH_ORDERINGS ( 2 )
typedef enum _CH_ORDER{
    CH_ACN = 1,
    CH_FUMA     /* first-order only */
}CH_ORDER;

#define BEAMFORMER_NUM_NORM_TYPES ( 3 )
typedef enum _NORM_TYPES{
    NORM_N3D = 1,
    NORM_SN3D,
    NORM_FUMA   /* first-order only */
}NORM_TYPES;
    
#define BEAMFORMER_MAX_NUM_BEAMS ( 64 )
    
    
/******************/
/* Main Functions */
/******************/
    
/* creates an instance of beamformer */
void beamformer_create(void** const phBeam);            /* address of beamformer handle */

/* destroys an instance of beamformer */
void beamformer_destroy(void** const phBeam);           /* address of beamformer handle */

/* initialises an instance of beamformer */
void beamformer_init(void* const hBeam,                 /* beamformer handle */
                     int samplerate);                   /* host sample rate */
    
/* decodes the input spherical harmonic signals to loudspeaker or headphone signals */
void beamformer_process(void* const hBeam,              /* beamformer handle */
                        float** const inputs,           /* input channels, [nInputs][nSampes] */
                        float** const outputs,          /* output channels, [nOutputs][nSampes] */
                        int nInputs,                    /* number of channels in 'inputs' matrix */
                        int nOutputs,                   /* number of channels in 'outputs' matrix */
                        int nSamples,                   /* number of samples in 'inputs' matrix */
                        int isPlaying);                 /* flag, 1: if there is signal in the buffers */

    
/*****************/
/* Set Functions */
/*****************/

/* Set reInit Flags to 1 */
void beamformer_refreshSettings(void* const hBeam);

/* Check if any reInit Flags are active, and reinitialise if they are. Only call when playback has stopped. */
void beamformer_checkReInit(void* const hBeam);
    
void beamformer_setBeamOrder(void* const hBeam,  int newValue);

void beamformer_setBeamAzi_deg(void* const hBeam, int index, float newAzi_deg);

void beamformer_setBeamElev_deg(void* const hBeam, int index, float newElev_deg);

void beamformer_setNumBeams(void* const hBeam, int new_nBeams);
 
void beamformer_setChOrder(void* const hBeam, int newOrder);

void beamformer_setNormType(void* const hBeam, int newType);
    
void beamformer_setBeamType(void* const hBeam, int newID);

    
/*****************/
/* Get Functions */
/*****************/
    
int beamformer_getBeamOrder(void* const hBeam);
    
int beamformer_getNumberOfBands(void);

float beamformer_getBeamAzi_deg(void* const hBeam, int index);

float beamformer_getBeamElev_deg(void* const hBeam, int index);

int beamformer_getNumBeams(void* const hBeam);

int beamformer_getMaxNumBeams(void);
    
int  beamformer_getNSHrequired(void* const hBeam);
   
int beamformer_getChOrder(void* const hBeam);

int beamformer_getNormType(void* const hBeam);
    
int beamformer_getBeamType(void* const hBeam); 
    
    
#ifdef __cplusplus
}
#endif


#endif /* __BEAMFORMER_H_INCLUDED__ */


