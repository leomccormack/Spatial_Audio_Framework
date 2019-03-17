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
 *     dirass.h (include header)
 * Description:
 *     A sound-field visualiser based on the directional re-assignment of beamformer energy,
 *     utilising the DoA estimates extracted from spatially-localised active-intensity
 *     (SLAI) vectors; which correspond to the scanning grid directions.
 *     For more information on the method, refer to:
 *         McCormack, L., Politis, A., and Pulkki, V. (2019). "Sharpening of angular
 *         spectra based on a directional re-assignment approach for ambisonic sound-field
 *         visualisation". IEEE International Conference on Acoustics, Speech and Signal
 *         Processing (ICASSP).
 *
 * Dependencies:
 *     saf_utilities, saf_vbap, saf_sh
 * Author, date created:
 *     Leo McCormack, 21.02.2019
 */

#ifndef __DIRASS_H_INCLUDED__
#define __DIRASS_H_INCLUDED__

#ifdef __cplusplus
extern "C" {
#endif
 
/***********/
/* Presets */
/***********/
 
typedef enum _INPUT_ORDERS{
    INPUT_ORDER_FIRST = 1,
    INPUT_ORDER_SECOND,
    INPUT_ORDER_THIRD,
    INPUT_ORDER_FOURTH,
    INPUT_ORDER_FIFTH,
    INPUT_ORDER_SIXTH,
    INPUT_ORDER_SEVENTH
    
}INPUT_ORDERS;
    
typedef enum _UPSCALE_ORDERS{
    UPSCALE_ORDER_FIRST = 1,
    UPSCALE_ORDER_SECOND,
    UPSCALE_ORDER_THIRD,
    UPSCALE_ORDER_FOURTH,
    UPSCALE_ORDER_FIFTH,
    UPSCALE_ORDER_SIXTH,
    UPSCALE_ORDER_SEVENTH,
    UPSCALE_ORDER_EIGHTH,
    UPSCALE_ORDER_NINTH,
    UPSCALE_ORDER_TENTH,
    UPSCALE_ORDER_ELEVENTH,
    UPSCALE_ORDER_TWELFTH,
    UPSCALE_ORDER_THIRTEENTH,
    UPSCALE_ORDER_FOURTEENTH,
    UPSCALE_ORDER_FIFTEENTH,
    UPSCALE_ORDER_SIXTHTEENTH,
    UPSCALE_ORDER_SEVENTEENTH,
    UPSCALE_ORDER_EIGHTEENTH,
    UPSCALE_ORDER_NINETEENTH,
    UPSCALE_ORDER_TWENTIETH
    
}UPSCALE_ORDERS;
    
typedef enum _GRID_OPTIONS{
    T_DESIGN_3 = 1,       /* 6 points */
    T_DESIGN_4,           /* 12 points */
    T_DESIGN_6,           /* 24 points */
    T_DESIGN_9,           /* 48 points */
    T_DESIGN_13,          /* 94 points */
    T_DESIGN_18,          /* 180 points */
    GRID_GEOSPHERE_6,     /* 362 points */
    T_DESIGN_30,          /* 480 points */
    GRID_GEOSPHERE_8,     /* 642 points */
    GRID_GEOSPHERE_9,     /* 812 points */
    GRID_GEOSPHERE_10,    /* 1002 points */
    GRID_GEOSPHERE_12     /* 1442 points */
    
}GRID_OPTIONS;

typedef enum _CH_ORDER{
    CH_ACN = 1
}CH_ORDER;

typedef enum _NORM_TYPES{
    NORM_N3D = 1,
    NORM_SN3D
}NORM_TYPES;

typedef enum _BEAM_TYPES {
    BEAM_TYPE_CARD = 1,
    BEAM_TYPE_HYPERCARD,
    BEAM_TYPE_MAX_EV
    
} BEAM_TYPES;
    
typedef enum _REASS_MODE {
    /* Re-assignment is disabled. i.e. just calculates (beamformer) energy-based maps */
    REASS_MODE_OFF = 1,
    /* Each sector beamformer energy is re-assigned to the nearest interpolation grid point,
     * based on the analysed DoA. More information can be found in:
     *     McCormack, L., Politis, A., and Pulkki, V. (2019). "Sharpening of angular
     *     spectra based on a directional re-assignment approach for ambisonic sound-field
     *     visualisation". IEEE International Conference on Acoustics, Speech and Signal
     *     Processing (ICASSP).*/
    REASS_NEAREST,
    /* Each sector beamformer is re-encoded into spherical harmonics of a higher order.
     * The map is then derived from the upscaled SHs as normal. More information can be found
     * in:
     *     McCormack, L., Politis, A., and Pulkki, V. (2019). "Sharpening of angular
     *     spectra based on a directional re-assignment approach for ambisonic sound-field
     *     visualisation". IEEE International Conference on Acoustics, Speech and Signal
     *     Processing (ICASSP).*/
    REASS_UPSCALE
    
} REASS_MODES;
    
typedef enum _HFOV_OPTIONS{
    HFOV_360 = 1,
    HFOV_180,
    HFOV_90,
    HFOV_60
    
}HFOV_OPTIONS;
    
typedef enum _ASPECT_RATIO_OPTIONS{
    ASPECT_RATIO_2_1 = 1,
    ASPECT_RATIO_16_9,
    ASPECT_RATIO_4_3
    
}ASPECT_RATIO_OPTIONS;
    

/******************/
/* Main Functions */
/******************/  

/* creates an instance of dirass */
void dirass_create(void** const phDir);                    /* address of dirass handle */

/* destroys an instance of accropac */
void dirass_destroy(void** const phDir);                   /* address of dirass handle */

/* initialises an instance of dirass */
void dirass_init(void* const hDir,                         /* dirass handle */
                 float  samplerate);                       /* host sample rate */
    
/* applies dirass analysis to input frame  */
void dirass_analysis(void* const hDir,                     /* dirass handle */
                     float** const inputs,                 /* input channels [NUM_SH_SIGNALS][FRAME_SIZE] */
                     int nInputs,                          /* number of channels in 'inputs' matrix */
                     int nSamples,                         /* number of samples in 'inputs' and 'outputs' matrices */
                     int isPlaying);                       /* flag, 0: no audio in buffer, 1: buffers have been filled */
    
   
/*****************/
/* Set Functions */
/*****************/

/* Set reInit Flags to 1 */
void dirass_refreshSettings(void* const hDir);

/* Check if any reInit Flags are active, and reinitialise if they are. Only call when playback has stopped. */
void dirass_checkReInit(void* const hDir);
    
void dirass_setBeamType(void* const hDir, int newMode);
    
void dirass_setInputOrder(void* const hDir,  int newValue);
    
void dirass_setDisplayGridOption(void* const hDir,  int newState);  // not thread safe!
    
void dirass_setDispWidth(void* const hDir,  int newValue); // not thread safe!
    
void dirass_setUpscaleOrder(void* const hDir,  int newState);
    
void dirass_setDiRAssMode(void* const hDir,  int newMode);
    
void dirass_setMinFreq(void* const hDir,  float newValue);

void dirass_setMaxFreq(void* const hDir,  float newValue);

void dirass_setChOrder(void* const hDir, int newOrder);

void dirass_setNormType(void* const hDir, int newType);
    
void dirass_setDispFOV(void* const hDir, int newOption);
    
void dirass_setAspectRatio(void* const hDir, int newOption);
    
void dirass_setMapAvgCoeff(void* const hDir, float newValue);
    
void dirass_requestPmapUpdate(void* const hDir);
    
    
/*****************/
/* Get Functions */
/*****************/
    
int dirass_getInputOrder(void* const hDir);
    
int dirass_getBeamType(void* const hDir);
    
int dirass_getDisplayGridOption(void* const hDir);
    
int dirass_getDispWidth(void* const hDir);

int dirass_getUpscaleOrder(void* const hDir);

int dirass_getDiRAssMode(void* const hDir); 

float dirass_getMinFreq(void* const hDir);

float dirass_getMaxFreq(void* const hDir);

int dirass_getSamplingRate(void* const hDir); 
    
int dirass_getNSHrequired(void* const hDir);

int dirass_getChOrder(void* const hDir);

int dirass_getNormType(void* const hDir);

int dirass_getDispFOV(void* const hDir);

int dirass_getAspectRatio(void* const hDir);
    
float dirass_getMapAvgCoeff(void* const hDir);
     
int dirass_getPmap(void* const hDir,
                   float** grid_dirs,
                   float** pmap,
                   int* nDirs,
                   int* pmapWidth,
                   int* hfov,
                   float* aspectRatio);


#ifdef __cplusplus
}
#endif


#endif /* __DIRASS_H_INCLUDED__ */





