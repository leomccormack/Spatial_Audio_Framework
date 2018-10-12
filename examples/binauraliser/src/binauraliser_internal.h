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
 *     binauraliser_internal.h
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

#ifndef __BINAURALISER_INTERNAL_H_INCLUDED__
#define __BINAURALISER_INTERNAL_H_INCLUDED__

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "binauraliser.h"
#include "binauraliser_database.h"
#define SAF_ENABLE_AFSTFT  /* for time-frequency transform */
#define SAF_ENABLE_HRIR    /* for HRIR->HRTF filterbank coefficients conversion etc. */
#define SAF_ENABLE_VBAP    /* for amplitude-normalised VBAP gains used for interpolating HRTFs */
#define SAF_ENABLE_SH
#include "saf.h"

#ifdef __cplusplus
extern "C" {
#endif

/***************/
/* Definitions */
/***************/
      
#define HOP_SIZE ( 128 )                                    /* STFT hop size = nBands */
#define HYBRID_BANDS ( HOP_SIZE + 5 )                       /* hybrid mode incurs an additional 5 bands  */
#define TIME_SLOTS ( FRAME_SIZE / HOP_SIZE )                /* 4/8/16 */
#define MAX_NUM_INPUTS ( 64 )                               /* Maximum permited channels for the VST standard */
#define MAX_NUM_OUTPUTS ( 64 )                              /* Maximum permited channels for the VST standard */
#define NUM_EARS ( 2 )                                      /* true for most humans */
 
#ifndef DEG2RAD
  #define DEG2RAD(x) (x * PI / 180.0f)
#endif
#ifndef RAD2DEG
  #define RAD2DEG(x) (x * 180.0f / PI)
#endif
    
/***********/
/* Structs */
/***********/

typedef struct _binauraliser
{
    /* audio buffers */
    float inputFrameTD[MAX_NUM_INPUTS][FRAME_SIZE];
    float outframeTD[NUM_EARS][FRAME_SIZE];
    float_complex inputframeTF[HYBRID_BANDS][MAX_NUM_INPUTS][TIME_SLOTS];
    float_complex outputframeTF[HYBRID_BANDS][NUM_EARS][TIME_SLOTS];
    complexVector** STFTInputFrameTF;
    complexVector** STFTOutputFrameTF;
    float** tempHopFrameTD;
    int fs;
    
    /* time-frequency transform */
    float freqVector[HYBRID_BANDS]; 
    void* hSTFT;
    
    /* sofa file info */
    char* sofa_filepath; 
    float* hrirs;
    float* hrir_dirs_deg;
    int N_hrir_dirs;
    int hrir_len;
    int hrir_fs;
    
    /* vbap gain table */
    int hrtf_vbapTableRes[2];
    int N_hrtf_vbap_gtable;
    int* hrtf_vbap_gtableIdx; /* N_hrtf_vbap_gtable x 3 */
    float* hrtf_vbap_gtableComp; /* N_hrtf_vbap_gtable x 3 */
    
    /* hrir filterbank coefficients */
    int useDefaultHRIRsFLAG; 
    float* itds_s; /* interaural-time differences for each HRIR (in seconds); nBands x 1 */
    float_complex* hrtf_fb; /* hrtf filterbank coefficients; nBands x nCH x N_hrirs */
    float* hrtf_fb_mag; /* magnitudes of the hrtf filterbank coefficients; nBands x nCH x N_hrirs */
    float_complex hrtf_interp[MAX_NUM_OUTPUTS][HYBRID_BANDS][NUM_EARS];
    
    /* flags */
    int recalc_hrtf_interpFLAG[MAX_NUM_OUTPUTS];
    int reInitHRTFsAndGainTables;
    int reInitTFT;
    
    /* misc. */
    float src_dirs_rot_deg[MAX_NUM_INPUTS][2];
    float src_dirs_rot_xyz[MAX_NUM_INPUTS][3];
    float src_dirs_xyz[MAX_NUM_INPUTS][3]; 
    int nTriangles;
    int input_nDims;  
    int output_nDims;
    
    /* user parameters */
    int nSources;
    int new_nSources;
    float src_dirs_deg[MAX_NUM_INPUTS][2];
    INTERP_MODES interpMode;
    int enableRotation;
    float yaw, roll, pitch;                  /* rotation angles in degrees */
    int bFlipYaw, bFlipPitch, bFlipRoll;     /* flag to flip the sign of the individual rotation angles */
    int useRollPitchYawFlag;                 /* rotation order flag, 1: r-p-y, 0: y-p-r */
    
} binauraliser_data;
     

/**********************/
/* Internal functions */
/**********************/
    
/* interpolates between 3 HRTFs via AN-VBAP gains. The HRTF magnitude responses and HRIR ITDs are interpolated seperately
 * before being re-combined */
void binauraliser_interpHRTFs(void* const hPan,                    /* pannerlib handle (includes VBAP gains, HRTFs and ITDs) */
                              float azimuth_deg,                   /* source azimuth in degrees */
                              float elevation_deg,                 /* source elevation in degrees */
                              float_complex h_intrp[HYBRID_BANDS][NUM_EARS]);
    
/* Initialise the HRTFs: either loading the default set or loading from a SOFA file, Then generate a VBAP gain table. */
void binauraliser_initHRTFsAndGainTables(void* const hBin);        /* binauraliser handle */
    
/* Initialise the filterbank used by binauraliser */
void binauraliser_initTFT(void* const hBin);                       /* binauraliser handle */
    
/* Loads directions from preset */
void binauraliser_loadPreset(PRESETS preset,                       /* PRESET enum */
                             float dirs_deg[MAX_NUM_INPUTS][2],    /* source/loudspeaker directions */
                             int* newNCH,                          /* & new number of channels */
                             int* nDims);                          /* & estimate of the number of dimensions (2 or 3) */

    
#ifdef __cplusplus
}
#endif

#endif /* __BINAURALISER_INTERNAL_H_INCLUDED__ */
