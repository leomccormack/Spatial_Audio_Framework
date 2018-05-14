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
 *     upmix_internal.h
 * Description:
 *     A (soon to be) collection of upmixing algorithms. However, currently, only stereo to
 *     5.x is supported, utilising a modified version of the direct-ambient decomposition
 *     approach described in: Faller, C. (2006). Multiple-loudspeaker playback of stereo
 *     signals. Journal of the Audio Engineering Society, 54(11), 1051-1064.
 * Dependencies:
 *     saf_utilities, afSTFTlib, saf_vbap, saf_sh
 * Author, date created:
 *     Leo McCormack, 04.04.2018
 */

#ifndef __UPMIX_INTERNAL_H_INCLUDED__
#define __UPMIX_INTERNAL_H_INCLUDED__

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "upmix.h"
#include "upmix_database.h"
#define SAF_ENABLE_VBAP     /* for VBAP gains */
#define SAF_ENABLE_SH       /* for cart2sph and sph2cart */
#define SAF_ENABLE_AFSTFT   /* for time-frequency transform */
#include "saf.h"

#ifdef __cplusplus
extern "C" {
#endif

/***************/
/* Definitions */
/***************/
    
#define HOP_SIZE ( 128 )                             /* number time-domain samples to be grouped into time-frequency slot */
#define HYBRID_BANDS ( HOP_SIZE + 5 )                /* number of frequency bands for processing */
#define TIME_SLOTS ( FRAME_SIZE / HOP_SIZE )         /* >=1 */
#define MAX_NUM_INPUT_CHANNELS ( 2 )                 /* Currently, only stereo to 5.x is supported. more in future! */
#define MAX_NUM_OUTPUT_CHANNELS ( 5 )                /* Currently, only stereo to 5.x is supported. more in future! */
#define MAX_GROUP_FREQ ( 18000 )                     /* Hz, past this point bands are grouped into 1 */
#define DIFFUSE_DELAY_MS ( 30.0f )                   /* diffuse stream delay in ms */
#define DIFFUSE_DELAY_TIME_SLOTS ( 11 )              /* diffuse stream delay in time slots (11==30ms approx). (bit lazy I know..) */
    
/***********/
/* Structs */
/***********/
    
typedef struct _codecPars
{
    /* 2D VBAP gain table */
    float* grid_vbap_gtable;                           /* 2D gain table to pan the source signal to estimated azimuth */
    int grid_nPairs;                                   /* number of loudspeaker pairs in vbap gain table */
    int grid_N_vbap_gtable;                            /* number of grid directions in vbap gain table */
    int vbap_azi_res;                                  /* azimuth step size in degrees; (min: 1) */
    
    /* band grouping */
    float maxGrpFreq;                                  /* maximum frequency in Hz, past this, all bands are grouped into one band */
    int* grp_idx;                                      /* the indices that define the band grouping; nGrpBands x 1 */
    float* grp_freqs;                                  /* the group frequencies; nGrpBands x 1 */
    int nGrpBands;                                     /* number of grouped bands */
    
    /* low-pass filter */
    float diff_lpf[HYBRID_BANDS];                      /* low-pass filter applied to the diffuse stream */
    
    /* for smoothing DoA over time */
    float* prev_est_dir;                               /* degrees, previous estimated source direction per band grouping; nGrpBands x 1 */
    
}codecPars;
    
typedef struct _upmix
{
    /* temporary audio buffers */
    float inputFrameTD[MAX_NUM_INPUT_CHANNELS][FRAME_SIZE];
    float outframeTD[MAX_NUM_OUTPUT_CHANNELS][FRAME_SIZE];
    float_complex inputframeTF[HYBRID_BANDS][MAX_NUM_INPUT_CHANNELS][TIME_SLOTS];
    float_complex outputframeTF[HYBRID_BANDS][MAX_NUM_OUTPUT_CHANNELS][TIME_SLOTS];
    complexVector** STFTInputFrameTF;
    complexVector** STFTOutputFrameTF;
    float** tempHopFrameTD;
    int fs;
    
    /* circular buffer for delaying input signal */
    float_complex inputframeTF_del[HYBRID_BANDS][MAX_NUM_INPUT_CHANNELS][TIME_SLOTS];
    float_complex inputframeTF_buffer[HYBRID_BANDS][MAX_NUM_INPUT_CHANNELS][DIFFUSE_DELAY_TIME_SLOTS];
    int buffer_rIdx, buffer_wIdx;       /* circular buffer read and write indices */
    
    /* time-frequency transform */
    void* hSTFT;                        /* handle for the afSTFT time.frequency transform */
    
    /* internal parameters */
    codecPars* pars;                    /* handle for codec parameters */
    float pValues[HYBRID_BANDS];        /* VBAP normalisation coefficients per band */
    float freqVector[HYBRID_BANDS];     /* frequency vector for processing */ 
    int reInitCodec;                    /* flag. 0: no init required, 1: init required, 2: init ongoing */
    float_complex Cx[HYBRID_BANDS][MAX_NUM_INPUT_CHANNELS][MAX_NUM_INPUT_CHANNELS];
    float_complex new_Ms[HYBRID_BANDS][MAX_NUM_OUTPUT_CHANNELS][MAX_NUM_INPUT_CHANNELS];
    float_complex new_Md[HYBRID_BANDS][MAX_NUM_OUTPUT_CHANNELS][MAX_NUM_INPUT_CHANNELS];
    float_complex directframeTF[HYBRID_BANDS][MAX_NUM_OUTPUT_CHANNELS][TIME_SLOTS];
    float_complex diffuseframeTF[HYBRID_BANDS][MAX_NUM_OUTPUT_CHANNELS][TIME_SLOTS];
    
    /* user parameters */
    float loudpkrs_dirs_deg[MAX_NUM_OUTPUT_CHANNELS][2]; /* currently only stereo to 5.x is supported */
    int nLoudspeakers;       /* number of loudspeakers in the target set-up */
    float pValueCoeff;       /* pValue coefficient; 0..1; 0: normal room, 0.5: listening room, 1: anechoic */
    float paramAvgCoeff;     /* coefficient for the one-pole filter that smooths the estimated parameters over time; 0..1 */
    float scaleDoAwidth;     /* influences the stage width. 0: only centre, 0.5: -90..90 azimuth, 1: -180..180 azimuth */
    float covAvg;            /* coefficient for the one-pole filter that smooths the covarience matrix over time; 0..1 */
    
} upmix_data;
     

/**********************/
/* Internal functions */
/**********************/
    
/* Intialises the codec parameters */
void upmix_initCodec(void* const hUpmx);           /* upmix handle */

    
#ifdef __cplusplus
}
#endif

#endif /* __UPMIX_INTERNAL_H_INCLUDED__ */






