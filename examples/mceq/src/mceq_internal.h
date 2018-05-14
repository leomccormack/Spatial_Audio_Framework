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
 *     mceq_internal.h
 * Description:
 *     Multi-channel equaliser.
 * Dependencies:
 *     saf_utilities, afSTFTlib,
 * Author, date created:
 *     Leo McCormack, 21.03.2018
 */

#ifndef __MCEQ_INTERNAL_H_INCLUDED__
#define __MCEQ_INTERNAL_H_INCLUDED__

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "mceq.h"
#define SAF_ENABLE_AFSTFT /* for time-frequency transform */
#include "saf.h"

#ifdef __cplusplus
extern "C" {
#endif

/***************
   Definitions
 ***************/
    
#define HOP_SIZE ( FRAME_SIZE )                             /* 256 */
#define NUM_BANDS ( HOP_SIZE + 1 )                          /* number of frequency bands for processing */
#define DISPLAY_FREQ_RES ( 2048 )
#define NUM_DISPLAY_FREQS ( DISPLAY_FREQ_RES + 1 )                          /* frequency resolution for  */
#define TIME_SLOTS ( FRAME_SIZE / HOP_SIZE )                /* >=1 */
#define MAX_NUM_CHANNELS ( 64 )                             /* Maximum permited channels for the VST standard */
#define MAX_NUM_FILTERS ( 10 )                              /* number of filters allowed */

typedef enum _FILTER_TYPES
{
    FILTER_LPF,
    FILTER_HPF,
    FILTER_HI_SHELF,
    FILTER_LO_SHELF,
    FILTER_PEAK
    
}FILTER_TYPES;
    
/***********
   Structs
 ***********/
    
typedef struct _filter
{
    float a[3];
    float b[3];
    FILTER_TYPES type;
    float fc, Q, G;
    float FBmag[NUM_BANDS];
    float disp_mags[NUM_DISPLAY_FREQS];
    int ID, bypassFLAG;
    
}filter;

typedef struct _mceq
{
    /* audio buffers */
    float inputFrameTD[MAX_NUM_CHANNELS][FRAME_SIZE];
    float outframeTD[MAX_NUM_CHANNELS][FRAME_SIZE];
    float_complex inputframeTF[MAX_NUM_CHANNELS][TIME_SLOTS][NUM_BANDS];
    float_complex outputframeTF[MAX_NUM_CHANNELS][TIME_SLOTS][NUM_BANDS];
    complexVector** STFTInputFrameTF;
    complexVector** STFTOutputFrameTF;
    float** tempHopFrameTD;
    int fs;
    
    /* time-frequency transform */
    void* hSTFT;
    
    /* internal parameters */ 
    float freqVector[NUM_BANDS];     /* frequency vector for processing */
    float freqVector_n[NUM_BANDS];   /* normalised frequency vector for processing */
    float disp_freqVector[NUM_DISPLAY_FREQS];
    float disp_freqVector_n[NUM_DISPLAY_FREQS];
    int reInitTFT; 
    int new_nChannels;
    
    /* user parameters */
    filter filters[MAX_NUM_FILTERS];
    int nChannels;
    int nFilters;
    
} mceq_data;
     

/**********************
   Internal functions
 **********************/
    
/* Initialise the filterbank used by mceq */
void mceq_initTFT(void* const hMEQ);                       /* mceq handle */
    
/*  */
void mceq_initFilter(filter* f,/* filter struct, with .type, .fc, .Q, and .G pre-defined */
                     float freqVector_n[NUM_BANDS],
                     float disp_freqVector_n[NUM_DISPLAY_FREQS],
                     float fs);
 
#ifdef __cplusplus
}
#endif

#endif /* __MCEQ_INTERNAL_H_INCLUDED__ */






