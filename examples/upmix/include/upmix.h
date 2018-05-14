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
 *     upmix.h (include header)
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

#ifndef __UPMIX_H_INCLUDED__
#define __UPMIX_H_INCLUDED__

#ifdef __cplusplus
extern "C" {
#endif

/******************/
/* Main Functions */
/******************/

/* creates an instance of the upmix */
void upmix_create(void** const phUpmx);                /* address of upmix handle */

/* destroys an instance of upmix */
void upmix_destroy(void** const phUpmx);               /* address of upmix handle */

/* initialises an instance of upmix */
void upmix_init(void* const hUpmx,                     /* upmix handle */
                int samplerate);                       /* host sample rate */
    
/* Apply upmixing */
void upmix_process(void* const hUpmx,                  /* upmix handle */
                   float** const inputs,               /* input channels; nInputs x nSamples */
                   float** const outputs,              /* output channels; nOutputs x nSamples */
                   int nInputs,                        /* number of channels in 'inputs' matrix */
                   int nOutputs,                       /* number of channels in 'outputs' matrix */
                   int nSamples,                       /* number of samples in 'inputs' and 'outputs' matrices */
                   int isPlaying);                     /* flag; set to 1 if there really is audio */
    
    
/*****************/
/* Set Functions */
/*****************/
    
void upmix_setPValueCoeff(void* const hUpmx, float newValue);
    
void upmix_setParamAvgCoeff(void* const hUpmx, float newValue);

void upmix_setScaleDoAwidth(void* const hUpmx, float newValue);

void upmix_setCovAvg(void* const hUpmx, float newValue);

    
/*****************/
/* Get Functions */
/*****************/
    
float upmix_getPValueCoeff(void* const hUpmx);

float upmix_getParamAvgCoeff(void* const hUpmx);

float upmix_getScaleDoAwidth(void* const hUpmx);

float upmix_getCovAvg(void* const hUpmx);
    

#ifdef __cplusplus
}
#endif


#endif /* __UPMIX_H_INCLUDED__ */





