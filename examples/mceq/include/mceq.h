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
 *     mceq.h (include header)
 * Description:
 *     Multi-channel equaliser.
 * Dependencies:
 *     saf_utilities, afSTFTlib,
 * Author, date created:
 *     Leo McCormack, 21.03.2018
 */

#ifndef __MCEQ_H_INCLUDED__
#define __MCEQ_H_INCLUDED__

#ifdef __cplusplus
extern "C" {
#endif

/******************/
/* Main Functions */
/******************/

/* creates an instance of the mceq */
void mceq_create(void** const phMEQ);               /* address of mceq handle */


/* destroys an instance of mceq */
void mceq_destroy(void** const phBp);               /* address of mceq handle */


/* initialises an instance of mceq */
void mceq_init(void* const hMEQ,                    /* mceq handle */
               int samplerate);                     /* host sample rate */
    
    
/* Equalise the multichannel audio */
void mceq_process(void* const hMEQ,                 /* mceq handle */
                  float** const inputs,             /* input channels [NUM_INPUTS][FRAME_SIZE] */
                  float** const outputs,            /* spherical harmonic signals */
                  int nInputs,                      /* number of channels in 'inputs' matrix */
                  int nOutputs,                     /* number of channels in 'outputs' matrix */
                  int nSamples,                     /* number of samples in 'inputs' and 'outputs' matrices */
                  int isPlaying);                   /* flag; set to 1 if there really is audio */
    
    
/*****************/
/* Set Functions */
/*****************/
    
void mceq_setNumChannels(void* const hMEQ, int newValue);
    
void mceq_setNumFilters(void* const hMEQ, int newValue);
    
void mceq_setFc(void* const hMEQ, float newValue, int filterIndex);
    
void mceq_addFilter(void* const hMEQ);
    
    
/*****************/
/* Get Functions */
/*****************/
    
int mceq_getNumChannels(void* const hMEQ);
    
int mceq_getNumFilters(void* const hMEQ);
 
float mceq_getFc(void* const hMEQ, int filterIndex);
    

#ifdef __cplusplus
}
#endif


#endif /* __MCEQ_H_INCLUDED__ */





