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
 *     rotator.h (include header)
 * Description:
 *     A simple spherical harmonic domain rotator.
 * Dependencies:
 *     saf_utilities, saf_sh
 * Author, date created:
 *     Leo McCormack, 02.11.2017
 */

#ifndef __ROTATOR_H_INCLUDED__
#define __ROTATOR_H_INCLUDED__

#ifdef __cplusplus
extern "C" {
#endif
    
/* creates an instance of rotator */
void rotator_create(void** const phRot);                /* address of rotator handle */

/* destroys an instance of rotator */
void rotator_destroy(void** const phRot);               /* address of rotator handle */

/* initialises an instance of rotator */
void rotator_init(void* const hRot,                     /* rotator handle */
                  int samplerate);                      /* host sample rate */
    
/* performs rotation of the SH signals */
void rotator_process(void* const hRot,                  /* rotator handle */
                     float** const inputs,              /* input channels, [nInputs][nSampes] */
                     float** const outputs,             /* output channels, [nOutputs][nSampes] */
                     int nInputs,                       /* number of channels in 'inputs' matrix */
                     int nOutputs,                      /* number of channels in 'outputs' matrix */
                     int nSamples,                      /* number of samples in 'inputs' matrix */
	                 int isPlaying);                    /* set to 1 if there is audio in the buffers */

/* SETS */

void rotator_setYaw(void* const hRot, float newYaw);
    
void rotator_setPitch(void* const hRot, float newPitch);
    
void rotator_setRoll(void* const hRot, float newRoll);
    
void rotator_setFlipYaw(void* const hRot, int newState);
    
void rotator_setFlipPitch(void* const hRot, int newState);
    
void rotator_setFlipRoll(void* const hRot, int newState);
    
void rotator_setChOrder(void* const hRot, int newOrder);

void rotator_setNormType(void* const hRot, int newType);
    
    
/* GETS */
    
float rotator_getYaw(void* const hRot);
    
float rotator_getPitch(void* const hRot);
    
float rotator_getRoll(void* const hRot);
    
int rotator_getFlipYaw(void* const hRot);
    
int rotator_getFlipPitch(void* const hRot);
    
int rotator_getFlipRoll(void* const hRot);
    
int rotator_getChOrder(void* const hRot);

int rotator_getNormType(void* const hRot);

    
#ifdef __cplusplus
}
#endif


#endif /* __ROTATOR_H_INCLUDED__ */


