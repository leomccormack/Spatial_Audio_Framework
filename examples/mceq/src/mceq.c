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
 *     mceq.c
 * Description:
 *     Multi-channel equaliser.
 * Dependencies:
 *     saf_utilities, afSTFTlib,
 * Author, date created:
 *     Leo McCormack, 21.03.2018
 */
 
#include "mceq_internal.h"

void mceq_create
(
    void ** const phMEQ
)
{
    mceq_data* pData = (mceq_data*)malloc(sizeof(mceq_data));
    if (pData == NULL) { return;/*error*/ }
    *phMEQ = (void*)pData;
    
    /* time-frequency transform + buffers */
    pData->hSTFT = NULL;
    pData->STFTInputFrameTF = NULL;
    pData->STFTOutputFrameTF = NULL;
    pData->tempHopFrameTD = NULL;
    
    /* internal parameters */
    pData->reInitTFT = 1;
     
    /* user parameters */
    pData->nFilters = 0;
    pData->nChannels = 2;
    pData->new_nChannels = pData->nChannels;
}


void mceq_destroy
(
    void ** const phMEQ
)
{
    mceq_data *pData = (mceq_data*)(*phMEQ);
	int t, ch;

	if (pData != NULL) {
        if(pData->hSTFT !=NULL)
            afSTFTfree(pData->hSTFT);
        for (t = 0; t<TIME_SLOTS; t++) {
            if(pData->STFTInputFrameTF!=NULL){
                for(ch=0; ch< pData->nChannels; ch++) {
                    free(pData->STFTInputFrameTF[t][ch].re);
                    free(pData->STFTInputFrameTF[t][ch].im);
                }
            }
            if(pData->STFTOutputFrameTF!=NULL){
                for (ch = 0; ch< pData->nChannels; ch++) {
                    free(pData->STFTOutputFrameTF[t][ch].re);
                    free(pData->STFTOutputFrameTF[t][ch].im);
                }
            }
        }
        if(pData->STFTInputFrameTF!=NULL)
            free2d((void**)pData->STFTInputFrameTF, TIME_SLOTS);
        if(pData->STFTOutputFrameTF!=NULL)
            free2d((void**)pData->STFTOutputFrameTF, TIME_SLOTS);
        if(pData->tempHopFrameTD!=NULL)
            free2d((void**)pData->tempHopFrameTD, pData->nChannels);
     
        free(pData);
        pData = NULL;
    }
}

void mceq_init
(
    void * const hMEQ,
    int          sampleRate
)
{
    mceq_data *pData = (mceq_data*)(hMEQ);
    int band;
    float freq_step;
    
    /* define frequency vectors */
    pData->fs = sampleRate;
    freq_step = (float)pData->fs/(2.0f*(float)HOP_SIZE);
    for(band=0; band <NUM_BANDS; band++)
        pData->freqVector[band] = (float)band * freq_step;
    freq_step = (float)pData->fs/(2.0f*(float)DISPLAY_FREQ_RES);
    for(band=0; band < NUM_DISPLAY_FREQS; band++)
        pData->freqVector[band] = (float)band * freq_step;
    
    mceq_addFilter(hMEQ);  
}

void mceq_process
(
    void  *  const hMEQ,
    float ** const inputs,
    float ** const outputs,
    int            nInputs,
    int            nOutputs,
    int            nSamples,
    int            isPlaying
)
{
    mceq_data *pData = (mceq_data*)(hMEQ);
    int t, sample, ch, i, band;
    float mag, arg;
    
    /* tmp parameters */
    
    /* reinitialise afSTFT if needed */
    if(pData->reInitTFT){
        pData->reInitTFT = 2;
        mceq_initTFT(hMEQ);
        pData->reInitTFT = 0;
    }
    if ((nSamples == FRAME_SIZE) && (isPlaying == 1) && (pData->reInitTFT == 0) ) {
        
        /* Load time-domain data */
        for(i=0; i < MIN(pData->nChannels,nInputs); i++)
            memcpy(pData->inputFrameTD[i], inputs[i], FRAME_SIZE * sizeof(float));
        for(; i<pData->nChannels; i++)
            memset(pData->inputFrameTD[i], 0, FRAME_SIZE * sizeof(float));
        
        /* Apply time-frequency transform (TFT) */
        for ( t=0; t< TIME_SLOTS; t++) {
            for( ch=0; ch < pData->nChannels; ch++)
                for ( sample=0; sample < HOP_SIZE; sample++)
                    pData->tempHopFrameTD[ch][sample] = pData->inputFrameTD[ch][sample + t*HOP_SIZE];
            afSTFTforward(pData->hSTFT, (float**)pData->tempHopFrameTD, (complexVector*)pData->STFTInputFrameTF[t]);
        }
        for( ch=0; ch < pData->nChannels; ch++)
            for ( t=0; t<TIME_SLOTS; t++)
                for(band=0; band<NUM_BANDS; band++)
                    pData->inputframeTF[ch][t][band] = cmplxf(pData->STFTInputFrameTF[t][ch].re[band], pData->STFTInputFrameTF[t][ch].im[band]);
   
        /* apply EQ */
        for( ch=0; ch < pData->nChannels; ch++){
            for ( t=0; t<TIME_SLOTS; t++){
                for(band=0; band<NUM_BANDS; band++){
                    mag = cabsf(pData->inputframeTF[ch][t][band]);
                    arg = atan2f(cimagf(pData->inputframeTF[ch][t][band]), crealf(pData->inputframeTF[ch][t][band])); 
                    pData->outputframeTF[ch][t][band] = ccmulf(cmplxf(pData->filters[0].FBmag[band] * mag,0.0f), cexpf(cmplxf(0.0f, arg)));
                }
            }
        }
        
        /* inverse-TFT */
        for (ch = 0; ch < pData->nChannels; ch++) {
            for (t = 0; t < TIME_SLOTS; t++) {
                for (band = 0; band < NUM_BANDS; band++) {
                    pData->STFTOutputFrameTF[t][ch].re[band] = crealf(pData->outputframeTF[ch][t][band]);
                    pData->STFTOutputFrameTF[t][ch].im[band] = cimagf(pData->outputframeTF[ch][t][band]);
                }
            }
        }
        for (t = 0; t < TIME_SLOTS; t++) {
            afSTFTinverse(pData->hSTFT, pData->STFTOutputFrameTF[t], pData->tempHopFrameTD);
            for (ch = 0; ch < MIN(pData->nChannels, nOutputs); ch++)
                for (sample = 0; sample < HOP_SIZE; sample++)
                    outputs[ch][sample + t* HOP_SIZE] = pData->tempHopFrameTD[ch][sample];
            for (; ch < nOutputs; ch++)
                for (sample = 0; sample < HOP_SIZE; sample++)
                    outputs[ch][sample + t* HOP_SIZE] = 0.0f;
        }
    }
    else{
        for (ch=0; ch < nOutputs; ch++)
            memset(outputs[ch],0, FRAME_SIZE*sizeof(float));
    } 
}


/* Set Functions */

void mceq_setNumChannels(void* const hMEQ, int newValue)
{
    mceq_data *pData = (mceq_data*)(hMEQ);
    pData->new_nChannels = newValue;
    if(pData->new_nChannels != pData->nChannels)
        pData->reInitTFT = 1;
}

void mceq_setNumFilters(void* const hMEQ, int newValue)
{
    mceq_data *pData = (mceq_data*)(hMEQ);
    
    
}

void mceq_setFc(void* const hMEQ, float newValue, int filterIndex)
{
    mceq_data *pData = (mceq_data*)(hMEQ);
}

void mceq_addFilter(void* const hMEQ)
{
    mceq_data *pData = (mceq_data*)(hMEQ);
    int fIdx;
    
    pData->nFilters++; 
    fIdx = pData->nFilters-1;
    pData->filters[fIdx].type = FILTER_PEAK;
    pData->filters[fIdx].fc = 1000.0f;
    pData->filters[fIdx].Q = 0.7071f;
    pData->filters[fIdx].G = 0.0f;
    
    mceq_initFilter(&(pData->filters[fIdx]), pData->freqVector_n, pData->disp_freqVector_n, (float)(pData->fs+0.5f));
}


/* Get Functions */
 
int mceq_getNumChannels(void* const hMEQ)
{
    mceq_data *pData = (mceq_data*)(hMEQ);
    return pData->nChannels;
}







    
    
