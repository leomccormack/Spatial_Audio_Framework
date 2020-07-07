/*
 * Copyright 2020 Leo McCormack
 *
 * Permission to use, copy, modify, and/or distribute this software for any
 * purpose with or without fee is hereby granted, provided that the above
 * copyright notice and this permission notice appear in all copies.
 *
 * THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES WITH
 * REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY
 * AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY SPECIAL, DIRECT,
 * INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM
 * LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR
 * OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR
 * PERFORMANCE OF THIS SOFTWARE.
 */

/**
 * @file decorrelator.c
 * @brief A multi-channel decorrelator
 *
 * @author Leo McCormack
 * @date 07.07.2020
 */

#include "decorrelator_internal.h"

void decorrelator_create
(
    void ** const phDecor
)
{
    decorrelator_data* pData = (decorrelator_data*)malloc1d(sizeof(decorrelator_data));
    *phDecor = (void*)pData;
    int ch;

    /* default user parameters */
    pData->nChannels = 1;
    
    /* afSTFT stuff */
    pData->hSTFT = NULL;
    pData->STFTOutputFrameTF = malloc1d(MAX_NUM_OUTPUTS * sizeof(complexVector));
    for(ch=0; ch< MAX_NUM_OUTPUTS; ch++) {
        pData->STFTOutputFrameTF[ch].re = (float*)calloc1d(HYBRID_BANDS, sizeof(float));
        pData->STFTOutputFrameTF[ch].im = (float*)calloc1d(HYBRID_BANDS, sizeof(float));
    }
    pData->tempHopFrameTD = (float**)malloc2d( MAX(MAX_NUM_INPUTS, MAX_NUM_OUTPUTS), HOP_SIZE, sizeof(float));
    pData->STFTInputFrameTF = malloc1d(MAX_NUM_INPUTS * sizeof(complexVector));
    for(ch=0; ch< MAX_NUM_INPUTS; ch++) {
        pData->STFTInputFrameTF[ch].re = (float*)calloc1d(HYBRID_BANDS, sizeof(float));
        pData->STFTInputFrameTF[ch].im = (float*)calloc1d(HYBRID_BANDS, sizeof(float));
    }
    pData->InputFrameTF = (float_complex***)malloc3d(HYBRID_BANDS, MAX_NUM_INPUTS, TIME_SLOTS, sizeof(float_complex));
    pData->OutputFrameTF = (float_complex***)malloc3d(HYBRID_BANDS, MAX_NUM_OUTPUTS, TIME_SLOTS, sizeof(float_complex));

    /* codec data */
    pData->hDec = NULL;
    pData->new_nChannels = pData->nChannels;
    pData->progressBar0_1 = 0.0f;
    pData->progressBarText = malloc1d(PROGRESSBARTEXT_CHAR_LENGTH*sizeof(char));
    strcpy(pData->progressBarText,"");
    
    /* flags */
    pData->procStatus = PROC_STATUS_NOT_ONGOING;
    pData->codecStatus = CODEC_STATUS_NOT_INITIALISED;
}

void decorrelator_destroy
(
    void ** const phDecor
)
{
    decorrelator_data *pData = (decorrelator_data*)(*phDecor);
    int ch;
    
    if (pData != NULL) {
        /* not safe to free memory during intialisation/processing loop */
        while (pData->codecStatus == CODEC_STATUS_INITIALISING ||
               pData->procStatus == PROC_STATUS_ONGOING){
            SAF_SLEEP(10);
        }
        
        /* free afSTFT and buffers */ 
        if(pData->hSTFT!=NULL)
            afSTFTfree(pData->hSTFT);
        if(pData->STFTInputFrameTF!=NULL){
            for (ch = 0; ch< MAX_NUM_INPUTS; ch++) {
                free(pData->STFTInputFrameTF[ch].re);
                free(pData->STFTInputFrameTF[ch].im);
            }
        }
        for (ch = 0; ch< MAX_NUM_OUTPUTS; ch++) {
            free(pData->STFTOutputFrameTF[ch].re);
            free(pData->STFTOutputFrameTF[ch].im);
        }
        free(pData->STFTInputFrameTF);
        free(pData->STFTOutputFrameTF);
        free(pData->tempHopFrameTD);
        free(pData->InputFrameTF);
        free(pData->OutputFrameTF);
        free(pData->progressBarText);

        if(pData->hDec!=NULL)
            latticeDecorrelator_destroy(&(pData->hDec));
        
        free(pData);
        pData = NULL;
        *phDecor = NULL;
    }
}

void decorrelator_init
(
    void * const hDecor,
    int          sampleRate
)
{
    decorrelator_data *pData = (decorrelator_data*)(hDecor);
    int band;
    
    /* define frequency vector */
    pData->fs = sampleRate;
    for(band=0; band <HYBRID_BANDS; band++){
        if(sampleRate == 44100)
            pData->freqVector[band] =  (float)__afCenterFreq44100[band];
        else /* Assume 48kHz */
            pData->freqVector[band] =  (float)__afCenterFreq48e3[band];
    }
}

void decorrelator_initCodec
(
    void* const hDecor
)
{
    decorrelator_data *pData = (decorrelator_data*)(hDecor);
    int nChannels;
    
    if (pData->codecStatus != CODEC_STATUS_NOT_INITIALISED)
        return; /* re-init not required, or already happening */
    while (pData->procStatus == PROC_STATUS_ONGOING){
        /* re-init required, but we need to wait for the current processing loop to end */
        pData->codecStatus = CODEC_STATUS_INITIALISING; /* indicate that we want to init */
        SAF_SLEEP(10);
    }
    
    /* for progress bar */
    pData->codecStatus = CODEC_STATUS_INITIALISING;
    strcpy(pData->progressBarText,"Preparing Decorrelators");
    pData->progressBar0_1 = 0.0f;
    
    /* (Re)Initialise afSTFT */
    nChannels = pData->new_nChannels; 
    if(pData->hSTFT==NULL)
        afSTFTinit(&(pData->hSTFT), HOP_SIZE, nChannels, nChannels, 0, 1);
    else if(pData->nChannels != nChannels) {/* Or change the number of channels */
        afSTFTchannelChange(pData->hSTFT, nChannels, nChannels);
        afSTFTclearBuffers(pData->hSTFT);
    }
    pData->nChannels = nChannels;

    /* Init decorrelators */
    //int orders[4] = {20, 15, 6, 3};
    int orders[4] = {15, 6, 3, 3};
    float freqCutoffs[4] = {700.0f, 2.8e3f, 4.5e3f, 12e3f};
    int fixedDelays[5] = {9, 8, 3, 2, 3};
    latticeDecorrelator_create(&(pData->hDec), nChannels, orders, freqCutoffs, fixedDelays, 4, pData->freqVector, HYBRID_BANDS);

    /* done! */
    strcpy(pData->progressBarText,"Done!");
    pData->progressBar0_1 = 1.0f;
    pData->codecStatus = CODEC_STATUS_INITIALISED;
}

void decorrelator_process
(
    void  *  const hDecor,
    float ** const inputs,
    float ** const outputs,
    int            nInputs,
    int            nOutputs,
    int            nSamples
)
{
    decorrelator_data *pData = (decorrelator_data*)(hDecor);
    int t, ch, i, band;
    
    /* local copies of user parameters */
    int nCH;
    nCH = pData->nChannels;

    /* Process frame */
    if (nSamples == FRAME_SIZE && (pData->codecStatus == CODEC_STATUS_INITIALISED) ) {
        pData->procStatus = PROC_STATUS_ONGOING;

        /* Load time-domain data */
        for(i=0; i < MIN(nCH, nInputs); i++)
            utility_svvcopy(inputs[i], FRAME_SIZE, pData->InputFrameTD[i]);
        for(; i<nCH; i++)
            memset(pData->InputFrameTD[i], 0, FRAME_SIZE * sizeof(float)); /* fill remaining channels with zeros */

        /* Apply time-frequency transform (TFT) */
        for(t=0; t< TIME_SLOTS; t++) {
            for(ch = 0; ch < nCH; ch++)
                utility_svvcopy(&(pData->InputFrameTD[ch][t*HOP_SIZE]), HOP_SIZE, pData->tempHopFrameTD[ch]);
            afSTFTforward(pData->hSTFT, pData->tempHopFrameTD, pData->STFTInputFrameTF);
            for(band=0; band<HYBRID_BANDS; band++)
                for(ch=0; ch < nCH; ch++)
                    pData->InputFrameTF[band][ch][t] = cmplxf(pData->STFTInputFrameTF[ch].re[band], pData->STFTInputFrameTF[ch].im[band]);
        }

        /* Main processing: */
        latticeDecorrelator_apply(pData->hDec, pData->InputFrameTF, TIME_SLOTS, pData->OutputFrameTF);

        /* inverse-TFT */
        for(t = 0; t < TIME_SLOTS; t++) {
            for(band = 0; band < HYBRID_BANDS; band++) {
                for(ch = 0; ch < nCH; ch++) {
                    if(ch==0){
                        pData->STFTOutputFrameTF[ch].re[band] = crealf(pData->InputFrameTF[band][ch][t]);
                        pData->STFTOutputFrameTF[ch].im[band] = cimagf(pData->InputFrameTF[band][ch][t]);
                    }
                    else{
                        pData->STFTOutputFrameTF[ch].re[band] = crealf(pData->OutputFrameTF[band][ch][t]);
                        pData->STFTOutputFrameTF[ch].im[band] = cimagf(pData->OutputFrameTF[band][ch][t]);
                    }
                }
            }
            afSTFTinverse(pData->hSTFT, pData->STFTOutputFrameTF, pData->tempHopFrameTD);
            for (ch = 0; ch < MIN(nCH, nOutputs); ch++)
                utility_svvcopy(pData->tempHopFrameTD[ch], HOP_SIZE, &(outputs[ch][t* HOP_SIZE]));
            for (; ch < nOutputs; ch++)
                memset(&(outputs[ch][t* HOP_SIZE]), 0, HOP_SIZE*sizeof(float));
        }
    }
    else
        for (ch=0; ch < nOutputs; ch++)
            memset(outputs[ch],0, FRAME_SIZE*sizeof(float));

    pData->procStatus = PROC_STATUS_NOT_ONGOING;
}


/* Set Functions */

void decorrelator_refreshParams(void* const hDecor)
{
    decorrelator_setCodecStatus(hDecor, CODEC_STATUS_NOT_INITIALISED);
}

void decorrelator_setNumberOfChannels(void* const hDecor, int newValue )
{
    decorrelator_data *pData = (decorrelator_data*)(hDecor);

    if(newValue!=pData->new_nChannels){
        pData->new_nChannels = newValue;
            decorrelator_setCodecStatus(hDecor, CODEC_STATUS_NOT_INITIALISED);
    }

}
 
/* Get Functions */

int decorrelator_getFrameSize(void)
{
    return FRAME_SIZE;
}

CODEC_STATUS decorrelator_getCodecStatus(void* const hDecor)
{
    decorrelator_data *pData = (decorrelator_data*)(hDecor);
    return pData->codecStatus;
}

float decorrelator_getProgressBar0_1(void* const hDecor)
{
    decorrelator_data *pData = (decorrelator_data*)(hDecor);
    return pData->progressBar0_1;
}

void decorrelator_getProgressBarText(void* const hDecor, char* text)
{
    decorrelator_data *pData = (decorrelator_data*)(hDecor);
    memcpy(text, pData->progressBarText, PROGRESSBARTEXT_CHAR_LENGTH*sizeof(char));
}

int decorrelator_getNumberOfChannels(void* const hDecor)
{
    decorrelator_data *pData = (decorrelator_data*)(hDecor);
    return pData->new_nChannels;
}

int decorrelator_getDAWsamplerate(void* const hDecor)
{
    decorrelator_data *pData = (decorrelator_data*)(hDecor);
    return pData->fs;
}

int decorrelator_getProcessingDelay()
{
    return 12*HOP_SIZE;
}
