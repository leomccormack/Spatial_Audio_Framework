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

    /* default user parameters */
    pData->nChannels = 1;
    pData->enableTransientDucker = 0;
    
    /* afSTFT stuff */
    pData->hSTFT = NULL;
    pData->InputFrameTD = (float**)malloc2d(MAX_NUM_INPUTS, FRAME_SIZE, sizeof(float));
    pData->OutputFrameTD = (float**)malloc2d(MAX_NUM_OUTPUTS, FRAME_SIZE, sizeof(float));
    pData->InputFrameTF = (float_complex***)malloc3d(HYBRID_BANDS, MAX_NUM_INPUTS, TIME_SLOTS, sizeof(float_complex));
    pData->OutputFrameTF = (float_complex***)malloc3d(HYBRID_BANDS, MAX_NUM_OUTPUTS, TIME_SLOTS, sizeof(float_complex));

    /* codec data */
    pData->hDec = NULL;
    pData->hDec2 = NULL;
    pData->hDucker = NULL;
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
    
    if (pData != NULL) {
        /* not safe to free memory during intialisation/processing loop */
        while (pData->codecStatus == CODEC_STATUS_INITIALISING ||
               pData->procStatus == PROC_STATUS_ONGOING){
            SAF_SLEEP(10);
        }
        
        /* free afSTFT and buffers */ 
        if(pData->hSTFT!=NULL)
            afSTFT_destroy(&(pData->hSTFT));
        free(pData->InputFrameTD);
        free(pData->OutputFrameTD);
        free(pData->InputFrameTF);
        free(pData->OutputFrameTF);
        free(pData->progressBarText);

        transientDucker_destroy(&(pData->hDucker));
        latticeDecorrelator_destroy(&(pData->hDec));
        latticeDecorrelator_destroy(&(pData->hDec2));


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
    
    /* define frequency vector */
    pData->fs = sampleRate;
    afSTFT_getCentreFreqs(pData->hSTFT, (float)sampleRate, HYBRID_BANDS, pData->freqVector);
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
        afSTFT_create(&(pData->hSTFT), nChannels, nChannels, HOP_SIZE, 0, 1, AFSTFT_BANDS_CH_TIME);
    else if(pData->nChannels != nChannels) {/* Or change the number of channels */
        afSTFT_channelChange(pData->hSTFT, nChannels, nChannels);
        afSTFT_clearBuffers(pData->hSTFT);
    }
    pData->nChannels = nChannels;

    /* Init transient ducker */
    transientDucker_destroy(&(pData->hDucker));
    transientDucker_create(&(pData->hDucker), nChannels, HYBRID_BANDS);

    /* Init decorrelator 1 (best to apply these decorrelators in 2 stages) */
    int orders[5] = {6, 6, 6, 3, 2};
    float freqCutoffs[5] = {700.0f, 2.4e3f, 4e3f, 12e3f, 20e3f};
    int fixedDelays[6] = {8, 8, 7, 2, 1, 2};
    latticeDecorrelator_destroy(&(pData->hDec));
    latticeDecorrelator_create(&(pData->hDec), nChannels, orders, freqCutoffs, fixedDelays, 5, pData->freqVector, 0, HYBRID_BANDS);

    /* Init decorrelator 2 */
    float freqCutoffs2[3] = {700.0f, 2.4e3f, 4e3f};
    int orders2[3] = {3, 3, 2};
    int fixedDelays2[4] = {2, 2, 1, 0};
    latticeDecorrelator_destroy(&(pData->hDec2));
    latticeDecorrelator_create(&(pData->hDec2), nChannels, orders2, freqCutoffs2, fixedDelays2, 3, pData->freqVector, nChannels, HYBRID_BANDS);

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
    int ch, i;
    
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
        afSTFT_forward(pData->hSTFT, pData->InputFrameTD, FRAME_SIZE, pData->InputFrameTF);

        /* Main processing: */
        if(pData->enableTransientDucker){
            transientDucker_apply(pData->hDucker, pData->InputFrameTF, TIME_SLOTS, pData->OutputFrameTF);
            latticeDecorrelator_apply(pData->hDec,  pData->OutputFrameTF,  TIME_SLOTS, pData->OutputFrameTF);
        }
        else
            latticeDecorrelator_apply(pData->hDec,  pData->InputFrameTF,  TIME_SLOTS, pData->OutputFrameTF);
        latticeDecorrelator_apply(pData->hDec2, pData->OutputFrameTF, TIME_SLOTS, pData->OutputFrameTF);

        /* inverse-TFT */
        afSTFT_backward(pData->hSTFT, pData->OutputFrameTF, FRAME_SIZE, pData->OutputFrameTD);

        /* Copy to output buffer */
        for (ch = 0; ch < MIN(nCH, nOutputs); ch++)
            utility_svvcopy(pData->OutputFrameTD[ch], FRAME_SIZE, outputs[ch]);
        for (; ch < nOutputs; ch++)
            memset(outputs[ch], 0, FRAME_SIZE*sizeof(float));
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
