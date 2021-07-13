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
 * @license ISC
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
    pData->decorAmount = 1.0f;
    pData->compensateLevel = 0;
    
    /* afSTFT stuff */
    pData->hSTFT = NULL;
    pData->InputFrameTD = (float**)malloc2d(MAX_NUM_CHANNELS, DECORRELATOR_FRAME_SIZE, sizeof(float));
    pData->OutputFrameTD = (float**)malloc2d(MAX_NUM_CHANNELS, DECORRELATOR_FRAME_SIZE, sizeof(float));
    pData->InputFrameTF = (float_complex***)malloc3d(HYBRID_BANDS, MAX_NUM_CHANNELS, TIME_SLOTS, sizeof(float_complex));
    pData->OutputFrameTF = (float_complex***)malloc3d(HYBRID_BANDS, MAX_NUM_CHANNELS, TIME_SLOTS, sizeof(float_complex));
    pData->transientFrameTF = (float_complex***)malloc3d(HYBRID_BANDS, MAX_NUM_CHANNELS, TIME_SLOTS, sizeof(float_complex));

    /* codec data */
    pData->hDecor = NULL;
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
        free(pData->transientFrameTF);
        free(pData->progressBarText);

        transientDucker_destroy(&(pData->hDucker));
        latticeDecorrelator_destroy(&(pData->hDecor));

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
    if(pData->hDecor!=NULL)
        latticeDecorrelator_reset(pData->hDecor);
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

    /* Init decorrelator  */
    const int orders[4] = {20, 15, 6, 3}; /* 20th order up to 700Hz, 15th->2.4kHz, 6th->4kHz, 3rd->12kHz, NONE(only delays)->Nyquist */
    const float freqCutoffs[4] = {600.0f, 2.4e3f, 4.0e3f, 12e3f};
    //const float freqCutoffs[4] = {900.0f, 6.8e3f, 12e3f, 16e3f};
    const int maxDelay = 8;
    latticeDecorrelator_destroy(&(pData->hDecor));
    latticeDecorrelator_create(&(pData->hDecor), pData->fs, HOP_SIZE, pData->freqVector, HYBRID_BANDS, pData->nChannels, (int*)orders, (float*)freqCutoffs, 4, maxDelay, 0, 0.75f);

    /* done! */
    strcpy(pData->progressBarText,"Done!");
    pData->progressBar0_1 = 1.0f;
    pData->codecStatus = CODEC_STATUS_INITIALISED;
}

void decorrelator_process
(
    void        *  const hDecor,
    const float *const * inputs,
    float       ** const outputs,
    int                  nInputs,
    int                  nOutputs,
    int                  nSamples
)
{
    decorrelator_data *pData = (decorrelator_data*)(hDecor);
    int ch, i, band, enableTransientDucker, compensateLevel;
    float decorAmount;
    
    /* local copies of user parameters */
    int nCH;
    nCH = pData->nChannels;
    decorAmount = pData->decorAmount;
    enableTransientDucker = pData->enableTransientDucker;
    compensateLevel = pData->compensateLevel;

    /* Process frame */
    if (nSamples == DECORRELATOR_FRAME_SIZE && (pData->codecStatus == CODEC_STATUS_INITIALISED) ) {
        pData->procStatus = PROC_STATUS_ONGOING;

        /* Load time-domain data */
        for(i=0; i < SAF_MIN(nCH, nInputs); i++)
            utility_svvcopy(inputs[i], DECORRELATOR_FRAME_SIZE, pData->InputFrameTD[i]);
        for(; i<nCH; i++)
            memset(pData->InputFrameTD[i], 0, DECORRELATOR_FRAME_SIZE * sizeof(float)); /* fill remaining channels with zeros */

        /* Apply time-frequency transform (TFT) */
        afSTFT_forward_knownDimensions(pData->hSTFT, pData->InputFrameTD, DECORRELATOR_FRAME_SIZE, MAX_NUM_CHANNELS, TIME_SLOTS, pData->InputFrameTF);

        /* Apply decorrelation */
        if(enableTransientDucker){
            /* remove transients */
            transientDucker_apply(pData->hDucker, pData->InputFrameTF, TIME_SLOTS, 0.95f, 0.995f, pData->OutputFrameTF, pData->transientFrameTF);
            /* decorrelate only the residual */
            latticeDecorrelator_apply(pData->hDecor,  pData->OutputFrameTF, TIME_SLOTS, pData->OutputFrameTF);
        }
        else
            latticeDecorrelator_apply(pData->hDecor,  pData->InputFrameTF, TIME_SLOTS, pData->OutputFrameTF);

        /* Optionally compensate for the level (as they channels wll no longer sum coherently) */
        if(compensateLevel){
            for(band=0; band<HYBRID_BANDS; band++)
                cblas_sscal(/*re+im*/2*nCH*TIME_SLOTS, 0.75f*(float)nCH/(sqrtf((float)nCH)), (float*)FLATTEN2D(pData->OutputFrameTF[band]), 1);
        }

        /* re-introduce the transient part */
        if(enableTransientDucker){
            //scalec =  cmplxf(1.0f, 0.0f);//!compensateLevel ? cmplxf(1.25f*(sqrtf((float)nCH)/(float)nCH), 0.0f) : cmplxf(1.0f, 0.0f);
            for(band=0; band<HYBRID_BANDS; band++)
                cblas_saxpy(/*re+im*/2*nCH*TIME_SLOTS, 1.0f, (float*)FLATTEN2D(pData->transientFrameTF[band]), 1, (float*)FLATTEN2D(pData->OutputFrameTF[band]), 1);
        }

        /* Mix  thedecorrelated audio with the input non-decorrelated audio */ 
        for(band=0; band<HYBRID_BANDS; band++){
            cblas_sscal(/*re+im*/2*nCH*TIME_SLOTS, decorAmount, (float*)FLATTEN2D(pData->OutputFrameTF[band]), 1);
            cblas_saxpy(/*re+im*/2*nCH*TIME_SLOTS, 1.0f-decorAmount, (float*)FLATTEN2D(pData->InputFrameTF[band]), 1, (float*)FLATTEN2D(pData->OutputFrameTF[band]), 1);
        }

        /* inverse-TFT */
        afSTFT_backward_knownDimensions(pData->hSTFT, pData->OutputFrameTF, DECORRELATOR_FRAME_SIZE, MAX_NUM_CHANNELS, TIME_SLOTS, pData->OutputFrameTD);

        /* Copy to output buffer */
        for (ch = 0; ch < SAF_MIN(nCH, nOutputs); ch++)
            utility_svvcopy(pData->OutputFrameTD[ch], DECORRELATOR_FRAME_SIZE, outputs[ch]);
        for (; ch < nOutputs; ch++)
            memset(outputs[ch], 0, DECORRELATOR_FRAME_SIZE*sizeof(float));
    }
    else
        for (ch=0; ch < nOutputs; ch++)
            memset(outputs[ch],0, DECORRELATOR_FRAME_SIZE*sizeof(float));

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

void decorrelator_setDecorrelationAmount(void* const hDecor, float newValue)
{
    decorrelator_data *pData = (decorrelator_data*)(hDecor);
    pData->decorAmount = SAF_CLAMP(newValue, 0.0f, 1.0f);
}

void decorrelator_setLevelCompensationFlag(void* const hDecor, int newValue)
{
    decorrelator_data *pData = (decorrelator_data*)(hDecor);
    saf_assert(newValue==0 || newValue==1, "newValue is a bool");
    pData->compensateLevel = newValue;
}

void decorrelator_setTransientBypassFlag(void* const hDecor, int newValue)
{
    decorrelator_data *pData = (decorrelator_data*)(hDecor);
    saf_assert(newValue==0 || newValue==1, "newValue is a bool");
    pData->enableTransientDucker = newValue;
}
 
/* Get Functions */

int decorrelator_getFrameSize(void)
{
    return DECORRELATOR_FRAME_SIZE;
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

float decorrelator_getDecorrelationAmount(void* const hDecor)
{
    decorrelator_data *pData = (decorrelator_data*)(hDecor);
    return pData->decorAmount;
}

int decorrelator_getLevelCompensationFlag(void* const hDecor)
{
    decorrelator_data *pData = (decorrelator_data*)(hDecor);
    return pData->compensateLevel;
}

int decorrelator_getTransientBypassFlag(void* const hDecor)
{
    decorrelator_data *pData = (decorrelator_data*)(hDecor);
    return pData->enableTransientDucker;
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
