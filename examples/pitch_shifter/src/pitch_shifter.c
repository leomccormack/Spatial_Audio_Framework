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
 * @file pitch_shifter.c
 * @brief A very basic multichannel pitch shifter
 *
 * @author Leo McCormack
 * @date 05.05.2020
 */

#include "pitch_shifter.h"
#include "pitch_shifter_internal.h"

void pitch_shifter_create
(
    void ** const phPS
)
{
    pitch_shifter_data* pData = (pitch_shifter_data*)malloc1d(sizeof(pitch_shifter_data));
    *phPS = (void*)pData;
 
    /* Default user parameters */
    pData->new_nChannels = pData->nChannels = 1;
    pData->pitchShift_factor = 0.5f; 
    pData->osamp_option = PITCH_SHIFTER_OSAMP_4;
    pData->fftsize_option = PITCH_SHIFTER_FFTSIZE_4096;

    /* internals */
    pData->hSmb = NULL;
    pData->progressBar0_1 = 0.0f;
    pData->progressBarText = malloc1d(PITCH_SHIFTER_PROGRESSBARTEXT_CHAR_LENGTH*sizeof(char));
    strcpy(pData->progressBarText,"");
    pData->sampleRate = 0.0f;

    /* flags */
    pData->procStatus = PROC_STATUS_NOT_ONGOING;
    pData->codecStatus = CODEC_STATUS_NOT_INITIALISED;
}

void pitch_shifter_destroy
(
    void ** const phPS
)
{
    pitch_shifter_data *pData = (pitch_shifter_data*)(*phPS);

    if (pData != NULL) {
        /* not safe to free memory during intialisation/processing loop */
        while (pData->codecStatus == CODEC_STATUS_INITIALISING ||
               pData->procStatus == PROC_STATUS_ONGOING){
            SAF_SLEEP(10);
        }

        if (pData->hSmb != NULL)
            smb_pitchShift_destroy(&(pData->hSmb));
        free(pData);
        pData = NULL;
    }
}

void pitch_shifter_init
(
    void * const hPS,
    int          sampleRate
)
{
    pitch_shifter_data *pData = (pitch_shifter_data*)(hPS);

    if(pData->sampleRate != sampleRate){
        pData->sampleRate = sampleRate;
        pitch_shifter_setCodecStatus(hPS, CODEC_STATUS_NOT_INITIALISED);
    } 
}

void pitch_shifter_initCodec
(
    void* const hPS
)
{
    pitch_shifter_data *pData = (pitch_shifter_data*)(hPS);
    int nChannels, fftSize, osamp;

    if (pData->codecStatus != CODEC_STATUS_NOT_INITIALISED)
        return; /* re-init not required, or already happening */
    while (pData->procStatus == PROC_STATUS_ONGOING){
        /* re-init required, but we need to wait for the current processing loop to end */
        pData->codecStatus = CODEC_STATUS_INITIALISING; /* indicate that we want to init */
        SAF_SLEEP(10);
    }

    /* for progress bar */
    pData->codecStatus = CODEC_STATUS_INITIALISING;
    strcpy(pData->progressBarText,"Initialising pitch shifter");
    pData->progressBar0_1 = 0.0f;

    nChannels = pData->new_nChannels;

    /* destroy current handle*/
    if (pData->hSmb != NULL)
        smb_pitchShift_destroy(&(pData->hSmb));

    /* Config */
    switch(pData->osamp_option){
        case PITCH_SHIFTER_OSAMP_2:  osamp = 2; break;
        case PITCH_SHIFTER_OSAMP_4:  osamp = 4; break;
        case PITCH_SHIFTER_OSAMP_8:  osamp = 8; break;
        case PITCH_SHIFTER_OSAMP_16: osamp = 16; break;
        case PITCH_SHIFTER_OSAMP_32: osamp = 32; break;
    }
    switch(pData->fftsize_option){
        case PITCH_SHIFTER_FFTSIZE_512:   fftSize = 512; break;
        case PITCH_SHIFTER_FFTSIZE_1024:  fftSize = 1024; break;
        case PITCH_SHIFTER_FFTSIZE_2048:  fftSize = 2048; break;
        case PITCH_SHIFTER_FFTSIZE_4096:  fftSize = 4096; break;
        case PITCH_SHIFTER_FFTSIZE_8192:  fftSize = 8192; break;
        case PITCH_SHIFTER_FFTSIZE_16384: fftSize = 16384; break;
    }

    /* Create new handle*/
    smb_pitchShift_create(&(pData->hSmb), nChannels, fftSize, osamp, pData->sampleRate);
    pData->nChannels = nChannels;

    /* done! */
    strcpy(pData->progressBarText,"Done!");
    pData->progressBar0_1 = 1.0f;
    pData->codecStatus = CODEC_STATUS_INITIALISED;
}

void pitch_shifter_process
(
    void  *  const hPS,
    float ** const inputs,
    float ** const outputs,
    int            nInputs,
    int            nOutputs,
    int            nSamples
)
{
    pitch_shifter_data *pData = (pitch_shifter_data*)(hPS);
    int ch, nChannels;
    nChannels = pData->nChannels;

    if ( (FRAME_SIZE==nSamples) && (pData->codecStatus == CODEC_STATUS_INITIALISED) ){
        pData->procStatus = PROC_STATUS_ONGOING;

        for(ch=0; ch<MIN(nChannels,nInputs); ch++)
            memcpy(pData->inputFrame[ch], inputs[ch], FRAME_SIZE*sizeof(float));
        for(; ch<nChannels; ch++)
            memset(pData->inputFrame[ch], 0, FRAME_SIZE*sizeof(float));

        smb_pitchShift_apply(pData->hSmb, pData->pitchShift_factor, nSamples, (float*)pData->inputFrame, (float*)pData->outputFrame);
        for(ch=0; ch<MIN(nChannels,nOutputs); ch++)
            memcpy(outputs[ch], pData->outputFrame[ch], FRAME_SIZE*sizeof(float));
    }

    pData->procStatus = PROC_STATUS_NOT_ONGOING;
}

/* sets */

void pitch_shifter_refreshParams(void* const hPS)
{
    pitch_shifter_setCodecStatus(hPS, CODEC_STATUS_NOT_INITIALISED);
}

void pitch_shifter_setPitchShiftFactor(void* const hPS, float newValue)
{
    pitch_shifter_data *pData = (pitch_shifter_data*)(hPS);
    pData->pitchShift_factor = newValue;
}

void pitch_shifter_setNumChannels(void* const hPS, int newValue)
{
    pitch_shifter_data *pData = (pitch_shifter_data*)(hPS);
    pData->new_nChannels = newValue;
    pitch_shifter_setCodecStatus(hPS, CODEC_STATUS_NOT_INITIALISED);
}

void pitch_shifter_setFFTSizeOption
(
    void* const hPS,
    PITCH_SHIFTER_FFTSIZE_OPTIONS newOption
)
{
    pitch_shifter_data *pData = (pitch_shifter_data*)(hPS);
    pData->fftsize_option = newOption;
    pitch_shifter_setCodecStatus(hPS, CODEC_STATUS_NOT_INITIALISED);
}

void pitch_shifter_setOSampOption
(
    void* const hPS,
    PITCH_SHIFTER_OSAMP_OPTIONS newOption
)
{
    pitch_shifter_data *pData = (pitch_shifter_data*)(hPS);
    pData->osamp_option = newOption;
    pitch_shifter_setCodecStatus(hPS, CODEC_STATUS_NOT_INITIALISED);
}


/* gets */

PITCH_SHIFTER_CODEC_STATUS pitch_shifter_getCodecStatus(void* const hBin)
{
    pitch_shifter_data *pData = (pitch_shifter_data*)(hBin);
    return pData->codecStatus;
}

float pitch_shifter_getProgressBar0_1(void* const hBin)
{
    pitch_shifter_data *pData = (pitch_shifter_data*)(hBin);
    return pData->progressBar0_1;
}

void pitch_shifter_getProgressBarText(void* const hBin, char* text)
{
    pitch_shifter_data *pData = (pitch_shifter_data*)(hBin);
    memcpy(text, pData->progressBarText, PITCH_SHIFTER_PROGRESSBARTEXT_CHAR_LENGTH*sizeof(char));
}

float pitch_shifter_getPitchShiftFactor(void* const hPS)
{
    pitch_shifter_data *pData = (pitch_shifter_data*)(hPS);
    return pData->pitchShift_factor;
}

PITCH_SHIFTER_FFTSIZE_OPTIONS pitch_shifter_getFFTSizeOption(void* const hPS)
{
    pitch_shifter_data *pData = (pitch_shifter_data*)(hPS);
    return pData->fftsize_option;
}

PITCH_SHIFTER_OSAMP_OPTIONS pitch_shifter_getOSampOption(void* const hPS)
{
    pitch_shifter_data *pData = (pitch_shifter_data*)(hPS);
    return pData->osamp_option;
}

int pitch_shifter_getNCHrequired(void* const hPS)
{
    pitch_shifter_data *pData = (pitch_shifter_data*)(hPS);
    return pData->new_nChannels;
}



