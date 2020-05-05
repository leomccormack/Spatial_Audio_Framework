//
//  main.c
//  test
//
//  Created by Leo McCormack on 27.4.2020. 
//

#include "unity.h"
#include "timer.h"
#include "saf.h"

tick_t start;

void setUp(void){
    // set stuff up here
}

void tearDown(void){
    // clean stuff up here
}

void test__saf_rfft(void){
    printf( "\n - test__saf_rfft\n");
    start = timer_current();
    int i, j, N;
    float* x_td, *test;
    float_complex* x_fd;
    void *hFFT;

    /* Config */
    const float acceptedTolerance = 0.000001f;
    const int fftSizesToTest[12] = {16,256,512,1024,2048,4096,8192,16384,32768,65536,1048576,67108864};

    /* Loop over the different FFT sizes */
    for (i=0; i<11; i++){
        N = fftSizesToTest[i];

        /* prep */
        x_td = malloc1d(N*sizeof(float));
        test = malloc1d(N*sizeof(float));
        x_fd = malloc1d((N/2+1)*sizeof(float_complex));
        rand_m1_1(x_td, N); /* populate with random numbers */
        saf_rfft_create(&hFFT, N);

        /* forward and backward transform */
        saf_rfft_forward(hFFT, x_td, x_fd);
        saf_rfft_backward(hFFT, x_fd, test);

        /* Check that, x_td==test */
        for(j=0; j<N; j++)
            TEST_ASSERT_FLOAT_WITHIN(acceptedTolerance, x_td[j], test[j]);

        /* clean-up */
        saf_rfft_destroy(&hFFT);
        free(x_fd);
        free(x_td);
        free(test);
    }

    printf( "Time elapsed: %lfs\n", (double)timer_elapsed( start ) );
}

void test__afSTFTMatrix(void){
    printf("\n - test__afSTFTMatrix\n");
    start = timer_current();
    int idx,frameIdx,c,t;
    int numChannels;
    float** inputTimeDomainData, **outputTimeDomainData, **tempFrame;
    float_complex*** frequencyDomainData;

    /* Config */
    const float acceptedTolerance = 0.01f; // Seems pretty high.. ?
    const int nTestFrames = 4;
    const int frameSize = 512;
    const int hopSize = 128;
    numChannels = 10;
    const int hybridMode = 1;

    /* prep */
    const int nTimeSlots = frameSize / hopSize;
    const int nBands = hopSize + (hybridMode ? 5 : 1);
    const int afSTFTdelay = hopSize * (hybridMode ? 12 : 9);
    const int lSig = nTestFrames*frameSize+afSTFTdelay;
    void* hSTFT;
    inputTimeDomainData = (float**) malloc2d(numChannels, lSig, sizeof(float));
    outputTimeDomainData = (float**) malloc2d(numChannels, lSig, sizeof(float));
    tempFrame = (float**) malloc2d(numChannels, frameSize, sizeof(float));
    frequencyDomainData = (float_complex***) malloc3d(nBands, numChannels, nTimeSlots, sizeof(float_complex));

    /* Initialise afSTFT and input data */
    afSTFTMatrixInit(&hSTFT, hopSize, numChannels, numChannels, 0, hybridMode, frameSize);
    rand_m1_1(ADR2D(inputTimeDomainData), numChannels*lSig); /* populate with random numbers */

    /* Pass input data through afSTFT */
    idx = 0;
    frameIdx = 0;
    while(idx<lSig){
        for(c=0; c<numChannels; c++)
            memcpy(tempFrame[c], &(inputTimeDomainData[c][frameIdx*frameSize]), frameSize*sizeof(float));

        /* forward and inverse */
        afSTFTMatrixForward(hSTFT, tempFrame, frequencyDomainData);
        afSTFTMatrixInverse(hSTFT, frequencyDomainData, tempFrame);

        for(c=0; c<numChannels; c++)
            memcpy(&(outputTimeDomainData[c][frameIdx*frameSize]), tempFrame[c], frameSize*sizeof(float)); 
        idx+=frameSize;
        frameIdx++;
    }

    /* Compensate for afSTFT delay, and check that input==output, given some numerical precision */
    for(c=0; c<numChannels; c++){
        memcpy(outputTimeDomainData[c], &(outputTimeDomainData[c][afSTFTdelay]), (lSig-afSTFTdelay) *sizeof(float));
        for(t=0; t<(lSig-afSTFTdelay); t++)
            TEST_ASSERT_FLOAT_WITHIN(acceptedTolerance, inputTimeDomainData[c][t], outputTimeDomainData[c][t]);
    }

    /* tidy-up */
    afSTFTMatrixFree(hSTFT);
    free(inputTimeDomainData);
    free(outputTimeDomainData);
    free(frequencyDomainData);
    printf( "Time elapsed: %lfs\n", (double)timer_elapsed( start ) );
}

void test__smb_pitchShifter(void){
    printf( "\n - test__smb_pitchShifter\n");
    start = timer_current();
    float* inputData, *outputData;
    void* hPS;
    float frequency;
    int i, nSamples,smbLatency;

    /* Config */
    const int signalLength_seconds = 10;
    const int sampleRate = 48e3;
    const int FFTsize = 8192;
    const int osfactor = 16;

    /* prep */
    smb_pitchShift_create(&hPS, 1, FFTsize, osfactor, sampleRate);
    nSamples = sampleRate*signalLength_seconds;
    inputData = malloc1d(nSamples*sizeof(float));
    outputData = calloc1d(nSamples,sizeof(float));
    frequency = (float)sampleRate/8.0f;
    for(i=0; i<nSamples; i++) /* sine tone at quarter Nyquist: */
        inputData[i] = sinf(2.0f * M_PI * (float)i * frequency/(float)sampleRate);
    smbLatency = FFTsize - (FFTsize/osfactor);

    /* Pitch shift down one octave */
    smb_pitchShift_apply(hPS, 0.5, nSamples, inputData, outputData);

    // TODO: Take FFT, the bin with the highest energy should correspond to 1/8 Nyquist

    /* clean-up */
    smb_pitchShift_destroy(&hPS);
    free(inputData);
    free(outputData);
    printf( "Time elapsed: %lfs\n", (double)timer_elapsed( start ) );
}

int main(void){
printf("*****************************************************************\n"
       "******** Spatial_Audio_Framework Unit Testing Program ***********\n"
       "*****************************************************************\n");

    /* initialise */
    timer_lib_initialize();
    UNITY_BEGIN();

    /* run each unit test */
    RUN_TEST(test__saf_rfft);
    RUN_TEST(test__afSTFTMatrix);
    RUN_TEST(test__smb_pitchShifter);

    /* close */
    timer_lib_shutdown();
    return UNITY_END();
}

