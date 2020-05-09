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
 * @file main.c
 * @brief testing program for the Spatial_Audio_Framework
 *
 * @author Leo McCormack
 * @date 27.04.2020
 */

#include "unity.h"
#include "timer.h"
#include "saf.h"

/* Prototypes for available unit tests */
void test__ims_shoebox(void);
void test__saf_rfft(void);
#ifdef AFSTFT_USE_FLOAT_COMPLEX
void test__afSTFTMatrix(void);
#endif
void test__afSTFT(void);
void test__smb_pitchShifter(void);
void test__sortf(void);
void test__realloc2d_r(void);
void test__getSHreal_recur(void);

/* ========================================================================== */
/*                                 Test Config                                */
/* ========================================================================== */

static tick_t start, start_test;
void setUp(void){ start_test = timer_current(); }
void tearDown(void){ }
static void timerResult(void) {
    printf( "    (Time elapsed: %lfs) \n", (double)timer_elapsed(start_test) );
}

#undef RUN_TEST
#define RUN_TEST(testfunc)  UNITY_NEW_TEST(#testfunc) \
    if (TEST_PROTECT()) {  setUp();  testfunc();  } \
    if (TEST_PROTECT() && (!TEST_IS_IGNORED))  {tearDown(); } \
    UnityConcludeTest(); timerResult();

int main(void){
printf("*****************************************************************\n"
       "********* Spatial_Audio_Framework Unit Testing Program **********\n"
       "*****************************************************************\n\n");

    /* initialise */
    timer_lib_initialize();
    start = timer_current();
    UNITY_BEGIN();

    /* run each unit test */
    RUN_TEST(test__ims_shoebox);
    RUN_TEST(test__saf_rfft);
#ifdef AFSTFT_USE_FLOAT_COMPLEX
    RUN_TEST(test__afSTFTMatrix);
#endif
    RUN_TEST(test__afSTFT);
    RUN_TEST(test__smb_pitchShifter);
    RUN_TEST(test__sortf);
    RUN_TEST(test__realloc2d_r);
    RUN_TEST(test__getSHreal_recur);

    /* close */
    timer_lib_shutdown();
    printf( "\nTotal time elapsed: %lfs", (double)timer_elapsed( start ) );
    return UNITY_END();
}


/* ========================================================================== */
/*                                 Unit Tests                                 */
/* ========================================================================== */

void test__ims_shoebox(void){
    void* hIms;
    float maxTime_s;
    float mov_src_pos[3], mov_rec_pos[3];
    long sourceID_1, sourceID_2, sourceID_3, sourceID_4, sourceID_5, receiverID;
    int i;

    /* Config */
    const int sh_order = 3;
    const int nBands = 7;
    const float abs_wall[nBands][6] =  /* Absorption Coefficients per Octave band, and per wall */
      { {0.180791250f, 0.207307300f, 0.134990800f, 0.229002250f, 0.212128400f, 0.241055000f},
        {0.225971250f, 0.259113700f, 0.168725200f, 0.286230250f, 0.265139600f, 0.301295000f},
        {0.258251250f, 0.296128100f, 0.192827600f, 0.327118250f, 0.303014800f, 0.344335000f},
        {0.301331250f, 0.345526500f, 0.224994001f, 0.381686250f, 0.353562000f, 0.401775000f},
        {0.361571250f, 0.414601700f, 0.269973200f, 0.457990250f, 0.424243600f, 0.482095000f},
        {0.451931250f, 0.518214500f, 0.337442000f, 0.572446250f, 0.530266000f, 0.602575000f},
        {0.602591250f, 0.690971300f, 0.449934800f, 0.763282250f, 0.707040400f, 0.803455000f} };
    const float src_pos[3] = {5.1f, 6.0f, 1.1f};
    const float src2_pos[3] = {2.1f, 1.0f, 1.3f};
    const float src3_pos[3] = {4.4f, 3.0f, 1.4f};
    const float src4_pos[3] = {6.4f, 4.0f, 1.3f};
    const float src5_pos[3] = {8.5f, 5.0f, 1.8f};
    const float rec_pos[3] = {8.8f, 5.5f, 0.9f};

    /* Set-up the shoebox room simulator, with two sources and one receiver */
    ims_shoebox_create(&hIms, 10, 7, 3, (float*)abs_wall, 125.0f, nBands, 343.0f, 48e3f);
    sourceID_1 = ims_shoebox_addSource(hIms, (float*)src_pos);
    sourceID_2 = ims_shoebox_addSource(hIms, (float*)src2_pos);
    receiverID = ims_shoebox_addReceiver(hIms, (float*)rec_pos);

    /* Moving source No.1 and the receiver */
    maxTime_s = 0.08f; /* 50ms */
    memcpy(mov_src_pos, src_pos, 3*sizeof(float));
    memcpy(mov_rec_pos, rec_pos, 3*sizeof(float));
    for(i=0; i<500; i++){
        mov_src_pos[1] = 2.0f * 500.0f/(float)i;
        mov_rec_pos[0] = 3.0f * 500.0f/(float)i;
        ims_shoebox_updateSource(hIms, sourceID_1, mov_src_pos);
        ims_shoebox_updateReceiver(hIms, receiverID, mov_rec_pos);
        ims_shoebox_renderEchogramSH(hIms, maxTime_s, sh_order);
    }

    /* Remove source No.1 */
    ims_shoebox_removeSource(hIms, 0);

    /* Add 3 more sources, then remove 2, and add one back again
     * (Just messing around, trying to trip up an IMS internal assertion) */
    sourceID_3 = ims_shoebox_addSource(hIms, (float*)src3_pos);
    sourceID_4 = ims_shoebox_addSource(hIms, (float*)src4_pos);
    sourceID_5 = ims_shoebox_addSource(hIms, (float*)src5_pos);
    ims_shoebox_removeSource(hIms, sourceID_3);
    ims_shoebox_removeSource(hIms, sourceID_4);
    sourceID_4 = ims_shoebox_addSource(hIms, (float*)src4_pos);

    /* Continue rendering */
    for(i=0; i<500; i++){
        mov_src_pos[1] = 2.0f * 500.0f/(float)i;
        mov_rec_pos[0] = 3.0f * 500.0f/(float)i;
        ims_shoebox_updateSource(hIms, sourceID_4, mov_src_pos);
        ims_shoebox_updateReceiver(hIms, receiverID, mov_rec_pos);
        ims_shoebox_renderEchogramSH(hIms, maxTime_s, sh_order);
    }

    /* clean-up */
    ims_shoebox_destroy(&hIms);
}

void test__saf_rfft(void){
    int i, j, N;
    float* x_td, *test;
    float_complex* x_fd;
    void *hFFT;

    /* Config */
    const float acceptedTolerance = 0.000001f;
    const int fftSizesToTest[12] =
        {16,256,512,1024,2048,4096,8192,16384,32768,65536,1048576,33554432};

    /* Loop over the different FFT sizes */
    for (i=0; i<12; i++){
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
}

#ifdef AFSTFT_USE_FLOAT_COMPLEX
void test__afSTFTMatrix(void){
    int idx,frameIdx,c,t;
    int numChannels;
    float** inputTimeDomainData, **outputTimeDomainData, **tempFrame;
    float_complex*** frequencyDomainData;

    /* Config */
    const float acceptedTolerance = 0.01f; // Seems pretty high.. ?
    const int nTestFrames = 1000;
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
}
#endif

void test__afSTFT(void){
    int idx,hopIdx,c,t;
    int numChannels;
    float** inputTimeDomainData, **outputTimeDomainData, **tempHop;
#ifdef AFSTFT_USE_FLOAT_COMPLEX
    float_complex** frequencyDomainData;
#else
    complexVector* frequencyDomainData;
#endif

    /* Config */
    const float acceptedTolerance = 0.01f; // Seems pretty high.. ? (same for AFSTFT_USE_FLOAT_COMPLEX defined/undefined)
    const int nTestHops = 5000;
    const int hopSize = 128;
    numChannels = 10;
    const int hybridMode = 1;

    /* prep */
    const int nBands = hopSize + (hybridMode ? 5 : 1);
    const int afSTFTdelay = hopSize * (hybridMode ? 12 : 9);
    const int lSig = nTestHops*hopSize+afSTFTdelay;
    void* hSTFT;
    inputTimeDomainData = (float**) malloc2d(numChannels, lSig, sizeof(float));
    outputTimeDomainData = (float**) malloc2d(numChannels, lSig, sizeof(float));
    tempHop = (float**) malloc2d(numChannels, hopSize, sizeof(float));
#ifdef AFSTFT_USE_FLOAT_COMPLEX
    frequencyDomainData = (float_complex**) malloc2d(numChannels, nBands, sizeof(float_complex));
#else
    frequencyDomainData = malloc1d(numChannels * sizeof(complexVector));
    for(c=0; c<numChannels; c++){
        frequencyDomainData[c].re = malloc1d(nBands*sizeof(float));
        frequencyDomainData[c].im = malloc1d(nBands*sizeof(float));
    }
#endif

    /* Initialise afSTFT and input data */
    afSTFTinit(&hSTFT, hopSize, numChannels, numChannels, 0, hybridMode);
    rand_m1_1(ADR2D(inputTimeDomainData), numChannels*lSig); /* populate with random numbers */

    /* Pass input data through afSTFT */
    idx = 0;
    hopIdx = 0;
    while(idx<lSig){
        for(c=0; c<numChannels; c++)
            memcpy(tempHop[c], &(inputTimeDomainData[c][hopIdx*hopSize]), hopSize*sizeof(float));

        /* forward and inverse */
        afSTFTforward(hSTFT, tempHop, frequencyDomainData);
        afSTFTinverse(hSTFT, frequencyDomainData, tempHop);

        for(c=0; c<numChannels; c++)
            memcpy(&(outputTimeDomainData[c][hopIdx*hopSize]), tempHop[c], hopSize*sizeof(float));
        idx+=hopSize;
        hopIdx++;
    }

    /* Compensate for afSTFT delay, and check that input==output, given some numerical precision */
    for(c=0; c<numChannels; c++){
        memcpy(outputTimeDomainData[c], &(outputTimeDomainData[c][afSTFTdelay]), (lSig-afSTFTdelay) *sizeof(float));
        for(t=0; t<(lSig-afSTFTdelay); t++)
            TEST_ASSERT_FLOAT_WITHIN(acceptedTolerance, inputTimeDomainData[c][t], outputTimeDomainData[c][t]);
    }

    /* tidy-up */
    afSTFTfree(hSTFT);
    free(inputTimeDomainData);
    free(outputTimeDomainData);
#ifdef AFSTFT_USE_FLOAT_COMPLEX
    free(frequencyDomainData);
#else
    for(c=0; c<numChannels; c++){
        free(frequencyDomainData[c].re);
        free(frequencyDomainData[c].im);
    }
    free(frequencyDomainData);
#endif
}

void test__smb_pitchShifter(void){
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
}

void test__sortf(void){
    float* values;
    int* sortedIdx;
    int i;

    /* Config */
    const int numValues = 10e5;

    /* Prep */
    sortedIdx = malloc1d(numValues*sizeof(int));
    values = malloc1d(numValues*sizeof(float));
    rand_m1_1(values, numValues); /* populate with random numbers */

    /* Sort in accending order */
    sortf(values, NULL, sortedIdx, numValues, 0);

    /* Check that the next value is either the same or greater than current value */
    for(i=0; i<numValues-1; i++)
        TEST_ASSERT_TRUE(values[sortedIdx[i]]<=values[sortedIdx[i+1]]);

    /* Sort in decending order */
    sortf(values, NULL, sortedIdx, numValues, 1);

    /* Check that the next value is either the same or less than current value */
    for(i=0; i<numValues-1; i++)
        TEST_ASSERT_TRUE(values[sortedIdx[i]]>=values[sortedIdx[i+1]]);

    /* clean-up */
    free(values);
    free(sortedIdx);
}

void test__realloc2d_r(void){
    int s, r, i, j, k;
    typedef struct _test_data{
        int ID;
        float val1, val2;
    }test_data;
    test_data** test;

    /* Configure reference data structures */
    test_data reference[6][6];
    for(i=0, k=0; i<6; i++){
        for(j=0; j<6; j++, k++){
            reference[i][j].ID = k;
            rand_m1_1(&reference[i][j].val1,1);
            rand_m1_1(&reference[i][j].val2,1);
        }
    }

    /* Starting size */
    test = (test_data**)malloc2d(1,3,sizeof(test_data));
    for(s=0; s<1; s++)
        for(r=0; r<3; r++)
            memcpy(&test[s][r], &reference[s][r], sizeof(test_data));

    /* Check that increasing the size of the array, still retains the previous data */
    test = (test_data**)realloc2d_r((void**)test, 4, 3, 1, 3, sizeof(test_data));
    for(i=0; i<1; i++){
        for(j=0; j<3; j++){
            TEST_ASSERT_TRUE(test[i][j].ID==reference[i][j].ID);
            TEST_ASSERT_TRUE(test[i][j].val1==reference[i][j].val1);
            TEST_ASSERT_TRUE(test[i][j].val2==reference[i][j].val2);
        }
    }

    /* Check that new data can then be added and indexed correctly */
    for(s=1; s<4; s++)
        for(r=0; r<3; r++)
             memcpy(&test[s][r], &reference[s][r], sizeof(test_data));
    for(i=0; i<4; i++){
        for(j=0; j<3; j++){
            TEST_ASSERT_TRUE(test[i][j].ID==reference[i][j].ID);
            TEST_ASSERT_TRUE(test[i][j].val1==reference[i][j].val1);
            TEST_ASSERT_TRUE(test[i][j].val2==reference[i][j].val2);
        }
    }

    /* Check that the array can be shrunk, but still retain the original data (except truncated) */
    test = (test_data**)realloc2d_r((void**)test, 4, 2, 4, 3, sizeof(test_data));
    for(i=0; i<4; i++){
        for(j=0; j<2; j++){
            TEST_ASSERT_TRUE(test[i][j].ID==reference[i][j].ID);
            TEST_ASSERT_TRUE(test[i][j].val1==reference[i][j].val1);
            TEST_ASSERT_TRUE(test[i][j].val2==reference[i][j].val2);
        }
    }

    /* clean-up */
    free(test);
}

void test__getSHreal_recur(void){
    int i, j;

    /* Config */
    /* In general, this recusive alternative is well below this tolerance value.
     * However, the error does get larger for higher-orders and when dir[1]
     * is near 0. */
    float acceptedTolerance = 0.005f;
    const int order = 7;
    const int nSH = ORDER2NSH(order);
    float dir[2];

    /* Check that the output of getSHreal_recur matches that of getSH_recur */
    float Yr[nSH];
    float Y[nSH];
    for(i=0; i<1e5; i++){
        rand_m1_1(&dir[0] , 1);
        rand_m1_1(&dir[1] , 1);
        dir[0] *= M_PI;
        dir[1] *= M_PI/2.0f;
        getSHreal_recur(order, (float*)dir, 1, (float*)Yr);
        getSHreal(order, (float*)dir, 1, (float*)Y);
        for(j=0; j<nSH; j++)
            TEST_ASSERT_FLOAT_WITHIN(acceptedTolerance, Yr[j], Y[j]);
            //assert(fabsf(Yr[j]-Y[j])<acceptedTolerance);
    }
}
