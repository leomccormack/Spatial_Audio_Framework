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
 * @file saf_test.c
 * @brief testing program for the Spatial_Audio_Framework
 *
 * @author Leo McCormack
 * @date 27.04.2020
 */

#include "saf_test.h"
#include "unity.h"   /* unit testing suite */
#include "timer.h"   /* for timing the individual tests */
#include "saf.h"     /* master framework include header */

#ifdef ENABLE_SAF_EXAMPLES_TESTS
/* SAF example headers: */
# include "ambi_bin.h"
# include "ambi_dec.h"
# include "ambi_drc.h"
# include "ambi_enc.h"
# include "array2sh.h"
# include "beamformer.h"
# include "binauraliser.h"
# include "dirass.h"
# include "matrixconv.h"
# include "multiconv.h"
# include "panner.h"
# include "pitch_shifter.h"
# include "powermap.h"
# include "rotator.h"
# include "sldoa.h"
#endif /* ENABLE_SAF_EXAMPLES_TESTS */

/* ========================================================================== */
/*                    Prototypes for available unit tests                     */
/* ========================================================================== */

/**
 * Testing for perfect recontruction of the saf_stft (when configured for 50%
 * window overlap) */
void test__saf_stft_50pc_overlap(void);
/**
 * Testing for perfect recontruction of the saf_stft (when configured for linear
 * time-invariant (LTI) filtering applications) */
void test__saf_stft_LTI(void);
/**
 * Testing the ims shoebox simulator, when applying the echograms in the time-
 * domain */
void test__ims_shoebox_TD(void);
/**
 * Testing the ims shoebox simulator, when generating room impulse respones
 * (RIRs) from the computed echograms */
void test__ims_shoebox_RIR(void);
/**
 * Testing the forward and backward real-(half)complex FFT (saf_rfft) */
void test__saf_rfft(void);
/**
 * Testing the saf_matrixConv */
void test__saf_matrixConv(void);
#ifdef AFSTFT_USE_FLOAT_COMPLEX
/**
 * Testing the alias-free STFT filterbank reconstruction */
void test__afSTFTMatrix(void);
#endif
/**
 * Testing the alias-free STFT filterbank reconstruction */
void test__afSTFT(void);
/**
 * Testing the smb_pitchShifter */
void test__smb_pitchShifter(void);
/**
 * Testing the sortf function (sorting real floating point numbers) */
void test__sortf(void);
/**
 * Testing the sortz function (sorting complex double-floating point numbers) */
void test__sortz(void);
/**
 * Testing the cmplxPairUp function (grouping up conjugate symmetric values) */
void test__cmplxPairUp(void);
/**
 * Testing the realloc2d_r function (reallocating 2-D array, while retaining the
 * previous data order; except truncated or extended) */
void test__realloc2d_r(void);
/**
 * Testing that the getSHreal_recur function is somewhat numerically simular
 * to the getSHreal function */
void test__getSHreal_recur(void);
/**
 * Testing that the coefficients computed with butterCoeffs are numerically
 * similar to the "butter" function in Matlab */
void test__butterCoeffs(void);
/**
 * Testing that the faf_IIRFilterbank can re-construct the original signal power
 */
void test__faf_IIRFilterbank(void);

#ifdef ENABLE_SAF_EXAMPLES_TESTS
/**
 * Testing the SAF ambi_bin example (this may also serve as a tutorial on how
 * to use it) */
void test__saf_example_ambi_bin(void);
/**
 * Testing the SAF ambi_dec example (this may also serve as a tutorial on how
 * to use it) */
void test__saf_example_ambi_dec(void);
/**
 * Testing the SAF array2sh example (this may also serve as a tutorial on how
 * to use it) */
void test__saf_example_array2sh(void);
#endif /* ENABLE_SAF_EXAMPLES_TESTS */


/* ========================================================================== */
/*                                 Test Config                                */
/* ========================================================================== */

static tick_t start;      /**< Start time for whole test program */
static tick_t start_test; /**< Start time for the current unit test */
/** Called before each unit test is executed */
void setUp(void) { start_test = timer_current(); }
/** Called after each unit test is executed */
void tearDown(void) { }
/** Displays the time taken to run the current unit test */
static void timerResult(void) {
    printf("    (Time elapsed: %lfs) \n", (double)timer_elapsed(start_test));
}

#undef RUN_TEST
/** A custom Unity RUN_TEST, which calls timerResult() upon exiting the test */
#define RUN_TEST(testfunc)  UNITY_NEW_TEST(#testfunc) \
if (TEST_PROTECT()) {  setUp();  testfunc();  } \
if (TEST_PROTECT() && (!TEST_IS_IGNORED))  {tearDown(); } \
UnityConcludeTest(); timerResult();

/** Main test program */
int main_test(void) {
    printf("%s", SAF_VERSION_BANNER);
    printf("Executing the Spatial_Audio_Framework unit testing program:\n\n");

    /* initialise */
    timer_lib_initialize();
    start = timer_current();
    UNITY_BEGIN();

    /* run each unit test */
    RUN_TEST(test__saf_stft_50pc_overlap);
    RUN_TEST(test__saf_stft_LTI);
    RUN_TEST(test__ims_shoebox_RIR);
    RUN_TEST(test__ims_shoebox_TD);
    RUN_TEST(test__saf_rfft);
    RUN_TEST(test__saf_matrixConv);
#ifdef AFSTFT_USE_FLOAT_COMPLEX
    RUN_TEST(test__afSTFTMatrix);
#endif
    RUN_TEST(test__afSTFT);
    RUN_TEST(test__smb_pitchShifter);
    RUN_TEST(test__sortf);
    RUN_TEST(test__sortz);
    RUN_TEST(test__cmplxPairUp);
    RUN_TEST(test__realloc2d_r);
    RUN_TEST(test__getSHreal_recur);
    RUN_TEST(test__butterCoeffs);
    RUN_TEST(test__faf_IIRFilterbank);
#ifdef ENABLE_SAF_EXAMPLES_TESTS
    RUN_TEST(test__saf_example_ambi_bin);
    RUN_TEST(test__saf_example_ambi_dec);
    RUN_TEST(test__saf_example_array2sh);
#endif /* ENABLE_SAF_EXAMPLES_TESTS */

    /* close */
    timer_lib_shutdown();
    printf("\nTotal time elapsed: %lfs", (double)timer_elapsed(start));
    return UNITY_END();
}


/* ========================================================================== */
/*                                 Unit Tests                                 */
/* ========================================================================== */

void test__saf_stft_50pc_overlap(void){
    int frame, winsize, hopsize, nFrames, ch, i, nBands, nTimesSlots, band;
    void* hSTFT;
    float** insig, **outsig, **inframe, **outframe;
    float_complex*** inspec, ***outspec;

    /* prep */
    const float acceptedTolerance = 0.000001f;
    const int fs = 48000;
    const int signalLength = 1*fs;
    const int framesize = 512;
    const int nCHin = 62;
    const int nCHout = 64;
    insig = (float**)malloc2d(nCHin,signalLength,sizeof(float)); /* One second long */
    outsig = (float**)malloc2d(nCHout,signalLength,sizeof(float));
    inframe = (float**)malloc2d(nCHin,framesize,sizeof(float));
    outframe = (float**)malloc2d(nCHout,framesize,sizeof(float));
    rand_m1_1(FLATTEN2D(insig), nCHin*signalLength); /* populate with random numbers */

    /* Set-up STFT for 50% overlapping */
    winsize = 128;
    hopsize = winsize/2;
    nBands = winsize+1;
    nTimesSlots = framesize/hopsize;
    inspec = (float_complex***)malloc3d(nBands, nCHin, nTimesSlots, sizeof(float_complex));
    outspec = (float_complex***)malloc3d(nBands, nCHout, nTimesSlots, sizeof(float_complex));
    saf_stft_create(&hSTFT, winsize, hopsize, nCHin, nCHout, SAF_STFT_BANDS_CH_TIME);
    saf_stft_channelChange(hSTFT, 123, 7);        /* messing about */
    saf_stft_flushBuffers(hSTFT);                 /* messing about */
    saf_stft_channelChange(hSTFT, nCHin, nCHout); /* change back */

    /* Pass insig through STFT, block-wise processing */
    nFrames = (int)((float)signalLength/(float)framesize);
    for(frame = 0; frame<nFrames; frame++){
        /* Forward */
        for(ch=0; ch<nCHin; ch++)
            memcpy(inframe[ch], &insig[ch][frame*framesize], framesize*sizeof(float));
        saf_stft_forward(hSTFT, inframe, framesize, inspec);

        /* Copy first channel of inspec to all outspec channels */
        for(band=0; band<nBands; band++)
            for(ch=0; ch<nCHout; ch++)
                memcpy(outspec[band][ch], inspec[band][0], nTimesSlots*sizeof(float_complex));

        /* Backward */
        saf_stft_backward(hSTFT, outspec, framesize, outframe);
        for(ch=0; ch<nCHout; ch++)
            memcpy(&outsig[ch][frame*framesize], outframe[ch], framesize*sizeof(float));
    }

    /* Check that input==output (given some numerical precision) */
    for(i=0; i<signalLength-framesize; i++)
        TEST_ASSERT_TRUE( fabsf(insig[0][i] - outsig[0][i+hopsize]) <= acceptedTolerance );

    /* Clean-up */
    saf_stft_destroy(&hSTFT);
    free(insig);
    free(outsig);
    free(inframe);
    free(outframe);
    free(inspec);
    free(outspec);
}

void test__saf_stft_LTI(void){
    int frame, winsize, hopsize, nFrames, ch, i, nBands, nTimesSlots, band;
    void* hSTFT;
    float** insig, **outsig, **inframe, **outframe;
    float_complex*** inspec, ***outspec;

    /* prep */
    const float acceptedTolerance = 0.000001f;
    const int fs = 48000;
    const int framesize = 128;
    const int nCHin = 62;
    const int nCHout = 64;
    insig = (float**)malloc2d(nCHin,fs,sizeof(float)); /* One second long */
    outsig = (float**)malloc2d(nCHout,fs,sizeof(float));
    inframe = (float**)malloc2d(nCHin,framesize,sizeof(float));
    outframe = (float**)malloc2d(nCHout,framesize,sizeof(float));
    rand_m1_1(FLATTEN2D(insig), nCHin*fs); /* populate with random numbers */

    /* Set-up STFT suitable for LTI filtering applications */
    winsize = hopsize = 128;
    nBands = winsize+1;
    nTimesSlots = framesize/hopsize;
    inspec = (float_complex***)malloc3d(nBands, nCHin, nTimesSlots, sizeof(float_complex));
    outspec = (float_complex***)malloc3d(nBands, nCHout, nTimesSlots, sizeof(float_complex));
    saf_stft_create(&hSTFT, winsize, hopsize, nCHin, nCHout, SAF_STFT_BANDS_CH_TIME);

    /* Pass insig through STFT, block-wise processing */
    nFrames = (int)((float)fs/(float)framesize);
    for(frame = 0; frame<nFrames; frame++){
        /* Forward */
        for(ch=0; ch<nCHin; ch++)
            memcpy(inframe[ch], &insig[ch][frame*framesize], framesize*sizeof(float));
        saf_stft_forward(hSTFT, inframe, framesize, inspec);

        /* Copy first channel of inspec to all outspec channels */
        for(band=0; band<nBands; band++)
            for(ch=0; ch<nCHout; ch++)
                memcpy(outspec[band][ch], inspec[band][0], nTimesSlots*sizeof(float_complex));

        /* Backward */
        saf_stft_backward(hSTFT, outspec, framesize, outframe);
        for(ch=0; ch<nCHout; ch++)
            memcpy(&outsig[ch][frame*framesize], outframe[ch], framesize*sizeof(float));
    }

    /* Check that input==output (given some numerical precision) */
    for(i=0; i<fs-framesize; i++)
        TEST_ASSERT_TRUE( fabsf(insig[0][i] - outsig[63][i]) <= acceptedTolerance );

    /* Clean-up */
    saf_stft_destroy(&hSTFT);
    free(insig);
    free(outsig);
    free(inframe);
    free(outframe);
    free(inspec);
    free(outspec);
}

void test__ims_shoebox_TD(void){
    void* hIms;
    float maxTime_s;
    float mov_src_pos[3], mov_rec_pos[3];
    float** src_sigs, ***rec_sh_outsigs;
    long sourceIDs[4], receiverIDs[1];
    int i;

    /* Config */
    const int signalLength = 10000;
    const int sh_order = 3;
    const int nBands = 5;
    const float abs_wall[5][6] =  /* Absorption Coefficients per Octave band, and per wall */
      { {0.180791250f, 0.207307300f, 0.134990800f, 0.229002250f, 0.212128400f, 0.241055000f},
        {0.225971250f, 0.259113700f, 0.168725200f, 0.286230250f, 0.265139600f, 0.301295000f},
        {0.258251250f, 0.296128100f, 0.192827600f, 0.327118250f, 0.303014800f, 0.344335000f},
        {0.301331250f, 0.345526500f, 0.224994001f, 0.381686250f, 0.353562000f, 0.401775000f},
        {0.361571250f, 0.414601700f, 0.269973200f, 0.457990250f, 0.424243600f, 0.482095000f} };
    const float src_pos[3]  = {5.1f, 6.0f, 1.1f};
    const float src2_pos[3] = {2.1f, 1.0f, 1.3f};
    const float src3_pos[3] = {3.1f, 5.0f, 2.3f};
    const float src4_pos[3] = {7.1f, 2.0f, 1.4f};
    const float rec_pos[3]  = {8.8f, 5.5f, 0.9f};

    /* Allocate memory for 4 sources and 1 spherical harmonic receiver */
    rec_sh_outsigs = (float***)malloc3d(1, ORDER2NSH(sh_order), signalLength, sizeof(float));
    src_sigs = (float**)malloc2d(4, signalLength, sizeof(float));
    rand_m1_1(FLATTEN2D(src_sigs), 4*signalLength);

    /* Set-up the shoebox room simulator for these four sources and SH receiver */
    ims_shoebox_create(&hIms, 10, 7, 3, (float*)abs_wall, 250.0f, nBands, 343.0f, 48e3f);
    sourceIDs[0] = ims_shoebox_addSource(hIms, (float*)src_pos, &src_sigs[0]);
    sourceIDs[1] = ims_shoebox_addSource(hIms, (float*)src2_pos, &src_sigs[1]);
    sourceIDs[2] = ims_shoebox_addSource(hIms, (float*)src3_pos, &src_sigs[2]);
    sourceIDs[3] = ims_shoebox_addSource(hIms, (float*)src4_pos, &src_sigs[3]);
    receiverIDs[0] = ims_shoebox_addReceiverSH(hIms, sh_order, (float*)rec_pos, &rec_sh_outsigs[0]);

    /* Moving source No.1 and the receiver */
    maxTime_s = 0.025f; /* 50ms */
    memcpy(mov_src_pos, src_pos, 3*sizeof(float));
    memcpy(mov_rec_pos, rec_pos, 3*sizeof(float));
    for(i=0; i<5; i++){
        mov_src_pos[1] = 2.0f + (float)i/100.0f;
        mov_rec_pos[0] = 3.0f + (float)i/100.0f;
        ims_shoebox_updateSource(hIms, sourceIDs[0], mov_src_pos);
        ims_shoebox_updateReceiver(hIms, receiverIDs[0], mov_rec_pos);
        ims_shoebox_computeEchograms(hIms, maxTime_s);
        ims_shoebox_applyEchogramTD(hIms, receiverIDs[0], signalLength, 0);
    }

    /* clean-up */
    free(src_sigs);
    free(rec_sh_outsigs);
    ims_shoebox_destroy(&hIms);
}

void test__ims_shoebox_RIR(void){
    void* hIms;
    float maxTime_s;
    float mov_src_pos[3], mov_rec_pos[3];
    long sourceID_1, sourceID_2, sourceID_3, sourceID_4, sourceID_5, receiverID;
    int i;

    /* Config */
    const int sh_order = 3;
    const int nBands = 7;
    const float abs_wall[7][6] =  /* Absorption Coefficients per Octave band, and per wall */
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

    /* Set-up the shoebox room simulator, with two sources and one spherical harmonic receiver */
    ims_shoebox_create(&hIms, 10, 7, 3, (float*)abs_wall, 125.0f, nBands, 343.0f, 48e3f);
    sourceID_1 = ims_shoebox_addSource(hIms, (float*)src_pos, NULL);
    sourceID_2 = ims_shoebox_addSource(hIms, (float*)src2_pos, NULL);
    receiverID = ims_shoebox_addReceiverSH(hIms, sh_order, (float*)rec_pos, NULL);

    /* Moving source No.1 and the receiver */
    maxTime_s = 0.05f; /* 50ms */
    memcpy(mov_src_pos, src_pos, 3*sizeof(float));
    memcpy(mov_rec_pos, rec_pos, 3*sizeof(float));
    for(i=0; i<0; i++){
        mov_src_pos[1] = 2.0f + (float)i/100.0f;
        mov_rec_pos[0] = 3.0f + (float)i/100.0f;
        ims_shoebox_updateSource(hIms, sourceID_1, mov_src_pos);
        ims_shoebox_updateReceiver(hIms, receiverID, mov_rec_pos);
        ims_shoebox_computeEchograms(hIms, maxTime_s);
        ims_shoebox_renderRIRs(hIms, 0);
    }

    /* Remove source No.1 */
    ims_shoebox_removeSource(hIms, 0);

    /* Add 3 more sources, then remove 2, and add one back again
     * (Just messing around, trying to trip up an IMS internal assertion) */
    sourceID_3 = ims_shoebox_addSource(hIms, (float*)src3_pos, NULL);
    sourceID_4 = ims_shoebox_addSource(hIms, (float*)src4_pos, NULL);
    sourceID_5 = ims_shoebox_addSource(hIms, (float*)src5_pos, NULL);
    ims_shoebox_removeSource(hIms, sourceID_3);
    ims_shoebox_removeSource(hIms, sourceID_4);
    sourceID_4 = ims_shoebox_addSource(hIms, (float*)src4_pos, NULL);

    /* Continue rendering */
    for(i=0; i<100; i++){
        mov_src_pos[1] = 2.0f + (float)i/1000.0f;
        mov_rec_pos[0] = 3.0f + (float)i/1000.0f;
        ims_shoebox_updateSource(hIms, sourceID_4, mov_src_pos);
        ims_shoebox_updateReceiver(hIms, receiverID, mov_rec_pos);
        ims_shoebox_computeEchograms(hIms, maxTime_s);
        ims_shoebox_renderRIRs(hIms, 0);
    }

    /* clean-up */
    ims_shoebox_destroy(&hIms);
}

void test__saf_matrixConv(void){
    int i, frame;
    float** inputTD, **outputTD, **inputFrameTD, **outputFrameTD;
    float*** filters;
    void* hMatrixConv;

    /* config */
    const int signalLength = 48000;
    const int hostBlockSize = 1024;
    const int filterLength = 512;
    const int nInputs = 64;
    const int nOutputs = 64;

    /* prep */
    inputTD = (float**)malloc2d(nInputs, signalLength, sizeof(float));
    outputTD = (float**)malloc2d(nOutputs, signalLength, sizeof(float));
    inputFrameTD = (float**)malloc2d(nInputs, hostBlockSize, sizeof(float));
    outputFrameTD = (float**)calloc2d(nOutputs, hostBlockSize, sizeof(float));
    filters = (float***)malloc3d(nOutputs, nInputs, filterLength, sizeof(float));
    rand_m1_1(FLATTEN3D(filters), nOutputs*nInputs*filterLength);
    rand_m1_1(FLATTEN2D(inputTD), nInputs*signalLength);
    saf_matrixConv_create(&hMatrixConv, hostBlockSize, FLATTEN3D(filters), filterLength,
                          nInputs, nOutputs, 0);

    /* Apply */
    for(frame = 0; frame<(int)signalLength/hostBlockSize; frame++){
        for(i = 0; i<nInputs; i++)
            memcpy(inputFrameTD[i], &inputTD[i][frame*hostBlockSize], hostBlockSize*sizeof(float));

         saf_matrixConv_apply(hMatrixConv, FLATTEN2D(inputFrameTD), FLATTEN2D(outputFrameTD));

        for(i = 0; i<nOutputs; i++)
            memcpy(&outputTD[i][frame*hostBlockSize], outputFrameTD[i], hostBlockSize*sizeof(float));
    }

    /* Clean-up */
    free(inputTD);
    free(outputTD);
    free(inputFrameTD);
    free(outputFrameTD);
    free(filters);
    saf_matrixConv_destroy(&hMatrixConv);
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
}

#ifdef AFSTFT_USE_FLOAT_COMPLEX
void test__afSTFTMatrix(void){
    int idx,frameIdx,c,t;
    int numChannels;
    float** inputTimeDomainData, **outputTimeDomainData, **tempFrame;
    float_complex*** frequencyDomainData;

    /* Config */
    const float acceptedTolerance_dB = -50.0f;
    const int nTestFrames = 250;
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
    rand_m1_1(FLATTEN2D(inputTimeDomainData), numChannels*lSig); /* populate with random numbers */

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
            TEST_ASSERT_TRUE( 20.0f*log10f(fabsf(inputTimeDomainData[c][t]-outputTimeDomainData[c][t]))<= acceptedTolerance_dB );
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
    const float acceptedTolerance_dB = -50.0f;
    const int nTestHops = 2000;
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
    rand_m1_1(FLATTEN2D(inputTimeDomainData), numChannels*lSig); /* populate with random numbers */

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
            TEST_ASSERT_TRUE( 20.0f*log10f(fabsf(inputTimeDomainData[c][t]-outputTimeDomainData[c][t]))<= acceptedTolerance_dB );
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
    void* hPS, *hFFT;
    float frequency;
    int i, smbLatency, ind;

    /* Config */
    const int sampleRate = 48000;
    const int FFTsize = 8192;
    const int osfactor = 16;
    const int nSamples = 8*FFTsize;

    /* prep */
    smb_pitchShift_create(&hPS, 1, FFTsize, osfactor, (float)sampleRate);
    inputData = malloc1d(nSamples*sizeof(float));
    outputData = calloc1d(nSamples,sizeof(float));
    frequency = (float)sampleRate/8.0f;
    for(i=0; i<nSamples; i++) /* sine tone at quarter Nyquist: */
        inputData[i] = sinf(2.0f * M_PI * (float)i * frequency/(float)sampleRate);
    smbLatency = FFTsize - (FFTsize/osfactor);

    /* Pitch shift down one octave */
    smb_pitchShift_apply(hPS, 0.5, nSamples, inputData, outputData);

    /* Take FFT, the bin with the highest energy should correspond to 1/8 Nyquist */
    float_complex* out_fft; // [nSamples / 2 + 1];
    out_fft = malloc1d((nSamples / 2 + 1) * sizeof(float_complex));
    saf_rfft_create(&hFFT, nSamples);
    saf_rfft_forward(hFFT, outputData, out_fft);
    utility_cimaxv(out_fft, nSamples/2+1, &ind);
    TEST_ASSERT_TRUE(ind == nSamples/16);
 
    /* clean-up */
    smb_pitchShift_destroy(&hPS);
    saf_rfft_destroy(&hFFT);
    free(inputData);
    free(outputData);
    free(out_fft);
}

void test__sortf(void){
    float* values;
    int* sortedIdx;
    int i;

    /* Config */
    const int numValues = 10000;

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

void test__sortz(void){
    int i;
    const int N = 36;
    double_complex vals[36] ={
        cmplx(1.0, 1.0), cmplx(7.0, 1.0),  cmplx(10.0, 5.0),
        cmplx(12.0, 4.0), cmplx(4.0, 4.0), cmplx(8.0, 0.0),
        cmplx(10.0, -1.0), cmplx(7.0, 5.0), cmplx(7.0, 2.0),
        cmplx(5.0, 1.0), cmplx(4.0, -1.0), cmplx(32.0, 3.0),
        cmplx(32.0, 32.5), cmplx(25.0, 0.0), cmplx(2.0, -2.0),
        cmplx(7.0, -2.0), cmplx(1.0, -1.0), cmplx(12.0, -1.0),
        cmplx(2.0, -1.0), cmplx(4.0, 2.0), cmplx(10.0, 6.0),
        cmplx(5.0, 2.0), cmplx(32.0, 1.5), cmplx(7.0, -10.0),
        cmplx(1.0, -1.5), cmplx(4.0, 25.0), cmplx(3.0, 2.0),
        cmplx(1.0, 4.5), cmplx(10.0, 5.0), cmplx(10.0, 2.0),
        cmplx(10.0, -3.5), cmplx(30.0, -10.0), cmplx(7.0, -12.0),
        cmplx(1.0, -13.5), cmplx(12.0, -12.0), cmplx(32.0, 23.0)
    };
    double_complex sorted_vals[36];

    /* Sort assending order */
    sortz(vals, sorted_vals, N, 0);

    /* Check that the next real(value) is either the same or greater than current real(value) */
    for(i=0; i<N-1; i++)
        TEST_ASSERT_TRUE(creal(sorted_vals[i])<=creal(sorted_vals[i+1]));

    /* Check that if the next real(value) is the same as the current real(value), then
     * the next imag(value) is greater that the current image(value)*/
    for(i=0; i<N-1; i++)
        if(fabs(creal(sorted_vals[i])-creal(sorted_vals[i+1])) < 0.00001 )
            TEST_ASSERT_TRUE(cimag(sorted_vals[i])<=cimag(sorted_vals[i+1]));

    /* Sort decending order */
    sortz(vals, sorted_vals, N, 1);

    /* Check that the next real(value) is either the same or smaller than current real(value) */
    for(i=0; i<N-1; i++)
        TEST_ASSERT_TRUE(creal(sorted_vals[i])>=creal(sorted_vals[i+1]));

    /* Check that if the next real(value) is the same as the current real(value), then
     * the next imag(value) is smaller that the current image(value)*/
    for(i=0; i<N-1; i++)
        if(fabs(creal(sorted_vals[i])-creal(sorted_vals[i+1])) < 0.00001 )
            TEST_ASSERT_TRUE(cimag(sorted_vals[i])>=cimag(sorted_vals[i+1]));
}

void test__cmplxPairUp(void){
    int i;
    const int N = 36;
    double_complex vals[36] ={
        cmplx(1.0, 1.0), cmplx(7.0, 1.0),  cmplx(10.0, 5.0),
        cmplx(12.0, 4.0), cmplx(4.0, 4.0), cmplx(8.0, 0.0),
        cmplx(10.0, -1.0), cmplx(7.0, 5.0), cmplx(7.0, 2.0),
        cmplx(5.0, 1.0), cmplx(4.0, -1.0), cmplx(32.0, 3.0),
        cmplx(32.0, 32.5), cmplx(25.0, 0.0), cmplx(2.0, -2.0),
        cmplx(7.0, -2.0), cmplx(1.0, -1.0), cmplx(12.0, -1.0),
        cmplx(2.0, -1.0), cmplx(4.0, 2.0), cmplx(10.0, 6.0),
        cmplx(5.0, 0.0), cmplx(32.0, 1.5), cmplx(7.0, -10.0),
        cmplx(1.0, -1.5), cmplx(4.0, 25.0), cmplx(3.0, 2.0),
        cmplx(1.0, 0.0), cmplx(10.0, 5.0), cmplx(10.0, 2.0),
        cmplx(10.0, -3.5), cmplx(30.0, -10.0), cmplx(7.0, -12.0),
        cmplx(1.0, -13.5), cmplx(12.0, -12.0), cmplx(32.0, 23.0)
    };
    double_complex sorted_vals[36];

    /* Sort assending order */
    cmplxPairUp(vals, sorted_vals, N);

    /* Check that the next real(value) is either the same or greater than current real(value),
     * Ignoring purely real numbers */
    for(i=0; i<N-1; i++)
        if( !(fabs(cimag(sorted_vals[i])) < 0.0001) && !(fabs(cimag(sorted_vals[i+1])) < 0.0001) )
            TEST_ASSERT_TRUE(creal(sorted_vals[i])<=creal(sorted_vals[i+1]));

    /* Check that the next real(value) is either the same or greater than current real(value),
     * Only considering purely real numbers */
    for(i=0; i<N-1; i++)
        if( (fabs(cimag(sorted_vals[i])) < 0.0001) && (fabs(cimag(sorted_vals[i+1])) < 0.0001) )
            TEST_ASSERT_TRUE(creal(sorted_vals[i])<=creal(sorted_vals[i+1]));

    /* Check that if the next real(value) is the same as the current real(value), then
     * the next imag(value) is greater that the current image(value)
     * Ignoring purely real numbers */
    for(i=0; i<N-1; i++)
        if(fabs(creal(sorted_vals[i])-creal(sorted_vals[i+1])) < 0.00001 )
            if( !(fabs(cimag(sorted_vals[i])) < 0.0001) && !(fabs(cimag(sorted_vals[i+1])) < 0.0001) )
                TEST_ASSERT_TRUE(cimag(sorted_vals[i])<=cimag(sorted_vals[i+1]));
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
    /* In general, the values from this recusive alternative are well below this
     * tolerance value. However, the error does get larger for higher-orders and
     * when dir[1] is near 0. */
    float acceptedTolerance = 0.005f;
    const int order = 15;
    const int nSH = ORDER2NSH(order);
    float dir[2];

    /* Check that the output of getSHreal_recur matches that of getSH_recur */
    float Yr[ORDER2NSH(15)];
    float Y[ORDER2NSH(15)];
    for(i=0; i<1e3; i++){
        rand_m1_1(&dir[0] , 1);
        rand_m1_1(&dir[1] , 1);
        dir[0] *= M_PI;
        dir[1] *= M_PI/2.0f;
        getSHreal_recur(order, (float*)dir, 1, (float*)Yr);
        getSHreal(order, (float*)dir, 1, (float*)Y);
        for(j=0; j<nSH; j++)
            TEST_ASSERT_FLOAT_WITHIN(acceptedTolerance, Yr[j], Y[j]);
    }
}

void test__butterCoeffs(void){
    int i;
    float fs, cutoff_freq, cutoff_freq2;
    int order;

    /* Config */
    const double acceptedTolerance = 0.00001f;

    /* 1st order Low-pass filter */
    fs = 48e3f;
    cutoff_freq = 3000.0f;
    order = 1;
    double a_test1[2], b_test1[2];
    butterCoeffs(BUTTER_FILTER_LPF, order, cutoff_freq, 0.0f, fs, (double*)b_test1, (double*)a_test1);
    const double a_ref1[2] = {1,-0.668178637919299};
    const double b_ref1[2] = {0.165910681040351,0.165910681040351};
    for(i=0; i<2; i++){ /* Compare with the values given by Matlab's butter function */
        TEST_ASSERT_FLOAT_WITHIN(acceptedTolerance, a_test1[i], a_ref1[i]);
        TEST_ASSERT_FLOAT_WITHIN(acceptedTolerance, b_test1[i], b_ref1[i]);
    }

    /* 2nd order Low-pass filter */
    fs = 48e3f;
    cutoff_freq = 12000.0f;
    order = 2;
    double a_test2[3], b_test2[3];
    butterCoeffs(BUTTER_FILTER_LPF, order, cutoff_freq, 0.0f, fs, (double*)b_test2, (double*)a_test2);
    const double a_ref2[3] = {1.0,-2.22044604925031e-16,0.171572875253810};
    const double b_ref2[3] = {0.292893218813452,0.585786437626905,0.292893218813452};
    for(i=0; i<3; i++){ /* Compare with the values given by Matlab's butter function */
        TEST_ASSERT_FLOAT_WITHIN(acceptedTolerance, a_test2[i], a_ref2[i]);
        TEST_ASSERT_FLOAT_WITHIN(acceptedTolerance, b_test2[i], b_ref2[i]);
    }

    /* 3rd order Low-pass filter */
    fs = 48e3f;
    cutoff_freq = 200.0f;
    order = 3;
    double a_test3[4], b_test3[4];
    butterCoeffs(BUTTER_FILTER_LPF, order, cutoff_freq, 0.0f, fs, (double*)b_test3, (double*)a_test3);
    const double a_ref3[4] = {1.0,-2.94764161678340,2.89664496645376,-0.948985866903327};
    const double b_ref3[4] = {2.18534587909103e-06,6.55603763727308e-06,6.55603763727308e-06,2.18534587909103e-06};
    for(i=0; i<4; i++){ /* Compare with the values given by Matlab's butter function */
        TEST_ASSERT_FLOAT_WITHIN(acceptedTolerance, a_test3[i], a_ref3[i]);
        TEST_ASSERT_FLOAT_WITHIN(acceptedTolerance, b_test3[i], b_ref3[i]);
    }

    /* 6th order Low-pass filter */
    fs = 48e3f;
    cutoff_freq = 1e3f;
    order = 6;
    double a_test4[7], b_test4[7];
    butterCoeffs(BUTTER_FILTER_LPF, order, cutoff_freq, 0.0f, fs, (double*)b_test4, (double*)a_test4);
    const double a_ref4[7] = {1,-5.49431292177096,12.5978414666894,-15.4285267903275,10.6436770055305,-3.92144696766748,0.602772146971300};
    const double b_ref4[7] = {6.15535184628202e-08,3.69321110776921e-07,9.23302776942303e-07,1.23107036925640e-06,9.23302776942303e-07,3.69321110776921e-07,6.15535184628202e-08};
    for(i=0; i<7; i++){ /* Compare with the values given by Matlab's butter function */
        TEST_ASSERT_FLOAT_WITHIN(acceptedTolerance, a_test4[i], a_ref4[i]);
        TEST_ASSERT_FLOAT_WITHIN(acceptedTolerance, b_test4[i], b_ref4[i]);
    }

    /* 3rd order High-pass filter */
    fs = 48e3f;
    cutoff_freq = 3000.0f;
    order = 3;
    double a_test5[4], b_test5[4];
    butterCoeffs(BUTTER_FILTER_HPF, order, cutoff_freq, 0.0f, fs, (double*)b_test5, (double*)a_test5);
    const double a_ref5[4] = {1,-2.21916861831167,1.71511783003340,-0.453545933365530};
    const double b_ref5[4] = {0.673479047713825,-2.02043714314147,2.02043714314147,-0.673479047713825};
    for(i=0; i<4; i++){ /* Compare with the values given by Matlab's butter function */
        TEST_ASSERT_FLOAT_WITHIN(acceptedTolerance, a_test5[i], a_ref5[i]);
        TEST_ASSERT_FLOAT_WITHIN(acceptedTolerance, b_test5[i], b_ref5[i]);
    }

    /* 4th order High-pass filter */
    fs = 48e3f;
    cutoff_freq = 100.0;
    order = 4;
    double a_test6[5], b_test6[5];
    butterCoeffs(BUTTER_FILTER_HPF, order, cutoff_freq, 0.0f, fs, (double*)b_test6, (double*)a_test6);
    const double a_ref6[5] = {1.0,-3.96579438007005,5.89796693861409,-3.89854491737242,0.966372387692057};
    const double b_ref6[5] = {0.983042413984288,-3.93216965593715,5.89825448390573,-3.93216965593715,0.983042413984288};
    for(i=0; i<5; i++){ /* Compare with the values given by Matlab's butter function */
        TEST_ASSERT_FLOAT_WITHIN(acceptedTolerance, a_test6[i], a_ref6[i]);
        TEST_ASSERT_FLOAT_WITHIN(acceptedTolerance, b_test6[i], b_ref6[i]);
    }

    /* 2nd order Band-pass filter */
    fs = 48e3f;
    cutoff_freq = 100.0;
    cutoff_freq2 = 400.0;
    order = 2;
    double a_test7[5], b_test7[5];
    butterCoeffs(BUTTER_FILTER_BPF, order, cutoff_freq, cutoff_freq2, fs, (double*)b_test7, (double*)a_test7);
    const double a_ref7[5] = {1.0,-3.94312581006024,5.83226704209421,-3.83511871130750,0.945977936232284};
    const double b_ref7[5] = {0.000375069616051004,0.0,-0.000750139232102008,0.0,0.000375069616051004};
    for(i=0; i<5; i++){ /* Compare with the values given by Matlab's butter function */
        TEST_ASSERT_FLOAT_WITHIN(acceptedTolerance, a_test7[i], a_ref7[i]);
        TEST_ASSERT_FLOAT_WITHIN(acceptedTolerance, b_test7[i], b_ref7[i]);
    }

    /* 3rd order Band-stop filter */
    fs = 48e3f;
    cutoff_freq = 240.0;
    cutoff_freq2 = 1600.0;
    order = 3;
    double a_test9[7], b_test9[7];
    butterCoeffs(BUTTER_FILTER_BSF, order, cutoff_freq, cutoff_freq2, fs, (double*)b_test9, (double*)a_test9);
    const double a_ref9[7] = {1,-5.62580309774365,13.2124846784594,-16.5822627287366,11.7304049556188,-4.43493124452282,0.700107676775329};
    const double b_ref9[7] = {0.836724592951539,-5.00379660039217,12.4847741945760,-16.6354041344203,12.4847741945760,-5.00379660039217,0.836724592951539};
    for(i=0; i<7; i++){ /* Compare with the values given by Matlab's butter function */
        TEST_ASSERT_FLOAT_WITHIN(acceptedTolerance, a_test9[i], a_ref9[i]);
        TEST_ASSERT_FLOAT_WITHIN(acceptedTolerance, b_test9[i], b_ref9[i]);
    }
}

void test__faf_IIRFilterbank(void){
    void* hFaF;
    int i, band;
    float* inSig, *outSig;
    float** outFrame, **outSig_bands;
    void* hFFT;
    float_complex* insig_fft, *outsig_fft;

    /* Config */
    const float acceptedTolerance_dB = 0.5f;
    const int signalLength = 256;
    const int frameSize = 16;
    float fs = 48e3;
    int order = 3;
    float fc[6] = {176.776695296637f, 353.553390593274f, 707.106781186547f, 1414.21356237309f, 2828.42712474619f, 5656.85424949238f};
    inSig = malloc1d(signalLength * sizeof(float));
    outSig_bands = (float**)malloc2d(7, signalLength, sizeof(float));
    outSig = calloc1d(signalLength, sizeof(float));

    insig_fft = malloc1d((signalLength / 2 + 1) * sizeof(float_complex));
    outsig_fft = malloc1d((signalLength / 2 + 1) * sizeof(float_complex));

    /* Impulse */
    memset(inSig, 0, signalLength*sizeof(float));
    inSig[0] = 1.0f;

    /* Pass impulse through filterbank */
    outFrame = (float**)malloc2d(7, frameSize, sizeof(float));
    faf_IIRFilterbank_create(&hFaF, order, (float*)fc, 6, fs, 512);
    for(i=0; i< signalLength/frameSize; i++){
        faf_IIRFilterbank_apply(hFaF, &inSig[i*frameSize], outFrame, frameSize);
        for(band=0; band<7; band++)
            memcpy(&outSig_bands[band][i*frameSize], outFrame[band], frameSize*sizeof(float));
    }
    faf_IIRFilterbank_destroy(&hFaF);

    /* Sum the individual bands */
    for(band=0; band<7; band++)
        utility_svvadd(outSig, outSig_bands[band], signalLength, outSig);

    /* Check that the magnitude difference between input and output is below 0.5dB */
    saf_rfft_create(&hFFT, signalLength);
    saf_rfft_forward(hFFT, inSig, insig_fft);
    saf_rfft_forward(hFFT, outSig, outsig_fft);
    for(i=0; i<signalLength/2+1; i++)
        TEST_ASSERT_FLOAT_WITHIN(acceptedTolerance_dB, 0.0f, 20.0f * log10f(cabsf( ccdivf(outsig_fft[i],insig_fft[i]) )));

    /* Now the same thing, but for 1st order */
    order = 1;
    faf_IIRFilterbank_create(&hFaF, order, (float*)fc, 6, fs, 512);
    for(i=0; i< signalLength/frameSize; i++){
        faf_IIRFilterbank_apply(hFaF, &inSig[i*frameSize], outFrame, frameSize);
        for(band=0; band<7; band++)
            memcpy(&outSig_bands[band][i*frameSize], outFrame[band], frameSize*sizeof(float));
    }
    faf_IIRFilterbank_destroy(&hFaF);
    memset(outSig, 0, signalLength*sizeof(float));
    for(band=0; band<7; band++)
        utility_svvadd(outSig, outSig_bands[band], signalLength, outSig);
    saf_rfft_forward(hFFT, outSig, outsig_fft);
    for(i=0; i<signalLength/2+1; i++)
        TEST_ASSERT_FLOAT_WITHIN(acceptedTolerance_dB, 0.0f, 20.0f * log10f(cabsf(ccdivf(outsig_fft[i], insig_fft[i]))));

    /* clean-up */
    saf_rfft_destroy(&hFFT);
    free(outFrame);
    free(inSig);
    free(outSig_bands);
    free(outSig);
    free(insig_fft);
    free(outsig_fft);
}

#ifdef ENABLE_SAF_EXAMPLES_TESTS
void test__saf_example_ambi_bin(void){
    int nSH, i;
    void* hAmbi;
    float leftEarEnergy, rightEarEnergy, direction_deg[2];
    float* inSig, *y;
    float** shSig, **binSig;

    /* Config */
    const int order = 4;
    const int fs = 48e3;
    const int signalLength = fs*2;

    /* Create and initialise an instance of ambi_bin */
    ambi_bin_create(&hAmbi);
    ambi_bin_init(hAmbi, fs); /* Cannot be called while "process" is on-going */

    /* Configure and initialise the ambi_bin codec */
    ambi_bin_setInputOrderPreset(hAmbi, (SH_ORDERS)order);
    ambi_bin_initCodec(hAmbi); /* Can be called whenever (thread-safe) */
    /* "initCodec" should be called after calling any of the "set" functions.
     * It should be noted that intialisations are only conducted if they are
     * needed, so calling this function periodically with a timer on a separate
     * thread is perfectly safe and viable. Also, if the intialisations take
     * longer than it takes to "process" the current block of samples, then the
     * output is simply muted/zeroed during this time. */

    /* Define input mono signal */
    nSH = ORDER2NSH(order);
    inSig = malloc1d(signalLength*sizeof(float));
    shSig = (float**)malloc2d(nSH,signalLength,sizeof(float));
    rand_m1_1(inSig, signalLength); /* Mono white-noise signal */

    /* Encode to get input spherical harmonic (Ambisonic) signal */
    direction_deg[0] = 90.0f; /* encode hard-left */
    direction_deg[1] = 0.0f;
    y = malloc1d(nSH*sizeof(float));
    getRSH(order, (float*)direction_deg, 1, y); /* SH plane-wave weights */
    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nSH, signalLength, 1, 1.0f,
                y, 1,
                inSig, signalLength, 0.0f,
                FLATTEN2D(shSig), signalLength);

    /* Decode to binaural */
    binSig = (float**)malloc2d(NUM_EARS,signalLength,sizeof(float));
    ambi_bin_process(hAmbi, shSig, binSig, nSH, NUM_EARS, signalLength);

    /* Assert that left ear energy is higher than the right ear */
    leftEarEnergy = rightEarEnergy = 0.0f;
    for(i=0; i<signalLength; i++){
        leftEarEnergy  += powf(fabsf(binSig[0][i]), 2.0f);
        rightEarEnergy += powf(fabsf(binSig[1][i]), 2.0f);
    }
    TEST_ASSERT_TRUE(leftEarEnergy>=rightEarEnergy);

    /* Clean-up */
    ambi_bin_destroy(&hAmbi);
    free(inSig);
    free(shSig);
    free(y);
    free(binSig);
}

void test__saf_example_ambi_dec(void){
    int nSH, i, j, max_ind;
    void* hAmbi;
    float loudspeakerEnergy[22], direction_deg[2];
    float* inSig, *y;
    float** shSig, **lsSig;

    /* Config */
    const int order = 4;
    const int fs = 48e3;
    const int signalLength = fs*2;

    /* Create and initialise an instance of ambi_dec */
    ambi_dec_create(&hAmbi);
    ambi_dec_init(hAmbi, fs); /* Cannot be called while "process" is on-going */

    /* Configure and initialise the ambi_dec codec */
    ambi_dec_setMasterDecOrder(hAmbi, (SH_ORDERS)order);
    /* 22.x loudspeaker layout, SAD decoder */
    ambi_dec_setOutputConfigPreset(hAmbi, LOUDSPEAKER_ARRAY_PRESET_22PX);
    ambi_dec_setDecMethod(hAmbi, DECODING_METHOD_SAD, 0/* low-freq decoder */);
    ambi_dec_setDecMethod(hAmbi, DECODING_METHOD_SAD, 1/* high-freq decoder */);
    ambi_dec_initCodec(hAmbi); /* Can be called whenever (thread-safe) */
    /* "initCodec" should be called after calling any of the "set" functions.
     * It should be noted that intialisations are only conducted if they are
     * needed, so calling this function periodically with a timer on a separate
     * thread is perfectly safe and viable. Also, if the intialisations take
     * longer than it takes to "process" the current block of samples, then the
     * output is simply muted/zeroed during this time. */

    /* Define input mono signal */
    nSH = ORDER2NSH(order);
    inSig = malloc1d(signalLength*sizeof(float));
    shSig = (float**)malloc2d(nSH,signalLength,sizeof(float));
    rand_m1_1(inSig, signalLength); /* Mono white-noise signal */

    /* Encode to get input spherical harmonic (Ambisonic) signal */
    direction_deg[0] = 90.0f; /* encode to loudspeaker direction: index 8 */
    direction_deg[1] = 0.0f;
    y = malloc1d(nSH*sizeof(float));
    getRSH(order, (float*)direction_deg, 1, y); /* SH plane-wave weights */
    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nSH, signalLength, 1, 1.0f,
                y, 1,
                inSig, signalLength, 0.0f,
                FLATTEN2D(shSig), signalLength);

    /* Decode to loudspeakers */
    lsSig = (float**)malloc2d(22,signalLength,sizeof(float));
    ambi_dec_process(hAmbi, shSig, lsSig, nSH, 22, signalLength);

    /* Assert that channel 8 (corresponding to the loudspeaker where the plane-
     * wave was encoded to) has the most energy */
    memset(loudspeakerEnergy, 0, 22*sizeof(float));
    for(i=0; i<signalLength; i++){
        for(j=0; j<22; j++)
            loudspeakerEnergy[j]  += powf(fabsf(lsSig[j][i]), 2.0f);
    }
    utility_simaxv(loudspeakerEnergy, 22, &max_ind);
    TEST_ASSERT_TRUE(max_ind==7);

    /* Clean-up */
    ambi_dec_destroy(&hAmbi);
    free(inSig);
    free(shSig);
    free(y);
    free(lsSig);
}

void test__saf_example_array2sh(void){
    int nSH, i, j;
    void* hA2sh, *safFFT, *hMC;
    float direction_deg[2], radius;
    float* inSig, *f;
    float** shSig, **inSig_32, **micSig, **h_array;
    double* kr;
    float_complex* tmp_H;
    float_complex*** H_array;

    /* Config */
    const int order = 4;
    const int fs = 48e3;
    const int signalLength = fs*2;
    const int nFFT = 1024;
    const int nBins = nFFT/2+1;

    /* Create and initialise an instance of array2sh for the Eigenmike32 */
    array2sh_create(&hA2sh);
    array2sh_init(hA2sh, fs); /* Cannot be called while "process" is on-going */
    array2sh_setPreset(hA2sh, MICROPHONE_ARRAY_PRESET_EIGENMIKE32);

    /* Define input mono signal */
    nSH = ORDER2NSH(order);
    inSig = malloc1d(signalLength*sizeof(float));
    rand_m1_1(inSig, signalLength); /* Mono white-noise signal */

    /* Simulate an Eigenmike in a free-field with a single plane-wave */
    f = malloc1d(nBins*sizeof(float));
    kr = malloc1d(nBins*sizeof(double));
    getUniformFreqVector(nFFT, fs, f);
    f[0] = f[1]/4.0f; /* To avoid NaNs at DC */
    radius = 0.042f;
    for(i=0; i<nBins; i++)
        kr[i] = 2.0*SAF_PId*(f[i])*(radius)/343.0f;
    direction_deg[0] = 90.0f;
    direction_deg[1] = 0.0f;
    H_array = (float_complex***)malloc3d(nBins, 32, 1, sizeof(float_complex));
    simulateSphArray(order, kr, kr, nBins, (float*)__Eigenmike32_coords_rad, 32,
                     (float*)direction_deg, 1, ARRAY_CONSTRUCTION_RIGID, 1.0f, FLATTEN3D(H_array));

    /* Inverse FFT to get the time-domain filters */
    tmp_H = malloc1d(nBins*sizeof(float_complex));
    h_array = (float**)malloc2d(32, nFFT, sizeof(float));
    saf_rfft_create(&safFFT, nFFT);
    for(i=0; i<32; i++){
        for(j=0; j<nBins; j++)
            tmp_H[j] = H_array[j][i][0];
        saf_rfft_backward(safFFT, tmp_H, h_array[i]);
    }

    /* Simulate the Eigenmike time-domain signals by convolving the mono signal
     * with each sensor transfer function */
    micSig = (float**)calloc2d(32, signalLength, sizeof(float));
    inSig_32 = (float**)malloc2d(32, signalLength, sizeof(float));
    for(i=0; i<32; i++) /* Replicate inSig for all 32 channels */
        memcpy(inSig_32[i], inSig, signalLength* sizeof(float));
    saf_multiConv_create(&hMC, 256, FLATTEN2D(h_array), nFFT, 32, 0);
    for(i=0; i<(int)((float)signalLength/256.0f); i++){
        saf_multiConv_apply(hMC, FLATTEN2D(inSig_32), FLATTEN2D(micSig));
    }

    /* Encode simulated Eigenmike signals into spherical harmonic signals */
    shSig = (float**)malloc2d(nSH,signalLength,sizeof(float));
    array2sh_process(hA2sh, micSig, shSig, 32, nSH, signalLength);

    /* Clean-up */
    array2sh_destroy(&hA2sh);
    saf_rfft_destroy(&safFFT);
    saf_multiConv_destroy(&hMC);
    free(inSig);
    free(shSig);
    free(inSig_32);
    free(f);
    free(kr);
    free(H_array);
    free(h_array);
    free(tmp_H);
}
#endif /* ENABLE_SAF_EXAMPLES_TESTS */
