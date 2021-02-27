/**
 * @file saf_test.c
 * @brief Unit testing program for the Spatial_Audio_Framework
 *
 * @author Leo McCormack
 * @date 27.04.2020
 */

#include "saf_test.h"
#include "unity.h"   /* unit testing suite */
#include "timer.h"   /* for timing the individual tests */
#include "saf.h"     /* master framework include header */
#include "saf_externals.h" /* to also include saf dependencies (cblas etc.) */

#if SAF_ENABLE_EXAMPLES_TESTS == 1
/* SAF example headers: */
# include "ambi_bin.h"
# include "ambi_dec.h"
# include "ambi_drc.h"
# include "ambi_enc.h"
# include "ambi_roomsim.h"
# include "array2sh.h"
# include "beamformer.h"
# include "binauraliser.h"
# include "decorrelator.h"
# include "dirass.h"
# include "matrixconv.h"
# include "multiconv.h"
# include "panner.h"
# include "pitch_shifter.h"
# include "powermap.h"
# include "rotator.h"
# include "sldoa.h"
#endif /* SAF_ENABLE_EXAMPLES_TESTS */

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
/** A custom Unity RUN_TEST, which calls timerResult() upon exiting each test */
#define RUN_TEST(testfunc)  UNITY_NEW_TEST(#testfunc) \
if (TEST_PROTECT()) {  setUp();  testfunc();  } \
if (TEST_PROTECT() && (!TEST_IS_IGNORED))  {tearDown(); } \
UnityConcludeTest(); timerResult();

/* Main test program */
int main_test(void) {
    printf("%s", SAF_VERSION_BANNER);
    printf("Executing the Spatial_Audio_Framework unit testing program");
#ifdef NDEBUG
    printf(" (Release):\n\n");
#else
    printf(" (Debug):\n\n");
#endif

    /* initialise */
    timer_lib_initialize();
    start = timer_current();
    UNITY_BEGIN();

    /* run each unit test */ 
    RUN_TEST(test__quaternion);
    RUN_TEST(test__dvf_dummyFunc);
    RUN_TEST(test__saf_stft_50pc_overlap);
    RUN_TEST(test__saf_stft_LTI);
    RUN_TEST(test__ims_shoebox_RIR);
    RUN_TEST(test__ims_shoebox_TD);
    RUN_TEST(test__saf_matrixConv);
    RUN_TEST(test__saf_rfft); 
    RUN_TEST(test__afSTFT);
    RUN_TEST(test__qmf);
    RUN_TEST(test__smb_pitchShifter);
    RUN_TEST(test__sortf);
    RUN_TEST(test__sortz);
    RUN_TEST(test__cmplxPairUp);
    RUN_TEST(test__getVoronoiWeights);
    RUN_TEST(test__unique_i);
    RUN_TEST(test__realloc2d_r);
    RUN_TEST(test__latticeDecorrelator);
    RUN_TEST(test__butterCoeffs);
    RUN_TEST(test__faf_IIRFilterbank);
    RUN_TEST(test__gexpm);
#if defined(SAF_ENABLE_SOFA_READER_MODULE)
    RUN_TEST(test__saf_sofa_open);
#endif
#ifdef SAF_ENABLE_TRACKER_MODULE
    RUN_TEST(test__tracker3d);
#endif
    RUN_TEST(test__formulate_M_and_Cr);
    RUN_TEST(test__formulate_M_and_Cr_cmplx);
    RUN_TEST(test__getLoudspeakerDecoderMtx);
    RUN_TEST(test__getSHreal);
    RUN_TEST(test__getSHreal_recur);
    RUN_TEST(test__getSHcomplex);
    RUN_TEST(test__getSHrotMtxReal);
    RUN_TEST(test__real2complexSHMtx);
    RUN_TEST(test__complex2realSHMtx);
    RUN_TEST(test__computeSectorCoeffsEP);
    RUN_TEST(test__checkCondNumberSHTReal);
    RUN_TEST(test__sphMUSIC);
    RUN_TEST(test__sphESPRIT);
    RUN_TEST(test__sphModalCoeffs);
    RUN_TEST(test__truncationEQ);
#if SAF_ENABLE_EXAMPLES_TESTS == 1
    RUN_TEST(test__saf_example_ambi_bin);
    RUN_TEST(test__saf_example_ambi_dec);
    RUN_TEST(test__saf_example_ambi_enc);
    RUN_TEST(test__saf_example_array2sh); 
    RUN_TEST(test__saf_example_rotator);
#endif /* SAF_ENABLE_EXAMPLES_TESTS */

    /* close */
    timer_lib_shutdown();
    printf("\nTotal time elapsed: %lfs", (double)timer_elapsed(start));
    return UNITY_END();
}


/* ========================================================================== */
/*                                 Unit Tests                                 */
/* ========================================================================== */

void test__quaternion(void){
    int i, j;
    float norm;
    float rot[3][3], rot2[3][3], residual[9];
    quaternion_data Q, Q1;

    for(i=0; i<1000; i++){
        /* Randomise the quaternion values */
        rand_m1_1(Q.Q, 4);

        /* Normalise to make it valid */
        norm = L2_norm(Q.Q, 4);
        Q.w /= norm;
        Q.x /= norm;
        Q.y /= norm;
        Q.z /= norm;
        /* Q.w = 0; Q.x = 0.0000563298236; Q.y = 0.947490811; Q.z = -0.319783032; // Problem case! */

        /* Convert to rotation matrix, then back, then to rotation matrix again */
        quaternion2rotationMatrix(&Q, rot);
        rotationMatrix2quaternion(rot, &Q1);
        quaternion2rotationMatrix(&Q1, rot2);

        /* Ensure that the difference between them is 0 */
        utility_svvsub((float*)rot, (float*)rot2, 9, residual);
        for(j=0; j<9; j++)
            TEST_ASSERT_TRUE(fabsf(residual[j])<1e-3f); 
    }
}

/* DVF tests **/
void test__dvf_dummyFunc(void){
    int a, b;
    a = 2;
    b = levelUp(a);
    
    TEST_ASSERT_EQUAL(a+1, b);
}

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
    for(i=0; i<10; i++){
        mov_src_pos[1] = 2.0f + (float)i/10.0f;
        mov_rec_pos[0] = 3.0f + (float)i/10.0f;
        ims_shoebox_updateSource(hIms, sourceID_4, mov_src_pos);
        ims_shoebox_updateReceiver(hIms, receiverID, mov_rec_pos);
        ims_shoebox_computeEchograms(hIms, maxTime_s);
        ims_shoebox_renderRIRs(hIms, 0);
    }

    /* clean-up */
    ims_shoebox_destroy(&hIms);
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
    for(i=0; i<1; i++){
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

void test__saf_matrixConv(void){
    int i, frame;
    float** inputTD, **outputTD, **inputFrameTD, **outputFrameTD;
    float*** filters;
    void* hMatrixConv;

    /* config */
    const int signalLength = 48000;
    const int hostBlockSize = 2048;
    const int filterLength = 512;
    const int nInputs = 32;
    const int nOutputs = 40;

    /* prep */
    inputTD = (float**)malloc2d(nInputs, signalLength, sizeof(float));
    outputTD = (float**)malloc2d(nOutputs, signalLength, sizeof(float));
    inputFrameTD = (float**)malloc2d(nInputs, hostBlockSize, sizeof(float));
    outputFrameTD = (float**)calloc2d(nOutputs, hostBlockSize, sizeof(float));
    filters = (float***)malloc3d(nOutputs, nInputs, filterLength, sizeof(float));
    rand_m1_1(FLATTEN3D(filters), nOutputs*nInputs*filterLength);
    rand_m1_1(FLATTEN2D(inputTD), nInputs*signalLength);
    saf_matrixConv_create(&hMatrixConv, hostBlockSize, FLATTEN3D(filters), filterLength,
                          nInputs, nOutputs, 1);

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

void test__afSTFT(void){
    int frame, nFrames, ch, i, nBands, procDelay, band, nHops;
    void* hSTFT;
    float* freqVector;
    float** insig, **outsig, **inframe, **outframe;
    float_complex*** inspec, ***outspec;

    /* prep */
    const float acceptedTolerance = 0.01f;
    const int fs = 48000;
    const int signalLength = 1*fs;
    const int framesize = 512;
    const int hopsize = 128;
    const int nCHin = 60;
    const int hybridMode = 1;
    const int nCHout = 64;
    insig = (float**)malloc2d(nCHin,signalLength,sizeof(float)); /* One second long */
    outsig = (float**)malloc2d(nCHout,signalLength,sizeof(float));
    inframe = (float**)malloc2d(nCHin,framesize,sizeof(float));
    outframe = (float**)malloc2d(nCHout,framesize,sizeof(float));
    rand_m1_1(FLATTEN2D(insig), nCHin*signalLength); /* populate with random numbers */

    /* Set-up */
    nHops = framesize/hopsize;
    afSTFT_create(&hSTFT, nCHin, nCHout, hopsize, 0, hybridMode, AFSTFT_BANDS_CH_TIME);
    procDelay = afSTFT_getProcDelay(hSTFT);
    nBands = afSTFT_getNBands(hSTFT);
    freqVector = malloc1d(nBands*sizeof(float));
    afSTFT_getCentreFreqs(hSTFT, (float)fs, nBands, freqVector);
    inspec = (float_complex***)malloc3d(nBands, nCHin, nHops, sizeof(float_complex));
    outspec = (float_complex***)malloc3d(nBands, nCHout, nHops, sizeof(float_complex));

    /* just some messing around... */
    afSTFT_channelChange(hSTFT, 100, 5);
    afSTFT_clearBuffers(hSTFT);
    afSTFT_channelChange(hSTFT, 39, 81);
    afSTFT_channelChange(hSTFT, nCHin, nCHout); /* back to original config */
    afSTFT_clearBuffers(hSTFT);

    /* Pass insig through the QMF filterbank, block-wise processing */
    nFrames = (int)((float)signalLength/(float)framesize);
    for(frame = 0; frame<nFrames; frame++){
        /* Forward transform */
        for(ch=0; ch<nCHin; ch++)
            memcpy(inframe[ch], &insig[ch][frame*framesize], framesize*sizeof(float));
        afSTFT_forward(hSTFT, inframe, framesize, inspec);

        /* Copy first channel of inspec to all outspec channels */
        for(band=0; band<nBands; band++)
            for(ch=0; ch<nCHout; ch++)
                memcpy(outspec[band][ch], inspec[band][0], nHops*sizeof(float_complex));

        /* Backwards transform */
        afSTFT_backward(hSTFT, outspec, framesize, outframe);
        for(ch=0; ch<nCHout; ch++)
            memcpy(&outsig[ch][frame*framesize], outframe[ch], framesize*sizeof(float));
    }

    /* Check that input==output (given some numerical precision) - channel 0 */
    for(i=0; i<signalLength-procDelay-framesize; i++)
        TEST_ASSERT_TRUE( fabsf(insig[0][i] - outsig[0][i+procDelay]) <= acceptedTolerance );

    /* Clean-up */
    afSTFT_destroy(&hSTFT);
    free(insig);
    free(outsig);
    free(inframe);
    free(outframe);
    free(inspec);
    free(outspec);
    free(freqVector);
}

void test__qmf(void){
    int frame, nFrames, ch, i, nBands, procDelay, band, nHops;
    void* hQMF;
    float* freqVector;
    float** insig, **outsig, **inframe, **outframe;
    float_complex*** inspec, ***outspec;

    /* prep */
    const float acceptedTolerance = 0.01f;
    const int fs = 48000;
    const int signalLength = 1*fs;
    const int framesize = 512;
    const int hopsize = 128;
    const int nCHin = 60;
    const int hybridMode = 1;
    const int nCHout = 64;
    insig = (float**)malloc2d(nCHin,signalLength,sizeof(float)); /* One second long */
    outsig = (float**)malloc2d(nCHout,signalLength,sizeof(float));
    inframe = (float**)malloc2d(nCHin,framesize,sizeof(float));
    outframe = (float**)malloc2d(nCHout,framesize,sizeof(float));
    rand_m1_1(FLATTEN2D(insig), nCHin*signalLength); /* populate with random numbers */

    /* Set-up */
    nHops = framesize/hopsize;
    qmf_create(&hQMF, nCHin, nCHout, hopsize, hybridMode, QMF_BANDS_CH_TIME);
    procDelay = qmf_getProcDelay(hQMF);
    nBands = qmf_getNBands(hQMF);
    freqVector = malloc1d(nBands*sizeof(float));
    qmf_getCentreFreqs(hQMF, (float)fs, nBands, freqVector);
    inspec = (float_complex***)malloc3d(nBands, nCHin, nHops, sizeof(float_complex));
    outspec = (float_complex***)malloc3d(nBands, nCHout, nHops, sizeof(float_complex));

    /* Pass insig through the QMF filterbank, block-wise processing */
    nFrames = (int)((float)signalLength/(float)framesize);
    for(frame = 0; frame<nFrames; frame++){
        /* QMF Analysis */
        for(ch=0; ch<nCHin; ch++)
            memcpy(inframe[ch], &insig[ch][frame*framesize], framesize*sizeof(float));
        qmf_analysis(hQMF, inframe, framesize, inspec);

        /* Copy first channel of inspec to all outspec channels */
        for(band=0; band<nBands; band++)
            for(ch=0; ch<nCHout; ch++)
                memcpy(outspec[band][ch], inspec[band][0], nHops*sizeof(float_complex));

        /* QMF Synthesis */
        qmf_synthesis(hQMF, outspec, framesize, outframe);
        for(ch=0; ch<nCHout; ch++)
            memcpy(&outsig[ch][frame*framesize], outframe[ch], framesize*sizeof(float));
    }

    /* Check that input==output (given some numerical precision) - channel 0 */
    for(i=0; i<signalLength-procDelay-framesize; i++)
        TEST_ASSERT_TRUE( fabsf(insig[0][i] - outsig[0][i+procDelay]) <= acceptedTolerance );

    /* Clean-up */
    qmf_destroy(&hQMF);
    free(insig);
    free(outsig);
    free(inframe);
    free(outframe);
    free(inspec);
    free(outspec);
    free(freqVector);
}

void test__smb_pitchShifter(void){
    float* inputData, *outputData;
    void* hPS, *hFFT;
    float frequency;
    int i, smbLatency, ind;

    /* Config */
    const int sampleRate = 48000;
    const int FFTsize = 8192;
    const int osfactor = 4;
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

void test__getVoronoiWeights(void){
    int i, it, td, nDirs;
    float* dirs_deg, *weights;
    float sum, tmp, scale;

    /* Config */
    const float acceptedTolerance = 0.01f;
    const int nIterations = 100;

    /* Loop over T-designs */
    for(td=2; td<21; td++){
        dirs_deg = (float*)__HANDLES_Tdesign_dirs_deg[td];
        nDirs = __Tdesign_nPoints_per_degree[td];

        /* Compute weights */
        weights = malloc1d(nDirs*sizeof(float));
        getVoronoiWeights(dirs_deg, nDirs, 0, weights);

        /* Assert that they sum to 4PI */
        sum = sumf(weights, nDirs);
        TEST_ASSERT_FLOAT_WITHIN(acceptedTolerance, 4.0f*SAF_PI, sum);

        /* Due to the uniform arrangement, all the weights should be the same */
        for(i=1; i<nDirs; i++)
            TEST_ASSERT_FLOAT_WITHIN(acceptedTolerance, weights[0], weights[i]);

        /* clean-up */
        free(weights);
    }

    /* Loop over some random arrangement of points */
    for(it=0; it<nIterations; it++){
        rand_0_1(&tmp, 1);
        nDirs = (int)(tmp*190.0f + 10.0f); /* random number between 10..200 */

        /* Random dirs (-180..180 azi, -180..180 elev) */
        dirs_deg = malloc1d(nDirs*2*sizeof(float));
        rand_m1_1(dirs_deg, nDirs*2);
        scale = 180.0f;
        utility_svsmul(dirs_deg, &scale, nDirs*2, dirs_deg);

        /* Compute weights */
        weights = malloc1d(nDirs*sizeof(float));
        getVoronoiWeights(dirs_deg, nDirs, 0, weights);

        /* Assert that they sum to 4PI */
        sum = sumf(weights, nDirs);
        TEST_ASSERT_FLOAT_WITHIN(acceptedTolerance, 4.0f*SAF_PI, sum);

        /* clean-up */
        free(dirs_deg);
        free(weights);
    }
}

void test__unique_i(void){
    int i, nUnique;
    int* uniqueVals;
    int* uniqueInds;

    /* test1 */
    int input[6] = {1, 2, 2, 10, 11, 12};
    int uniqueVals_ref[5] = {1, 2, 10, 11, 12};
    int uniqueInds_ref[5] = {0, 2, 3, 4, 5};
    unique_i(input, 6, &uniqueVals, &uniqueInds, &nUnique);
    TEST_ASSERT_EQUAL(5, nUnique);
    for(i=0; i<nUnique; i++){
        TEST_ASSERT_EQUAL(uniqueVals_ref[i], uniqueVals[i]);
        TEST_ASSERT_EQUAL(uniqueInds_ref[i], uniqueInds[i]);
    }
    free(uniqueVals);
    free(uniqueInds);

    /* test2 */
    int input2[12] = {1, 10, 1, 3, 1, 3, 4, 7, 8, 10, 10, 2};
    int uniqueVals_ref2[7] = {1, 3, 4, 7, 8, 10, 2};
    int uniqueInds_ref2[7] = {4, 5, 6, 7, 8, 10, 11};
    unique_i(input2, 12, &uniqueVals, &uniqueInds, &nUnique);
    TEST_ASSERT_EQUAL(7, nUnique);
    for(i=0; i<nUnique; i++){
        TEST_ASSERT_EQUAL(uniqueVals_ref2[i], uniqueVals[i]);
        TEST_ASSERT_EQUAL(uniqueInds_ref2[i], uniqueInds[i]);
    }
    free(uniqueVals);
    free(uniqueInds);
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

void test__latticeDecorrelator(void){
    int c, band, nBands, idx, hopIdx, i;
    void* hDecor, *hDecor2, *hSTFT;
    float icc, tmp, tmp2;
    float* freqVector;
    float** inputTimeDomainData, **outputTimeDomainData, **tempHop;
    float_complex*** inTFframe, ***outTFframe;

    /* config */
    const float acceptedICC = 0.02f;
    const int nCH = 24;
    const int nTestHops = 2000;
    const int hopSize = 128;
    const int afSTFTdelay = hopSize*12;
    const int lSig = nTestHops*hopSize+afSTFTdelay;
    const float fs = 48e3f;
    nBands = hopSize+5;

    /* audio buffers */
    inputTimeDomainData = (float**) malloc2d(1, lSig, sizeof(float));
    outputTimeDomainData = (float**) malloc2d(nCH, lSig, sizeof(float));
    inTFframe = (float_complex***)malloc3d(nBands, nCH, 1, sizeof(float_complex));
    outTFframe = (float_complex***)malloc3d(nBands, nCH, 1, sizeof(float_complex));
    tempHop = (float**) malloc2d(nCH, hopSize, sizeof(float));

    /* Initialise afSTFT and input data */
    afSTFT_create(&hSTFT, 1, nCH, hopSize, 0, 1, AFSTFT_BANDS_CH_TIME);
    rand_m1_1(FLATTEN2D(inputTimeDomainData), 1*lSig); /* populate with random numbers */
    freqVector = malloc1d(nBands*sizeof(float));
    afSTFT_getCentreFreqs(hSTFT, fs, nBands, freqVector);

    /* setup decorrelator 1 */
    int orders[5] = {6, 6, 6, 3, 2};
    float freqCutoffs[5] = {700.0f, 2.4e3f, 4e3f, 12e3f, 20e3f};
    int fixedDelays[6] = {8, 8, 7, 2, 1, 2};
    latticeDecorrelator_create(&hDecor, nCH, orders, freqCutoffs, fixedDelays, 5, freqVector, 0, 133);

    /* setup decorrelator 2 */
    float freqCutoffs2[3] = {700.0f, 2.4e3f, 4e3f};
    int orders2[3] = {2, 3, 2};
    int fixedDelays2[4] = {2, 2, 1, 0};
    latticeDecorrelator_create(&hDecor2, nCH, orders2, freqCutoffs2, fixedDelays2, 3, freqVector, nCH, 133);

    /* Pass input data through afSTFT */
    idx = 0;
    hopIdx = 0;
    while(idx<lSig){
        for(c=0; c<1; c++)
            memcpy(tempHop[c], &(inputTimeDomainData[c][hopIdx*hopSize]), hopSize*sizeof(float));

        /* forward TF transform, and replicate to all channels */
        afSTFT_forward(hSTFT, tempHop, hopSize, inTFframe);
        for(band=0; band<nBands; band++)
            for(i=1; i<nCH;i++)
                inTFframe[band][i][0] = inTFframe[band][0][0];

        /* decorrelate */
        latticeDecorrelator_apply(hDecor, inTFframe, 1, outTFframe);
        latticeDecorrelator_apply(hDecor2, outTFframe, 1, outTFframe);

        /*  backward TF transform */
        afSTFT_backward(hSTFT, outTFframe, hopSize, tempHop);

        /* Copy frame to output TD buffer */
        for(c=0; c<nCH; c++)
            memcpy(&(outputTimeDomainData[c][hopIdx*hopSize]), tempHop[c], hopSize*sizeof(float));
        idx+=hopSize;
        hopIdx++;
    }

    /* Compensate for afSTFT delay, and check that the inter-channel correlation
     * coefficient is below the accepted threshold (ideally 0, if fully
     * decorrelated...) */
    for(c=0; c<nCH; c++){
        utility_svvdot(inputTimeDomainData[0], &outputTimeDomainData[c][afSTFTdelay], (lSig-afSTFTdelay), &icc);
        utility_svvdot(inputTimeDomainData[0], inputTimeDomainData[0], (lSig-afSTFTdelay), &tmp);
        utility_svvdot(&outputTimeDomainData[c][afSTFTdelay], &outputTimeDomainData[c][afSTFTdelay], (lSig-afSTFTdelay), &tmp2);

        icc = icc/sqrtf(tmp*tmp2); /* normalise */
        TEST_ASSERT_TRUE(fabsf(icc)<acceptedICC);
    }
#if 0
    /* Check for mutually decorrelated channels... */
    int c2;
    for(c=0; c<nCH; c++){
        for(c2=0; c2<nCH; c2++){
            utility_svvdot(&outputTimeDomainData[c][afSTFTdelay], &outputTimeDomainData[c2][afSTFTdelay], (lSig-afSTFTdelay), &icc);
            utility_svvdot(&outputTimeDomainData[c2][afSTFTdelay], &outputTimeDomainData[c2][afSTFTdelay], (lSig-afSTFTdelay), &tmp);
            utility_svvdot(&outputTimeDomainData[c][afSTFTdelay], &outputTimeDomainData[c][afSTFTdelay], (lSig-afSTFTdelay), &tmp2);

            icc = icc/sqrtf(tmp*tmp2); /* normalise */
           // TEST_ASSERT_TRUE(fabsf(icc)<acceptedICC);
        }
    }
#endif

    /* Clean-up */
    latticeDecorrelator_destroy(&hDecor);
    latticeDecorrelator_destroy(&hDecor2);
    free(inTFframe);
    free(outTFframe);
    free(tempHop);
    afSTFT_destroy(&hSTFT);
    free(inputTimeDomainData);
    free(outputTimeDomainData);
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

void test__gexpm(void){
    int i, j;
    float outM[6][6];

    /* prep */
    const float acceptedTolerance = 0.0001f;
    const float inM[6][6] = {
        {-0.376858200853762,0.656790634216694,0.124479178614046,-0.334752428307223,1.50745241578235,0.0290651989052969},
        {0.608382058262806,0.581930485432986,3.23135406998058,-0.712003744668929,-1.33872571354702,-0.334742482743222},
        {-0.795741418256672,0.690709474622409,0.620971281129248,1.38749471231620,0.897245329198841,-0.0693670166113321},
        {0.179789913109994,-1.06135084902804,-1.10032635271188,0.612441344250358,-2.43213807790664,-0.479265889956047},
        {-0.277441781278754,-0.0732116130293688,-0.572551795688137,1.02024767389969,0.167385894565923,1.45210312619277},
        {-0.205305770089918,-1.59783032780633,1.08539265129120,0.460057585947626,-1.02420974042838,1.04117461500218}
    };
    const float outM_ref[6][6] = {
        {0.385163650730121,0.0865151585709784,0.898406722231524,0.877640791713973,0.435244824708340,0.888866982998854},
        {-0.664938511314777,5.02943129352875,8.24444951891833,2.23840978101979,-0.942669833528886,-2.38535530623266},
        {-0.388189314743059,0.429308537172675,1.13870842882926,1.60875776611798,-1.44249911796405,-1.51822150286392},
        {1.05630187656688,0.256606570814868,-2.42701873560847,-1.42372526577009,-0.335273289873574,-1.94362909671742},
        {0.0261470437116839,-3.03329326250434,-3.50207776203591,0.412043775125377,-0.536000387729306,1.61801775548557},
        {-0.292024827617294,-4.31537192033477,-3.99160103133879,0.312499067924889,-1.46924802440347,1.98522802303672}
    };

    /* Compute matrix exponential */
    gexpm((float*)inM, 6, 0, (float*)outM);

    /* Check that output of SAF's gexpm, is similar to Matlab's expm: */
    for(i=0; i<6; i++)
        for(j=0; j<6; j++)
            TEST_ASSERT_TRUE( fabsf(outM[i][j] - outM_ref[i][j]) <= acceptedTolerance );
}

#if defined(SAF_ENABLE_SOFA_READER_MODULE)
void test__saf_sofa_open(void){
    SAF_SOFA_ERROR_CODES error;
    saf_sofa_container sofa;
    error = saf_sofa_open(&sofa, "/Users/mccorml1/Documents/FABIAN_HRTF_DATABASE_V1/1 HRIRs/SOFA/FABIAN_HRIR_measured_HATO_20.sofa");
    saf_sofa_close(&sofa);
}
#endif

#ifdef SAF_ENABLE_TRACKER_MODULE
void test__tracker3d(void){
    int hop, i, j, k, nSH, nGrid, rand_idx, dropouts;
    int inds[2];
    int* target_IDs;
    void* hT3d, *hMUSIC;
    float measNoiseSD_deg, noiseSpecDen_deg, scale, rand01;
    float est_dirs_deg[2][2], est_dirs_xyz[2][3];
    float* grid_dirs_deg;
    float** insigs, **inputSH, **inputSH_noise, **inputSH_hop, **Y, **Cx, **V, **Vn;
    float_complex** Vn_cmplx;
    int nTargets;
    float *target_dirs_xyz;

    /* Test configuration */
    const float acceptedTolerance = 0.001f;
    const int order = 2;
    const float fs = 48e3;
    const int hopsize = 128;
    const float sigLen = fs*5;
    const int nSources = 2; /* cannot be changed, hard-coded for 2 */
    const float src_dirs_deg[2][2] = { {-35.0f, 30.0f}, {120.0f, 0.0f} };

    /* Configure the tracker */
    tracker3d_config tpars;
    /* Number of Monte Carlo samples/particles. The more complex the
     * distribution is, the more particles required (but also, the more
     * computationally expensive the tracker becomes). */
    tpars.Np = 20;
    tpars.ARE_UNIT_VECTORS = 1;
    tpars.maxNactiveTargets = 4; /* about 2 higher than expected is good */
    /* Likelihood of an estimate being noise/clutter */
    tpars.noiseLikelihood = 0.2f; /* between [0..1] */
    /* Measurement noise - e.g. to assume that estimates within the range +/-20
     * degrees belong to the same target, set SDmnoise_deg = 20 */
    measNoiseSD_deg = 20.0f;
    tpars.measNoiseSD = 1.0f-cosf(measNoiseSD_deg*SAF_PI/180.0f); /* Measurement noise standard deviation */
    /* Noise spectral density - not fully understood. But it influences the
     * smoothness of the target tracks */
    noiseSpecDen_deg = 1.0f;
    tpars.noiseSpecDen = 1.0f-cosf(noiseSpecDen_deg*SAF_PI/180.0f);  /* Noise spectral density */
    /* FLAG - whether to allow for multiple target deaths in the same tracker
     * prediction step */
    tpars.ALLOW_MULTI_DEATH = 1;
    /* Probability of birth and death */
    tpars.init_birth = 0.5f; /* value between [0 1] - Prior probability of birth */
    tpars.alpha_death = 20.0f; /* always >= 1; 1 is good */      /* 20-> means death is very unlikely... */
    tpars.beta_death = 1.0f; /* always >= 1; 1 is good */
    /* Elapsed time (in seconds) between observations */
    tpars.dt = 1.0f/(fs/(float)hopsize); /* Hop length of frames */
    /* Whether or not to allow multiple active sources for each update */
    /* Real-time tracking is based on the particle with highest weight. A
     * one-pole averaging filter is used to smooth the weights over time. */
    tpars.W_avg_coeff = 0.5f;
    /* Force kill targets that are close to another target. In these cases, the
     * target that has been 'alive' for the least amount of time, is killed */
    tpars.FORCE_KILL_TARGETS = 1;
    tpars.forceKillDistance = 0.2f;
    /* Mean position priors x,y,z (assuming directly in-front) */
    tpars.M0[0] = 1.0f; tpars.M0[1] = 0.0f; tpars.M0[2] = 0.0f;
    /* Mean Velocity priors x,y,z (assuming stationary) */
    tpars.M0[3] = 0.0f; tpars.M0[4] = 0.0f; tpars.M0[5] = 0.0f;
    /* Target velocity - e.g. to assume that a target can move 20 degrees in
    * two seconds along the horizontal, set V_azi = 20/2 */
    const float Vazi_deg = 3.0f;  /* Velocity of target on azimuthal plane */
    const float Vele_deg = 3.0f;  /* Velocity of target on median plane */
    memset(tpars.P0, 0, 6*6*sizeof(float));
    /* Variance PRIORs of estimates along the x,y,z axes, respectively. Assuming
     * coordinates will lay on the unit sphere +/- x,y,z, so a range of 2, and
     * hence a variance of 2^2: */
    tpars.P0[0][0] = 4.0f; tpars.P0[1][1] = 4.0f; tpars.P0[2][2] = 4.0f;
    /* Velocity PRIORs of estimates the x,y,z axes */
    tpars.P0[3][3] = 1.0f-cosf(Vazi_deg*SAF_PI/180.0f); /* x */
    tpars.P0[4][4] = tpars.P0[3][3];                    /* y */
    tpars.P0[5][5] = 1.0f-cosf(Vele_deg*SAF_PI/180.0f); /* z */
    /* PRIOR probabilities of noise. (Assuming the noise is uniformly
     * distributed in the entire spatial grid). */
    tpars.cd = 1.0f/(4.0f*SAF_PI);

    /* Create tracker */
    tracker3d_create(&hT3d, tpars);

    /* Create spherical harmonic input signals */
    insigs = (float**)malloc2d(nSources, sigLen, sizeof(float));
    rand_m1_1(FLATTEN2D(insigs), nSources*sigLen);
    nSH = ORDER2NSH(order);
    Y = (float**)malloc2d(nSH, nSources, sizeof(float));
    getRSH(order, (float*)src_dirs_deg, nSources, FLATTEN2D(Y));
    scale = 1.0f/(float)nSources;
    utility_svsmul(FLATTEN2D(Y), &scale, nSH*nSources, NULL);
    inputSH = (float**)malloc2d(nSH, sigLen, sizeof(float));
    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nSH, sigLen, nSources, 1.0f,
                FLATTEN2D(Y), nSources,
                FLATTEN2D(insigs), sigLen, 0.0f,
                FLATTEN2D(inputSH), sigLen);

    /* Add some noise */
    inputSH_noise = (float**)malloc2d(nSH, sigLen, sizeof(float));
    rand_m1_1(FLATTEN2D(inputSH_noise), nSH*sigLen);
    scale = 0.05f;
    utility_svsmul(FLATTEN2D(inputSH_noise), &scale, nSH*sigLen, NULL);
    utility_svvadd(FLATTEN2D(inputSH), FLATTEN2D(inputSH_noise), nSH*sigLen, FLATTEN2D(inputSH));

    /* Create DoA estimator */
    nGrid = 240; /* number of points in a t-design of degree 21 */
    grid_dirs_deg = (float*)__Tdesign_degree_21_dirs_deg;
    sphMUSIC_create(&hMUSIC, order, grid_dirs_deg, nGrid);

    /* Memory allocations */
    inputSH_hop = (float**)malloc2d(nSH, hopsize, sizeof(float));
    Cx = (float**)malloc2d(nSH, nSH, sizeof(float));
    V = (float**)malloc2d(nSH, nSH, sizeof(float));
    Vn = (float**)malloc2d(nSH, (nSH-nSources), sizeof(float)); /* noise subspace */
    Vn_cmplx = (float_complex**)malloc2d(nSH, (nSH-nSources), sizeof(float_complex));
    target_dirs_xyz = NULL;
    target_IDs = NULL;

    /* Loop over hops */
    dropouts = 0;
    for(hop=0; hop<(int)((float)sigLen/(float)hopsize); hop++){
        /* Grab current hop */
        for(i=0; i<nSH; i++)
            memcpy(inputSH_hop[i], &inputSH[i][hop*hopsize], hopsize*sizeof(float));

        /* Eigenvalue decomposition and truncation of eigen vectors to obtain
         * noise subspace (based on source number) */
        cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasTrans, nSH, nSH, hopsize, 1.0f,
                    FLATTEN2D(inputSH_hop), hopsize,
                    FLATTEN2D(inputSH_hop), hopsize, 0.0f,
                    FLATTEN2D(Cx), nSH);
        utility_sseig(FLATTEN2D(Cx), nSH, 1, FLATTEN2D(V), NULL, NULL);
        for(i=0; i<nSH; i++)
            for(j=0, k=nSources; j<nSH-nSources; j++, k++)
                Vn[i][j] = V[i][k];
        for(i=0; i<nSH; i++)
            for(j=0; j<nSH-nSources; j++)
                Vn_cmplx[i][j] = cmplxf(Vn[i][j], 0.0f);

        /* DoA estimation */
        sphMUSIC_compute(hMUSIC, FLATTEN2D(Vn_cmplx), nSources, NULL, (int*)inds);
        est_dirs_deg[0][0] = grid_dirs_deg[inds[0]*2+0];
        est_dirs_deg[0][1] = grid_dirs_deg[inds[0]*2+1];
        est_dirs_deg[1][0] = grid_dirs_deg[inds[1]*2+0];
        est_dirs_deg[1][1] = grid_dirs_deg[inds[1]*2+1];
        unitSph2cart((float*)est_dirs_deg, nSources, 1, (float*)est_dirs_xyz);

        /* Pick an estimate at random */
        rand_0_1(&rand01, 1);
        rand_idx = (int)(rand01*(float)nSources);

        /* Feed tracker */
        tracker3d_step(hT3d, (float*)&est_dirs_xyz[rand_idx], 1, &target_dirs_xyz, &target_IDs, &nTargets);

        /* Give the tracker a couple of steps, and then assert that it is keeping track of these two targets */
        if(hop>10){
            TEST_ASSERT_TRUE( nTargets <= nSources );
            if(nTargets==nSources){
                TEST_ASSERT_TRUE( fabsf(est_dirs_xyz[0][0] - target_dirs_xyz[0*3+0]) <= acceptedTolerance ||
                                  fabsf(est_dirs_xyz[0][0] - target_dirs_xyz[1*3+0]) <= acceptedTolerance);
                TEST_ASSERT_TRUE( fabsf(est_dirs_xyz[0][1] - target_dirs_xyz[0*3+1]) <= acceptedTolerance ||
                                  fabsf(est_dirs_xyz[0][1] - target_dirs_xyz[1*3+1]) <= acceptedTolerance);
                TEST_ASSERT_TRUE( fabsf(est_dirs_xyz[0][2] - target_dirs_xyz[0*3+2]) <= acceptedTolerance ||
                                  fabsf(est_dirs_xyz[0][2] - target_dirs_xyz[1*3+2]) <= acceptedTolerance);
                TEST_ASSERT_TRUE( fabsf(est_dirs_xyz[1][0] - target_dirs_xyz[0*3+0]) <= acceptedTolerance ||
                                  fabsf(est_dirs_xyz[1][0] - target_dirs_xyz[1*3+0]) <= acceptedTolerance);
                TEST_ASSERT_TRUE( fabsf(est_dirs_xyz[1][1] - target_dirs_xyz[0*3+1]) <= acceptedTolerance ||
                                  fabsf(est_dirs_xyz[1][1] - target_dirs_xyz[1*3+1]) <= acceptedTolerance);
                TEST_ASSERT_TRUE( fabsf(est_dirs_xyz[1][2] - target_dirs_xyz[0*3+2]) <= acceptedTolerance ||
                                  fabsf(est_dirs_xyz[1][2] - target_dirs_xyz[1*3+2]) <= acceptedTolerance);
            }
            else
                dropouts++; /* Should be very unlikely, (as the probably of death set to be so low), but it can still happen... */
        }
    }
    TEST_ASSERT_TRUE(dropouts<5);

    /* Clean-up */
    tracker3d_destroy(&hT3d);
    sphMUSIC_destroy(&hMUSIC);
    free(target_dirs_xyz);
    free(target_IDs);
    free(insigs);
    free(inputSH);
    free(inputSH_noise);
    free(inputSH_hop);
    free(Y);
    free(Cx);
    free(V);
    free(Vn);
    free(Vn_cmplx);
}
#endif

void test__formulate_M_and_Cr(void){
    int i, j, it, nCHin, nCHout, lenSig;
    float reg, tmp;
    float** Q, **x, **y, **z, **Cx, **Cy, **Cz, **M, **Cr, **Mr, **Q_Cx, **Cp;
    float** decor, **z_r, **eye_nCHout;
    void* hCdf, *hCdf_res;

    /* Config */
    const float acceptedTolerance = 0.1f; /* Due to regularisation, the result will never be exact */
    /* However, this is a very generous tolerance value. If the number of input
     * and output channels are similar, then this tolerance can be much lower
     * (0.00001). The error is only ever high when there is a large discrepency
     * between the number of input and output channels. */
    const int nIterations = 1000;

    /* Loop through iterations */
    for(it=0; it<nIterations; it++){
        rand_0_1(&tmp, 1);
        nCHin = (int)(tmp*12.0f + 4.1f); /* random number between 4 and 16 */
        rand_0_1(&tmp, 1);
        nCHout = (int)(tmp*12.0f + 4.1f); /* random number between 4 and 16 */
        rand_0_1(&tmp, 1);
        lenSig = (int)(tmp*384.0f + 128.1f); /* random number between 128 and 512 */

        /* Define prototype decoder and compute input signal covariance matrix */
        Q = (float**)calloc2d(nCHout, nCHin, sizeof(float));
        for(i=0; i<MIN(nCHin, nCHout); i++)
            Q[i][i] = 1.0f; /* Identity */
        x = (float**)malloc2d(nCHin, lenSig, sizeof(float));
        rand_m1_1(FLATTEN2D(x), nCHin*lenSig);
        Cx = (float**)malloc2d(nCHin, nCHin, sizeof(float));
        cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasTrans, nCHin, nCHin, lenSig, 1.0f,
                    FLATTEN2D(x), lenSig,
                    FLATTEN2D(x), lenSig, 0.0f,
                    FLATTEN2D(Cx), nCHin);

        /* Compute target covariance matrix */
        y = (float**)malloc2d(nCHout, lenSig, sizeof(float));
        rand_m1_1(FLATTEN2D(y), nCHout*lenSig);
        Cy = (float**)malloc2d(nCHout, nCHout, sizeof(float));
        cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasTrans, nCHout, nCHout, lenSig, 1.0f,
                    FLATTEN2D(y), lenSig,
                    FLATTEN2D(y), lenSig, 0.0f,
                    FLATTEN2D(Cy), nCHout);

        /* Compute optimal mixing matrix - with energy compensation enabled */
        M = (float**)malloc2d(nCHout, nCHin, sizeof(float));
        reg = 0.2f;
        cdf4sap_create(&hCdf, nCHin, nCHout);
        formulate_M_and_Cr(hCdf, FLATTEN2D(Cx), FLATTEN2D(Cy), FLATTEN2D(Q), 1, reg, FLATTEN2D(M), NULL);

        /* Apply mixing matrix to 'x' and assert that it's covariance matrix matches
         * the target covariance matrix */
        z = (float**)malloc2d(nCHout, lenSig, sizeof(float));
        cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nCHout, lenSig, nCHin, 1.0f,
                    FLATTEN2D(M), nCHin,
                    FLATTEN2D(x), lenSig, 0.0f,
                    FLATTEN2D(z), lenSig);
        Cz = (float**)malloc2d(nCHout, nCHout, sizeof(float));
        cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasTrans, nCHout, nCHout, lenSig, 1.0f,
                    FLATTEN2D(z), lenSig,
                    FLATTEN2D(z), lenSig, 0.0f,
                    FLATTEN2D(Cz), nCHout);
        if(nCHin>=nCHout){
            for(i=0; i<nCHout; i++)
                for(j=0; j<nCHout; j++)
                    TEST_ASSERT_FLOAT_WITHIN(acceptedTolerance, Cy[i][j], Cz[i][j]);
        }
        else{ /* if nCHin<nCHout, then only the diagonal elements will match */
            for(i=0; i<nCHout; i++)
                TEST_ASSERT_FLOAT_WITHIN(acceptedTolerance, Cy[i][i], Cz[i][i]);
        }

        /* Determine prototype covariance matrix */
        Q_Cx = (float**)malloc2d(nCHout, nCHin, sizeof(float));
        cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nCHout, nCHin, nCHin, 1.0f,
                    FLATTEN2D(Q), nCHin,
                    FLATTEN2D(Cx), nCHin, 0.0f,
                    FLATTEN2D(Q_Cx), nCHin);
        Cp = (float**)malloc2d(nCHout, nCHout, sizeof(float));
        cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasTrans, nCHout, nCHout, nCHin, 1.0f,
                    FLATTEN2D(Q_Cx), nCHin,
                    FLATTEN2D(Q), nCHin, 0.0f,
                    FLATTEN2D(Cp), nCHout);
        for(i=0; i<nCHout; i++)
            for(j=0; j<nCHout; j++)
                if(i!=j)
                    Cp[i][j] = 0.0f; /* Zero non-diagonal elements */

        /* Create perfectly incoherent frame. Note, in practice this would instead
         * be a decorrelated version of the prototype signals, [i.e.
         * decorrelate(Q*x) ]*/
        decor = (float**)malloc2d(nCHout, lenSig, sizeof(float));
        rand_m1_1(FLATTEN2D(decor), nCHout*lenSig);

        /* Now compute optimal mixing matrix, but this time also including the
         * residual mixing matrix */
        M = (float**)malloc2d(nCHout, nCHin, sizeof(float));
        reg = 0.2f;
        Cr = (float**)malloc2d(nCHout, nCHout, sizeof(float));
        formulate_M_and_Cr(hCdf, FLATTEN2D(Cx), FLATTEN2D(Cy), FLATTEN2D(Q), 0, reg, FLATTEN2D(M), FLATTEN2D(Cr));
        cdf4sap_create(&hCdf_res, nCHout, nCHout);
        Mr = (float**)calloc2d(nCHout, nCHout, sizeof(float));
        eye_nCHout = (float**)calloc2d(nCHout, nCHout, sizeof(float));
        for(i=0; i<nCHout; i++)
            eye_nCHout[i][i] = 1.0f;
        formulate_M_and_Cr(hCdf_res, FLATTEN2D(Cp), FLATTEN2D(Cr), FLATTEN2D(eye_nCHout), 0, reg, FLATTEN2D(Mr), NULL);

        /* Apply mixing matrix to x, and residual mixing matrix to the decorrelated
         * prototype signals, and sum */
        cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nCHout, lenSig, nCHin, 1.0f,
                    FLATTEN2D(M), nCHin,
                    FLATTEN2D(x), lenSig, 0.0f,
                    FLATTEN2D(z), lenSig);
        z_r = (float**)malloc2d(nCHout, lenSig, sizeof(float));
        cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nCHout, lenSig, nCHout, 1.0f,
                    FLATTEN2D(Mr), nCHout,
                    FLATTEN2D(decor), lenSig, 0.0f,
                    FLATTEN2D(z_r), lenSig);
        utility_svvadd(FLATTEN2D(z), FLATTEN2D(z_r), nCHout*lenSig, FLATTEN2D(z));

        /* Assert that the covariance matrix of 'z' matches the target covariance
         * matrix */
        if(nCHin>=nCHout){
            for(i=0; i<nCHout; i++)
                for(j=0; j<nCHout; j++)
                    TEST_ASSERT_FLOAT_WITHIN(acceptedTolerance, Cy[i][j], Cz[i][j]);
        }
        else{ /* if nCHin<nCHout, then only the diagonal elements will match */
            for(i=0; i<nCHout; i++)
                TEST_ASSERT_FLOAT_WITHIN(acceptedTolerance, Cy[i][i], Cz[i][i]);
        }

        /* Clean-up */
        cdf4sap_destroy(&hCdf);
        cdf4sap_destroy(&hCdf_res);
        free(Q);
        free(x);
        free(y);
        free(z);
        free(Cx);
        free(Cy);
        free(Cz);
        free(M);
        free(Cr);
        free(Mr);
        free(Q_Cx);
        free(Cp);
        free(decor);
        free(z_r);
        free(eye_nCHout);
    }
}

void test__formulate_M_and_Cr_cmplx(void){
    int i, j, it, nCHin, nCHout, lenSig;
    float reg, tmp;
    float_complex** Q, **x, **y, **z, **Cx, **Cy, **Cz, **M, **Cr, **Mr, **Q_Cx, **Cp;
    float_complex** decor, **z_r, **eye_nCHout;
    void* hCdf, *hCdf_res;
    const float_complex calpha = cmplxf(1.0f, 0.0f), cbeta = cmplxf(0.0f, 0.0f);

    /* Config */
    const float acceptedTolerance = 0.1f; /* Due to regularisation, the result will never be exact */
    /* However, this is a very generous tolerance value. If the number of input
     * and output channels are similar, then this tolerance can be much lower
     * (0.00001). The error is only ever high when there is a large discrepency
     * between the number of input and output channels. */
    const int nIterations = 300;

    /* Loop through iterations */
    for(it=0; it<nIterations; it++){
        rand_0_1(&tmp, 1);
        nCHin = (int)(tmp*12.0f + 4.1f); /* random number between 4 and 16 */
        rand_0_1(&tmp, 1);
        nCHout = (int)(tmp*12.0f + 4.1f); /* random number between 4 and 16 */
        rand_0_1(&tmp, 1);
        lenSig = (int)(tmp*384.0f + 128.1f); /* random number between 128 and 512 */

        /* Define prototype decoder and compute input signal covariance matrix */
        Q = (float_complex**)calloc2d(nCHout, nCHin, sizeof(float_complex));
        for(i=0; i<MIN(nCHin, nCHout); i++)
            Q[i][i] = cmplxf(1.0f, 0.0f); /* Identity */
        x = (float_complex**)malloc2d(nCHin, lenSig, sizeof(float_complex));
        rand_cmplx_m1_1(FLATTEN2D(x), nCHin*lenSig);
        Cx = (float_complex**)malloc2d(nCHin, nCHin, sizeof(float_complex));
        cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasConjTrans, nCHin, nCHin, lenSig, &calpha,
                    FLATTEN2D(x), lenSig,
                    FLATTEN2D(x), lenSig, &cbeta,
                    FLATTEN2D(Cx), nCHin);

        /* Compute target covariance matrix */
        y = (float_complex**)malloc2d(nCHout, lenSig, sizeof(float_complex));
        rand_cmplx_m1_1(FLATTEN2D(y), nCHout*lenSig);
        Cy = (float_complex**)malloc2d(nCHout, nCHout, sizeof(float_complex));
        cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasConjTrans, nCHout, nCHout, lenSig, &calpha,
                    FLATTEN2D(y), lenSig,
                    FLATTEN2D(y), lenSig, &cbeta,
                    FLATTEN2D(Cy), nCHout);

        /* Compute optimal mixing matrix - with energy compensation enabled */
        M = (float_complex**)malloc2d(nCHout, nCHin, sizeof(float_complex));
        reg = 0.2f;
        cdf4sap_cmplx_create(&hCdf, nCHin, nCHout);
        formulate_M_and_Cr_cmplx(hCdf, FLATTEN2D(Cx), FLATTEN2D(Cy), FLATTEN2D(Q), 1, reg, FLATTEN2D(M), NULL);

        /* Apply mixing matrix to 'x' and assert that it's covariance matrix matches
         * the target covariance matrix */
        z = (float_complex**)malloc2d(nCHout, lenSig, sizeof(float_complex));
        cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nCHout, lenSig, nCHin, &calpha,
                    FLATTEN2D(M), nCHin,
                    FLATTEN2D(x), lenSig, &cbeta,
                    FLATTEN2D(z), lenSig);
        Cz = (float_complex**)malloc2d(nCHout, nCHout, sizeof(float_complex));
        cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasConjTrans, nCHout, nCHout, lenSig, &calpha,
                    FLATTEN2D(z), lenSig,
                    FLATTEN2D(z), lenSig, &cbeta,
                    FLATTEN2D(Cz), nCHout);
        if(nCHin>=nCHout){
            for(i=0; i<nCHout; i++){
                for(j=0; j<nCHout; j++){
                    TEST_ASSERT_FLOAT_WITHIN(acceptedTolerance, crealf(Cy[i][j]), crealf(Cz[i][j]));
                    TEST_ASSERT_FLOAT_WITHIN(acceptedTolerance, cimagf(Cy[i][j]), cimagf(Cz[i][j]));
                }
            }
        }
        else{ /* if nCHin<nCHout, then only the diagonal elements will match */
            for(i=0; i<nCHout; i++){
                TEST_ASSERT_FLOAT_WITHIN(acceptedTolerance, crealf(Cy[i][i]), crealf(Cz[i][i]));
                TEST_ASSERT_FLOAT_WITHIN(acceptedTolerance, cimagf(Cy[i][i]), cimagf(Cz[i][i]));
            }
        }

        /* Determine prototype covariance matrix */
        Q_Cx = (float_complex**)malloc2d(nCHout, nCHin, sizeof(float_complex));
        cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nCHout, nCHin, nCHin, &calpha,
                    FLATTEN2D(Q), nCHin,
                    FLATTEN2D(Cx), nCHin, &cbeta,
                    FLATTEN2D(Q_Cx), nCHin);
        Cp = (float_complex**)malloc2d(nCHout, nCHout, sizeof(float_complex));
        cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasConjTrans, nCHout, nCHout, nCHin, &calpha,
                    FLATTEN2D(Q_Cx), nCHin,
                    FLATTEN2D(Q), nCHin, &cbeta,
                    FLATTEN2D(Cp), nCHout);
        for(i=0; i<nCHout; i++)
            for(j=0; j<nCHout; j++)
                if(i!=j)
                    Cp[i][j] = cmplxf(0.0f, 0.0f); /* Zero non-diagonal elements */

        /* Create perfectly incoherent frame. Note, in practice this would instead
         * be a decorrelated version of the prototype signals, [i.e.
         * decorrelate(Q*x) ]*/
        decor = (float_complex**)malloc2d(nCHout, lenSig, sizeof(float_complex));
        rand_cmplx_m1_1(FLATTEN2D(decor), nCHout*lenSig);

        /* Now compute optimal mixing matrix, but this time also including the
         * residual mixing matrix */
        M = (float_complex**)malloc2d(nCHout, nCHin, sizeof(float_complex));
        reg = 0.2f;
        Cr = (float_complex**)malloc2d(nCHout, nCHout, sizeof(float_complex));
        formulate_M_and_Cr_cmplx(hCdf, FLATTEN2D(Cx), FLATTEN2D(Cy), FLATTEN2D(Q), 0, reg, FLATTEN2D(M), FLATTEN2D(Cr));
        cdf4sap_cmplx_create(&hCdf_res, nCHout, nCHout);
        Mr = (float_complex**)calloc2d(nCHout, nCHout, sizeof(float_complex));
        eye_nCHout = (float_complex**)calloc2d(nCHout, nCHout, sizeof(float_complex));
        for(i=0; i<nCHout; i++)
            eye_nCHout[i][i] = cmplxf(1.0f, 0.0f);
        formulate_M_and_Cr_cmplx(hCdf_res, FLATTEN2D(Cp), FLATTEN2D(Cr), FLATTEN2D(eye_nCHout), 0, reg, FLATTEN2D(Mr), NULL);

        /* Apply mixing matrix to x, and residual mixing matrix to the decorrelated
         * prototype signals, and sum */
        cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nCHout, lenSig, nCHin, &calpha,
                    FLATTEN2D(M), nCHin,
                    FLATTEN2D(x), lenSig, &cbeta,
                    FLATTEN2D(z), lenSig);
        z_r = (float_complex**)malloc2d(nCHout, lenSig, sizeof(float_complex));
        cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nCHout, lenSig, nCHout, &calpha,
                    FLATTEN2D(Mr), nCHout,
                    FLATTEN2D(decor), lenSig, &cbeta,
                    FLATTEN2D(z_r), lenSig);
        utility_cvvadd(FLATTEN2D(z), FLATTEN2D(z_r), nCHout*lenSig, FLATTEN2D(z));

        /* Assert that the covariance matrix of 'z' matches the target covariance
         * matrix */
        if(nCHin>=nCHout){
            for(i=0; i<nCHout; i++){
                for(j=0; j<nCHout; j++){
                    TEST_ASSERT_FLOAT_WITHIN(acceptedTolerance, crealf(Cy[i][j]), crealf(Cz[i][j]));
                    TEST_ASSERT_FLOAT_WITHIN(acceptedTolerance, cimagf(Cy[i][j]), cimagf(Cz[i][j]));
                }
            }
        }
        else{ /* if nCHin<nCHout, then only the diagonal elements will match */
            for(i=0; i<nCHout; i++){
                TEST_ASSERT_FLOAT_WITHIN(acceptedTolerance, crealf(Cy[i][i]), crealf(Cz[i][i]));
                TEST_ASSERT_FLOAT_WITHIN(acceptedTolerance, cimagf(Cy[i][i]), cimagf(Cz[i][i]));
            }
        }

        /* Clean-up */
        cdf4sap_cmplx_destroy(&hCdf);
        cdf4sap_cmplx_destroy(&hCdf_res);
        free(Q);
        free(x);
        free(y);
        free(z);
        free(Cx);
        free(Cy);
        free(Cz);
        free(M);
        free(Cr);
        free(Mr);
        free(Q_Cx);
        free(Cp);
        free(decor);
        free(z_r);
        free(eye_nCHout);
    }
}

void test__getLoudspeakerDecoderMtx(void){
    int i, j, k, nLS, order, nSH;
    float scale;
    float* ls_dirs_deg, *amp, *en;
    float** ls_dirs_rad, **decMtx_SAD, **decMtx_MMD, **decMtx_EPAD, **decMtx_AllRAD, **Ysrc, **LSout;

    /* Config */
    const float acceptedTolerance = 0.00001f;
    int nTestOrders = 10;
    int testOrders[10] = {1,2,3,4,5,6,7,8,9,10};

    /* Loop over orders */
    for(i=0; i<nTestOrders; i++) {
        order = testOrders[i];
        nSH = ORDER2NSH(order);

        /* Pull an appropriate t-design for this order */
        ls_dirs_deg = (float*)__HANDLES_Tdesign_dirs_deg[2 * order-1];
        nLS = __Tdesign_nPoints_per_degree[2 * order-1];
        ls_dirs_rad = (float**)malloc2d(nLS, 2, sizeof(float));
        for(j=0; j<nLS; j++){
            ls_dirs_rad[j][0] = ls_dirs_deg[j*2] * M_PI/180.0f;
            ls_dirs_rad[j][1] = M_PI/2.0f - ls_dirs_deg[j*2+1] * M_PI/180.0f; /* elevation->inclination */
        }

        /* Compute decoders */
        decMtx_SAD = (float**)malloc2d(nLS, nSH, sizeof(float));
        decMtx_MMD = (float**)malloc2d(nLS, nSH, sizeof(float));
        decMtx_EPAD = (float**)malloc2d(nLS, nSH, sizeof(float));
        decMtx_AllRAD = (float**)malloc2d(nLS, nSH, sizeof(float));
        getLoudspeakerDecoderMtx(ls_dirs_deg, nLS, LOUDSPEAKER_DECODER_SAD, order, 0, FLATTEN2D(decMtx_SAD));
        getLoudspeakerDecoderMtx(ls_dirs_deg, nLS, LOUDSPEAKER_DECODER_MMD, order, 0, FLATTEN2D(decMtx_MMD));
        getLoudspeakerDecoderMtx(ls_dirs_deg, nLS, LOUDSPEAKER_DECODER_EPAD, order, 0, FLATTEN2D(decMtx_EPAD));
        getLoudspeakerDecoderMtx(ls_dirs_deg, nLS, LOUDSPEAKER_DECODER_ALLRAD, order, 0, FLATTEN2D(decMtx_AllRAD));

        /* SAD/MMD/EPAD should all be equivalent in this special/uniform case */
        for(j=0; j<nLS; j++)
            for(k=0; k<nSH; k++)
                TEST_ASSERT_FLOAT_WITHIN(acceptedTolerance, decMtx_SAD[j][k], decMtx_MMD[j][k]);
        for(j=0; j<nLS; j++)
            for(k=0; k<nSH; k++)
                TEST_ASSERT_FLOAT_WITHIN(acceptedTolerance, decMtx_SAD[j][k], decMtx_EPAD[j][k]);

        /* Compute output for PWs in direction of Loudspeakers: */
        Ysrc = (float**)malloc2d(nSH, nLS, sizeof(float));
        getSHreal(order, FLATTEN2D(ls_dirs_rad), nLS, FLATTEN2D(Ysrc));
        LSout = (float**)malloc2d(nLS, nLS, sizeof(float));
        cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nLS, nLS, nSH, 1.0f,
                    FLATTEN2D(decMtx_EPAD), nSH,
                    FLATTEN2D(Ysrc), nLS, 0.0f,
                    FLATTEN2D(LSout), nLS);
        
        /* Compute amplitude and energy for each source */
        amp = (float*)calloc1d(nLS, sizeof(float));
        en = (float*)calloc1d(nLS, sizeof(float));
        for (int idxSrc=0; idxSrc<nLS; idxSrc++)
        {
            for (int idxLS=0; idxLS<nLS; idxLS++)
            {
                amp[idxSrc] += LSout[idxLS][idxSrc];
                en[idxSrc] += LSout[idxLS][idxSrc] * LSout[idxLS][idxSrc];
            }
        }
        /* Check output amplitude and Energy */
        for (int idxSrc=0; idxSrc<nLS; idxSrc++)
        {
            TEST_ASSERT_FLOAT_WITHIN(acceptedTolerance, amp[idxSrc], 1.0);
            TEST_ASSERT_FLOAT_WITHIN(acceptedTolerance, en[idxSrc], (float)nSH / (float)nLS);
        }

        /* Clean-up */
        free(decMtx_SAD);
        free(decMtx_MMD);
        free(decMtx_EPAD);
        free(decMtx_AllRAD);
        free(ls_dirs_rad);
        free(amp);
        free(en);
        free(Ysrc);
        free(LSout);
    }
}

void test__getSHreal(void){
    int i, j, k, order, nDirs, nSH;
    float scale;
    float* t_dirs_deg;
    float** t_dirs_rad, **Y, **YYT;

    /* Config */
    const float acceptedTolerance = 0.00001f;
    int nTestOrders = 10;
    int testOrders[10] = {1,2,3,4,5,6,7,8,9,10};

    /* Loop over orders */
    for(i=0; i<nTestOrders; i++){
        order = testOrders[i];
        nSH = ORDER2NSH(order);

        /* Pull an appropriate t-design */
        t_dirs_deg = (float*)__HANDLES_Tdesign_dirs_deg[2*order];
        nDirs = __Tdesign_nPoints_per_degree[2*order];
        t_dirs_rad = (float**)malloc2d(nDirs, 2, sizeof(float));
        for(j=0; j<nDirs; j++){
            t_dirs_rad[j][0] = t_dirs_deg[j*2] * M_PI/180.0f;
            t_dirs_rad[j][1] = M_PI/2.0f - t_dirs_deg[j*2+1] * M_PI/180.0f; /* elevation->inclination */
        }

        /* Compute spherical harmonic coefficients */
        Y = (float**)malloc2d(nSH, nDirs, sizeof(float));
        getSHreal(order, FLATTEN2D(t_dirs_rad), nDirs, FLATTEN2D(Y));
        scale = SQRT4PI;
        utility_svsmul(FLATTEN2D(Y), &scale, nSH*nDirs, FLATTEN2D(Y));

        /* Check Y is orthogonal: */
        YYT = (float**)malloc2d(nSH, nSH, sizeof(float));
        cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasTrans, nSH, nSH, nDirs, 1.0f,
                    FLATTEN2D(Y), nDirs,
                    FLATTEN2D(Y), nDirs, 0.0f,
                    FLATTEN2D(YYT), nSH);

        /* Should be Identity: */
        scale = 1.0f/(float)nDirs;
        utility_svsmul(FLATTEN2D(YYT), &scale, nSH*nSH, FLATTEN2D(YYT));
        for(j=0; j<nSH; j++){
            for(k=0; k<nSH; k++) {
                if(j==k)
                    TEST_ASSERT_FLOAT_WITHIN(acceptedTolerance, 1.0f, YYT[j][k]);
                else
                    TEST_ASSERT_FLOAT_WITHIN(acceptedTolerance, 0.0f, YYT[j][k]);
            }
        }

        /* clean-up */
        free(t_dirs_rad);
        free(Y);
        free(YYT);
    }
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

void test__getSHcomplex(void){
    int i, j, k, order, nDirs, nSH;
    float_complex scale;
    float* t_dirs_deg;
    float** t_dirs_rad;
    float_complex **Y, **YYH;
    const float_complex calpha = cmplxf(1.0f, 0.0f), cbeta = cmplxf(0.0f, 0.0f);

    /* Config */
    const float acceptedTolerance = 0.00001f;
    int nTestOrders = 10;
    int testOrders[10] = {1,2,3,4,5,6,7,8,9,10};

    /* Loop over orders */
    for(i=0; i<nTestOrders; i++){
        order = testOrders[i];
        nSH = ORDER2NSH(order);

        /* Pull an appropriate t-design */
        t_dirs_deg = (float*)__HANDLES_Tdesign_dirs_deg[2*order];
        nDirs = __Tdesign_nPoints_per_degree[2*order];
        t_dirs_rad = (float**)malloc2d(nDirs, 2, sizeof(float));
        for(j=0; j<nDirs; j++){
            t_dirs_rad[j][0] = t_dirs_deg[j*2] * M_PI/180.0f;
            t_dirs_rad[j][1] = M_PI/2.0f - t_dirs_deg[j*2+1] * M_PI/180.0f; /* elevation->inclination */
        }

        /* Compute spherical harmonic coefficients */
        Y = (float_complex**)malloc2d(nSH, nDirs, sizeof(float_complex));
        getSHcomplex(order, FLATTEN2D(t_dirs_rad), nDirs, FLATTEN2D(Y));
        scale = cmplxf(SQRT4PI, 0.0f);
        utility_cvsmul(FLATTEN2D(Y), &scale, nSH*nDirs, FLATTEN2D(Y));

        /* Check Y is orthogonal: */
        YYH = (float_complex**)malloc2d(nSH, nSH, sizeof(float_complex));
        cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasConjTrans, nSH, nSH, nDirs, &calpha,
                    FLATTEN2D(Y), nDirs,
                    FLATTEN2D(Y), nDirs, &cbeta,
                    FLATTEN2D(YYH), nSH);

        /* Should be Identity: */
        scale = cmplxf(1.0f/(float)nDirs, 0.0f);
        utility_cvsmul(FLATTEN2D(YYH), &scale, nSH*nSH, FLATTEN2D(YYH));
        for(j=0; j<nSH; j++){
            for(k=0; k<nSH; k++) {
                if(j==k)
                    TEST_ASSERT_FLOAT_WITHIN(acceptedTolerance, 1.0f, crealf(YYH[j][k]));
                else
                    TEST_ASSERT_FLOAT_WITHIN(acceptedTolerance, 0.0f, crealf(YYH[j][k]));
            }
        }

        /* clean-up */
        free(t_dirs_rad);
        free(Y);
        free(YYH);
    }
}

void test__getSHrotMtxReal(void){
    int i,j,nSH;
    float Rzyx[3][3];
    float** Mrot;

    /* Config */
    const float acceptedTolerance = 0.00001f;
    int order = 22;

    /* Rotation matrix for 0,0,0 should be identity */
    yawPitchRoll2Rzyx(0.0f, 0.0f, 0.0f, 0, Rzyx);
    nSH = ORDER2NSH(order);
    Mrot = (float**)malloc2d(nSH, nSH, sizeof(float));
    getSHrotMtxReal(Rzyx, FLATTEN2D(Mrot), order);
    for(i=0; i<nSH; i++){
        for(j=0; j<nSH; j++) {
            if(j==i)
                TEST_ASSERT_FLOAT_WITHIN(acceptedTolerance, 1.0f, Mrot[i][j]);
            else
                TEST_ASSERT_FLOAT_WITHIN(acceptedTolerance, 0.0f, Mrot[i][j]);
        }
    }
    free(Mrot);

    /* Compare to the getSHrotMtx() Matlab function  */
    order = 4;
    nSH = ORDER2NSH(order);
    Mrot = (float**)malloc2d(nSH, nSH, sizeof(float));
    yawPitchRoll2Rzyx(0.04f, 0.54f,-0.4f, 0, Rzyx);
    getSHrotMtxReal(Rzyx, FLATTEN2D(Mrot), order);
    double Mrot_ref[25][25] = {
            {1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},{0,0.912317819470322,-0.334007492880439,-0.236886451652771,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
            {0,0.408043822669133,0.790002010621868,0.457599237319041,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
            {0,0.0342991990938353,-0.514135991653113,0.857022605902780,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
            {0,0,0,0,0.773751979486127,-0.480511616313319,0.297436898769771,-0.164460121209763,-0.234308814625387,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
            {0,0,0,0,0.320815885111266,0.584443217512645,-0.457030341925157,-0.339982347095703,-0.480664710153360,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
            {0,0,0,0,0.323409465640717,0.558336000748573,0.436154765179890,0.626143845136656,0.0371501522262563,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
            {0,0,0,0,0.365398067572425,-0.182693579159072,-0.703504421517165,0.441781344152855,0.378177314513551,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
            {0,0,0,0,0.245510920021695,0.287086534852415,0.132306868781138,-0.519748017168846,0.754759962358177,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
            {0,0,0,0,0,0,0,0,0,0.642754542747763,-0.587652464622319,0.146359326676735,-0.179940097166632,0.249957116297551,-0.161211805496773,-0.315061710316419,0,0,0,0,0,0,0,0,0},
            {0,0,0,0,0,0,0,0,0,0.316547622267400,0.324276933833715,-0.489415761677808,0.525421745728824,-0.0811795764406443,-0.0642914639380568,-0.517998801533831,0,0,0,0,0,0,0,0,0},
            {0,0,0,0,0,0,0,0,0,-0.0477608186606479,0.302122638638019,0.214473275742620,-0.433723919089070,-0.427443247772927,-0.611726955971008,-0.339717518973177,0,0,0,0,0,0,0,0,0},
            {0,0,0,0,0,0,0,0,0,0.148935636035543,0.571302238306694,0.529863460253249,0.0476038953094580,0.594213419796629,0.0656256769672685,-0.104948528910382,0,0,0,0,0,0,0,0,0},
            {0,0,0,0,0,0,0,0,0,0.311309233760352,0.304630835298635,-0.396153335826512,-0.667628966408715,-0.0103234397880398,0.454946318162605,0.0231945482299087,0,0,0,0,0,0,0,0,0},
            {0,0,0,0,0,0,0,0,0,0.514785682894208,0.113244732089517,0.407883773582348,0.233719845299723,-0.593950310633879,0.241281704427283,0.300305444687571,0,0,0,0,0,0,0,0,0},
            {0,0,0,0,0,0,0,0,0,0.316675769196523,0.161927142796105,-0.298312669792114,0.0285933354722383,0.205549150173188,-0.571110978701303,0.644414328446904,0,0,0,0,0,0,0,0,0},
            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.526471642263643,-0.616929911516989,0.267922897453092,0.0235630456100945,0.0776050535864247,-0.190481327947399,0.295565129451190,-0.0753134473777231,-0.366811472459093},
            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.234144273956354,0.0978665390875757,-0.545910447747527,0.175528558261790,-0.376101588123769,0.335795191612168,-0.141736252789070,-0.0455702308901721,-0.574798644029333},
            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.0718436126062899,0.305262278899232,-0.0197737560173443,-0.298299395229287,0.646776790379034,0.111401675977437,0.0997398996043224,-0.463839920427382,-0.395542458465569},
            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.155033529872408,-0.118946002867737,0.138228495430813,-0.0977208017941514,-0.285522105871139,-0.450196541284017,-0.600496309285322,-0.520682311298467,-0.131355606942160},
            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.0236933293789157,0.311297649179989,0.703254159219873,0.348811131545197,-0.261303521121084,0.391172954707122,0.0807830377413570,-0.219358047572331,-0.101769931423874},
            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.146767948839247,0.439950893376704,0.0598087344890290,-0.520771343866458,-0.439502688322895,-0.362741803354952,0.407296904607327,0.0826968395396408,-0.112466610956744},
            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.386795790652846,0.451176951621299,0.0223488932476933,0.463808781391941,0.287701399151563,-0.482347736946315,-0.226762742725175,0.241251512069808,-0.0784553883303562},
            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.576800968786616,0.0555128465726625,0.144555412279657,-0.473213285269062,0.0597643274078365,0.343735767588532,-0.480720100388111,0.108090832343090,0.234286982126144},
            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.366598721881537,0.0733558553140817,-0.301930038675134,0.195400170636906,-0.0699710544219968,-0.0214401526687090,0.258994980191915,-0.617374325026823,0.526589247038282}};
    for(i=0; i<nSH; i++)
        for(j=0; j<nSH; j++)
            TEST_ASSERT_FLOAT_WITHIN(acceptedTolerance, (float)Mrot_ref[i][j], Mrot[i][j]);
    free(Mrot);
}

void test__real2complexSHMtx(void){
    int o, it, j, nSH, order;
    float* Y_real_ref;
    float_complex* Y_complex_ref, * tmp, *Y_complex_test;
    float_complex** T_r2c;
    float dir[2];
    const float_complex calpha = cmplxf(1.0f, 0.0f), cbeta = cmplxf(0.0f, 0.0f);

    /* Config */
    const float acceptedTolerance = 0.0000001f;
    int nTestOrders = 10;
    int testOrders[10] = {1,2,3,4,5,6,7,8,9,10};
    int nIter = 400;

    /* Loop over orders */
    for(o=0; o<nTestOrders; o++){
        order = testOrders[o];
        nSH = ORDER2NSH(order);
        Y_real_ref = malloc1d(nSH * sizeof(float));
        tmp = malloc1d(nSH * sizeof(float_complex));
        Y_complex_ref = malloc1d(nSH * sizeof(float_complex));
        Y_complex_test = malloc1d(nSH * sizeof(float_complex));
        T_r2c = (float_complex**)malloc2d(nSH, nSH, sizeof(float_complex));

        /* Loop over iterations */
        for(it=0; it<nIter; it++) {
            /* Random direction */
            rand_m1_1(&dir[0] , 1);
            rand_m1_1(&dir[1] , 1);
            dir[0] *= M_PI;
            dir[1] *= M_PI/2.0f;

            /* Compute reference spherical harmonic weights */
            getSHcomplex(order, (float*)dir, 1, Y_complex_ref);
            getSHreal(order, (float*)dir, 1, Y_real_ref);

            /* Convert to complex weights */
            real2complexSHMtx(order, FLATTEN2D(T_r2c));
            for(j=0; j<nSH; j++)
                tmp[j] = cmplxf(Y_real_ref[j], 0.0f);
            cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasTrans, 1, nSH, nSH, &calpha,
                        tmp, nSH,
                        FLATTEN2D(T_r2c), nSH, &cbeta, /* Had to transpose it! */
                        Y_complex_test, nSH);

            /* Should be equal to the reference */
            for (j = 0; j < nSH; j++) {
                TEST_ASSERT_FLOAT_WITHIN(acceptedTolerance, crealf(Y_complex_ref[j]), crealf(Y_complex_test[j]));
                TEST_ASSERT_FLOAT_WITHIN(acceptedTolerance, cimagf(Y_complex_ref[j]), cimagf(Y_complex_test[j]));
            }
        }

        /* Clean-up */
        free(tmp);
        free(Y_real_ref);
        free(Y_complex_ref);
        free(Y_complex_test);
        free(T_r2c);
    }
}

void test__complex2realSHMtx(void){
    int o, it, j, nSH, order;
    float* Y_real_ref;
    float_complex* Y_complex_ref, *Y_real_test;
    float_complex** T_c2r;
    float dir[2];
    const float_complex calpha = cmplxf(1.0f, 0.0f), cbeta = cmplxf(0.0f, 0.0f);

    /* Config */
    const float acceptedTolerance = 0.000001f;
    int nTestOrders = 10;
    int testOrders[10] = {1,2,3,4,5,6,7,8,9,10};
    int nIter = 400;

    /* Loop over orders */
    for(o=0; o<nTestOrders; o++){
        order = testOrders[o];
        nSH = ORDER2NSH(order);
        Y_real_ref = malloc1d(nSH * sizeof(float));
        Y_complex_ref = malloc1d(nSH * sizeof(float_complex));
        Y_real_test = malloc1d(nSH * sizeof(float_complex));
        T_c2r = (float_complex**)malloc2d(nSH, nSH, sizeof(float_complex));

        /* Loop over iterations */
        for(it=0; it<nIter; it++) {
            /* Random direction */
            rand_m1_1(&dir[0] , 1);
            rand_m1_1(&dir[1] , 1);
            dir[0] *= M_PI;
            dir[1] *= M_PI/2.0f;

            /* Compute reference spherical harmonic weights */
            getSHcomplex(order, (float*)dir, 1, Y_complex_ref);
            getSHreal(order, (float*)dir, 1, Y_real_ref);

            /* Convert to complex weights */
            complex2realSHMtx(order, FLATTEN2D(T_c2r));
            cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasTrans, 1, nSH, nSH, &calpha,
                        Y_complex_ref, nSH,
                        FLATTEN2D(T_c2r), nSH, &cbeta, /* Had to transpose it! */
                        Y_real_test, nSH);

            /* Should be equal to the reference */
            for (j = 0; j < nSH; j++) {
                TEST_ASSERT_FLOAT_WITHIN(acceptedTolerance, Y_real_ref[j], crealf(Y_real_test[j]));
            }
        }

        /* Clean-up */
        free(Y_real_ref);
        free(Y_complex_ref);
        free(Y_real_test);
        free(T_c2r);
    }
}

void test__computeSectorCoeffsEP(void){
    int i, j, numSec, order_sec, nSH_sec, nSH;
    float* sec_dirs_deg;
    float** sectorCoeffs;
    float_complex*** A_xyz;

    /* Config */
    const float acceptedTolerance = 0.000001f;
    const int order = 2;

    /* Sector design and compute coefficients */
    order_sec = order-1;
    numSec = __Tdesign_nPoints_per_degree[2*order_sec-1];
    sec_dirs_deg = (float*)__HANDLES_Tdesign_dirs_deg[2*order_sec-1];
    nSH = ORDER2NSH(order);
    nSH_sec = ORDER2NSH(order_sec);
    A_xyz = (float_complex***)malloc3d(nSH, nSH_sec, 3, sizeof(float_complex));
    computeVelCoeffsMtx(order_sec, FLATTEN3D(A_xyz));
    sectorCoeffs = (float**)malloc2d((numSec*4),nSH,sizeof(float));
    computeSectorCoeffsEP(order_sec, FLATTEN3D(A_xyz), SECTOR_PATTERN_PWD, sec_dirs_deg, numSec, FLATTEN2D(sectorCoeffs));

    /* Check with Matlab reference */
    double sectorCoeffs_ref[9][16]= {
            {0.886226925452758,0.511663353973244,0.511663353973244,0.511663353973244,0.886226925452758,0.511663353973244,-0.511663353973244,-0.511663353973244,0.886226925452758,-0.511663353973244,0.511663353973244,-0.511663353973244,0.886226925452758,-0.511663353973244,-0.511663353973244,0.511663353973244},
            {0.886226925452758,0,0.511663353973244,0,-0.886226925452758,0,0.511663353973244,0,0.886226925452758,0,0.511663353973244,0,-0.886226925452758,0,0.511663353973244,0},
            {0.886226925452758,0,0,0.511663353973244,-0.886226925452758,0,0,0.511663353973244,-0.886226925452758,0,0,0.511663353973244,0.886226925452758,0,0,0.511663353973244},
            {0.886226925452758,0.511663353973244,0,0,0.886226925452758,0.511663353973244,0,0,-0.886226925452758,0.511663353973244,0,0,-0.886226925452758,0.511663353973244,0,0},
            {0,0.396332729760601,0.396332729760601,0,0,-0.396332729760601,0.396332729760601,0,0,0.396332729760601,-0.396332729760601,0,0,-0.396332729760601,-0.396332729760601,0},
            {0,0,0.396332729760601,0.396332729760601,0,0,-0.396332729760601,-0.396332729760601,0,0,-0.396332729760601,0.396332729760601,0,0,0.396332729760601,-0.396332729760601},
            {0,-0.228822808215942,-0.228822808215942,0.457645616431885,0,-0.228822808215942,0.228822808215942,-0.457645616431885,0,0.228822808215942,-0.228822808215942,-0.457645616431885,0,0.228822808215942,0.228822808215942,0.457645616431885},
            {0,0.396332729760601,0,0.396332729760601,0,-0.396332729760601,0,0.396332729760601,0,-0.396332729760601,0,-0.396332729760601,0,0.396332729760601,0,-0.396332729760601},
            {0,0.396332729760601,-0.396332729760601,0,0,0.396332729760601,0.396332729760601,0,0,-0.396332729760601,-0.396332729760601,0,0,-0.396332729760601,0.396332729760601,0}
    };
    for(i=0; i<9; i++)
        for(j=0; j<16; j++)
            TEST_ASSERT_FLOAT_WITHIN(acceptedTolerance, (float)sectorCoeffs_ref[i][j], sectorCoeffs[j][i]);

    free(sectorCoeffs);
    free(A_xyz);
}

void test__checkCondNumberSHTReal(void){
    int i, j, nDirs, order, nSH;
    float* t_dirs_deg, *cond_N;
    float** t_dirs_rad;

    /* Config */
    const float acceptedTolerance = 0.00001f;
    int nTestOrders = 10;
    int testOrders[10] = {1,2,3,4,5,6,7,8,9,10};

    /* Loop over orders */
    for(i=0; i<nTestOrders; i++) {
        order = testOrders[i];
        nSH = ORDER2NSH(order);

        /* Pull an appropriate t-design */
        t_dirs_deg = (float*)__HANDLES_Tdesign_dirs_deg[2 * order];
        nDirs = __Tdesign_nPoints_per_degree[2 * order];
        t_dirs_rad = (float **) malloc2d(nDirs, 2, sizeof(float));
        for (j = 0; j < nDirs; j++) {
            t_dirs_rad[j][0] = t_dirs_deg[j * 2] * M_PI / 180.0f;
            t_dirs_rad[j][1] = M_PI / 2.0f - t_dirs_deg[j * 2 + 1] * M_PI /
                                             180.0f; /* elevation->inclination */
        }

        /* Condition numbers for an appropriate t-design should be 1 */
        cond_N = malloc1d((order+1)*sizeof(float));
        checkCondNumberSHTReal(order, FLATTEN2D(t_dirs_rad), nDirs, NULL, cond_N);
        for(j=0; j<order+1; j++)
            TEST_ASSERT_FLOAT_WITHIN(acceptedTolerance, 1.0f, cond_N[j]);

        /* Clean-up */
        free(t_dirs_rad);
        free(cond_N);
    }
}

void test__sphMUSIC(void){
    int i, j, k, nGrid, nSH, nSrcs, srcInd_1, srcInd_2;
    float test_dirs_deg[2][2];
    float* grid_dirs_deg;
    float** Y_src, **src_sigs, **src_sigs_sh, **Cx, **V, **Vn;
    float_complex** Vn_cmplx;
    void* hMUSIC;

    /* config */
    const int order = 3;
    const int lsig = 48000;

    /* define scanning grid directions */
    nGrid = 240;
    grid_dirs_deg = (float*)__Tdesign_degree_21_dirs_deg;

    /* test scenario and signals */
    nSrcs = 2;
    srcInd_1 = 139;
    srcInd_2 = 204;
    test_dirs_deg[0][0] = grid_dirs_deg[srcInd_1*2];
    test_dirs_deg[0][1] = grid_dirs_deg[srcInd_1*2+1];
    test_dirs_deg[1][0] = grid_dirs_deg[srcInd_2*2];
    test_dirs_deg[1][1] = grid_dirs_deg[srcInd_2*2+1];
    nSH = ORDER2NSH(order);
    Y_src = (float**)malloc2d(nSH, nSrcs, sizeof(float));
    getRSH(order, (float*)test_dirs_deg, nSrcs, FLATTEN2D(Y_src));
    src_sigs = (float**)malloc2d(nSrcs, lsig, sizeof(float));
    rand_m1_1(FLATTEN2D(src_sigs), nSrcs*lsig); /* uncorrelated noise sources */

    /* encode to SH and compute spatial covariance matrix */
    src_sigs_sh = (float**)malloc2d(nSH, lsig, sizeof(float));
    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nSH, lsig, nSrcs, 1.0f,
                FLATTEN2D(Y_src), nSrcs,
                FLATTEN2D(src_sigs), lsig, 0.0f,
                FLATTEN2D(src_sigs_sh), lsig);
    Cx = (float**)malloc2d(nSH, nSH, sizeof(float));
    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasTrans, nSH, nSH, lsig, 1.0f,
                FLATTEN2D(src_sigs_sh), lsig,
                FLATTEN2D(src_sigs_sh), lsig, 0.0f,
                FLATTEN2D(Cx), nSH);

    /* Eigenvalue decomposition and truncation of eigen vectors to obtain
     * noise subspace (based on source number) */
    V = (float**)malloc2d(nSH, nSH, sizeof(float));
    utility_sseig(FLATTEN2D(Cx), nSH, 1, FLATTEN2D(V), NULL, NULL);
    Vn = (float**)malloc2d(nSH, (nSH-nSrcs), sizeof(float)); /* noise subspace */
    for(i=0; i<nSH; i++)
        for(j=0, k=nSrcs; j<nSH-nSrcs; j++, k++)
            Vn[i][j] = V[i][k];
    Vn_cmplx = (float_complex**)malloc2d(nSH, (nSH-nSrcs), sizeof(float_complex)); /* noise subspace (complex) */
    for(i=0; i<nSH; i++)
        for(j=0; j<nSH-nSrcs; j++)
            Vn_cmplx[i][j] = cmplxf(Vn[i][j], 0.0f);

    /* compute sphMUSIC, returning "peak-find" indices */
    int inds[2];
    sphMUSIC_create(&hMUSIC, order, grid_dirs_deg, nGrid);
    sphMUSIC_compute(hMUSIC, FLATTEN2D(Vn_cmplx), nSrcs, NULL, (int*)inds);

    /* Assert that the true source indices were found (note that the order can flip) */
    TEST_ASSERT_TRUE(inds[0] == srcInd_1 || inds[0] == srcInd_2);
    TEST_ASSERT_TRUE(inds[1] == srcInd_1 || inds[1] == srcInd_2);

    /* clean-up */
    sphMUSIC_destroy(&hMUSIC);
    free(Y_src);
    free(src_sigs);
    free(src_sigs_sh);
    free(Cx);
    free(V);
    free(Vn);
    free(Vn_cmplx);
}

void test__sphESPRIT(void){
    int i,j,nSH, nSrcs;
    void* hESPRIT;
    float test_dirs_deg[2][2], estdirs_deg[2][2];
    float** Y_src, **src_sigs, **src_sigs_sh, **tmpCx;
    float_complex** Cx, **C_Cx, **T_r2c, **Cx_R, **U, **Us;
    const float_complex calpha = cmplxf(1.0f, 0.0f), cbeta = cmplxf(0.0f, 0.0f);

    /* config */
    const float acceptedTolerance = 0.01f; /* degrees */
    const int order = 3;
    const int lsig = 48000;

    /* test scenario and signals */
    nSrcs = 2;
    test_dirs_deg[0][0] = -90.0f;
    test_dirs_deg[0][1] = 10.0f;
    test_dirs_deg[1][0] = 20.0f;
    test_dirs_deg[1][1] = -40.0f;
    nSH = ORDER2NSH(order);
    Y_src = (float**)malloc2d(nSH, nSrcs, sizeof(float));
    getRSH(order, (float*)test_dirs_deg, nSrcs, FLATTEN2D(Y_src));
    src_sigs = (float**)malloc2d(nSrcs, lsig, sizeof(float));
    rand_m1_1(FLATTEN2D(src_sigs), nSrcs*lsig); /* uncorrelated noise sources */

    /* encode to SH and compute spatial covariance matrix */
    src_sigs_sh = (float**)malloc2d(nSH, lsig, sizeof(float));
    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nSH, lsig, nSrcs, 1.0f,
                FLATTEN2D(Y_src), nSrcs,
                FLATTEN2D(src_sigs), lsig, 0.0f,
                FLATTEN2D(src_sigs_sh), lsig);
    tmpCx = (float**)malloc2d(nSH, nSH, sizeof(float));
    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasTrans, nSH, nSH, lsig, 1.0f,
                FLATTEN2D(src_sigs_sh), lsig,
                FLATTEN2D(src_sigs_sh), lsig, 0.0f,
                FLATTEN2D(tmpCx), nSH);
    Cx = (float_complex**)malloc2d(nSH, nSH, sizeof(float_complex));
    for(i=0; i<nSH; i++)
        for(j=0; j<nSH; j++)
            Cx[i][j] = cmplxf(tmpCx[i][j], 0.0f); /* real->complex data-type */

    /* Convert to complex basis */
    T_r2c = (float_complex**)malloc2d(nSH, nSH, sizeof(float_complex));
    real2complexSHMtx(order, FLATTEN2D(T_r2c));
    for(i=0; i<nSH; i++)
        for(j=0; j<nSH; j++)
            T_r2c[i][j] = conjf(T_r2c[i][j]);
    Cx_R = (float_complex**)malloc2d(nSH, nSH, sizeof(float_complex));
    cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasConjTrans, nSH, nSH, nSH, &calpha,
                FLATTEN2D(Cx), nSH,
                FLATTEN2D(T_r2c), nSH, &cbeta,
                FLATTEN2D(Cx_R), nSH);
    C_Cx = (float_complex**)malloc2d(nSH, nSH, sizeof(float_complex));
    cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nSH, nSH, nSH, &calpha,
                FLATTEN2D(T_r2c), nSH,
                FLATTEN2D(Cx_R), nSH, &cbeta,
                FLATTEN2D(C_Cx), nSH);

    /* Eigenvalue decomposition and truncation of eigen vectors to obtain
     * signal subspace (based on source number) */
    U = (float_complex**)malloc2d(nSH, nSH, sizeof(float_complex));
    utility_cseig(FLATTEN2D(C_Cx), nSH, 1, FLATTEN2D(U), NULL, NULL);
    Us = (float_complex**)malloc2d(nSH, nSrcs, sizeof(float_complex)); /* signal subspace */
    for(i=0; i<nSH; i++)
        for(j=0; j<nSrcs; j++)
            Us[i][j] = U[i][j];

    /* use sphESPRIT to estimate source directions... */
    sphESPRIT_create(&hESPRIT, order);
    sphESPRIT_estimateDirs(hESPRIT, FLATTEN2D(Us), nSrcs, (float*)estdirs_deg);
    for(i=0; i<nSrcs; i++){
        estdirs_deg[i][0]*=180.0f/M_PI; /* rad->deg */
        estdirs_deg[i][1]*=180.0f/M_PI;
    }

    /* Assert that the true source directions were found (note that the order can flip) */
    TEST_ASSERT_TRUE(fabsf(estdirs_deg[0][0]-test_dirs_deg[0][0])<acceptedTolerance ||
                     fabsf(estdirs_deg[1][0]-test_dirs_deg[0][0])<acceptedTolerance);
    TEST_ASSERT_TRUE(fabsf(estdirs_deg[0][1]-test_dirs_deg[0][1])<acceptedTolerance ||
                     fabsf(estdirs_deg[1][1]-test_dirs_deg[0][1])<acceptedTolerance);
    TEST_ASSERT_TRUE(fabsf(estdirs_deg[0][0]-test_dirs_deg[1][0])<acceptedTolerance ||
                     fabsf(estdirs_deg[1][0]-test_dirs_deg[1][0])<acceptedTolerance);
    TEST_ASSERT_TRUE(fabsf(estdirs_deg[0][1]-test_dirs_deg[1][1])<acceptedTolerance ||
                     fabsf(estdirs_deg[1][1]-test_dirs_deg[1][1])<acceptedTolerance);

    /* clean-up */
    sphESPRIT_destroy(&hESPRIT);
    free(Y_src);
    free(src_sigs);
    free(src_sigs_sh);
    free(tmpCx);
    free(Cx);
    free(C_Cx);
    free(T_r2c);
    free(Cx_R);
    free(U);
    free(Us);
}

void test__sphModalCoeffs(void){
    int i, j;
    float* freqVector;
    double* kr;
    double_complex** b_N_dipole, **b_N_card, **b_N_omni, **b_N_omni_test;

    /* Config */
    const double acceptedTolerance = 0.000001f;
    const int order = 4;
    const int N = 16;
    const float fs = 48000;
    const double radius = 0.04;
    const double c = 343.0;

    /* prep */
    freqVector = malloc1d((N/2+1)*sizeof(float));
    getUniformFreqVector(N, fs, freqVector);
    kr = malloc1d((N/2+1)*sizeof(double));
    for(i=0; i<N/2+1; i++)
        kr[i] = 2.0*SAF_PId* (double)freqVector[i] * radius/c;
    b_N_dipole = (double_complex**)malloc2d((N/2+1), (order+1), sizeof(double_complex));
    b_N_card = (double_complex**)malloc2d((N/2+1), (order+1), sizeof(double_complex));
    b_N_omni = (double_complex**)malloc2d((N/2+1), (order+1), sizeof(double_complex));
    b_N_omni_test = (double_complex**)malloc2d((N/2+1), (order+1), sizeof(double_complex));

    /* Compute modal coefficients */
    sphModalCoeffs(order, kr, (N/2+1), ARRAY_CONSTRUCTION_OPEN_DIRECTIONAL, 0.0, FLATTEN2D(b_N_dipole));
    sphModalCoeffs(order, kr, (N/2+1), ARRAY_CONSTRUCTION_OPEN_DIRECTIONAL, 0.5, FLATTEN2D(b_N_card));
    sphModalCoeffs(order, kr, (N/2+1), ARRAY_CONSTRUCTION_OPEN_DIRECTIONAL, 1.0, FLATTEN2D(b_N_omni));
    sphModalCoeffs(order, kr, (N/2+1), ARRAY_CONSTRUCTION_OPEN, 666.0 /* not used */, FLATTEN2D(b_N_omni_test));

    /* Check that "open directional", with "dirCoeff=1" is identical to just "open" */
    for(i=0; i<N/2+1; i++){
        for(j=0; j< order+1; j++){
            TEST_ASSERT_TRUE( fabs(creal(b_N_omni[i][j]) - creal(b_N_omni_test[i][j])) <= acceptedTolerance );
            TEST_ASSERT_TRUE( fabs(cimag(b_N_omni[i][j]) - cimag(b_N_omni_test[i][j])) <= acceptedTolerance );
        }
    }

    /* clean-up */
    free(b_N_dipole);
    free(b_N_card);
    free(b_N_omni);
    free(b_N_omni_test);
}

void test__truncationEQ(void)
{
    double *kr;
    float *w_n, *gain, *gainDB;

    /* Config */
    const int order_truncated = 4;
    const int order_target = 42;
    const float softThreshold = 12.0f;
    const int enableMaxRE = 1;
    const double fs = 48000;
    const int nBands = 128;
    kr = malloc1d(nBands * sizeof(double));
    const double r = 0.085;
    const double c = 343.;
    w_n = calloc1d((order_truncated+1), sizeof(float));
    gain = malloc1d(nBands * sizeof(float));

    /* Prep */
    double* freqVector = malloc1d(nBands*sizeof(double));
    for (int k=0; k<nBands; k++)
    {
        freqVector[k] = (double)k * fs/(2.0*((double)nBands-1));
        kr[k] = 2.0*SAF_PId / c * freqVector[k] * r;
    }
    if (enableMaxRE) {
        // maxRE as order weighting
        float *maxRECoeffs = malloc1d((order_truncated+1) * sizeof(float));
        beamWeightsMaxEV(order_truncated, maxRECoeffs);
        for (int idx_n=0; idx_n<order_truncated+1; idx_n++) {
            w_n[idx_n] = maxRECoeffs[idx_n];
            w_n[idx_n] /= sqrtf((float)(2*idx_n+1) / (4.0*SAF_PI));
        }
        float w_0 = w_n[0];
        for (int idx_n=0; idx_n<order_truncated+1; idx_n++)
            w_n[idx_n] /= w_0;
        free(maxRECoeffs);
    }
    else {
        // just truncation, no tapering
        for (int idx_n=0; idx_n<order_truncated+1; idx_n++)
            w_n[idx_n] = 1.0f;
    }

    truncationEQ(w_n, order_truncated, order_target, kr, nBands, softThreshold, gain);

    /* Asserting gain offset */
    TEST_ASSERT_TRUE(gain[0]-1.0 < 2.0e-6);

    /* Asserting that gain within 0 and 12 (+6db soft clip) */
    gainDB = malloc1d(nBands * sizeof(double));
    for (int idxBand=0; idxBand<nBands; idxBand++){
        gainDB[idxBand] = 20.0*log10(gain[idxBand]);
        TEST_ASSERT_TRUE(gainDB[idxBand] > 0-2.0e-6);
        TEST_ASSERT_TRUE(gainDB[idxBand] < softThreshold + 6.0 + 0-2.0e-6);
    }

    /* clean-up */
    free(kr);
    free(w_n);
    free(freqVector);
    free(gain);
    free(gainDB);
}

#if SAF_ENABLE_EXAMPLES_TESTS == 1
void test__saf_example_ambi_bin(void){
    int nSH, i, ch, framesize;
    void* hAmbi;
    float leftEarEnergy, rightEarEnergy, direction_deg[2];
    float* inSig, *y;
    float** shSig, **binSig, **shSig_frame, **binSig_frame;

    /* Config */
    const int order = 4;
    const int fs = 48000;
    const int signalLength = fs*2;

    /* Create and initialise an instance of ambi_bin */
    ambi_bin_create(&hAmbi);

    /* Configure and initialise the ambi_bin codec */
    ambi_bin_setNormType(hAmbi, NORM_N3D);
    ambi_bin_setInputOrderPreset(hAmbi, (SH_ORDERS)order);
    ambi_bin_initCodec(hAmbi); /* Can be called whenever (thread-safe) */
    /* "initCodec" should be called after calling any of the "set" functions.
     * It should be noted that intialisations are only conducted if they are
     * needed, so calling this function periodically with a timer on a separate
     * thread is perfectly safe and viable. Also, if the intialisations take
     * longer than it takes to "process" the current block of samples, then the
     * output is simply muted/zeroed during this time. */

    ambi_bin_init(hAmbi, fs); /* Should be called before calling "process"
                               * Cannot be called while "process" is on-going */
    ambi_bin_initCodec(hAmbi); /* Can be called whenever (thread-safe) */

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
    framesize = ambi_bin_getFrameSize();
    binSig = (float**)calloc2d(NUM_EARS,signalLength,sizeof(float));
    shSig_frame = (float**)malloc1d(nSH*sizeof(float*));
    binSig_frame = (float**)malloc1d(NUM_EARS*sizeof(float*));
    for(i=0; i<(int)((float)signalLength/(float)framesize); i++){
        for(ch=0; ch<nSH; ch++)
            shSig_frame[ch] = &shSig[ch][i*framesize];
        for(ch=0; ch<NUM_EARS; ch++)
            binSig_frame[ch] = &binSig[ch][i*framesize];

        ambi_bin_process(hAmbi, shSig_frame, binSig_frame, nSH, NUM_EARS, framesize);
    }

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
    free(shSig_frame);
    free(binSig_frame);
}

void test__saf_example_ambi_dec(void){
    int nSH, i, j, ch, max_ind, framesize;
    void* hAmbi;
    float loudspeakerEnergy[22], direction_deg[2];
    float* inSig, *y;
    float** shSig, **lsSig, **shSig_frame, **lsSig_frame;

    /* Config */
    const int order = 4;
    const int fs = 48000;
    const int signalLength = fs*2;

    /* Create and initialise an instance of ambi_dec */
    ambi_dec_create(&hAmbi);

    /* Configure and initialise the ambi_dec codec */
    ambi_dec_setNormType(hAmbi, NORM_N3D);
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

    ambi_dec_init(hAmbi, fs); /* Should be called before calling "process"
                               * Cannot be called while "process" is on-going */

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
    framesize = ambi_dec_getFrameSize();
    lsSig = (float**)calloc2d(22,signalLength,sizeof(float));
    shSig_frame = (float**)malloc1d(nSH*sizeof(float*));
    lsSig_frame = (float**)malloc1d(22*sizeof(float*));
    for(i=0; i<(int)((float)signalLength/(float)framesize); i++){
        for(ch=0; ch<nSH; ch++)
            shSig_frame[ch] = &shSig[ch][i*framesize];
        for(ch=0; ch<22; ch++)
            lsSig_frame[ch] = &lsSig[ch][i*framesize];

        ambi_dec_process(hAmbi, shSig_frame, lsSig_frame, nSH, 22, framesize);
    }

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
    free(shSig_frame);
    free(lsSig_frame);
}

void test__saf_example_ambi_enc(void){
    int nSH, i, ch, framesize, j, delay;
    void* hAmbi;
    float direction_deg[2][2];
    float** inSig, *y;
    float** shSig, **shSig_ref, **inSig_frame, **shSig_frame;

    /* Config */
    const float acceptedTolerance = 0.000001f;
    const int order = 4;
    const int fs = 48000;
    const int signalLength = fs*2;
    direction_deg[0][0] = 90.0f; /* encode to loudspeaker direction: index 8 */
    direction_deg[0][1] = 0.0f;
    direction_deg[1][0] = 20.0f; /* encode to loudspeaker direction: index 8 */
    direction_deg[1][1] = -45.0f;
    delay = ambi_enc_getProcessingDelay();

    /* Create and initialise an instance of ambi_enc */
    ambi_enc_create(&hAmbi);
    ambi_enc_init(hAmbi, fs); /* Cannot be called while "process" is on-going */

    /* Configure ambi_enc */
    ambi_enc_setOutputOrder(hAmbi, (SH_ORDERS)order);
    ambi_enc_setNormType(hAmbi, NORM_N3D); /* (The default for all SH-related examples is SN3D) */
    ambi_enc_setEnablePostScaling(hAmbi, 0); /* Disable scaling output by number of input channels */
    ambi_enc_setNumSources(hAmbi, 2);
    ambi_enc_setSourceAzi_deg(hAmbi, 0, direction_deg[0][0]);
    ambi_enc_setSourceElev_deg(hAmbi, 0, direction_deg[0][1]);
    ambi_enc_setSourceAzi_deg(hAmbi, 1, direction_deg[1][0]);
    ambi_enc_setSourceElev_deg(hAmbi, 1, direction_deg[1][1]);

    /* Define input mono signal */
    nSH = ORDER2NSH(order);
    inSig = (float**)malloc2d(2,signalLength,sizeof(float));
    shSig_ref = (float**)malloc2d(nSH,signalLength,sizeof(float));
    rand_m1_1(FLATTEN2D(inSig), 2*signalLength); /* Mono white-noise signal */

    /* Encode reference */
    y = malloc1d(nSH*2*sizeof(float));
    getRSH(order, (float*)direction_deg, 2, y); /* SH plane-wave weights */
    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nSH, signalLength, 2, 1.0f,
                y, 2,
                FLATTEN2D(inSig), signalLength, 0.0f,
                FLATTEN2D(shSig_ref), signalLength);

    /* Encode via ambi_enc */
    framesize = ambi_enc_getFrameSize();
    shSig = (float**)calloc2d(nSH,signalLength,sizeof(float));
    inSig_frame = (float**)malloc1d(2*sizeof(float*));
    shSig_frame = (float**)malloc1d(nSH*sizeof(float*));
    for(i=0; i<(int)((float)signalLength/(float)framesize); i++){
        for(ch=0; ch<2; ch++)
            inSig_frame[ch] = &inSig[ch][i*framesize];
        for(ch=0; ch<nSH; ch++)
            shSig_frame[ch] = &shSig[ch][i*framesize];

        ambi_enc_process(hAmbi, inSig_frame, shSig_frame, 2, nSH, framesize);
    }

    /* ambi_enc should be equivalent to the reference, except delayed due to the
     * temporal interpolation employed in ambi_enc */
    for(i=0; i<nSH; i++)
        for(j=0; j<signalLength-delay-framesize; j++)
            TEST_ASSERT_FLOAT_WITHIN(acceptedTolerance, shSig_ref[i][j], shSig[i][j+delay]);

    /* Clean-up */
    ambi_enc_destroy(&hAmbi);
    free(inSig);
    free(shSig);
    free(shSig_ref);
    free(y);
}

void test__saf_example_array2sh(void){
    int nSH, i, j, framesize, ch;
    void* hA2sh, *safFFT, *hMC;
    float direction_deg[2], radius;
    float* inSig, *f;
    float** shSig, **inSig_32, **micSig, **h_array, **micSig_frame, **shSig_frame;
    double* kr;
    float_complex* tmp_H;
    float_complex*** H_array;

    /* Config */
    const int order = 4;
    const int fs = 48000;
    const int signalLength = fs*2;
    const int nFFT = 1024;
    const int nBins = nFFT/2+1;

    /* Create and initialise an instance of array2sh for the Eigenmike32 */
    array2sh_create(&hA2sh);
    array2sh_init(hA2sh, fs); /* Cannot be called while "process" is on-going */
    array2sh_setPreset(hA2sh, MICROPHONE_ARRAY_PRESET_EIGENMIKE32);
    array2sh_setNormType(hA2sh, NORM_N3D);

    /* Define input mono signal */
    nSH = ORDER2NSH(order);
    inSig = malloc1d(signalLength*sizeof(float));
    rand_m1_1(inSig, signalLength); /* Mono white-noise signal */

    /* Simulate an Eigenmike in a free-field with a single plane-wave */
    f = malloc1d(nBins*sizeof(float));
    kr = malloc1d(nBins*sizeof(double));
    getUniformFreqVector(nFFT, (float)fs, f);
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
    for(i=0; i<(int)((float)signalLength/256.0f); i++)
        saf_multiConv_apply(hMC, FLATTEN2D(inSig_32), FLATTEN2D(micSig));

    /* Encode simulated Eigenmike signals into spherical harmonic signals */
    framesize = array2sh_getFrameSize();
    shSig = (float**)malloc2d(nSH,signalLength,sizeof(float));
    micSig_frame = (float**)malloc1d(32*sizeof(float*));
    shSig_frame = (float**)malloc1d(nSH*sizeof(float*));
    for(i=0; i<(int)((float)signalLength/(float)framesize); i++){
        for(ch=0; ch<32; ch++)
            micSig_frame[ch] = &micSig[ch][i*framesize];
        for(ch=0; ch<nSH; ch++)
            shSig_frame[ch] = &shSig[ch][i*framesize];

        array2sh_process(hA2sh, micSig_frame, shSig_frame, 32, nSH, framesize);
    }

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
    free(micSig_frame);
    free(shSig_frame);
}

void test__saf_example_rotator(void){
    int ch, nSH, i, j, delay, framesize;
    void* hRot;
    float direction_deg[2], ypr[3], Rzyx[3][3];
    float** inSig, *y, **shSig_frame, **shSig_rot_frame;
    float** shSig, **shSig_rot, **shSig_rot_ref, **Mrot;

    /* Config */
    const float acceptedTolerance = 0.000001f;
    const int order = 4;
    const int fs = 48000;
    const int signalLength = fs*2;
    direction_deg[0] = 90.0f; /* encode to loudspeaker direction: index 8 */
    direction_deg[1] = 0.0f;
    ypr[0] = -0.4f;
    ypr[1] = -1.4f;
    ypr[2] = 2.1f;
    delay = rotator_getProcessingDelay();

    /* Create and initialise an instance of rotator */
    rotator_create(&hRot);
    rotator_init(hRot, fs); /* Cannot be called while "process" is on-going */

    /* Configure rotator codec */
    rotator_setOrder(hRot, (SH_ORDERS)order);
    rotator_setNormType(hRot, NORM_N3D);
    rotator_setYaw(hRot, ypr[0]*180.0f/M_PI); /* rad->degrees */
    rotator_setPitch(hRot, ypr[1]*180.0f/M_PI);
    rotator_setRoll(hRot, ypr[2]*180.0f/M_PI);

    /* Define input mono signal */
    nSH = ORDER2NSH(order);
    inSig = (float**)malloc2d(1,signalLength,sizeof(float));
    shSig = (float**)malloc2d(nSH,signalLength,sizeof(float));
    rand_m1_1(FLATTEN2D(inSig), signalLength); /* Mono white-noise signal */

    /* Encode */
    y = malloc1d(nSH*sizeof(float));
    getRSH(order, (float*)direction_deg, 1, y); /* SH plane-wave weights */
    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nSH, signalLength, 1, 1.0f,
                y, 1,
                FLATTEN2D(inSig), signalLength, 0.0f,
                FLATTEN2D(shSig), signalLength);

    /* Rotated version reference */
    Mrot = (float**)malloc2d(nSH, nSH, sizeof(float));
    yawPitchRoll2Rzyx(ypr[0], ypr[1], ypr[2], 0, Rzyx);
    getSHrotMtxReal(Rzyx, FLATTEN2D(Mrot), order);
    shSig_rot_ref = (float**)malloc2d(nSH,signalLength,sizeof(float));
    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nSH, signalLength, nSH, 1.0f,
                FLATTEN2D(Mrot), nSH,
                FLATTEN2D(shSig), signalLength, 0.0f,
                FLATTEN2D(shSig_rot_ref), signalLength);

    /* Rotate with rotator */
    framesize = rotator_getFrameSize();
    shSig_rot = (float**)malloc2d(nSH,signalLength,sizeof(float));
    shSig_frame = (float**)malloc1d(nSH*sizeof(float*));
    shSig_rot_frame = (float**)malloc1d(nSH*sizeof(float*));
    for(i=0; i<(int)((float)signalLength/(float)framesize); i++){
        for(ch=0; ch<nSH; ch++)
            shSig_frame[ch] = &shSig[ch][i*framesize];
        for(ch=0; ch<nSH; ch++)
            shSig_rot_frame[ch] = &shSig_rot[ch][i*framesize];

        rotator_process(hRot, shSig_frame, shSig_rot_frame, nSH, nSH, framesize);
    }

    /* ambi_enc should be equivalent to the reference, except delayed due to the
     * temporal interpolation employed in ambi_enc */
    for(i=0; i<nSH; i++)
        for(j=0; j<signalLength-delay; j++)
            TEST_ASSERT_FLOAT_WITHIN(acceptedTolerance, shSig_rot_ref[i][j], shSig_rot[i][j+delay]);

    /* Clean-up */
    rotator_destroy(&hRot);
    free(inSig);
    free(shSig);
    free(shSig_rot_ref);
    free(Mrot);
    free(y);
    free(shSig_frame);
    free(shSig_rot_frame);
}
#endif /* SAF_ENABLE_EXAMPLES_TESTS */
