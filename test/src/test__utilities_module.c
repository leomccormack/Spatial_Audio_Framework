/*
 * Copyright 2020-2021 Leo McCormack
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
 * @file test__utilities_module.c
 * @brief Unit tests for the SAF utilities module
 * @author Leo McCormack, Michael McCrea
 * @date 27.04.2020
 * @license ISC
 */

#include "saf_test.h"

void test__cylindricalBesselFunctions(void){ // TODO: may as well check the derivatives too...
    int i;
    double J_n[10], Y_n[10];

    /* Config */
    const float acceptedTolerance = 0.00001f;
    int testOrder = 7; /* note, REF values hardcoded for order 7 */
    double z[10] = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0}; /* note, REF values hardcoded for these values */

    /* Reference values computed in MATLAB with:
     * J_n = besselj(N, z);
     * y_n = bessely(N, z);
     */
    double J_nREF[10] = {0.0, 1.50232581743681e-06, 0.000174944074868274, 0.00254729445180469, 0.0151760694220584, 0.0533764101558907, 0.129586651841481, 0.233583569505696, 0.320589077979826, 0.327460879242453};
    double Y_nREF[10] = {0.0, -30588.9570521240, -271.548025367994, -19.8399354089864, -3.70622393164077, -1.26289883576932, -0.656590825719075, -0.405371018606768, -0.200063904600409, 0.0172445799076681};

    /* test bessel_Jn */
    bessel_Jn(testOrder, z, 10, J_n, NULL);
    for(i=0; i<10; i++)
        TEST_ASSERT_TRUE(fabs(J_n[i]-J_nREF[i])<acceptedTolerance);

    /* test bessel_Yn */
    bessel_Yn(testOrder, z, 10, Y_n, NULL);
    for(i=0; i<10; i++)
        TEST_ASSERT_TRUE(fabs(Y_n[i]-Y_nREF[i])<acceptedTolerance);
}

void test__sphericalBesselFunctions(void){ // TODO: may as well check the derivatives too...
    int i, successFlag;
    double j_n[10], i_n[10], y_n[10], k_n[10];

    /* Config */
    const float acceptedTolerance = 0.00001f;
    int testOrder = 7; /* note, REF values hardcoded for order 7 */
    double z[10] = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0}; /* note, REF values hardcoded for these values */

    /* Reference values computed in MATLAB with:
     * j_n = sqrt(pi./(2*z)).*besselj(N+0.5, z);
     * i_n = sqrt(pi./(2*z)).*besseli(N+0.5, z);
     * y_n = sqrt(pi./(2*z)).*bessely(N+0.5, z);
     * k_n = sqrt(pi./(2*z)).*besselk(N+0.5, z);
     */
    double j_nREF[10] = {0.0, 4.79013419873948e-07, 5.60965570334894e-05, 0.000824843253217635, 0.00498650846172602, 0.0179027781779895, 0.0447223808293482, 0.0839226228445072, 0.122272711565833, 0.137946585027486};
    double i_nREF[10] = {0.0, 5.08036087257580e-07, 7.09794452304064e-05, 0.00140087680258227, 0.0127983365433790, 0.0783315436379810, 0.377879458299915, 1.56419501808402, 5.83626393050750, 20.2384754394417};
    double y_nREF[10] = {0.0, -140452.852366906, -617.054329642527, -29.4761692244538, -3.98778927238432, -1.02739463881260, -0.425887203702750, -0.237025274765842, -0.132622247946352, -0.0402143438632017};
    double k_nREF[10] = {0.0, 204287.522076393, 712.406907885478, 23.1153112578315, 1.80293583642309, 0.222213613092395, 0.0360276414091966, 0.00698538879470478, 0.00153285534574965, 0.000367847412220325};

    /* test bessel_jn */
    successFlag = bessel_jn(testOrder, z, 10, j_n, NULL);
    TEST_ASSERT_TRUE(successFlag);
    for(i=0; i<10; i++)
        TEST_ASSERT_TRUE(fabs(j_n[i]-j_nREF[i])<acceptedTolerance);

    /* test bessel_in */
    successFlag = bessel_in(testOrder, z, 10, i_n, NULL);
    TEST_ASSERT_TRUE(successFlag);
    for(i=0; i<10; i++)
        TEST_ASSERT_TRUE(fabs(i_n[i]-i_nREF[i])<acceptedTolerance);

    /* test bessel_yn */
    successFlag = bessel_yn(testOrder, z, 10, y_n, NULL);
    TEST_ASSERT_TRUE(successFlag);
    for(i=0; i<10; i++)
        TEST_ASSERT_TRUE(fabs(y_n[i]-y_nREF[i])<acceptedTolerance);

    /* test bessel_kn */
    successFlag = bessel_kn(testOrder, z, 10, k_n, NULL);
    TEST_ASSERT_TRUE(successFlag);
    for(i=0; i<10; i++)
        TEST_ASSERT_TRUE(fabs(k_n[i]-k_nREF[i])<acceptedTolerance);
}

void test__cart2sph(void){
    const float acceptedTolerance = 0.00001f;
    float cordCar [100][3];
    float cordSph [100][3];
    float cordCarTest [100][3];

    /* Generate some random Cartesian coordinates */
    rand_m1_1((float*) cordCar, 100*3);

    /* rad */
    cart2sph((float*) cordCar, 100, SAF_FALSE, (float*) cordSph);
    sph2cart((float*) cordSph, 100, SAF_FALSE, (float*) cordCarTest);
    for (int i=0; i<100; i++)
        for (int j=0; j<3; j++)
            TEST_ASSERT_FLOAT_WITHIN(acceptedTolerance, cordCar[i][j], cordCarTest[i][j]);

    /* deg */
    cart2sph((float*) cordCar, 100, SAF_TRUE, (float*) cordSph);
    sph2cart((float*) cordSph, 100, SAF_TRUE, (float*) cordCarTest);
    for (int i=0; i<100; i++)
        for (int j=0; j<3; j++)
            TEST_ASSERT_FLOAT_WITHIN(acceptedTolerance, cordCar[i][j], cordCarTest[i][j]);
}

void test__delaunaynd(void){
    int nMesh;
    int* mesh;

    /* Not really a unit test... You have to copy the mesh indices into e.g. Matlab, plot, and see... */

    /* 2D 3 points */
    float three_xy[3][2] =
      { {7.0f,7.0f},{2.0f,7.0f},{2.0f,1.0f} };
    mesh = NULL;
    delaunaynd((float*)three_xy, 3, 2/*nDims*/, &mesh, &nMesh);
    free(mesh);

    /* 2D 4 points */
    float four_xy[4][2] =
      { {7.0f,7.0f},{2.0f,7.0f},{2.0f,1.0f},{7.0f,1.0f} };
    mesh = NULL;
    delaunaynd((float*)four_xy, 4, 2/*nDims*/, &mesh, &nMesh);
    free(mesh);

    /* 2D Square */
    float square_xy[26][2] =
      { {-1.0f,-1.0f},{-1.0f,-0.5f},{-1.0f,0.0f},{-1.0f,0.5f},{-1.0f,1.0f},{-0.5f,-1.0f},{-0.5f,-0.5f},{-0.5f,0.0f},{-0.5f,0.5f},
        {-0.5f,1.0f},{0.0f,-1.0f},{0.0f,-0.5f},{0.0f,0.0f},{0.0f,0.5f},{0.0f,1.0f},{0.5f,-1.0f},
        {0.5f,-0.5f},{0.5f,0.0f},{0.5f,0.5f},{0.5f,1.0f},{1.0f,-1.0f},{1.0f,-0.5f},
        {1.0f,0.0f},{1.0f,0.5f},{1.0f,1.0f},{0.0f,0.0f} };
    mesh = NULL;
    delaunaynd((float*)square_xy, 26, 2/*nDims*/, &mesh, &nMesh);
    free(mesh);

    /* 3D Cube */
    float cube_xyz[8][3] =
      { {-1.0f,-1.0f,-1.0f},{-1.0f,1.0f,-1.0f},{1.0f,-1.0f,-1.0f},{1.0f,1.0f,-1.0f},
        {-1.0f,-1.0f,1.0f}, {-1.0f,1.0f,1.0f}, {1.0f,-1.0f,1.0f}, {1.0f,1.0f,1.0f} };
    mesh = NULL;
    delaunaynd((float*)cube_xyz, 8, 3/*nDims*/, &mesh, &nMesh);
    free(mesh);

    /* 3D Cube with a point in the centre */
    float cube_xyz2[9][3] =
      { {-1.0f,-1.0f,-1.0f},{-1.0f,1.0f,-1.0f},{1.0f,-1.0f,-1.0f},{1.0f,1.0f,-1.0f},
        {-1.0f,-1.0f,1.0f}, {-1.0f,1.0f,1.0f}, {1.0f,-1.0f,1.0f}, {1.0f,1.0f,1.0f}, {0.0f,0.0f,0.0f} };
    mesh = NULL;
    delaunaynd((float*)cube_xyz2, 9, 3/*nDims*/, &mesh, &nMesh);
    free(mesh);
}

void test__quaternion(void){
    int i, j;
    float norm;
    float rot[3][3], rot2[3][3], residual[9], ypr[3], test_ypr[3];
    quaternion_data Q, Q1, Q2;

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

        /* Testing that quaternion2euler() and euler2Quaternion() are invertable */
        quaternion2euler(&Q1, 1, EULER_ROTATION_YAW_PITCH_ROLL, &ypr[0], &ypr[1], &ypr[2]);
        euler2Quaternion(ypr[0], ypr[1], ypr[2], 1, EULER_ROTATION_YAW_PITCH_ROLL, &Q2);
        quaternion2euler(&Q2, 1, EULER_ROTATION_YAW_PITCH_ROLL, &test_ypr[0], &test_ypr[1], &test_ypr[2]);
        for(j=0; j<3; j++)
            TEST_ASSERT_TRUE(fabsf(test_ypr[j]-ypr[j])<1e-2f);
    }
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
                          nInputs, nOutputs, SAF_TRUE);

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
    const float acceptedTolerance = 0.00001f;
    const int fftSizesToTest[24] =
        {16,256,512,1024,2048,4096,8192,16384,32768,65536,1048576,     /*     2^x */
         80,160,320,640,1280,240,480,960,1920,3840,7680,15360,30720 }; /* non-2^x, (but still supported by vDSP) */

    /* Loop over the different FFT sizes */
    for (i=0; i<24; i++){
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

void test__saf_fft(void){
    int i, j, N;
    float_complex* x_td, *test;
    float_complex* x_fd;
    void *hFFT;

    /* Config */
    const float acceptedTolerance = 0.00001f;
    const int fftSizesToTest[24] =
        {16,256,512,1024,2048,4096,8192,16384,32768,65536,1048576,     /*     2^x */
         80,160,320,640,1280,240,480,960,1920,3840,7680,15360,30720 }; /* non-2^x, (but still supported by vDSP) */

    /* Loop over the different FFT sizes */
    for (i=0; i<24; i++){
        N = fftSizesToTest[i];

        /* prep */
        x_td = malloc1d(N*sizeof(float_complex));
        test = malloc1d(N*sizeof(float_complex));
        x_fd = malloc1d(N*sizeof(float_complex));
        rand_m1_1((float*)x_td, N*2); /* populate with random numbers */
        saf_fft_create(&hFFT, N);

        /* forward and backward transform */
        saf_fft_forward(hFFT, x_td, x_fd);
        saf_fft_backward(hFFT, x_fd, test);

        /* Check that, x_td==test */
        for(j=0; j<N; j++){
            TEST_ASSERT_FLOAT_WITHIN(acceptedTolerance, crealf(x_td[j]), crealf(test[j]));
            TEST_ASSERT_FLOAT_WITHIN(acceptedTolerance, cimagf(x_td[j]), cimagf(test[j]));
        }

        /* clean-up */
        saf_fft_destroy(&hFFT);
        free(x_fd);
        free(x_td);
        free(test);
    }
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
        inputData[i] = sinf(2.0f * SAF_PI * (float)i * frequency/(float)sampleRate);
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

void test__latticeDecorrelator(void){
    int c, band, nBands, idx, hopIdx, i;
    void* hDecor, *hSTFT;
    float icc, tmp, tmp2;
    float* freqVector;
    float** inputTimeDomainData, **outputTimeDomainData, **tempHop;
    float_complex*** inTFframe, ***outTFframe;

    /* config */
    const float acceptedICC = 0.05f;
    const int nCH = 24;
    const int nTestHops = 800;
    const int hopSize = 128;
    const int procDelay = hopSize*12 + 12;
    const int lSig = nTestHops*hopSize+procDelay;
    const float fs = 48e3f;
    nBands = hopSize+5;

    /* audio buffers */
    inputTimeDomainData = (float**) calloc2d(1, lSig, sizeof(float));
    outputTimeDomainData = (float**) calloc2d(nCH, lSig, sizeof(float));
    inTFframe = (float_complex***)malloc3d(nBands, nCH, 1, sizeof(float_complex));
    outTFframe = (float_complex***)malloc3d(nBands, nCH, 1, sizeof(float_complex));
    tempHop = (float**) malloc2d(nCH, hopSize, sizeof(float));

    /* Initialise afSTFT and input data */
    afSTFT_create(&hSTFT, 1, nCH, hopSize, 0, 1, AFSTFT_BANDS_CH_TIME);
    rand_m1_1(FLATTEN2D(inputTimeDomainData), 1*lSig); /* populate with random numbers */
    freqVector = malloc1d(nBands*sizeof(float));
    afSTFT_getCentreFreqs(hSTFT, fs, nBands, freqVector);

    /* setup decorrelator */
    int orders[4] = {20, 15, 6, 6}; /* 20th order up to 700Hz, 15th->2.4kHz, 6th->4kHz, 3rd->12kHz, NONE(only delays)->Nyquist */
    //float freqCutoffs[4] = {600.0f, 2.6e3f, 4.5e3f, 12e3f};
    float freqCutoffs[4] = {900.0f, 6.8e3f, 12e3f, 24e3f};
    const int maxDelay = 12;
    latticeDecorrelator_create(&hDecor, fs, hopSize, freqVector, nBands, nCH, orders, freqCutoffs, 4, maxDelay, 0, 0.75f);

    /* Processing loop */
    idx = 0;
    hopIdx = 0;
    while(idx<lSig-hopSize*2){
        for(c=0; c<1; c++)
            memcpy(tempHop[c], &(inputTimeDomainData[c][hopIdx*hopSize]), hopSize*sizeof(float));

        /* forward TF transform, and replicate to all channels */
        afSTFT_forward(hSTFT, tempHop, hopSize, inTFframe);
        for(band=0; band<nBands; band++)
            for(i=1; i<nCH;i++)
                inTFframe[band][i][0] = inTFframe[band][0][0];

        /* decorrelate */
        latticeDecorrelator_apply(hDecor, inTFframe, 1, outTFframe);

        /*  backward TF transform */
        afSTFT_backward(hSTFT, outTFframe, hopSize, tempHop);

        /* Copy frame to output TD buffer */
        for(c=0; c<nCH; c++)
            memcpy(&(outputTimeDomainData[c][hopIdx*hopSize]), tempHop[c], hopSize*sizeof(float));
        idx+=hopSize;
        hopIdx++;
    }

    /* Compensate for processing delay, and check that the inter-channel correlation
     * coefficient is below the accepted threshold (ideally 0, if fully
     * decorrelated...) */
    for(c=0; c<nCH; c++){
        utility_svvdot(inputTimeDomainData[0], &outputTimeDomainData[c][procDelay], (lSig-procDelay), &icc);
        utility_svvdot(inputTimeDomainData[0], inputTimeDomainData[0], (lSig-procDelay), &tmp);
        utility_svvdot(&outputTimeDomainData[c][procDelay], &outputTimeDomainData[c][procDelay], (lSig-procDelay), &tmp2);

        icc = icc/sqrtf(tmp*tmp2); /* normalise */
        TEST_ASSERT_TRUE(fabsf(icc)<acceptedICC);
    }
#if 0
    /* Check for mutually decorrelated channels... */
    int c2;
    for(c=0; c<nCH; c++){
        for(c2=0; c2<nCH; c2++){
            utility_svvdot(&outputTimeDomainData[c][procDelay], &outputTimeDomainData[c2][procDelay], (lSig-procDelay), &icc);
            utility_svvdot(&outputTimeDomainData[c2][procDelay], &outputTimeDomainData[c2][procDelay], (lSig-procDelay), &tmp);
            utility_svvdot(&outputTimeDomainData[c][procDelay], &outputTimeDomainData[c][procDelay], (lSig-procDelay), &tmp2);

            if (c!=c2){
                icc = icc/sqrtf(tmp*tmp2); /* normalise */
                TEST_ASSERT_TRUE(fabsf(icc)<acceptedICC);
            }
        }
    }
#endif

    /* Clean-up */
    latticeDecorrelator_destroy(&hDecor);
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

/* Test evalIIRTransferFunction()
 * Coefficients for the first 7 tests below are taken from the butterworth tests
 * above, so can be compared against the mag/phase results given by MATLAB's
 * freqz function. The 8th test loop evaluates results for 1st order shelving
 * filters with coefficients generated by the DVF filter functions. The last
 * test loop runs those same tests, but with the floating point version of
 * evalIIRTransferFunctionf() (valid for low order filters)
 */
void test__evalIIRTransferFunction(void) {
    int i, nCoeffs;
    float phase_ref, mag_ref, mag_ref_db;
    
    /* Config */
    const double phaseTolerance = 0.0174533f * 5.f; /* ~ 1 degree * mul */
    const double magToleranceDb = 0.1f; /* tolerance in dB, for a target magnitude of 0dB */
    const float errScale = 2.f / 120.f; /* tolerance grows for lower dB target: toleranceLevel/atLevel.
                                         * e.g. 2/120 = 2dB tolerance for -120 target dB */
    const int nFreqs = 10;
    const float freqs[10] = {147.21423,270.49564,411.40091,687.90202,1395.3072,2024.3936,3696.9416,6784.4745,9798.67,17594.058};
    float mag[10], phase[10];
    float fs = 44.1e3f;
    
    /* evalIIRTransferFunction(): coeffs of type double */
    
    /* Test 1 * 1st order Low-pass filter */
    nCoeffs = 2;
    const double a_t1[2] = {1,-0.6681786};
    const double b_t1[2] = {0.1659107,0.1659107};
    const double mag_ref1[10] = {0.99861294,0.99533929,0.98931332,0.97092312,0.89393904,0.80765661,0.59366549,0.35440473,0.23070415,0.06521195};
    const double phase_ref1[10] = {-0.052676048,-0.096585083,-0.14632679,-0.24173907,-0.46473796,-0.63062928,-0.93519006,-1.2085189,-1.337995,-1.5055381};
    evalIIRTransferFunction((double*)b_t1, (double*)a_t1, nCoeffs, (float*)freqs, nFreqs, fs, 0, (float*)mag, (float*)phase);
    for(i = 0; i < nFreqs; i++) {
        phase_ref = phase_ref1[i];
        mag_ref = mag_ref1[i];
        mag_ref_db = 20*log10(mag_ref);
        TEST_ASSERT_FLOAT_WITHIN(magToleranceDb + errScale*fabsf(mag_ref_db), mag_ref_db, 20*log10(mag[i]));
        TEST_ASSERT_FLOAT_WITHIN(phaseTolerance, phase_ref, phase[i]);
    }
    
    /* Test 2 * 2nd order Low-pass filter */
    nCoeffs = 3;
    const double a_t2[3] = {1,-0,0.1715729};
    const double b_t2[3] = {0.2928932,0.5857864,0.2928932};
    const float mag_ref2[10] = {0.99999991,0.99999985,0.99999955,0.99999702,0.99995046,0.99977761,0.99736787,0.96409579,0.81776268,0.1073164};
    const float phase_ref2[10] = {-0.014832279,-0.027258003,-0.041470589,-0.069414188,-0.14150066,-0.20679977,-0.39012494,-0.7974393,-1.3261562,-2.6614056};
    evalIIRTransferFunction((double*)b_t2, (double*)a_t2, nCoeffs, (float*)freqs, nFreqs, fs, 0, (float*)mag, (float*)phase);
    for(i = 0; i < nFreqs; i++) {
        phase_ref = phase_ref2[i];
        mag_ref = mag_ref2[i];
        mag_ref_db = 20*log10(mag_ref);
        TEST_ASSERT_FLOAT_WITHIN(magToleranceDb + errScale*fabsf(mag_ref_db), mag_ref_db, 20*log10(mag[i]));
        TEST_ASSERT_FLOAT_WITHIN(phaseTolerance, phase_ref, phase[i]);
    }

    /* Test 3 * 3rd order Low-pass filter */
    nCoeffs = 4;
    const double a_t3[4] = {1,-2.9476416,2.896645,-0.9489859};
    const double b_t3[4] = {2.2e-06,6.6e-06,6.6e-06,2.2e-06};
    const float mag_ref3[10] = {0.8954383,0.3011618,0.0892913,0.0191409,0.0022769,0.0007374,0.0001152,1.56e-05,3.8e-06,1e-07};
    const float phase_ref3[10] = {-1.8249735,3.0678618,2.4995092,2.1114704,1.8340934,1.751328,1.6679375,1.6206872,1.6020054,1.579398};
    evalIIRTransferFunction((double*)b_t3, (double*)a_t3, nCoeffs, (float*)freqs, nFreqs, fs, 0, (float*)mag, (float*)phase);
    for(i = 0; i < nFreqs; i++) {
        phase_ref = phase_ref3[i];
        mag_ref = mag_ref3[i];
        mag_ref_db = 20*log10(mag_ref);
        TEST_ASSERT_FLOAT_WITHIN(magToleranceDb + errScale*fabsf(mag_ref_db), mag_ref_db, 20*log10(mag[i]));
        TEST_ASSERT_FLOAT_WITHIN(phaseTolerance, phase_ref, phase[i]);
    }
    
    /* Test 4 * 6th order Low-pass filter */
    nCoeffs = 7;
    const double a_t4[7] = {1,-5.49431292177096,12.5978414666894,-15.4285267903275,10.6436770055305,-3.92144696766748,0.6027721469713};
    const double b_t4[7] = {6.15535184628202e-08,3.69321110776921e-07,9.23302776942303e-07,1.2310703692564e-06,9.23302776942303e-07,3.69321110776921e-07,6.15535184628202e-08};
    const float mag_ref4[10] = {0.9999999834907868,0.9999997831836054,0.9999679556869572,0.9849426248859378,0.08033081621985783,0.008452216914022819,0.0002063542729228268,3.793812554381118e-06,2.274031694371124e-07,9.970589432354785e-11};
    const float phase_ref4[10] = {-0.6201852189230334,-1.148525513374147,-1.774695369143539,3.109543344373707,-0.4296773811384472,-1.349824316530828,-2.195405632723407,-2.65814688739603,-2.839508904295157,-3.058387834019209};
    evalIIRTransferFunction((double*)b_t4, (double*)a_t4, nCoeffs, (float*)freqs, nFreqs, fs, 0, (float*)mag, (float*)phase);
    for(i = 0; i < nFreqs; i++) {
        float phase_ref = phase_ref4[i];
        float mag_ref = mag_ref4[i];
        float mag_ref_db = 20*log10(mag_ref);
        TEST_ASSERT_FLOAT_WITHIN(magToleranceDb + errScale*fabsf(mag_ref_db), mag_ref_db, 20*log10(mag[i]));
        TEST_ASSERT_FLOAT_WITHIN(phaseTolerance, phase_ref, phase[i]);
    }
    
    /* Test 5 * 3rd order High-pass filter */
    nCoeffs = 4;
    const double a_t5[4] = {1,-2.2191686,1.7151178,-0.4535459};
    const double b_t5[4] = {0.673479,-2.0204371,2.0204371,-0.673479};
    const float mag_ref5[10] = {0.0001466,0.0009096,0.0032014,0.0149875,0.125037,0.362653,0.927991,0.9985214,0.9999112,0.9999999};
    const float phase_ref5[10] = {-1.6762949,-1.7648759,-1.866651,-2.0692621,-2.6256366,3.0800183,1.6530258,0.7789431,0.4789307,0.1307956};
    evalIIRTransferFunction((double*)b_t5, (double*)a_t5, nCoeffs, (float*)freqs, nFreqs, fs, 0, (float*)mag, (float*)phase);
    for(i = 0; i < nFreqs; i++) {
        float phase_ref = phase_ref5[i];
        float mag_ref = mag_ref5[i];
        float mag_ref_db = 20*log10(mag_ref);
        TEST_ASSERT_FLOAT_WITHIN(magToleranceDb + errScale*fabsf(mag_ref_db), mag_ref_db, 20*log10(mag[i]));
        TEST_ASSERT_FLOAT_WITHIN(phaseTolerance, phase_ref, phase[i]);
    }
    
    /* Test 6 * 4th order High-pass filter
     * 400 Hz cut (differs from butterworth test above) */
    nCoeffs = 5;
    const double a_t6[5] = {1,-3.863184622426,5.598835456747838,-3.607752453919942,0.872108645089876};
    const double b_t6[5] = {4.3909323578772e-07,1.75637294315089e-06,2.63455941472633e-06,1.75637294315089e-06,4.3909323578772e-07};
    const float mag_ref6[10] = {0.9996691528983467,0.9595570109649983,0.5370184819357747,0.08100263003740536,0.004753436194609436,0.001057169058757887,8.896712774518116e-05,6.197328265811134e-06,9.491865964914827e-07,5.478157027512644e-09};
    const float phase_ref6[10] = {-1.072517623166929,-2.13344694428915,2.732267641095127,1.462991201859678,0.6929733816699927,0.4733493046075806,0.2541184532330854,0.130425028023503,0.08157492611996242,0.02248140228360206};
    evalIIRTransferFunction((double*)b_t6, (double*)a_t6, nCoeffs, (float*)freqs, nFreqs, fs, 0, (float*)mag, (float*)phase);
    for(i = 0; i < nFreqs; i++) {
        float phase_ref = phase_ref6[i];
        float mag_ref = mag_ref6[i];
        float mag_ref_db = 20*log10(mag_ref);
        TEST_ASSERT_FLOAT_WITHIN(magToleranceDb + errScale*fabsf(mag_ref_db), mag_ref_db, 20*log10(mag[i]));
        TEST_ASSERT_FLOAT_WITHIN(phaseTolerance, phase_ref, phase[i]);
    }
    
    /* Test 7 * 2nd order Band-pass filter */
    nCoeffs = 5;
    const double a_t7[5] = {1,-3.9431258,5.832267,-3.8351187,0.9459779};
    const double b_t7[5] = {0.0003751,0,-0.0007501,0,0.0003751};
    const float mag_ref7[10] = {0.7829909,0.9051549,0.5636772,0.1816557,0.0400635,0.0185759,0.0053305,0.0014022,0.0005484,4.16e-05};
    const float phase_ref7[10] = {0.4017825,-0.7852502,-1.8127451,-2.4983166,-2.8544848,-2.9475768,-3.0381483,-3.0886103,-3.1084696,-3.1324667};
    evalIIRTransferFunction((double*)b_t7, (double*)a_t7, nCoeffs, (float*)freqs, nFreqs, fs, 0, (float*)mag, (float*)phase);    for(i=0; i<nFreqs; i++){
        float phase_ref = phase_ref7[i];
        float mag_ref = mag_ref7[i];
        float mag_ref_db = 20*log10(mag_ref);
        TEST_ASSERT_FLOAT_WITHIN(magToleranceDb + errScale*fabsf(mag_ref_db), mag_ref_db, 20*log10(mag[i]));
        TEST_ASSERT_FLOAT_WITHIN(phaseTolerance, phase_ref, phase[i]);
    }
    
    /* Test the response for 12 settings used by the DVF filters
     * (1st order sheving filter) */
    const int nTest = 12;
    nCoeffs = 2;
    double as_dvf[12][2] = { /* nTest x nCoeffs */
        {1,-0.95864619},{1,-0.96599375},{1,-0.9648155},{1,-0.84969467},{1,-0.93327999},{1,-0.95974372},
        {1,-0.83460338},{1,-0.74744027},{1,-0.67445272},{1,-0.76911048},{1,-0.64266857},{1,-0.54043336}
    };
    double bs_dvf[12][2] = { /* nTest x nCoeffs */
        {8.0841171,-7.5217374},{2.1404179,-2.0439339},{1.3371466,-1.2837154},{0.73353449,-0.52570938},
        {1.1312827,-1.0324739},{1.1334695,-1.0826754},{0.18784397,-0.098191093},{0.43493823,-0.24012268},
        {0.72850398,-0.42469436},{0.1577158,-0.075691788},{0.34545179,-0.16259809},{0.60618525,-0.27852873}
    };
    const float mags_dvf[12][10] = { /* nTest x nFreqs */
        {12.68472,11.390262,10.245287,9.081507,8.2881923,8.1239126,8.0139907,7.9799796,7.9724969,7.9680408},
        {2.6652868,2.469902,2.3321513,2.2179604,2.1523794,2.139901,2.1317514,2.1292619,2.1287161,2.1283915},
        {1.473629,1.4224625,1.3865239,1.35693,1.3400517,1.3368519,1.3347643,1.3341269,1.3339872,1.3339041},
        {1.3740782,1.3545086,1.3209933,1.2349497,1.0204791,0.89934391,0.76424296,0.70524874,0.69060103,0.68154193},
        {1.4538524,1.4034324,1.3412611,1.2506542,1.1632856,1.1414561,1.125968,1.1210229,1.1199249,1.1192693},
        {1.2358328,1.2022487,1.1755623,1.151321,1.1364575,1.1335478,1.1316331,1.1310458,1.130917,1.1308403},
        {0.53871826,0.53107297,0.51772931,0.48194542,0.3814555,0.3150888,0.22673701,0.17897735,0.16548424,0.15666687},
        {0.76984932,0.76629998,0.75986062,0.74092839,0.6717326,0.60914565,0.49873472,0.42504312,0.40260097,0.38760994},
        {0.93261062,0.93115748,0.92849149,0.92042692,0.88786149,0.85374491,0.78101307,0.72261208,0.7032242,0.68987066},
        {0.3542684,0.3519694,0.34781956,0.33577027,0.29342253,0.25693916,0.1950882,0.15408756,0.14134236,0.13268952},
        {0.51134341,0.51045389,0.50881513,0.50380352,0.48269293,0.45893406,0.4014785,0.34644275,0.32577048,0.31064565},
        {0.71281398,0.71244818,0.7117705,0.70966741,0.70027544,0.68857561,0.65428371,0.61109598,0.59151358,0.57580457}
    };
    const float phases_dvf[12][10] = { /* nTest x nFreqs */
        {-0.17782001,-0.24874011,-0.2637155,-0.22716735,-0.13811864,-0.098858451,-0.054719219,-0.028348151,-0.01776687,-0.0049023852},
        {-0.11818797,-0.14315719,-0.13342746,-0.10042039,-0.055483624,-0.038914731,-0.021247104,-0.010960723,-0.006863078,-0.001892663},
        {-0.054685172,-0.064793725,-0.059294582,-0.043867626,-0.023978567,-0.016782288,-0.0091501707,-0.0047182622,-0.0029540703,-0.00081461202},
        {-0.064892969,-0.11661773,-0.17043929,-0.25417001,-0.34353751,-0.33903683,-0.25655254,-0.15106232,-0.097685198,-0.027478557},
        {-0.069273584,-0.10993957,-0.13346817,-0.13654806,-0.096248577,-0.071345378,-0.040468966,-0.0211288,-0.013264997,-0.0036639447},
        {-0.042920762,-0.054369797,-0.052367865,-0.040534678,-0.022767208,-0.01601877,-0.0087641883,-0.0045240297,-0.0028331222,-0.00078136769},
        {-0.082361693,-0.14918816,-0.22112994,-0.34301241,-0.52772513,-0.5813414,-0.53771111,-0.3682489,-0.25049777,-0.073003569},
        {-0.036111073,-0.06587829,-0.098882791,-0.15880356,-0.2712657,-0.32156729,-0.32729224,-0.23402028,-0.16071005,-0.047082247},
        {-0.014103031,-0.025780252,-0.03883755,-0.063048578,-0.11207771,-0.13775677,-0.14908283,-0.11045335,-0.07655027,-0.022550556},
        {-0.050350716,-0.09181969,-0.13772611,-0.22079247,-0.37596653,-0.44682619,-0.46554397,-0.34629266,-0.24235026,-0.072094835},
        {-0.019043671,-0.03486822,-0.052685779,-0.086317621,-0.15955383,-0.2051296,-0.24897791,-0.20818672,-0.15155625,-0.046353162},
        {-0.006828925,-0.012518705,-0.018958244,-0.031275577,-0.059563883,-0.079317458,-0.10560439,-0.097624085,-0.0740598,-0.023377524}
    };
    for(int t = 0; t < nTest; t++) {
        float phase_ref, mag_ref, mag_ref_db;
        double* pAs = &as_dvf[t][0];
        double* pBs = &bs_dvf[t][0];
        evalIIRTransferFunction(pBs, pAs, nCoeffs, (float*)freqs, nFreqs, fs, 0, (float*)mag, (float*)phase);
        for(i = 0; i < nFreqs; i++){
            phase_ref = phases_dvf[t][i];
            mag_ref = mags_dvf[t][i];
            mag_ref_db = 20*log10(mag_ref);
            TEST_ASSERT_FLOAT_WITHIN(magToleranceDb + errScale*fabsf(mag_ref_db), mag_ref_db, 20*log10(mag[i]));
            TEST_ASSERT_FLOAT_WITHIN(phaseTolerance, phase_ref, phase[i]);
        }
    }
    /* using the same parameters as above, test evalIIRTransferFunctionf():
     * coeffs of type float */
    float as_dvf_f[12][2] = { /* nTest x nCoeffs */
        {1,-0.95864619},{1,-0.96599375},{1,-0.9648155},{1,-0.84969467},{1,-0.93327999},{1,-0.95974372},
        {1,-0.83460338},{1,-0.74744027},{1,-0.67445272},{1,-0.76911048},{1,-0.64266857},{1,-0.54043336}
    };
    float bs_dvf_f[12][2] = { /* nTest x nCoeffs */
        {8.0841171,-7.5217374},{2.1404179,-2.0439339},{1.3371466,-1.2837154},{0.73353449,-0.52570938},
        {1.1312827,-1.0324739},{1.1334695,-1.0826754},{0.18784397,-0.098191093},{0.43493823,-0.24012268},
        {0.72850398,-0.42469436},{0.1577158,-0.075691788},{0.34545179,-0.16259809},{0.60618525,-0.27852873}
    };
    for(int t = 0; t < nTest; t++) {
        float* pAs = &as_dvf_f[t][0];
        float* pBs = &bs_dvf_f[t][0];
        evalIIRTransferFunctionf(pBs, pAs, nCoeffs, (float*)freqs, nFreqs, fs, 0, (float*)mag, (float*)phase);
        for(i = 0; i < nFreqs; i++){
            phase_ref = phases_dvf[t][i];
            mag_ref = mags_dvf[t][i];
            mag_ref_db = 20*log10(mag_ref);
            TEST_ASSERT_FLOAT_WITHIN(magToleranceDb + errScale*fabsf(mag_ref_db), mag_ref_db, 20*log10(mag[i]));
            TEST_ASSERT_FLOAT_WITHIN(phaseTolerance, phase_ref, phase[i]);
        }
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
    const int frameSize = 256;//16;
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
    faf_IIRFilterbank_create(&hFaF, order, (float*)fc, 6, fs, frameSize);
    for(i=0; i< signalLength/frameSize; i++){
        faf_IIRFilterbank_apply(hFaF, &inSig[i*frameSize], outFrame, frameSize);
        for(band=0; band<7; band++)
            memcpy(&outSig_bands[band][i*frameSize], outFrame[band], frameSize*sizeof(float));
    }
    faf_IIRFilterbank_destroy(&hFaF);

    /* Sum the individual bands */
    for(band=0; band<7; band++)
        cblas_saxpy(signalLength, 1.0f, outSig_bands[band], 1, outSig, 1);

    /* Check that the magnitude difference between input and output is below 0.5dB */
    saf_rfft_create(&hFFT, signalLength);
    saf_rfft_forward(hFFT, inSig, insig_fft);
    saf_rfft_forward(hFFT, outSig, outsig_fft);
    for(i=0; i<signalLength/2+1; i++)
        TEST_ASSERT_FLOAT_WITHIN(acceptedTolerance_dB, 0.0f, 20.0f * log10f(cabsf( ccdivf(outsig_fft[i],insig_fft[i]) )));

    /* Now the same thing, but for 1st order */
    order = 1;
    faf_IIRFilterbank_create(&hFaF, order, (float*)fc, 6, fs, frameSize);
    for(i=0; i< signalLength/frameSize; i++){
        faf_IIRFilterbank_apply(hFaF, &inSig[i*frameSize], outFrame, frameSize);
        for(band=0; band<7; band++)
            memcpy(&outSig_bands[band][i*frameSize], outFrame[band], frameSize*sizeof(float));
    }
    faf_IIRFilterbank_destroy(&hFaF);
    memset(outSig, 0, signalLength*sizeof(float));
    for(band=0; band<7; band++)
        cblas_saxpy(signalLength, 1.0f, outSig_bands[band], 1, outSig, 1);
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
        {-0.376858200853762f,0.656790634216694f,0.124479178614046f,-0.334752428307223f,1.50745241578235f,0.0290651989052969f},
        {0.608382058262806f,0.581930485432986f,3.23135406998058f,-0.712003744668929f,-1.33872571354702f,-0.334742482743222f},
        {-0.795741418256672f,0.690709474622409f,0.620971281129248f,1.38749471231620f,0.897245329198841f,-0.0693670166113321f},
        {0.179789913109994f,-1.06135084902804f,-1.10032635271188f,0.612441344250358f,-2.43213807790664f,-0.479265889956047f},
        {-0.277441781278754f,-0.0732116130293688f,-0.572551795688137f,1.02024767389969f,0.167385894565923f,1.45210312619277f},
        {-0.205305770089918f,-1.59783032780633f,1.08539265129120f,0.460057585947626f,-1.02420974042838f,1.04117461500218f}
    };
    const float outM_ref[6][6] = {
        {0.385163650730121f,0.0865151585709784f,0.898406722231524f,0.877640791713973f,0.435244824708340f,0.888866982998854f},
        {-0.664938511314777f,5.02943129352875f,8.24444951891833f,2.23840978101979f,-0.942669833528886f,-2.38535530623266f},
        {-0.388189314743059f,0.429308537172675f,1.13870842882926f,1.60875776611798f,-1.44249911796405f,-1.51822150286392f},
        {1.05630187656688f,0.256606570814868f,-2.42701873560847f,-1.42372526577009f,-0.335273289873574f,-1.94362909671742f},
        {0.0261470437116839f,-3.03329326250434f,-3.50207776203591f,0.412043775125377f,-0.536000387729306f,1.61801775548557f},
        {-0.292024827617294f,-4.31537192033477f,-3.99160103133879f,0.312499067924889f,-1.46924802440347f,1.98522802303672f}
    };

    /* Compute matrix exponential */
    gexpm((float*)inM, 6, 0, (float*)outM);

    /* Check that output of SAF's gexpm, is similar to Matlab's expm: */
    for(i=0; i<6; i++)
        for(j=0; j<6; j++)
            TEST_ASSERT_TRUE( fabsf(outM[i][j] - outM_ref[i][j]) <= acceptedTolerance );
}

/* Target values are generated by MATLAB functions in `generate_coeffs_for_plugin_tests.m`
 * (corresponding function names are noted above each data set), which is not included in this
 * repository but are in the author's development repository `nearfield_rangeextrapolation`.
 */
void test__dvf_calcDVFShelfParams(void){

    /* setup */
    const float acceptedTolerance = 0.00001f;
    const float acceptedTolerance_fc = 0.1f;
    const int nTheta = 19;
    const int nRho = 5;
    const float rho[5] = {1.150000f,1.250000f,1.570000f,2.381000f,3.990000f};
    const float g0_ref[5][19] = { /* testRhoCoeffs_g_0 */
        {22.670282f,17.717752f,11.902597f,7.906282f,4.720884f,2.217061f,0.134088f,-1.613730f,-3.095960f,-5.279052f,-5.433653f,-6.342905f,-7.107677f,-7.744796f,-8.236084f,-8.613662f,-8.864276f,-9.065870f,-9.089385f},
        {18.295933f,15.510441f,11.452312f,7.951083f,4.997235f,2.609075f,0.592613f,-1.107964f,-2.557504f,-4.547256f,-4.853912f,-5.750024f,-6.504702f,-7.133244f,-7.621092f,-7.993574f,-8.244015f,-8.438287f,-8.467470f},
        {11.937032f,11.093339f,9.245757f,7.118216f,4.990070f,3.083402f,1.371444f,-0.121838f,-1.427296f,-2.979742f,-3.542803f,-4.381065f,-5.091220f,-5.683427f,-6.149122f,-6.508598f,-6.748356f,-6.923465f,-6.961620f},
        {6.676990f,6.451424f,5.818377f,4.924700f,3.861979f,2.760683f,1.662668f,0.629080f,-0.327831f,-1.328149f,-1.970549f,-2.649238f,-3.234743f,-3.727775f,-4.122829f,-4.433178f,-4.640902f,-4.783351f,-4.823625f},
        {3.628860f,3.534311f,3.298166f,2.922799f,2.438587f,1.888286f,1.296135f,0.698518f,0.112899f,-0.473626f,-0.960644f,-1.428032f,-1.841763f,-2.196404f,-2.487131f,-2.717121f,-2.873915f,-2.978235f,-3.010937f}
    };
    const float gInf_ref[5][19] = { /* testRhoCoeffs_g_inf */
        {-4.643651f,-4.225287f,-4.134752f,-4.386332f,-5.244711f,-6.439307f,-7.659091f,-8.887172f,-10.004796f,-10.694171f,-11.190476f,-10.876569f,-10.140292f,-9.913242f,-9.411469f,-8.981807f,-8.723677f,-8.529900f,-8.574359f},
        {-4.128221f,-3.834507f,-3.575000f,-3.637788f,-4.278932f,-5.310000f,-6.609705f,-7.815000f,-8.925450f,-9.646588f,-10.000000f,-9.784733f,-9.301643f,-8.862963f,-8.370815f,-7.953778f,-7.693305f,-7.500645f,-7.518260f},
        {-3.094135f,-2.963709f,-2.721834f,-2.573043f,-2.793627f,-3.414652f,-4.403297f,-5.518539f,-6.578461f,-7.332562f,-7.702192f,-7.582977f,-7.376856f,-6.745349f,-6.279895f,-5.891862f,-5.636418f,-5.456323f,-5.437006f},
        {-1.937289f,-1.889079f,-1.765709f,-1.620800f,-1.598110f,-1.815613f,-2.314443f,-3.041183f,-3.857777f,-4.533446f,-4.931544f,-4.962571f,-4.717069f,-4.357935f,-3.971281f,-3.646312f,-3.422461f,-3.269044f,-3.231471f},
        {-1.126412f,-1.103440f,-1.049199f,-0.969714f,-0.917898f,-0.962176f,-1.182409f,-1.566237f,-2.065834f,-2.552771f,-2.884909f,-2.977707f,-2.811758f,-2.629199f,-2.355800f,-2.118920f,-1.949860f,-1.834291f,-1.800638f}
    };
    const float fc_ref[5][19] = { /* testRhoCoeffs_f_c */
        {525.636204f,409.078426f,427.552571f,936.671283f,1635.128987f,2622.394035f,3167.199181f,3899.649293f,4176.703569f,4361.226917f,4634.448076f,4516.401848f,4567.834168f,4685.234222f,4908.786495f,4966.258562f,4936.982049f,4927.963688f,5210.861482f},
        {410.072475f,389.319631f,398.844102f,613.394238f,1116.223303f,2095.651724f,2847.557763f,3726.141143f,4080.406901f,4304.960791f,4463.911798f,4449.375495f,4501.166349f,4623.582375f,4757.884246f,4911.093999f,4935.074404f,4940.266143f,5155.085794f},
        {358.247441f,352.931439f,352.752741f,402.566754f,590.733021f,1127.131294f,2007.589994f,3160.896502f,3808.131027f,4155.246718f,4336.853155f,4375.553567f,4406.656373f,4543.636509f,4649.878140f,4849.374583f,4974.986343f,5006.214905f,5164.504029f},
        {318.842699f,318.199637f,315.776327f,326.423309f,364.498469f,500.548368f,980.626776f,2174.301881f,3296.777215f,3904.656864f,4203.152454f,4329.347194f,4338.652755f,4492.976051f,4579.879128f,4849.678327f,5052.801340f,5116.753611f,5267.402018f},
        {297.454930f,297.570719f,296.701047f,300.362959f,308.255747f,342.596563f,509.934390f,1379.970914f,2702.827191f,3646.599635f,4078.866661f,4301.570222f,4303.807248f,4472.223890f,4535.654099f,4855.399825f,5119.558700f,5210.329993f,5380.750972f}
    };
    /* params to test */
    float g0, gInf, fc;

    for(int ri = 0; ri < nRho; ri++){
        float r = rho[ri];
        for(int ti = 0; ti < nTheta; ti++){
            calcDVFShelfParams(ti, r, &g0, &gInf, &fc);
            TEST_ASSERT_FLOAT_WITHIN(acceptedTolerance, g0_ref[ri][ti], g0);
            TEST_ASSERT_FLOAT_WITHIN(acceptedTolerance, gInf_ref[ri][ti], gInf);
            TEST_ASSERT_FLOAT_WITHIN(acceptedTolerance_fc, fc_ref[ri][ti], fc);
        }
    }
}

/* Parameter interpolation is implicitly checked in test__dvf_dvfShelfCoeffs. */
void test__dvf_interpDVFShelfParams(void){

    /* interpDVFShelfParams() calls calcDVFShelfParams() twice to generate the
     * high shelf parameters for the nearests angles in the lookup table. Those
     * parameters are subsequently interpolated. So the success of
     * interpHighShelfCoeffs() relies first on the calcDVFShelfParams(),
     * so that should be tested first. */

    /* setup */
    const float acceptedTolerance = 0.0001f;
    const float acceptedTolerance_frq = 0.01f;
    const int   nTheta = 6;
    const float theta[6] = {0.000000f,2.300000f,47.614000f,98.600000f,166.200000f,180.000000f};
    const int   nRho = 5;
    const float rho[5] = {1.150000f,1.250000f,1.570000f,2.381000f,3.990000f};
    const float iG0_ref[5][6] = { /* testShelfParamsInterp_iG_0 */
        {22.670282f,21.531200f,2.814473f,-5.412009f,-8.989264f,-9.089385f},
        {18.295933f,17.655270f,3.178890f,-4.810981f,-8.364464f,-8.467470f},
        {11.937032f,11.742982f,3.538333f,-3.463974f,-6.856924f,-6.961620f},
        {6.676990f,6.625110f,3.023452f,-1.880613f,-4.729220f,-4.823625f},
        {3.628860f,3.607114f,2.019588f,-0.892461f,-2.938593f,-3.010937f}
    };
    const float iGInf_ref[5][6] = { /* testShelfParamsInterp_iG_inf */
        {-4.643651f,-4.547427f,-6.154277f,-11.120993f,-8.603536f,-8.574359f},
        {-4.128221f,-4.060667f,-5.063987f,-9.950522f,-7.573856f,-7.518260f},
        {-3.094135f,-3.064137f,-3.266476f,-7.650444f,-5.524759f,-5.437006f},
        {-1.937289f,-1.926201f,-1.763717f,-4.875811f,-3.327342f,-3.231471f},
        {-1.126412f,-1.121129f,-0.951611f,-2.838410f,-1.878207f,-1.800638f}
    };
    const float iFc_ref[5][6] = { /* testShelfParamsInterp_if_c */
        {525.636204f,498.827915f,2386.832594f,4596.197114f,4931.390665f,5210.861482f},
        {410.072475f,405.299321f,1861.960103f,4441.658657f,4938.293282f,5155.085794f},
        {358.247441f,357.024760f,999.146666f,4311.428254f,4994.348051f,5164.504029f},
        {318.842699f,318.694795f,468.086862f,4161.363071f,5092.451748f,5267.402018f},
        {297.454930f,297.481562f,334.402844f,4018.349277f,5175.836901f,5380.750972f}
    };
    /* High shelf parameters to be tested */
    float iG0, iGInf, iFc;

    for(int ri = 0; ri < nRho; ri++){
        float r = rho[ri];
        for(int ti = 0; ti < nTheta; ti++){
            interpDVFShelfParams(theta[ti], r, &iG0, &iGInf, &iFc);
            TEST_ASSERT_FLOAT_WITHIN(acceptedTolerance, iG0, iG0_ref[ri][ti]);
            TEST_ASSERT_FLOAT_WITHIN(acceptedTolerance, iGInf, iGInf_ref[ri][ti]);
            TEST_ASSERT_FLOAT_WITHIN(acceptedTolerance_frq, iFc, iFc_ref[ri][ti]);
        }
    }
}

void test__dvf_dvfShelfCoeffs(void){
    /* setup */
    const float acceptedTolerance = 0.00001f;
    const int   nTheta = 6;
    const float theta[6] = {0.000000f,2.300000f,47.614000f,98.600000f,166.200000f,180.000000f};
    const int   nRho = 5;
    const float rho[5] = {1.150000f,1.250000f,1.570000f,2.381000f,3.990000f};
    const float fs = 44100;
    const float b0_ref[5][6] = { /* testIIRCoefs_b0 */
        {8.084117f,7.162779f,0.733534f,0.181211f,0.157716f,0.157753f},
        {5.162983f,4.832104f,0.847478f,0.218379f,0.188114f,0.188100f},
        {2.787888f,2.735502f,1.052975f,0.322169f,0.274301f,0.274359f},
        {1.733188f,1.725027f,1.162720f,0.508950f,0.432220f,0.432405f},
        {1.337147f,1.334601f,1.133469f,0.693396f,0.606185f,0.606313f}
    };
    const float b1_ref[5][6] = { /* testIIRCoefs_b1 */
        {-7.521737f,-6.689086f,-0.525709f,-0.092147f,-0.075692f,-0.072026f},
        {-4.880667f,-4.570874f,-0.654751f,-0.113974f,-0.090171f,-0.086752f},
        {-2.654257f,-2.604818f,-0.917824f,-0.171818f,-0.130191f,-0.126320f},
		{-1.659057f,-1.651278f,-1.090421f,-0.278191f,-0.201595f,-0.195404f},
        {-1.283715f,-1.281267f,-1.082675f,-0.387883f,-0.278529f,-0.268341f}
    };
    const float a1_ref[5][6] = { /* testIIRCoefs_a1 */
        {-0.958646f,-0.960287f,-0.849695f,-0.833925f,-0.769110f,-0.755889f},
        {-0.965649f,-0.965782f,-0.866341f,-0.818335f,-0.743436f,-0.731349f},
        {-0.966189f,-0.966188f,-0.910070f,-0.775971f,-0.682649f,-0.670043f},
        {-0.965632f,-0.965605f,-0.948954f,-0.713458f,-0.602472f,-0.587018f},
        {-0.964816f,-0.964791f,-0.959744f,-0.661426f,-0.540433f,-0.522001f}
    };

    float b0, b1, a1;       /* IIR coeffs to be tested */
    float iG0, iGInf, iFc;

    for(int ri = 0; ri < nRho; ri++){
        float r = rho[ri];
        for(int ti = 0; ti < nTheta; ti++){
            interpDVFShelfParams(theta[ti], r, &iG0, &iGInf, &iFc);
            dvfShelfCoeffs(iG0, iGInf, iFc, fs, &b0, &b1, &a1);
            TEST_ASSERT_FLOAT_WITHIN(acceptedTolerance, b0, b0_ref[ri][ti]);
            TEST_ASSERT_FLOAT_WITHIN(acceptedTolerance, b1, b1_ref[ri][ti]);
            TEST_ASSERT_FLOAT_WITHIN(acceptedTolerance, a1, a1_ref[ri][ti]);
        }
    }
}
