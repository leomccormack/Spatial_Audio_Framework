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
 * @file test__examples.c
 * @brief Unit tests for the SAF examples
 * @author Leo McCormack
 * @date 27.04.2020
 * @license ISC
 */

#include "saf_test.h"

#ifdef SAF_ENABLE_EXAMPLES_TESTS

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
    ambi_bin_setEnableRotation(hAmbi, 1);
    ambi_bin_setYaw(hAmbi, 180.0f); /* turn the listener around */

    /* Define input mono signal */
    nSH = ORDER2NSH(order);
    inSig = malloc1d(signalLength*sizeof(float));
    shSig = (float**)malloc2d(nSH,signalLength,sizeof(float));
    rand_m1_1(inSig, signalLength); /* Mono white-noise signal */

    /* Encode to get input spherical harmonic (Ambisonic) signal */
    direction_deg[0] = -90.0f; /* encode hard-right */
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

        ambi_bin_process(hAmbi, (const float* const*)shSig_frame, binSig_frame, nSH, NUM_EARS, framesize);
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

        ambi_dec_process(hAmbi, (const float* const*)shSig_frame, lsSig_frame, nSH, 22, framesize);
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

        ambi_enc_process(hAmbi, (const float* const*)inSig_frame, shSig_frame, 2, nSH, framesize);
    }

    /* ambi_enc should be equivalent to the reference */
    for(i=0; i<nSH; i++)
        for(j=0; j<signalLength-delay; j++)
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

        array2sh_process(hA2sh, (const float* const*)micSig_frame, shSig_frame, 32, nSH, framesize);
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
    rotator_setYaw(hRot, ypr[0]*180.0f/SAF_PI); /* rad->degrees */
    rotator_setPitch(hRot, ypr[1]*180.0f/SAF_PI);
    rotator_setRoll(hRot, ypr[2]*180.0f/SAF_PI);

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

        rotator_process(hRot, (const float* const*)shSig_frame, shSig_rot_frame, nSH, nSH, framesize);
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

void test__saf_example_spreader(void){
    int i, ch, framesize, nOutputs;
    void* hSpr;
    float** inSigs, **outSigs;
    float** inSig_frame, **outSig_frame;

    /* Config */
    const int fs = 48000;
    const int nInputs = 1;
    const int signalLength = fs*2;

    /* Create and initialise an instance of spreader */
    spreader_create(&hSpr);

    /* Configure and initialise the spreader codec */
    spreader_setUseDefaultHRIRsflag(hSpr, 1);
    nOutputs = NUM_EARS; /* the default is binaural operation */
    spreader_setNumSources(hSpr, nInputs);
    spreader_init(hSpr, fs); /* Should be called before calling "process"
                               * Cannot be called while "process" is on-going */
    spreader_initCodec(hSpr); /* Can be called whenever (thread-safe) */

    /* Define input mono signal */
    inSigs = (float**)malloc2d(nInputs,signalLength,sizeof(float));
    outSigs = (float**)malloc2d(nOutputs,signalLength,sizeof(float));
    rand_m1_1(FLATTEN2D(inSigs), nInputs*signalLength); /* white-noise signals */

    /* Apply spreader */
    framesize = spreader_getFrameSize();
    inSig_frame = (float**)malloc1d(nInputs*sizeof(float*));
    outSig_frame = (float**)malloc1d(nOutputs*sizeof(float*));
    for(i=0; i<(int)((float)signalLength/(float)framesize); i++){
        for(ch=0; ch<nInputs; ch++)
            inSig_frame[ch] = &inSigs[ch][i*framesize];
        for(ch=0; ch<nOutputs; ch++)
            outSig_frame[ch] = &outSigs[ch][i*framesize];

        spreader_process(hSpr, (const float* const*)inSig_frame, outSig_frame, nInputs, nOutputs, framesize);
    }

    /* Clean-up */
    spreader_destroy(&hSpr);
    free(inSigs);
    free(outSigs);
    free(inSig_frame);
    free(outSig_frame);
}

#endif /* SAF_ENABLE_EXAMPLES_TESTS */
