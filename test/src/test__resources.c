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
 * @file test__resources.c
 * @brief Unit tests for the SAF third-party resources
 * @author Leo McCormack
 * @date 27.04.2020
 * @license ISC
 */

#include "saf_test.h"

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
    //inspec = (float_complex***)malloc3d(nBands, nCHin+12, nHops+10, sizeof(float_complex));
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
        //afSTFT_forward_knownSize(hSTFT, inframe, framesize, nCHin+12, nHops+10, inspec);

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

void test__malloc4d(void){
    int i,j,k,l;
    int REF[3][4][2][5];
    int CPY[3][4][2][5];
    int**** test_malloc_4d;
    test_malloc_4d = (int****)malloc4d(3, 4, 2, 5, sizeof(int));

    /* Fill the reference static 4D array, and the dynamically allocated 4D array with the same values */
    for(i=0; i<3; i++){
        for(j=0; j<4; j++){
            for(k=0; k<2; k++){
                for(l=0; l<5; l++){
                    test_malloc_4d[i][j][k][l] = i*4*2*5 + j*2*5 + k*5 + l;
                    REF[i][j][k][l] = i*4*2*5 + j*2*5 + k*5 + l;
                }
            }
        }
    }

    /* Copy the dynamically allocated array to a static copy (to check that the data has actually been contiguously allocated) */
    memcpy(CPY, FLATTEN4D(test_malloc_4d), 3*4*2*5*sizeof(int));
    for(i=0; i<3; i++)
        for(j=0; j<4; j++)
            for(k=0; k<2; k++)
                for(l=0; l<5; l++) /* The copy should be identical to the reference */
                    TEST_ASSERT_TRUE(CPY[i][j][k][l] == REF[i][j][k][l]);

    /* Clean-up */
    free(test_malloc_4d);
}

void test__malloc5d(void){
    int i,j,k,l,p;
    int REF[2][4][3][5][2];
    int CPY[2][4][3][5][2];
    int***** test_malloc_5d;
    test_malloc_5d = (int*****)malloc5d(2, 4, 3, 5, 2, sizeof(int));

    /* Fill the reference static 5D array, and the dynamically allocated 5D array with the same values */
    for(i=0; i<2; i++){
        for(j=0; j<4; j++){
            for(k=0; k<3; k++){
                for(l=0; l<5; l++){
                    for(p=0; p<2; p++){
                        test_malloc_5d[i][j][k][l][p] = i*4*3*5*2 + j*3*5*2 + k*5*2 + l*2 + p;
                        REF[i][j][k][l][p] = i*4*3*5*2 + j*3*5*2 + k*5*2 + l*2 + p;
                    }
                }
            }
        }
    }

    /* Copy the dynamically allocated array to a static copy (to check that the data has actually been contiguously allocated) */
    memcpy(CPY, FLATTEN5D(test_malloc_5d), 2*4*3*5*2*sizeof(int));
    for(i=0; i<2; i++)
        for(j=0; j<4; j++)
            for(k=0; k<3; k++)
                for(l=0; l<5; l++)
                    for(p=0; p<2; p++)/* The copy should be identical to the reference */
                        TEST_ASSERT_TRUE(CPY[i][j][k][l][p] == REF[i][j][k][l][p]);

    /* Clean-up */
    free(test_malloc_5d);
}

void test__malloc6d(void){
    int i,j,k,l,p,o;
    int REF[2][3][2][4][2][3];
    int CPY[2][3][2][4][2][3];
    int****** test_malloc_6d;
    test_malloc_6d = (int******)malloc6d(2, 3, 2, 4, 2, 3, sizeof(int));

    /* Fill the reference static 5D array, and the dynamically allocated 5D array with the same values */
    for(i=0; i<2; i++){
        for(j=0; j<3; j++){
            for(k=0; k<2; k++){
                for(l=0; l<4; l++){
                    for(p=0; p<2; p++){
                        for(o=0; o<3; o++){
                            test_malloc_6d[i][j][k][l][p][o] = i*3*2*4*2*3 + j*2*4*2*3 + k*4*2*3 + l*2*3 + p*3 + o;
                            REF[i][j][k][l][p][o] = i*3*2*4*2*3 + j*2*4*2*3 + k*4*2*3 + l*2*3 + p*3 + o;
                        }
                    }
                }
            }
        }
    }

    /* Copy the dynamically allocated array to a static copy (to check that the data has actually been contiguously allocated) */
    memcpy(CPY, FLATTEN6D(test_malloc_6d), 2*3*2*4*2*3*sizeof(int));
    for(i=0; i<2; i++)
        for(j=0; j<3; j++)
            for(k=0; k<2; k++)
                for(l=0; l<4; l++)
                    for(p=0; p<2; p++)
                        for(o=0; o<3; o++)/* The copy should be identical to the reference */
                            TEST_ASSERT_TRUE(CPY[i][j][k][l][p][o] == REF[i][j][k][l][p][o]);

    /* Clean-up */
    free(test_malloc_6d);
}
