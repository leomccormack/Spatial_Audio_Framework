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
 * @file test__reverb_module.c
 * @brief Unit tests for the SAF reverb module
 * @author Leo McCormack
 * @date 27.04.2020
 * @license ISC
 */

#include "saf_test.h"

void test__ims_shoebox_RIR(void){
    void* hIms;
    float maxTime_s;
    float mov_src_pos[3], mov_rec_pos[3];
    int sourceID_1, sourceID_2, sourceID_3, sourceID_4, sourceID_5, receiverID;
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
    const float roomdims[3] = {10.0f, 7.0f, 3.0f};

    /* Set-up the shoebox room simulator, with two sources and one spherical harmonic receiver */
    ims_shoebox_create(&hIms, (float*)roomdims, (float*)abs_wall, 125.0f, nBands, 343.0f, 48e3f);
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
        ims_shoebox_computeEchograms(hIms, -1, maxTime_s);
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
        ims_shoebox_computeEchograms(hIms, -1, maxTime_s);
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
    int sourceIDs[4], receiverIDs[1];
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
    const float roomdims[3] = {10.0f, 7.0f, 3.0f};

    /* Allocate memory for 4 sources and 1 spherical harmonic receiver */
    rec_sh_outsigs = (float***)malloc3d(1, ORDER2NSH(sh_order), signalLength, sizeof(float));
    src_sigs = (float**)malloc2d(4, signalLength, sizeof(float));
    rand_m1_1(FLATTEN2D(src_sigs), 4*signalLength);

    /* Set-up the shoebox room simulator for these four sources and SH receiver */
    ims_shoebox_create(&hIms, (float*)roomdims, (float*)abs_wall, 250.0f, nBands, 343.0f, 48e3f);
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
        ims_shoebox_computeEchograms(hIms, -1, maxTime_s);
        ims_shoebox_applyEchogramTD(hIms, receiverIDs[0], signalLength, 0);
    }

    /* clean-up */
    free(src_sigs);
    free(rec_sh_outsigs);
    ims_shoebox_destroy(&hIms);
}
