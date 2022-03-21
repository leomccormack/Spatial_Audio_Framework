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
 * @file test__sh_module.c
 * @brief Unit tests for the SAF sh module
 * @author Leo McCormack
 * @date 27.04.2020
 * @license ISC
 */

#include "saf_test.h"      /* the SAF unit test declarations */

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
            t_dirs_rad[j][0] = t_dirs_deg[j*2] * SAF_PI/180.0f;
            t_dirs_rad[j][1] = SAF_PI/2.0f - t_dirs_deg[j*2+1] * SAF_PI/180.0f; /* elevation->inclination */
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
        dir[0] *= SAF_PI;
        dir[1] *= SAF_PI/2.0f;
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
            t_dirs_rad[j][0] = t_dirs_deg[j*2] * SAF_PI/180.0f;
            t_dirs_rad[j][1] = SAF_PI/2.0f - t_dirs_deg[j*2+1] * SAF_PI/180.0f; /* elevation->inclination */
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
    int i,j,nSH,order;
    float Rzyx[3][3];
    float** Mrot;

    /* Config */
    const float acceptedTolerance = 0.00001f;

    /* Rotation matrix for 0,0,0 should be identity */
    order = 22;
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
            dir[0] *= SAF_PI;
            dir[1] *= SAF_PI/2.0f;

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
            dir[0] *= SAF_PI;
            dir[1] *= SAF_PI/2.0f;

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
            t_dirs_rad[j][0] = t_dirs_deg[j * 2] * SAF_PI / 180.0f;
            t_dirs_rad[j][1] = SAF_PI / 2.0f - t_dirs_deg[j * 2 + 1] * SAF_PI /
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

void test__calculateGridWeights(void){
    int i, nDirs, order;
    float* t_dirs_deg, *t_dirs_rad, *w;

    /* Config */
    const float acceptedTolerance = 0.00001f;
    int testOrder = 3;
    
    /* Pull an appropriate t-design */
    t_dirs_deg = (float*)__HANDLES_Tdesign_dirs_deg[2 * testOrder -1];
    nDirs = __Tdesign_nPoints_per_degree[2 * testOrder - 1];
    t_dirs_rad = (float *) malloc1d(nDirs*2 *sizeof(float));
    for(i=0; i<2*nDirs; i++)
        t_dirs_rad[i] = DEG2RAD(t_dirs_deg[i]);
    for(i=0; i<nDirs; i++)
        t_dirs_rad[i*2+1] = SAF_PI/2.0f - t_dirs_rad[i*2+1];
    w = malloc1d(nDirs * sizeof(float));
    order =  calculateGridWeights(t_dirs_rad,nDirs,-1,w);

    TEST_ASSERT_EQUAL(testOrder, order);
    for (i=0; i<nDirs; i++)
        TEST_ASSERT_FLOAT_WITHIN(acceptedTolerance,FOURPI/nDirs, w[i]);

    /* clean-up */
    free(t_dirs_rad);
    free(w);
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
    utility_sseig(NULL, FLATTEN2D(Cx), nSH, 1, FLATTEN2D(V), NULL, NULL);
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

void test__sphPWD(void){
    int nGrid, nSH, nSrcs, srcInd_1, srcInd_2;
    float test_dirs_deg[2][2];
    float* grid_dirs_deg;
    float** Y_src, **src_sigs, **src_sigs_sh, **Cx;
    float_complex** Cx_cmplx;
    void* hPWD;

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
    Cx_cmplx = (float_complex**)calloc2d(nSH, nSH, sizeof(float_complex));
    cblas_scopy(nSH*nSH, FLATTEN2D(Cx), 1, (float*)FLATTEN2D(Cx_cmplx), 2);

    /* compute sphMUSIC, returning "peak-find" indices */
    int inds[2];
    sphPWD_create(&hPWD, order, grid_dirs_deg, nGrid);
    sphPWD_compute(hPWD, FLATTEN2D(Cx_cmplx), nSrcs, NULL, (int*)inds);

    /* Assert that the true source indices were found (note that the order can flip) */
    TEST_ASSERT_TRUE(inds[0] == srcInd_1 || inds[0] == srcInd_2);
    TEST_ASSERT_TRUE(inds[1] == srcInd_1 || inds[1] == srcInd_2);

    /* clean-up */
    sphPWD_destroy(&hPWD);
    free(Y_src);
    free(src_sigs);
    free(src_sigs_sh);
    free(Cx);
    free(Cx_cmplx);
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
    utility_cseig(NULL, FLATTEN2D(C_Cx), nSH, 1, FLATTEN2D(U), NULL, NULL);
    Us = (float_complex**)malloc2d(nSH, nSrcs, sizeof(float_complex)); /* signal subspace */
    for(i=0; i<nSH; i++)
        for(j=0; j<nSrcs; j++)
            Us[i][j] = U[i][j];

    /* use sphESPRIT to estimate source directions... */
    sphESPRIT_create(&hESPRIT, order);
    sphESPRIT_estimateDirs(hESPRIT, FLATTEN2D(Us), nSrcs, (float*)estdirs_deg);
    for(i=0; i<nSrcs; i++){
        estdirs_deg[i][0]*=180.0f/SAF_PI; /* rad->deg */
        estdirs_deg[i][1]*=180.0f/SAF_PI;
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
