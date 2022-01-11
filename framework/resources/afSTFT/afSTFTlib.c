/*
 Copyright (c) 2015 Juha Vilkamo
 
 Permission is hereby granted, free of charge, to any person obtaining a copy
 of this software and associated documentation files (the "Software"), to deal
 in the Software without restriction, including without limitation the rights
 to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:
 
 The above copyright notice and this permission notice shall be included in
 all copies or substantial portions of the Software.
 
 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 THE SOFTWARE.
 */

/**
 * @file afSTFTlib.c
 * @brief A modified version of afSTFTlib
 *
 * The original afSTFT code (by Juha Vilkamo) can be found here:
 * https://github.com/jvilkamo/afSTFT
 * This version is slightly modified to be more in-line with how the rest of SAF
 * is structured.
 * The files afSTFTlib.h/.c act as the interface to afSTFT, which is then
 * implemented in afSTFT_internal.h/.c.
 *
 * This version also adds functionality to change the number of channels on the
 * fly, flush the run-time buffers with zeros, return the current frequency
 * vector and the current processing delay.
 * It also incorporates SAF utilities (for the vectorisation and FFT).
 *
 * The afSTFT design is also described in more detail in [1]
 *
 * @see [1] Vilkamo, J., & Ba"ckstro"m, T. (2018). Time--Frequency Processing:
 *          Methods and Tools. In Parametric Time--Frequency Domain Spatial
 *          Audio. John Wiley & Sons.
 *
 * @author Juha Vilkamo
 * @date 08.04.2015
 * @license MIT
 */

#include "afSTFTlib.h"
#include "afSTFT_internal.h"

/** afSTFT centre frequencies for 128 hop size, hybrid-mode enabled, 48kHz */
static const double __afCenterFreq48e3[133] =
{    0.000000000, 140.644316361, 234.355478108, 328.144332285, 421.855497937, 515.644326841, 609.355515147, 703.144330614, 796.855543885, 937.500032020, 1125.000017338, 1312.500035449, 1500.000075751, 1687.500031782, 1875.000024239, 2062.499975101, 2250.000053703, 2437.500044271, 2625.000002315, 2812.500019782, 3000.000041692, 3187.499983930, 3374.999995137, 3562.499994173, 3750.000018557, 3937.500021643, 4125.000009859, 4312.500011528, 4500.000010423, 4687.500014446, 4875.000013588, 5062.500013570, 5250.000007575, 5437.500010288, 5625.000004178, 5812.500003421, 6000.000005158, 6187.500003404, 6375.000003488, 6562.500007191, 6750.000005972, 6937.500008499, 7125.000006936, 7312.500008549, 7500.000005032, 7687.500004875, 7875.000004878, 8062.500007586, 8250.000006218, 8437.499999805, 8625.000000113, 8812.499997984, 9000.000008860, 9187.500004401, 9375.000001529, 9562.500006565, 9750.000006335, 9937.499999557, 10125.000002928, 10312.500002384, 10500.000004406, 10687.500002820, 10875.000001403, 11062.500002219, 11250.000001097, 11437.500001292, 11625.000000815, 11812.500000140, 12000.000000000, 12187.499999584, 12374.999999473, 12562.499999294, 12749.999998799, 12937.499997514, 13124.999998543, 13312.499997602, 13499.999995904, 13687.499996961, 13874.999996550, 14062.500000495, 14249.999993960, 14437.499993440, 14624.999997861, 14812.499995461, 14999.999991137, 15187.500001756, 15374.999999428, 15562.500000999, 15749.999993809, 15937.499992382, 16124.999995683, 16312.499995240, 16499.999994365, 16687.499991354, 16874.999992234, 17062.499991361, 17249.999994298, 17437.499992410, 17624.999995960, 17812.499995945, 17999.999994836, 18187.499996913, 18374.999996125, 18562.499990092, 18749.999991865, 18937.499986965, 19124.999985762, 19312.499985261, 19499.999989766, 19687.499988292, 19874.999989851, 20062.499978542, 20249.999981602, 20437.500005879, 20625.000004853, 20812.500015815, 20999.999958305, 21187.499980259, 21374.999997733, 21562.499955794, 21749.999946298, 21937.500025004, 22124.999975461, 22312.499968567, 22499.999924162, 22687.499964503, 22874.999982475, 23062.499968048, 23249.999976609, 23437.499982579, 23624.999922020, 23812.499893152, 24000.000000000    };

/** afSTFT centre frequencies for 128 hop size, hybrid-mode enabled, 44.1kHz */
static const double __afCenterFreq44100[133] =
{    0.000000000, 129.216965656, 215.314095512, 301.482605287, 387.579738729, 473.748225285, 559.845379541, 646.013853751, 732.111030944, 861.328154418, 1033.593765929, 1205.859407569, 1378.125069596, 1550.390654200, 1722.656272269, 1894.921852124, 2067.187549340, 2239.453165674, 2411.718752127, 2583.984393174, 2756.250038304, 2928.515610236, 3100.781245532, 3273.046869646, 3445.312517049, 3617.578144885, 3789.843759058, 3962.109385592, 4134.375009576, 4306.640638272, 4478.906262484, 4651.171887467, 4823.437506959, 4995.703134452, 5167.968753839, 5340.234378143, 5512.500004739, 5684.765628127, 5857.031253205, 6029.296881607, 6201.562505487, 6373.828132809, 6546.093756373, 6718.359382855, 6890.625004623, 7062.890629479, 7235.156254481, 7407.421881970, 7579.687505713, 7751.953124821, 7924.218750103, 8096.484373148, 8268.750008140, 8441.015629043, 8613.281251405, 8785.546881031, 8957.812505821, 9130.078124593, 9302.343752690, 9474.609377190, 9646.875004048, 9819.140627591, 9991.406251289, 10163.671877038, 10335.937501008, 10508.203126187, 10680.468750748, 10852.734375129, 11025.000000000, 11197.265624618, 11369.531249516, 11541.796874351, 11714.062498897, 11886.328122716, 12058.593748662, 12230.859372797, 12403.124996237, 12575.390622207, 12747.656246830, 12919.921875455, 13092.187494451, 13264.453118973, 13436.718748035, 13608.984370830, 13781.249991857, 13953.515626614, 14125.781249475, 14298.046875918, 14470.312494312, 14642.578118001, 14814.843746034, 14987.109370627, 15159.374994823, 15331.640617057, 15503.906242865, 15676.171867063, 15848.437494762, 16020.703118027, 16192.968746288, 16365.234371274, 16537.499995256, 16709.765622164, 16882.031246440, 17054.296865897, 17226.562492526, 17398.828113024, 17571.093736919, 17743.359361459, 17915.624990597, 18087.890614243, 18260.156240676, 18432.421855285, 18604.687483097, 18776.953130401, 18949.218754459, 19121.484389530, 19293.749961692, 19466.015606863, 19638.281247918, 19810.546834386, 19982.812450661, 20155.078147972, 20327.343727455, 20499.609346121, 20671.874930324, 20844.140592387, 21016.406233899, 21188.671845644, 21360.937478510, 21533.203108994, 21705.468678356, 21877.734276834, 22050.000000000    };

/**
 * Matrix for converting the centre frequencies of the first 5 stft bins into
 * the centre frequencies for the 9 hybrid bands.
 */
static const float __stft2hybCentreFreq[9][5] =
{ { 1.0f, 0.0f, 0.0f, 0.0f, 0.0f },
  { 0.0f, 0.7501f, 0.0f, 0.0f, 0.0f },
  { 0.0f, 1.2499f, 0.0f, 0.0f, 0.0f },
  { 0.0f, 0.0f, 0.8751f, 0.0f, 0.0f },
  { 0.0f, 0.0f, 1.1249f, 0.0f, 0.0f },
  { 0.0f, 0.0f, 0.0f, 0.9167f, 0.0f },
  { 0.0f, 0.0f, 0.0f, 1.0833f, 0.0f },
  { 0.0f, 0.0f, 0.0f, 0.0f, 0.9375f },
  { 0.0f, 0.0f, 0.0f, 0.0f, 1.0625f } };

/** Passes input time-domain data through the afSTFT filterbank.
 *  Copyright (c) 2015 Juha Vilkamo, MIT license */
static void afAnalyse
(
    float* inTD/* nSamplesTD x nCH */,
    int nSamplesTD,
    int nCH,
    int hopSize,
    int LDmode,
    int hybridmode, 
    float_complex* outTF /* out_nBands x nTimeslots x nCH */
)
{
    int t, ch, sample, band;
    void* hSTFT;
    float_complex*** FrameTF;
    float** tempFrameTD;
    int nTimeSlots, nBands;

    nBands = hopSize + (hybridmode ? 5 : 1);
    nTimeSlots = nSamplesTD/hopSize;

    /* allocate memory */
    afSTFT_create(&(hSTFT), nCH, 1, hopSize, LDmode, hybridmode, AFSTFT_TIME_CH_BANDS);
    FrameTF = (float_complex***)malloc3d(nTimeSlots, nCH, nBands, sizeof(float_complex));
    tempFrameTD = (float**)calloc2d(nCH, nTimeSlots*hopSize, sizeof(float));

    /* perform TF transform */
    for(ch=0; ch<nCH; ch++)
        for(sample=0; sample<nSamplesTD; sample++)
            tempFrameTD[ch][sample] = inTD[sample* nCH + ch];
    afSTFT_forward(hSTFT, tempFrameTD, nTimeSlots*hopSize, FrameTF);

    /* save result to output */
    for (band = 0; band < nBands; band++)
        for (t = 0; t < nTimeSlots; t++)
            for (ch = 0; ch < nCH; ch++)
                outTF[band * nTimeSlots * nCH + t * nCH + ch] = FrameTF[t][ch][band];

    /* clean-up */
    afSTFT_destroy(&hSTFT);
    free(FrameTF);
    free(tempFrameTD);
}

/** Data structure for the afSTFT filterbank */
typedef struct _afSTFT_data {
    int hopsize;                      /**< Hop size in samples */
    int hybridmode;                   /**< 1: hybrid filtering enabled; 0: disabled */
    int nCHin;                        /**< Number of input channels */
    int nCHout;                       /**< Number of output channels*/
    int nBands;                       /**< Number of frequency bands */
    AFSTFT_FDDATA_FORMAT format;      /**< see #AFSTFT_FDDATA_FORMAT */
    void* hInt;                       /**< Internal handle for afSTFT */
    complexVector* STFTInputFrameTF;  /**< Internal input complex buffer */
    complexVector* STFTOutputFrameTF; /**< Internal output complex buffer */
    int afSTFTdelay;                  /**< Processing delay in samples */
    float** tempHopFrameTD;           /**< temporary multi-channel time-domain buffer of size "HOP_SIZE". */

}afSTFT_data;


/* ========================================================================== */
/*                                Main functions                              */
/* ========================================================================== */

void afSTFT_create
(
    void ** const phSTFT,
    int nCHin,
    int nCHout,
    int hopsize,
    int lowDelayMode,
    int hybridmode,
    AFSTFT_FDDATA_FORMAT format
)
{
    *phSTFT = malloc1d(sizeof(afSTFT_data));
    afSTFT_data *h = (afSTFT_data*)(*phSTFT);
    int ch;

    if(hybridmode)
        assert(hopsize==64 || hopsize==128 || hopsize==256);
    assert(1024 % hopsize == 0 );

    h->nCHin = nCHin;
    h->nCHout = nCHout;
    h->hopsize = hopsize;
    h->hybridmode = hybridmode;
    h->nBands = hybridmode ? hopsize+5 : hopsize+1; /* hybrid mode incurs an additional 4 bands */
    if(lowDelayMode)
        h->afSTFTdelay = hybridmode ? 7*hopsize  : 4*hopsize; /* hybrid mode incurs 3 additional hops of latency */
    else
        h->afSTFTdelay = hybridmode ? 12*hopsize : 9*hopsize; /* hybrid mode incurs 3 additional hops of latency */
    h->format = format;

    /* init afSTFT core */
    afSTFTlib_init(&(h->hInt), hopsize, nCHin, nCHout, lowDelayMode, hybridmode);

    /* temp buffers */
    if(nCHout>0){
        h->STFTOutputFrameTF = malloc1d(nCHout * sizeof(complexVector));
        for(ch=0; ch < nCHout; ch++) {
            h->STFTOutputFrameTF[ch].re = (float*)calloc1d(h->nBands, sizeof(float));
            h->STFTOutputFrameTF[ch].im = (float*)calloc1d(h->nBands, sizeof(float));
        }
    }
    else
        h->STFTOutputFrameTF = NULL;
    if(nCHout > 0 || nCHin > 0)
        h->tempHopFrameTD = (float**)malloc2d( SAF_MAX(nCHin, nCHout), hopsize, sizeof(float));
    if(nCHin>0){
        h->STFTInputFrameTF = malloc1d(nCHin * sizeof(complexVector));
        for(ch=0; ch < nCHin; ch++) {
            h->STFTInputFrameTF[ch].re = (float*)calloc1d(h->nBands, sizeof(float));
            h->STFTInputFrameTF[ch].im = (float*)calloc1d(h->nBands, sizeof(float));
        }
    }
    else
        h->STFTInputFrameTF = NULL;
}

void afSTFT_destroy
(
    void ** const phSTFT
)
{
    afSTFT_data *h = (afSTFT_data*)(*phSTFT);
    int ch;

    if(h!=NULL){
        /* For run-time */
        afSTFTlib_free(h->hInt);
        if(h->STFTInputFrameTF!=NULL){
            for (ch = 0; ch< h->nCHin; ch++) {
                free(h->STFTInputFrameTF[ch].re);
                free(h->STFTInputFrameTF[ch].im);
            }
        }
        for (ch = 0; ch< h->nCHout; ch++) {
            free(h->STFTOutputFrameTF[ch].re);
            free(h->STFTOutputFrameTF[ch].im);
        }
        free(h->STFTInputFrameTF);
        free(h->STFTOutputFrameTF);
        free(h->tempHopFrameTD);

        free(h);
        h=NULL;
        *phSTFT = NULL;
    }
}

void afSTFT_forward
(
    void * const hSTFT,
    float** dataTD,
    int framesize,
    float_complex*** dataFD
)
{
    afSTFT_data *h = (afSTFT_data*)(hSTFT);
    int ch, t, nHops, band;

    assert(framesize % h->hopsize == 0); /* framesize must be multiple of hopsize */
    nHops = framesize/h->hopsize;

    /* Loop over hops */
    for(t=0; t < nHops; t++) {
        /* forward transform */
        for(ch = 0; ch < h->nCHin; ch++)
            utility_svvcopy(&(dataTD[ch][t*(h->hopsize)]), (h->hopsize), h->tempHopFrameTD[ch]);
        afSTFTlib_forward(h->hInt, h->tempHopFrameTD, h->STFTInputFrameTF);

        /* store */
        switch(h->format){
            case AFSTFT_BANDS_CH_TIME:
                for(band=0; band<h->nBands; band++)
                    for(ch=0; ch < h->nCHin; ch++)
                        dataFD[band][ch][t] = cmplxf(h->STFTInputFrameTF[ch].re[band], h->STFTInputFrameTF[ch].im[band]);
                break;
            case AFSTFT_TIME_CH_BANDS:
                for(ch=0; ch < h->nCHin; ch++){
                    cblas_scopy(h->nBands, h->STFTInputFrameTF[ch].re, 1, (float*)dataFD[t][ch], 2);
                    cblas_scopy(h->nBands, h->STFTInputFrameTF[ch].im, 1, &((float*)dataFD[t][ch])[1], 2);
                }
                break;
        }
    }
}

void afSTFT_forward_knownDimensions
(
    void * const hSTFT,
    float** dataTD,
    int framesize,
    int dataFD_nCH,
    int dataFD_nHops,
    float_complex*** dataFD
)
{
    afSTFT_data *h = (afSTFT_data*)(hSTFT);
    int ch, t, nHops;
    float_complex* pDataFD;

    assert(framesize % h->hopsize == 0); /* framesize must be multiple of hopsize */
    nHops = framesize/h->hopsize;
    pDataFD = &dataFD[0][0][0];

    /* Loop over hops */
    for(t=0; t < nHops; t++) {
        /* forward transform */
        for(ch = 0; ch < h->nCHin; ch++)
            utility_svvcopy(&(dataTD[ch][t*(h->hopsize)]), (h->hopsize), h->tempHopFrameTD[ch]);
        afSTFTlib_forward(h->hInt, h->tempHopFrameTD, h->STFTInputFrameTF);

        /* store */
        switch(h->format){
            case AFSTFT_BANDS_CH_TIME:
                for(ch=0; ch < h->nCHin; ch++){
                    cblas_scopy(h->nBands, h->STFTInputFrameTF[ch].re, 1, (float*)&pDataFD[0*dataFD_nCH*nHops + ch*dataFD_nHops + t], dataFD_nCH*dataFD_nHops*2);
                    cblas_scopy(h->nBands, h->STFTInputFrameTF[ch].im, 1, &((float*)&pDataFD[0*dataFD_nCH*nHops + ch*dataFD_nHops + t])[1], dataFD_nCH*dataFD_nHops*2);
                }
                break;
            case AFSTFT_TIME_CH_BANDS:
                for(ch=0; ch < h->nCHin; ch++){
                    cblas_scopy(h->nBands, h->STFTInputFrameTF[ch].re, 1, (float*)dataFD[t][ch], 2);
                    cblas_scopy(h->nBands, h->STFTInputFrameTF[ch].im, 1, &((float*)dataFD[t][ch])[1], 2);
                }
                break;
        }
    }
}

void afSTFT_forward_flat
(
    void * const hSTFT,
    float* dataTD,
    int framesize,
    float_complex* dataFD
)
{
    afSTFT_data *h = (afSTFT_data*)(hSTFT);
    int ch, t, nHops, band;

    assert(framesize % h->hopsize == 0); /* framesize must be multiple of hopsize */
    nHops = framesize/h->hopsize;

    /* Loop over hops */
    for(t=0; t < nHops; t++) {
        /* forward transform */
        for(ch = 0; ch < h->nCHin; ch++)
            utility_svvcopy(&(dataTD[ch * framesize + t*(h->hopsize)]), (h->hopsize), h->tempHopFrameTD[ch]);
        afSTFTlib_forward(h->hInt, h->tempHopFrameTD, h->STFTInputFrameTF);

        /* store */
        switch(h->format){
            case AFSTFT_BANDS_CH_TIME:
                for(band=0; band<h->nBands; band++)
                    for(ch=0; ch < h->nCHin; ch++)
                        dataFD[band * (h->nCHin) * nHops + ch * nHops + t] = cmplxf(h->STFTInputFrameTF[ch].re[band], h->STFTInputFrameTF[ch].im[band]);
                break;
            case AFSTFT_TIME_CH_BANDS:
                for(band=0; band<h->nBands; band++)
                    for(ch=0; ch < h->nCHin; ch++)
                        dataFD[t * (h->nCHin) * (band<h->nBands) + ch * (band<h->nBands) + band] = cmplxf(h->STFTInputFrameTF[ch].re[band], h->STFTInputFrameTF[ch].im[band]);
                break;
        }
    }
}

void afSTFT_backward
(
    void * const hSTFT,
    float_complex*** dataFD,
    int framesize,
    float** dataTD
)
{
    afSTFT_data *h = (afSTFT_data*)(hSTFT);
    int ch, t, nHops, band;

    assert(framesize % h->hopsize == 0); /* framesize must be multiple of hopsize */
    nHops = framesize/h->hopsize;

    /* Loop over hops */
    for(t = 0; t < nHops; t++) {
        /* backward transform */
        switch(h->format){
            case AFSTFT_BANDS_CH_TIME:
                for(band = 0; band < h->nBands; band++) {
                    for(ch = 0; ch < h->nCHout; ch++) {
                        h->STFTOutputFrameTF[ch].re[band] = crealf(dataFD[band][ch][t]);
                        h->STFTOutputFrameTF[ch].im[band] = cimagf(dataFD[band][ch][t]);
                    }
                }
                break;
            case AFSTFT_TIME_CH_BANDS:
                for(band = 0; band < h->nBands; band++) {
                    for(ch = 0; ch < h->nCHout; ch++) {
                        h->STFTOutputFrameTF[ch].re[band] = crealf(dataFD[t][ch][band]);
                        h->STFTOutputFrameTF[ch].im[band] = cimagf(dataFD[t][ch][band]);
                    }
                }
                break;
        }
        afSTFTlib_inverse(h->hInt, h->STFTOutputFrameTF, h->tempHopFrameTD);

        /* store */
        for (ch = 0; ch <  h->nCHout; ch++)
            memcpy(&(dataTD[ch][t*(h->hopsize)]), h->tempHopFrameTD[ch], h->hopsize*sizeof(float));
    }
}

void afSTFT_backward_knownDimensions
(
    void * const hSTFT,
    float_complex*** dataFD,
    int framesize,
    int dataFD_nCH,
    int dataFD_nHops,
    float** dataTD
)
{
    afSTFT_data *h = (afSTFT_data*)(hSTFT);
    int ch, t, nHops;
    float_complex* pDataFD;

    assert(framesize % h->hopsize == 0); /* framesize must be multiple of hopsize */
    nHops = framesize/h->hopsize;
    pDataFD = &dataFD[0][0][0];

    /* Loop over hops */
    for(t = 0; t < nHops; t++) {
        /* backward transform */
        switch(h->format){
            case AFSTFT_BANDS_CH_TIME:
                for(ch=0; ch < h->nCHout; ch++){
                    cblas_scopy(h->nBands, (float*)&pDataFD[0*dataFD_nCH*nHops + ch*dataFD_nHops + t], dataFD_nCH*dataFD_nHops*2, h->STFTOutputFrameTF[ch].re, 1);
                    cblas_scopy(h->nBands, &((float*)&pDataFD[0*dataFD_nCH*nHops + ch*dataFD_nHops + t])[1], dataFD_nCH*dataFD_nHops*2, h->STFTOutputFrameTF[ch].im, 1);
                }
                break;
            case AFSTFT_TIME_CH_BANDS:
                for(ch = 0; ch < h->nCHout; ch++) {
                    cblas_scopy(h->nBands, (float*)dataFD[t][ch], 2, h->STFTOutputFrameTF[ch].re, 1);
                    cblas_scopy(h->nBands, &((float*)dataFD[t][ch])[1], 2, h->STFTOutputFrameTF[ch].im, 1);
                }
                break;
        }
        afSTFTlib_inverse(h->hInt, h->STFTOutputFrameTF, h->tempHopFrameTD);

        /* store */
        for (ch = 0; ch <  h->nCHout; ch++)
            memcpy(&(dataTD[ch][t*(h->hopsize)]), h->tempHopFrameTD[ch], h->hopsize*sizeof(float));
    }
}

void afSTFT_backward_flat
(
    void * const hSTFT,
    float_complex* dataFD,
    int framesize,
    float* dataTD
)
{
    afSTFT_data *h = (afSTFT_data*)(hSTFT);
    int ch, t, nHops, band;

    assert(framesize % h->hopsize == 0); /* framesize must be multiple of hopsize */
    nHops = framesize/h->hopsize;

    /* Loop over hops */
    for(t = 0; t < nHops; t++) {
        /* backward transform */
        switch(h->format){
            case AFSTFT_BANDS_CH_TIME:
                for(band = 0; band < h->nBands; band++) {
                    for(ch = 0; ch < h->nCHout; ch++) {
                        h->STFTOutputFrameTF[ch].re[band] = crealf(dataFD[band * (h->nCHout) * nHops + ch * nHops + t]);
                        h->STFTOutputFrameTF[ch].im[band] = cimagf(dataFD[band * (h->nCHout) * nHops + ch * nHops + t]);
                    }
                }
                break;
            case AFSTFT_TIME_CH_BANDS:
                for(band = 0; band < h->nBands; band++) {
                    for(ch = 0; ch < h->nCHout; ch++) {
                        h->STFTOutputFrameTF[ch].re[band] = crealf(dataFD[t * (h->nCHout) * (h->nBands) + ch * (h->nBands) + band]);
                        h->STFTOutputFrameTF[ch].im[band] = cimagf(dataFD[t * (h->nCHout) * (h->nBands) + ch * (h->nBands) + band]);
                    }
                }
                break;
        }
        afSTFTlib_inverse(h->hInt, h->STFTOutputFrameTF, h->tempHopFrameTD);

        /* store */
        for (ch = 0; ch <  h->nCHout; ch++)
            memcpy(&(dataTD[ch * framesize + t*(h->hopsize)]), h->tempHopFrameTD[ch], h->hopsize*sizeof(float));
    }
}

void afSTFT_channelChange
(
    void * const hSTFT,
    int new_nCHin,
    int new_nCHout
)
{
    afSTFT_data *h = (afSTFT_data*)(hSTFT);
    int i;

    afSTFTlib_channelChange(h->hInt, new_nCHin, new_nCHout);

    /* resize buffers */
    if(h->nCHin!=new_nCHin){
        for(i=new_nCHin; i<h->nCHin; i++){
            free(h->STFTInputFrameTF[i].re);
            free(h->STFTInputFrameTF[i].im);
        }
        h->STFTInputFrameTF = realloc1d(h->STFTInputFrameTF, sizeof(complexVector)*new_nCHin);
        for(i=h->nCHin; i<new_nCHin; i++){
            h->STFTInputFrameTF[i].re = (float*)calloc1d(h->nBands, sizeof(float));
            h->STFTInputFrameTF[i].im = (float*)calloc1d(h->nBands, sizeof(float));
        }
    }
    if(h->nCHout!=new_nCHout){
        for(i=new_nCHout; i<h->nCHout; i++){
            free(h->STFTOutputFrameTF[i].re);
            free(h->STFTOutputFrameTF[i].im);
        }
        h->STFTOutputFrameTF = realloc1d(h->STFTOutputFrameTF, sizeof(complexVector)*new_nCHout);
        for(i=h->nCHout; i<new_nCHout; i++){
            h->STFTOutputFrameTF[i].re = (float*)calloc1d(h->nBands, sizeof(float));
            h->STFTOutputFrameTF[i].im = (float*)calloc1d(h->nBands, sizeof(float));
        }
    }
    if( SAF_MAX(h->nCHin, h->nCHout) != SAF_MAX(new_nCHin, new_nCHout))
        h->tempHopFrameTD = (float**)realloc2d((void**)h->tempHopFrameTD, SAF_MAX(new_nCHin, new_nCHout), h->hopsize, sizeof(float));

    h->nCHin = new_nCHin;
    h->nCHout = new_nCHout;
}

void afSTFT_clearBuffers
(
    void * const hSTFT
)
{
    afSTFT_data *h = (afSTFT_data*)(hSTFT);
    afSTFTlib_clearBuffers(h->hInt);
}

int afSTFT_getNBands
(
    void * const hSTFT
)
{
    afSTFT_data *h = (afSTFT_data*)(hSTFT);
    return h->nBands;
}

int afSTFT_getProcDelay
(
    void * const hSTFT
)
{
    afSTFT_data *h = (afSTFT_data*)(hSTFT);
    return h->afSTFTdelay;
}

void afSTFT_getCentreFreqs
(
    void * const hSTFT,
    float fs,
    int nBands,
    float* freqVector
)
{
    int band;
    if(hSTFT==NULL){
        assert(nBands>=133);
        for(band=0; band < nBands; band++){
            if(fs==44100.0f)
                freqVector[band] =  (float)__afCenterFreq44100[band];
            else /* assume 48e3 */
                freqVector[band] =  (float)__afCenterFreq48e3[band];
        }
        return;
    }

    afSTFT_data *h = (afSTFT_data*)(hSTFT);
    float* centreFreq_tmp;
    int i, j;

    assert(nBands >= h->nBands); /* just to check that "centreFreq" is of sufficient length */

    if(h->hybridmode){
        /* uniform frequency vector: */
        centreFreq_tmp = malloc1d((h->hopsize+1)*sizeof(float));
        getUniformFreqVector(h->hopsize*2, fs, centreFreq_tmp);

        /* convert first 5 to hybrid centre frequencies */
        cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 9, 1, 5, 1.0f,
                   (const float*)__stft2hybCentreFreq, 5,
                   centreFreq_tmp, 1, 0.0f,
                   freqVector, 1);

        /* Remaining centre frequencies are then the uniform  centre frequencies 6:end  */
        for(i=9, j=5; i<h->nBands; i++, j++)
            freqVector[i] = centreFreq_tmp[j];

        free(centreFreq_tmp);
    }
    else
        getUniformFreqVector(h->hopsize*2, fs, freqVector);
}

void afSTFT_FIRtoFilterbankCoeffs
(
    float* hIR /*N_dirs x nCH x ir_len*/,
    int N_dirs,
    int nCH,
    int ir_len,
    int hopSize,
    int LDmode,
    int hybridmode,
    float_complex* hFB /* nBands x nCH x N_dirs */
)
{
    int i, j, t, nd, nm, nTimeSlots, ir_pad, nBands;
    int* maxIdx;
    float maxVal, idxDel, irFB_energy, irFB_gain, phase;
    float* centerImpulse, *centerImpulseFB_energy, *ir;
    float_complex cross;
    float_complex* centerImpulseFB, *irFB;

    nBands = hopSize + (hybridmode ? 5 : 1);
    ir_pad = 1024;//+512;
    nTimeSlots = (SAF_MAX(ir_len,hopSize)+ir_pad)/hopSize;
    maxIdx = calloc1d(nCH,sizeof(int));
    centerImpulse = calloc1d(SAF_MAX(ir_len,hopSize)+ir_pad, sizeof(float));

    /* pick a direction to estimate the center of FIR delays */
    for(j=0; j<nCH; j++){
        maxVal = 2.23e-13f;
        for(i=0; i<ir_len; i++){
            if(hIR[j*ir_len + i] > maxVal){
                maxVal = hIR[j*ir_len + i];
                maxIdx[j] = i;
            }
        }
    }
    idxDel = 0.0f;
    for(j=0; j<nCH; j++)
        idxDel += (float)maxIdx[j];
    idxDel /= (float)nCH;
    idxDel = (idxDel + 1.5f);

    /* ideal impulse at mean delay */
    centerImpulse[(int)idxDel] = 1.0f;

    /* analyse impulse with the filterbank */
    centerImpulseFB = malloc1d(nBands*nTimeSlots*nCH*sizeof(float_complex));
    afAnalyse(centerImpulse, SAF_MAX(ir_len,hopSize)+ir_pad, 1, hopSize, LDmode, hybridmode, centerImpulseFB);
    centerImpulseFB_energy = calloc1d(nBands, sizeof(float));
    for(i=0; i<nBands; i++)
        for(t=0; t<nTimeSlots; t++)
            centerImpulseFB_energy[i] += powf(cabsf(centerImpulseFB[i*nTimeSlots + t]), 2.0f);

    /* initialise FB coefficients */
    ir = calloc1d( (SAF_MAX(ir_len,hopSize)+ir_pad) * nCH, sizeof(float));
    irFB = malloc1d(nBands*nCH*nTimeSlots*sizeof(float_complex));
    for(nd=0; nd<N_dirs; nd++){
        for(j=0; j<ir_len; j++)
            for(i=0; i<nCH; i++)
                ir[j*nCH+i] = hIR[nd*nCH*ir_len + i*ir_len + j];
        afAnalyse(ir, SAF_MAX(ir_len,hopSize)+ir_pad, nCH, hopSize, LDmode, hybridmode, irFB);
        for(nm=0; nm<nCH; nm++){
            for(i=0; i<nBands; i++){
                irFB_energy = 0;
                for(t=0; t<nTimeSlots; t++)
                    irFB_energy += powf(cabsf(irFB[i*nTimeSlots*nCH + t*nCH + nm]), 2.0f); /* out_nBands x nTimeslots x nCH */
                irFB_gain = sqrtf(irFB_energy/SAF_MAX(centerImpulseFB_energy[i], 2.23e-8f));
                cross = cmplxf(0.0f,0.0f);
                for(t=0; t<nTimeSlots; t++)
                    cross = ccaddf(cross, ccmulf(irFB[i*nTimeSlots*nCH + t*nCH + nm], conjf(centerImpulseFB[i*nTimeSlots + t])));
                phase = atan2f(cimagf(cross), crealf(cross));
                hFB[i*nCH*N_dirs + nm*N_dirs + nd] = crmulf( cexpf(cmplxf(0.0f, phase)), irFB_gain);
            }
        }
    }

    /* clean-up */
    free(maxIdx);
    free(centerImpulse);
    free(centerImpulseFB_energy);
    free(centerImpulseFB);
    free(ir);
    free(irFB);
}

