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
 * @brief Slightly modified version of afSTFTlib
 *
 * The original afSTFT code, written by Juha Vilkamo, can be found here:
 * https://github.com/jvilkamo/afSTFT
 * This version is slightly modified. It adds a function to change the number of
 * channels on the fly and includes vectors for the hybrid mode centre
 * frequencies @44.1kHz/48kHz with 128 hop size for convenience.
 * It also supports the use of SAF utilities (for the vectorisation and FFT).
 */

#include "afSTFTlib.h"
#include "afSTFT_protoFilter.h"
#ifdef AFSTFT_USE_SAF_UTILITIES
# include "../../modules/saf_utilities/saf_utilities.h"
#else
/* the vecTools.h/c and fft4g.h/c files, used by the original afSTFTlib
 * may be found here: https://github.com/jvilkamo/afSTFT */
# include "vecTools.h"
#endif

#ifdef AFSTFT_USE_FLOAT_COMPLEX
# ifndef AFSTFT_USE_SAF_UTILITIES
# error You must use SAF utilities to make afSTFT output float_complex
# endif
#endif

/* Internal function prototypes */

void afHybridInit(void** handle, int hopSize, int inChannels, int outChannels);

#ifdef AFSTFT_USE_FLOAT_COMPLEX
void afHybridForward(void* handle, float_complex**);
#else
void afHybridForward(void* handle, complexVector* FD);
#endif

#ifdef AFSTFT_USE_FLOAT_COMPLEX
void afHybridInverse(void* handle, float_complex**);
#else
void afHybridInverse(void* handle, complexVector* FD);
#endif

void afHybridFree(void* handle);

/* Coefficients for a half-band filter, i.e., the "hybrid filter" applied optionally at the bands 1--4. */
#define COEFF1 0.031273141818515176604f
#define COEFF2 0.28127313041521179171f
#define COEFF3 0.5f

/* centre frequencies */
const double __afCenterFreq48e3[133] =
{    0.000000000, 140.644316361, 234.355478108, 328.144332285, 421.855497937, 515.644326841, 609.355515147, 703.144330614, 796.855543885, 937.500032020, 1125.000017338, 1312.500035449, 1500.000075751, 1687.500031782, 1875.000024239, 2062.499975101, 2250.000053703, 2437.500044271, 2625.000002315, 2812.500019782, 3000.000041692, 3187.499983930, 3374.999995137, 3562.499994173, 3750.000018557, 3937.500021643, 4125.000009859, 4312.500011528, 4500.000010423, 4687.500014446, 4875.000013588, 5062.500013570, 5250.000007575, 5437.500010288, 5625.000004178, 5812.500003421, 6000.000005158, 6187.500003404, 6375.000003488, 6562.500007191, 6750.000005972, 6937.500008499, 7125.000006936, 7312.500008549, 7500.000005032, 7687.500004875, 7875.000004878, 8062.500007586, 8250.000006218, 8437.499999805, 8625.000000113, 8812.499997984, 9000.000008860, 9187.500004401, 9375.000001529, 9562.500006565, 9750.000006335, 9937.499999557, 10125.000002928, 10312.500002384, 10500.000004406, 10687.500002820, 10875.000001403, 11062.500002219, 11250.000001097, 11437.500001292, 11625.000000815, 11812.500000140, 12000.000000000, 12187.499999584, 12374.999999473, 12562.499999294, 12749.999998799, 12937.499997514, 13124.999998543, 13312.499997602, 13499.999995904, 13687.499996961, 13874.999996550, 14062.500000495, 14249.999993960, 14437.499993440, 14624.999997861, 14812.499995461, 14999.999991137, 15187.500001756, 15374.999999428, 15562.500000999, 15749.999993809, 15937.499992382, 16124.999995683, 16312.499995240, 16499.999994365, 16687.499991354, 16874.999992234, 17062.499991361, 17249.999994298, 17437.499992410, 17624.999995960, 17812.499995945, 17999.999994836, 18187.499996913, 18374.999996125, 18562.499990092, 18749.999991865, 18937.499986965, 19124.999985762, 19312.499985261, 19499.999989766, 19687.499988292, 19874.999989851, 20062.499978542, 20249.999981602, 20437.500005879, 20625.000004853, 20812.500015815, 20999.999958305, 21187.499980259, 21374.999997733, 21562.499955794, 21749.999946298, 21937.500025004, 22124.999975461, 22312.499968567, 22499.999924162, 22687.499964503, 22874.999982475, 23062.499968048, 23249.999976609, 23437.499982579, 23624.999922020, 23812.499893152, 24000.000000000    };

const double __afCenterFreq44100[133] =
{    0.000000000, 129.216965656, 215.314095512, 301.482605287, 387.579738729, 473.748225285, 559.845379541, 646.013853751, 732.111030944, 861.328154418, 1033.593765929, 1205.859407569, 1378.125069596, 1550.390654200, 1722.656272269, 1894.921852124, 2067.187549340, 2239.453165674, 2411.718752127, 2583.984393174, 2756.250038304, 2928.515610236, 3100.781245532, 3273.046869646, 3445.312517049, 3617.578144885, 3789.843759058, 3962.109385592, 4134.375009576, 4306.640638272, 4478.906262484, 4651.171887467, 4823.437506959, 4995.703134452, 5167.968753839, 5340.234378143, 5512.500004739, 5684.765628127, 5857.031253205, 6029.296881607, 6201.562505487, 6373.828132809, 6546.093756373, 6718.359382855, 6890.625004623, 7062.890629479, 7235.156254481, 7407.421881970, 7579.687505713, 7751.953124821, 7924.218750103, 8096.484373148, 8268.750008140, 8441.015629043, 8613.281251405, 8785.546881031, 8957.812505821, 9130.078124593, 9302.343752690, 9474.609377190, 9646.875004048, 9819.140627591, 9991.406251289, 10163.671877038, 10335.937501008, 10508.203126187, 10680.468750748, 10852.734375129, 11025.000000000, 11197.265624618, 11369.531249516, 11541.796874351, 11714.062498897, 11886.328122716, 12058.593748662, 12230.859372797, 12403.124996237, 12575.390622207, 12747.656246830, 12919.921875455, 13092.187494451, 13264.453118973, 13436.718748035, 13608.984370830, 13781.249991857, 13953.515626614, 14125.781249475, 14298.046875918, 14470.312494312, 14642.578118001, 14814.843746034, 14987.109370627, 15159.374994823, 15331.640617057, 15503.906242865, 15676.171867063, 15848.437494762, 16020.703118027, 16192.968746288, 16365.234371274, 16537.499995256, 16709.765622164, 16882.031246440, 17054.296865897, 17226.562492526, 17398.828113024, 17571.093736919, 17743.359361459, 17915.624990597, 18087.890614243, 18260.156240676, 18432.421855285, 18604.687483097, 18776.953130401, 18949.218754459, 19121.484389530, 19293.749961692, 19466.015606863, 19638.281247918, 19810.546834386, 19982.812450661, 20155.078147972, 20327.343727455, 20499.609346121, 20671.874930324, 20844.140592387, 21016.406233899, 21188.671845644, 21360.937478510, 21533.203108994, 21705.468678356, 21877.734276834, 22050.000000000    };

/**
 * Main data structure for afSTFTlib
 */
typedef struct{
    int inChannels;
    int outChannels; 
    int hopSize;
    int hLen; 
    int LDmode;
    int hopIndexIn;
    int hopIndexOut;
    int totalHops;
    float *protoFilter;
    float *protoFilterI;
    float **inBuffer;
    float *fftProcessFrameTD;
    float **outBuffer;
#ifdef AFSTFT_USE_SAF_UTILITIES
    void* hSafFFT;
    float_complex *fftProcessFrameFD;
    float* tempHopBuffer;
#else
    int pr;
    int log2n;
    void *vtFFT;
    float *fftProcessFrameFD;
#endif
    void *h_afHybrid;
    int hybridMode;
} afSTFT;

/**
 * Data structure for the hybrid filtering employed by afSTFTlib.
 *
 * The purpose of this filtering is to further divide the 4 lowest FFT-bins,
 * to improve the frequency resolution at low-frequencies. For example, 129
 * bins would become 133 hybrid-bins.
 */
typedef struct{
    int inChannels;
    int outChannels;
    int hopSize;
    float hybridCoeffs[3];
    complexVector **analysisBuffer;
    int loopPointer;
#ifdef AFSTFT_USE_FLOAT_COMPLEX
    complexVector *FDtmp;
#endif
} afHybrid;

typedef struct {
    int nBands;
    int nSamples;
    int nHops;
    int hopSize;
    int maxCh;
    float **TDptrs;
    float_complex** FDtmp;
    void* a;
} afMatrix;

void afSTFTinit(void** handle, int hopSize, int inChannels, int outChannels, int LDmode, int hybridMode)
{
    int k, ch, dsFactor;
    float eq;
    
    *handle = malloc(sizeof(afSTFT));
    
    /* Basic initializations */
    afSTFT *h = (afSTFT*)(*handle);
    
    h->inChannels = inChannels;
    h->outChannels = outChannels;
    h->hopSize = hopSize;
    dsFactor=1024/hopSize;
    h->hLen = 10240/dsFactor;
    h->totalHops=10;
    h->hopIndexIn=0;
    h->hopIndexOut=0;
    h->LDmode = LDmode;
    h->protoFilter = (float*)malloc(sizeof(float)*h->hLen);
    h->protoFilterI = (float*)malloc(sizeof(float)*h->hLen);
    h->inBuffer = (float**)malloc(sizeof(float*)*h->inChannels);
    h->outBuffer = (float**)malloc(sizeof(float*)*h->outChannels);
    h->fftProcessFrameTD = (float*)calloc(sizeof(float),h->hopSize*2);
#ifdef AFSTFT_USE_SAF_UTILITIES
    saf_rfft_create(&(h->hSafFFT), h->hopSize*2);
    h->fftProcessFrameFD  = calloc((h->hopSize+1), sizeof(float_complex));
    h->tempHopBuffer = malloc(h->hopSize*sizeof(float));
#else
    switch (hopSize) {
        case 32:
            h->log2n=6;
            break;
        case 64:
            h->log2n=7;
            break;
        case 128:
            h->log2n=8;
            break;
        case 256:
            h->log2n=9;
            break;
        case 512:
            h->log2n=10;
            break;
        case 1024:
            h->log2n=11;
            break;
        default:
            // No other modes defined
            return;
            break;
    }
    h->fftProcessFrameFD  = (float*)calloc(sizeof(float),(h->hopSize+1)*2);
    vtInitFFT(&(h->vtFFT), h->fftProcessFrameTD, h->fftProcessFrameFD, h->log2n);
#endif
    
    /* Normalization to ensure 0dB gain */
    if (h->LDmode==0) {
#ifdef AFSTFT_USE_SAF_UTILITIES
        eq = 2.0f/sqrtf(5.487604141f);
#else
        eq = 1.0f/sqrtf((float)h->hopSize*5.487604141f);
#endif
        
        for (k=0; k<h->hLen; k++) {
            h->protoFilter[h->hLen-k-1] = protoFilter1024[k*dsFactor]*eq;
            h->protoFilterI[h->hLen-k-1] = protoFilter1024[k*dsFactor]*eq;
        }
    }
    else
    {
#ifdef AFSTFT_USE_SAF_UTILITIES
        eq = 2.0f/sqrtf(4.544559956f);
#else
        eq = 1.0f/sqrtf((float)h->hopSize*4.544559956f);
#endif
        for (k=0; k<h->hLen; k++) {
            h->protoFilter[h->hLen-k-1] = protoFilter1024LD[k*dsFactor]*eq;
            h->protoFilterI[k]=protoFilter1024LD[k*dsFactor]*eq; 
        }
    }
    for(ch=0;ch<h->inChannels;ch++)
        h->inBuffer[ch] = (float*)calloc(h->hLen,sizeof(float));
    
    for(ch=0;ch<h->outChannels;ch++)
        h->outBuffer[ch] = (float*)calloc(h->hLen,sizeof(float));
    
    /* Initialize the hybrid filter memory etc. */
    h->hybridMode=hybridMode;
    if (h->hybridMode)
        afHybridInit(&(h->h_afHybrid), h->hopSize, h->inChannels,h->outChannels);
}

void afSTFTchannelChange(void* handle, int new_inChannels, int new_outChannels)
{
    afSTFT *h = (afSTFT*)(handle);
    afHybrid *hyb_h = h->h_afHybrid;
    
    int i, ch, sample;
    if(h->inChannels!=new_inChannels){
        for(i=new_inChannels; i<h->inChannels; i++)
            free(h->inBuffer[i]);
        h->inBuffer = (float**)realloc(h->inBuffer, sizeof(float*)*new_inChannels);
        for(i=h->inChannels; i<new_inChannels; i++)
            h->inBuffer[i] = (float*)calloc(h->hLen,sizeof(float));
    }
    
    if(h->outChannels!=new_outChannels){
        for(i=new_outChannels; i<h->outChannels; i++)
            free(h->outBuffer[i]);
        h->outBuffer = (float**)realloc(h->outBuffer, sizeof(float*)*new_outChannels);
        for(i=h->outChannels; i<new_outChannels; i++)
            h->outBuffer[i] = (float*)calloc(h->hLen,sizeof(float));
    }
    
    if (h->hybridMode) {
        hyb_h = h->h_afHybrid;
        if (hyb_h->inChannels != new_inChannels) {
            for (ch = new_inChannels; ch < hyb_h->inChannels; ch++) {
                for (sample = 0; sample < 7; sample++) {
                    free(hyb_h->analysisBuffer[ch][sample].re);
                    free(hyb_h->analysisBuffer[ch][sample].im);
                }
                free(hyb_h->analysisBuffer[ch]);
            }
            hyb_h->analysisBuffer = (complexVector**) realloc(hyb_h->analysisBuffer, sizeof(complexVector*) * new_inChannels);
            for (ch = hyb_h->inChannels; ch < new_inChannels; ch++) {
                hyb_h->analysisBuffer[ch] = (complexVector*) malloc(sizeof(complexVector) * 7);
                for (sample = 0; sample < 7; sample++) {
                    hyb_h->analysisBuffer[ch][sample].re = (float*) calloc(sizeof(float), h->hopSize + 1);
                    hyb_h->analysisBuffer[ch][sample].im = (float*) calloc(sizeof(float), h->hopSize + 1);
                }
            }
        }
#ifdef AFSTFT_USE_FLOAT_COMPLEX
        if (hyb_h->inChannels != new_inChannels || hyb_h->outChannels != new_outChannels) {
            int old_nCh = hyb_h->inChannels < hyb_h->outChannels ? hyb_h->outChannels : hyb_h->inChannels;
            int new_nCh = new_inChannels < new_outChannels ? new_outChannels : new_inChannels;
            for (ch = new_nCh; ch < old_nCh; ch++) {
                free(hyb_h->FDtmp[ch].re);
                free(hyb_h->FDtmp[ch].im);
            }
            hyb_h->FDtmp = (complexVector*) realloc(hyb_h->FDtmp,
                sizeof(complexVector) * new_nCh);
            for (ch = old_nCh; ch < new_nCh; ch++) {
                hyb_h->FDtmp[ch].re = (float*) calloc(h->hopSize + 5, sizeof(float));
                hyb_h->FDtmp[ch].im = (float*) calloc(h->hopSize + 5, sizeof(float));
            }
        }
#endif
    }
    h->inChannels = new_inChannels;
    h->outChannels = new_outChannels;
    if (h->hybridMode){
        hyb_h->inChannels = new_inChannels;
        hyb_h->outChannels = new_outChannels;
    }
}

void afSTFTclearBuffers(void* handle)
{
    afSTFT *h = (afSTFT*)(handle);
    afHybrid *hyb_h = h->h_afHybrid;
    int i, ch, sample;
    
    for(i=0; i<h->inChannels; i++)
        memset(h->inBuffer[i], 0, h->hLen*sizeof(float));
    for(i=0; i<h->outChannels; i++)
        memset(h->outBuffer[i], 0, h->hLen*sizeof(float));
    if (h->hybridMode){
        for(ch=0; ch<hyb_h->inChannels; ch++) {
            for (sample=0;sample<7;sample++) {
                memset(hyb_h->analysisBuffer[ch][sample].re, 0, sizeof(float)*(h->hopSize+1));
                memset(hyb_h->analysisBuffer[ch][sample].im, 0, sizeof(float)*(h->hopSize+1));

            }
        }
#ifdef AFSTFT_USE_FLOAT_COMPLEX
        // TODO is this actually necessary?
        int nCh = h->inChannels < h->outChannels ? h->outChannels : h->inChannels;
        for (ch = 0; ch < nCh; ch++) {
            memset(hyb_h->FDtmp[ch].re, 0, sizeof(float) * (h->hopSize + 5));
            memset(hyb_h->FDtmp[ch].im, 0, sizeof(float) * (h->hopSize + 5));
        }
#endif
    }
}

#ifdef AFSTFT_USE_FLOAT_COMPLEX
void afSTFTforward(void* handle, float** inTD, float_complex** outFD)
#else
void afSTFTforward(void* handle, float** inTD, complexVector* outFD)
#endif
{
    afSTFT *h = (afSTFT*)(handle);
    int ch,k,j,hopIndex_this,hopIndex_this2;
    float *p1,*p2,*p3;
#ifndef AFSTFT_USE_SAF_UTILITIES
    float *p4;
#endif
    int lr;
    
    for (ch=0;ch<h->inChannels;ch++)
    {
        /* Copy the input frame into the memory buffer */
        hopIndex_this2 = h->hopIndexIn;
        p1=&(h->inBuffer[ch][hopIndex_this2*h->hopSize]);
        p2=inTD[ch];
        memcpy((void*)p1,(void*)p2,sizeof(float)*(h->hopSize));
        
        hopIndex_this2++;
        if (hopIndex_this2 >= h->totalHops)
        {
            hopIndex_this2 = 0;
        }
        
        /* Apply prototype filter to the collected data in the memory buffer, and fold the result (for the FFT operation). */
        p1 = h->fftProcessFrameTD;
#ifdef AFSTFT_USE_SAF_UTILITIES
        memset(p1, 0, h->hopSize*2*sizeof(float));
#else
        vtClr(p1, h->hopSize*2);
#endif
        lr=0; /* Left or right part of the frame */
        hopIndex_this = hopIndex_this2;
        for (k=0;k<h->totalHops;k++)
        {
            p1=&(h->inBuffer[ch][h->hopSize*hopIndex_this]);
            p2=&(h->protoFilter[k*h->hopSize]);
            if (lr==1)
            {
                p3=&(h->fftProcessFrameTD[h->hopSize]);
                lr=0;
            }
            else
            {
                p3=&(h->fftProcessFrameTD[0]);
                lr=1;
            }
#ifdef AFSTFT_USE_SAF_UTILITIES
//            utility_svvmul(p1, p2, h->hopSize, h->tempHopBuffer);
//            utility_svvadd(p3, h->tempHopBuffer, h->hopSize, p3);
            for (j=0;j<h->hopSize;j++)
                p3[j] += (p1[j])*(p2[j]);
#else
            vtVma(p1, p2, p3, h->hopSize);  /* Vector multiply-add */
#endif
            hopIndex_this++;
            if (hopIndex_this >= h->totalHops)
            {
                hopIndex_this = 0;
            }
        }
        
        /* Apply FFT and copy the data to the output vector */
#ifdef AFSTFT_USE_SAF_UTILITIES
# ifdef AFSTFT_USE_FLOAT_COMPLEX
        saf_rfft_forward(h->hSafFFT, h->fftProcessFrameTD, outFD[ch]);
# else
        saf_rfft_forward(h->hSafFFT, h->fftProcessFrameTD, h->fftProcessFrameFD);
        for(k = 0; k<h->hopSize+1; k++){
            outFD[ch].re[k] = crealf(h->fftProcessFrameFD[k]);
            outFD[ch].im[k] = cimagf(h->fftProcessFrameFD[k]);
        }
# endif // AFSTFT_USE_FLOAT_COMPLEX
#else
        vtRunFFT(h->vtFFT,1);
        outFD[ch].re[0]=h->fftProcessFrameFD[0];
        outFD[ch].im[0]=0.0f; /* DC im = 0 */
        outFD[ch].re[h->hopSize]=h->fftProcessFrameFD[h->hopSize];
        outFD[ch].im[h->hopSize]=0.0f; /* Nyquist im = 0 */
        p1 = outFD[ch].re + 1;
        p2 = outFD[ch].im + 1;
        p3 = h->fftProcessFrameFD + 1;
        p4 = h->fftProcessFrameFD + 1 + h->hopSize;
        memcpy((void*)p1,(void*)p3,sizeof(float)*(h->hopSize - 1));
        memcpy((void*)p2,(void*)p4,sizeof(float)*(h->hopSize - 1));
#endif
    }
    h->hopIndexIn++;
    if (h->hopIndexIn >= h->totalHops)
    {
        h->hopIndexIn = 0;
    }
    
    /* Subdivide lowest bands with half-band filters if hybrid mode is enabled */
    if (h->hybridMode)
    {
        afHybridForward(h->h_afHybrid, outFD);
    }
}

#ifdef AFSTFT_USE_FLOAT_COMPLEX
void afSTFTinverse(void* handle, float_complex** inFD, float** outTD)
#else
void afSTFTinverse(void* handle, complexVector* inFD, float** outTD)
#endif
{
    afSTFT *h = (afSTFT*)(handle);
    int ch,k,j,hopIndex_this,hopIndex_this2;
    float *p1,*p2,*p3;
#ifndef AFSTFT_USE_SAF_UTILITIES
    float *p4;
#endif
    int lr;
    
    /* Combine subdivided lowest bands if hybrid mode is enabled */
    if (h->hybridMode)
    {
        afHybridInverse(h->h_afHybrid, inFD);
    }
    
    for (ch=0;ch<h->outChannels;ch++)
    {
        /* Copy data from input to internal memory */
        hopIndex_this2 = h->hopIndexOut;
        
        /* Inverse FFT */
#ifdef AFSTFT_USE_SAF_UTILITIES
# ifdef AFSTFT_USE_FLOAT_COMPLEX
        utility_cvvcopy(inFD[ch], h->hopSize + 1, h->fftProcessFrameFD);
# else
        for(k = 0; k<h->hopSize+1; k++)
            h->fftProcessFrameFD[k] = cmplxf(inFD[ch].re[k], inFD[ch].im[k]);
# endif
        
        /* The low delay mode requires this procedure corresponding to the circular shift of the data in the time domain */
        if (h->LDmode == 1)
            for (k=1; k<h->hopSize; k+=2)
                h->fftProcessFrameFD[k] = crmulf(h->fftProcessFrameFD[k], -1.0f);
        
        saf_rfft_backward(h->hSafFFT, h->fftProcessFrameFD, h->fftProcessFrameTD);
#else
        h->fftProcessFrameFD[0] = inFD[ch].re[0]; /* DC */
        h->fftProcessFrameFD[h->hopSize] = inFD[ch].re[h->hopSize]; /* Nyquist */
        p1 = inFD[ch].re + 1;
        p2 = inFD[ch].im + 1;
        p3 = h->fftProcessFrameFD + 1;
        p4 = h->fftProcessFrameFD + 1 + h->hopSize;
        memcpy((void*)p3,(void*)p1,sizeof(float)*(h->hopSize - 1));
        memcpy((void*)p4,(void*)p2,sizeof(float)*(h->hopSize - 1));
        
        /* The low delay mode requires this procedure corresponding to the circular shift of the data in the time domain */
        if (h->LDmode == 1) {
            for (k=1;k<h->hopSize;k+=2) {
                *p3 = -*p3;
                *p4 = -*p4;
                p3+=2;
                p4+=2;
            }
        }
        
        vtRunFFT(h->vtFFT, -1);
#endif
        
        /* Clear buffer at the pointer location and increment the pointer */
        p1 = &(h->outBuffer[ch][hopIndex_this2*h->hopSize]);
#ifdef AFSTFT_USE_SAF_UTILITIES
        memset(p1, 0, h->hopSize*sizeof(float));
#else
        vtClr(p1,h->hopSize);
#endif
        hopIndex_this2++;
        if (hopIndex_this2 >= h->totalHops)
        {
            hopIndex_this2=0;
        }
        hopIndex_this = hopIndex_this2;
        
        lr=0; /* Left or right part of the frame */
        for (k=0;k<h->totalHops;k++)
        {
            /* Apply the prototype filter to the repeated version of the IFFT'd data. */
            p1=&(h->outBuffer[ch][h->hopSize*hopIndex_this]);
            p2=&(h->protoFilterI[k*h->hopSize]);
            
            if (lr==1)
            {
                p3=&(h->fftProcessFrameTD[h->hopSize]);
                lr=0;
            }
            else
            {
                p3=&(h->fftProcessFrameTD[0]);
                lr=1;
            }
 
            /* Overlap-add to the existing data in the memory buffer (from previous frames). */
#ifdef AFSTFT_USE_SAF_UTILITIES
//            utility_svvmul(p2, p3, h->hopSize, h->tempHopBuffer);
//            utility_svvadd(p1, h->tempHopBuffer, h->hopSize, p1);
            for (j=0;j<h->hopSize;j++)
                p1[j] += (p2[j])*(p3[j]); 
#else
            vtVma(p2, p3, p1, h->hopSize); /* Vector multiply-add */
#endif
            
            hopIndex_this++;
            if (hopIndex_this >= h->totalHops)
            {
                hopIndex_this = 0;
            }
        }
        
        /* Copy a frame from work memory to the output */
        p1 = outTD[ch];
        p2 = &(h->outBuffer[ch][h->hopSize*hopIndex_this]);
        memcpy((void*)p1,(void*)p2,sizeof(float)*(h->hopSize));
        
    }
    h->hopIndexOut++;
    if (h->hopIndexOut >= h->totalHops)
    {
        h->hopIndexOut=0;
    }
    
}

void afSTFTfree(void* handle)
{
    afSTFT *h = (afSTFT*)(handle);
    int ch;
    if (h->hybridMode)
    {
        afHybridFree(h->h_afHybrid);
    }
    for(ch=0;ch<h->inChannels;ch++)
    {
        free(h->inBuffer[ch]);
    }
    
    for(ch=0;ch<h->outChannels;ch++)
    {
        free(h->outBuffer[ch]);
    }
    free(h->protoFilter);
    free(h->protoFilterI);
    free(h->inBuffer);
    free(h->outBuffer);
    free(h->fftProcessFrameTD);
    free(h->fftProcessFrameFD);
#ifdef AFSTFT_USE_SAF_UTILITIES
    saf_rfft_destroy(&(h->hSafFFT));
    free(h->tempHopBuffer);
#else
    vtFreeFFT(h->vtFFT);
#endif
    free(h);
}

void afHybridInit(void** handle, int hopSize, int inChannels, int outChannels)
{
    /* Allocates 7 samples of memory for FIR filtering at lowest bands, and for delays at other bands. */
    int ch,sample;
    *handle = malloc(sizeof(afHybrid));
    afHybrid *h = (afHybrid*)(*handle);
    h->inChannels = inChannels;
    h->hopSize = hopSize;
    h->outChannels = outChannels;
    h->analysisBuffer = (complexVector**)malloc(sizeof(complexVector*)*h->inChannels);
    h->loopPointer=0;
    for (ch=0;ch<h->inChannels;ch++)
    {
        h->analysisBuffer[ch] = (complexVector*)malloc(sizeof(complexVector)*7);
        for (sample=0;sample<7;sample++)
        {
            h->analysisBuffer[ch][sample].re=(float*)calloc(sizeof(float),h->hopSize+1);
            h->analysisBuffer[ch][sample].im=(float*)calloc(sizeof(float),h->hopSize+1);
        }
    }
#ifdef AFSTFT_USE_FLOAT_COMPLEX
    int nCh = inChannels < outChannels ? outChannels : inChannels;
    h->FDtmp = (complexVector*) malloc(nCh * sizeof(complexVector));
    for (ch = 0; ch < nCh; ch++) {
        h->FDtmp[ch].re = (float*) calloc((h->hopSize + 5), sizeof(float));
        h->FDtmp[ch].im = (float*) calloc((h->hopSize + 5), sizeof(float));
    }
#endif
}

#ifdef AFSTFT_USE_FLOAT_COMPLEX
void afHybridForward(void* handle, float_complex** FDin)
#else
void afHybridForward(void* handle, complexVector* FD)
#endif
{
    afHybrid *h = (afHybrid*)(handle);
    int ch,band,sample,realImag;
    float *pr1, *pr2, *pi1, *pi2;
    float re,im;
    int sampleIndices[7];
    int loopPointerThis;
    h->loopPointer++;
    if( h->loopPointer == 7)
    {
        h->loopPointer = 0;
    }
#ifdef AFSTFT_USE_FLOAT_COMPLEX
    complexVector *FD = h->FDtmp;
    int b;
    for (ch = 0; ch < h->inChannels; ch++) {
        for (b = 0; b < h->hopSize + 5; b++) {
            FD[ch].re[b] = crealf(FDin[ch][b]);
            FD[ch].im[b] = cimagf(FDin[ch][b]);
        }
    }
#endif
    for (ch=0;ch<h->inChannels;ch++)
    {
        /* Copy data from input to the memory buffer */
        pr1 = FD[ch].re;
        pi1 = FD[ch].im;
        pr2 = h->analysisBuffer[ch][h->loopPointer].re;
        pi2 = h->analysisBuffer[ch][h->loopPointer].im;
        memcpy((void*)pr2,(void*)pr1,sizeof(float)*(h->hopSize+1));
        memcpy((void*)pi2,(void*)pi1,sizeof(float)*(h->hopSize+1));
        
        /* Get the pointer to a position corresponding to the group delay of the linear-phase half-band filter. */
        loopPointerThis = h->loopPointer - 3;
        if( loopPointerThis < 0)
        {
            loopPointerThis += 7;
        }
        pr1 = FD[ch].re;
        pr2 = h->analysisBuffer[ch][loopPointerThis].re;
        for (realImag=0;realImag<2;realImag++)
        {
            /* The 0.5 multipliers are the center coefficients of the half-band FIR filters. Data is duplicated for the half-bands. */
            *pr1 = *pr2;
            *(pr1+1) = *(pr2+1)*0.5f;
            *(pr1+2) = *(pr1+1);
            *(pr1+3) = *(pr2+2)*0.5f;
            *(pr1+4) = *(pr1+3);
            *(pr1+5) = *(pr2+3)*0.5f;
            *(pr1+6) = *(pr1+5);
            *(pr1+7) = *(pr2+4)*0.5f;
            *(pr1+8) = *(pr1+7);
            
            /* The rest of the bands are shifted upwards in the frequency indices, and delayed by the group delay of the half-band filters */
            memcpy((void*)(pr1+9),(void*)(pr2+5),sizeof(float)*(h->hopSize-4));
            
            /* Repeat process for the imaginary part, at next iteration. */
            pr1 = FD[ch].im;
            pr2 = h->analysisBuffer[ch][loopPointerThis].im;
        }
    
        for (sample=0;sample<7;sample++)
        {
            sampleIndices[sample]=h->loopPointer+1+sample;
            if(sampleIndices[sample] > 6)
            {
                sampleIndices[sample]-=7;
            }
            
        }
        for (band=1; band<5; band++)
        {
            /* The rest of the half-band FIR filtering is implemented below. The real<->imaginary shifts are for shifting the half-band filter spectra. */
            re = -COEFF1*h->analysisBuffer[ch][sampleIndices[6]].im[band];
            im =  COEFF1*h->analysisBuffer[ch][sampleIndices[6]].re[band];
            re -= COEFF2*h->analysisBuffer[ch][sampleIndices[4]].im[band];
            im += COEFF2*h->analysisBuffer[ch][sampleIndices[4]].re[band];
            re += COEFF2*h->analysisBuffer[ch][sampleIndices[2]].im[band];
            im -= COEFF2*h->analysisBuffer[ch][sampleIndices[2]].re[band];
            re += COEFF1*h->analysisBuffer[ch][sampleIndices[0]].im[band];
            im -= COEFF1*h->analysisBuffer[ch][sampleIndices[0]].re[band];
            
            /* The addition or subtraction process below provides the upper and lower half-band spectra (the coefficient 0.5 had the same sign for both bands).
               The half-band orders are switched for bands=1,3 with respect to band=2,4, because of the organization of the spectral data at the downsampled frequency band signals. As the result of the order switching, the bands are organized by the ascending spectral position. */
            if (band == 1 || band== 3)
            {
                FD[ch].re[band*2-1] -= re;
                FD[ch].im[band*2-1] -= im;
                FD[ch].re[band*2] += re;
                FD[ch].im[band*2] += im;
            }
            else
            {
                FD[ch].re[band*2-1] += re;
                FD[ch].im[band*2-1] += im;
                FD[ch].re[band*2] -= re;
                FD[ch].im[band*2] -= im;
            }
            
        }
    }
#ifdef AFSTFT_USE_FLOAT_COMPLEX
    for (ch = 0; ch < h->inChannels; ch++) {
        for (b = 0; b < h->hopSize + 5; b++) {
            FDin[ch][b] = cmplxf(FD[ch].re[b], FD[ch].im[b]);
        }
    }
#endif
}

#ifdef AFSTFT_USE_FLOAT_COMPLEX
void afHybridInverse(void* handle, float_complex** FDin)
#else
void afHybridInverse(void* handle, complexVector* FD)
#endif
{
    afHybrid *h = (afHybrid*)(handle);
    int ch,realImag;
    float *pr;
#ifdef AFSTFT_USE_FLOAT_COMPLEX
    int b;
    complexVector *FD = h->FDtmp;

    for (ch = 0; ch < h->outChannels; ch++) {
        for (b = 0; b < h->hopSize + 5; b++) {
            FD[ch].re[b] = crealf(FDin[ch][b]);
            FD[ch].im[b] = cimagf(FDin[ch][b]);
        }
    }
#endif

    for (ch=0;ch<h->outChannels;ch++)
    {
        pr = FD[ch].re;
        for (realImag=0;realImag<2;realImag++)
        {
            /* Since no downsampling was applied, the inverse hybrid filtering is just sum of the bands */
            *(pr+1) = *(pr+1) + *(pr+2);
            *(pr+2) = *(pr+3) + *(pr+4);
            *(pr+3) = *(pr+5) + *(pr+6);
            *(pr+4) = *(pr+7) + *(pr+8);
            
            /* The rest of the bands are shifted to their original positions */
            memmove((void*)(pr+5),(void*)(pr+9),sizeof(float)*(h->hopSize-4));
            
            /* Repeat process for the imaginary part, at next iteration. */
            pr = FD[ch].im;
        }
    }
#ifdef AFSTFT_USE_FLOAT_COMPLEX
    for (ch = 0; ch < h->outChannels; ch++) {
        for (b = 0; b < h->hopSize + 5; b++) {
            FDin[ch][b] = cmplxf(FD[ch].re[b], FD[ch].im[b]);
        }
    }
#endif
}

void afHybridFree(void* handle)
{
    int ch,sample;
    afHybrid *h = (afHybrid*)(handle);
    for (ch=0;ch<h->inChannels;ch++)
    {
        
        for (sample=0;sample<7;sample++)
        {
            free(h->analysisBuffer[ch][sample].re);
            free(h->analysisBuffer[ch][sample].im);
        }
        free(h->analysisBuffer[ch]);
          
    }
    free(h->analysisBuffer);
#ifdef AFSTFT_USE_FLOAT_COMPLEX
    int nCh = h->inChannels < h->outChannels ? h->outChannels : h->inChannels;
    for (ch = 0; ch < nCh; ch++) {
        free(h->FDtmp[ch].re);
        free(h->FDtmp[ch].im);
    }
    free(h->FDtmp);
#endif
    free(handle);
}

#ifdef AFSTFT_USE_FLOAT_COMPLEX
void afSTFTMatrixInit(void** handle,
                      const int hopSize,
                      const int inChannels,
                      const int outChannels,
                      const int LDmode,
                      const int hybridMode,
                      const int nSamples)
{
    // we only work with nicely-behaved ratios
    if (nSamples % hopSize != 0)
        return;

    *handle = malloc(sizeof(afMatrix));
    afMatrix *h = (afMatrix*) (*handle);
    h->nBands = hopSize + (hybridMode ? 5 : 1);
    h->nSamples = nSamples;
    h->hopSize = hopSize;
    h->nHops = nSamples / hopSize;
    h->maxCh = inChannels > outChannels ? inChannels : outChannels;
    h->TDptrs = malloc(h->maxCh * sizeof(float*));
    h->FDtmp = (float_complex**) malloc2d(h->maxCh, h->nBands, sizeof(float_complex));
    afSTFTinit(&(h->a), h->hopSize, inChannels, outChannels, LDmode, hybridMode);
}

void afSTFTMatrixChannelChange(void* handle, int new_inChannels, int new_outChannels)
{ 
    afMatrix *h = (afMatrix*) handle;
    afSTFTchannelChange(h->a, new_inChannels, new_outChannels);

    h->maxCh = new_inChannels > new_outChannels ? new_inChannels : new_outChannels;
    h->TDptrs = realloc(h->TDptrs, h->maxCh * sizeof(float*));
    free(h->FDtmp);
    h->FDtmp = (float_complex**) malloc2d(h->maxCh, h->nBands, sizeof(float_complex));
}

void afSTFTMatrixFree(void *handle) {
    afMatrix *h = (afMatrix*) handle;
    free(h->FDtmp);
    free(h->TDptrs);
    afSTFTfree(h->a);
    free(handle);
}

void afSTFTMatrixForward(void* handle, float** inTD, float_complex*** outFD)
{
    afMatrix *h = (afMatrix*) handle;
    afSTFT *a = (afSTFT*)(h->a);

    int nCh = a->inChannels;
    int hop, ch, b;

    for (hop = 0; hop < h->nHops; hop++) {
        for (ch = 0; ch < nCh; ch++) {
            h->TDptrs[ch] = inTD[ch] + hop * h->hopSize;
        }
        afSTFTforward(h->a, h->TDptrs, h->FDtmp);
        for (ch = 0; ch < nCh; ch++)
            for (b = 0; b < h->nBands; b++)
                outFD[b][ch][hop] = h->FDtmp[ch][b];
     }
}

void afSTFTMatrixInverse(void* handle, float_complex*** inFD, float** outTD)
{
    afMatrix *h = (afMatrix*) handle;
    afSTFT *a = (afSTFT*)(h->a);

    int nCh = a->outChannels;
    int hop, ch, b;

    for (hop = 0; hop < h->nHops; hop++){
        for (ch = 0; ch < nCh; ch++)
            for (b = 0; b < h->nBands; b++)
                h->FDtmp[ch][b] = inFD[b][ch][hop];
        for (ch = 0; ch < nCh; ch++) {
            h->TDptrs[ch] = outTD[ch] + hop * h->hopSize;
        }
        afSTFTinverse(h->a, h->FDtmp, h->TDptrs);
    }
    //memset(ADR2D(outTD), 0, h->hopSize * h->nHops * nCh*sizeof(float));
}
#endif
