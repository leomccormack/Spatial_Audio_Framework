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

#include "afSTFTlib.h"
#include "vecTools.h"
#include "afSTFT_protoFilter.h"

/* Coefficients for a half-band filter, i.e., the "hybrid filter" applied optionally at the bands 1--4. */
#define COEFF1 0.031273141818515176604f
#define COEFF2 0.28127313041521179171f
#define COEFF3 0.5f

void afSTFTinit(void** handle, int hopSize, int inChannels, int outChannels, int LDmode, int hybridMode)
{
    int k,ch,log2n, dsFactor;
    float eq;
    
    *handle = malloc(sizeof(afSTFT));
    switch (hopSize) {
        case 32:
            log2n=6;
            break;
        case 64:
            log2n=7;
            break;
        case 128:
            log2n=8;
            break;
        case 256:
            log2n=9;
            break;
        case 512:
            log2n=10;
            break;
        case 1024:
            log2n=11;
            break;
        default:
            // No other modes defined
            return;
            break;
    }
    
    /* Basic initializations */
    afSTFT *h = (afSTFT*)(*handle);
    h->log2n=log2n;
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
    h->inBuffer = (float**)malloc(sizeof(float*)* h->inChannels);
    h->outBuffer = (float**)malloc(sizeof(float*)* h->outChannels);
    h->fftProcessFrameTD = (float*)calloc(sizeof(float),h->hopSize*2);
    h->fftProcessFrameFD  = (float*)calloc(sizeof(float),(h->hopSize+1)*2);
	//h->fftProcessFrameTD_mkl = (float_complex*)calloc(sizeof(float_complex), h->hopSize * 2);
	//h->fftProcessFrameFD_mkl = (float_complex*)calloc(sizeof(float_complex), 2*(h->hopSize + 1));
    //vtInitFFT(&(h->vtFFT),h->fftProcessFrameTD, h->fftProcessFrameFD, h->fftProcessFrameTD_mkl, h->fftProcessFrameFD_mkl, h->log2n);
    vtInitFFT(&(h->vtFFT),h->fftProcessFrameTD, h->fftProcessFrameFD, h->log2n);
    
    /* Normalization to ensure 0dB gain */
    if (h->LDmode==0){
        eq = 1.0f/sqrtf((float)h->hopSize*5.487604141f);
        for (k=0; k<h->hLen; k++){
            h->protoFilter[h->hLen-k-1] = protoFilter1024[k*dsFactor]*eq;
            h->protoFilterI[h->hLen-k-1] = protoFilter1024[k*dsFactor]*eq;
        }
    }
    else{
        eq = 1.0f/sqrtf((float)h->hopSize*4.544559956f);
        for (k=0; k<h->hLen; k++){
            h->protoFilter[h->hLen-k-1] = protoFilter1024LD[k*dsFactor]*eq;
            h->protoFilterI[k]=protoFilter1024LD[k*dsFactor]*eq;
            
        }
    }
    for(ch=0;ch< h->inChannels;ch++)
        h->inBuffer[ch] = (float*)calloc(h->hLen,sizeof(float));
    
    for(ch=0;ch< h->outChannels;ch++)
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
    
    if (h->hybridMode){
        hyb_h = h->h_afHybrid;
        if(hyb_h->inChannels!=new_inChannels){
            for (ch=new_inChannels; ch<hyb_h->inChannels; ch++) {
                for (sample=0;sample<7;sample++) {
                    free(hyb_h->analysisBuffer[ch][sample].re);
                    free(hyb_h->analysisBuffer[ch][sample].im);
                }
                free(hyb_h->analysisBuffer[ch]);
            }
            hyb_h->analysisBuffer = (complexVector**)realloc(hyb_h->analysisBuffer, sizeof(complexVector*)*new_inChannels);
            for (ch=hyb_h->inChannels; ch<new_inChannels; ch++) {
                hyb_h->analysisBuffer[ch] = (complexVector*)malloc(sizeof(complexVector)*7);
                for (sample=0;sample<7;sample++) {
                    hyb_h->analysisBuffer[ch][sample].re=(float*)calloc(sizeof(float),h->hopSize+1);
                    hyb_h->analysisBuffer[ch][sample].im=(float*)calloc(sizeof(float),h->hopSize+1);
                }
            }
        }
    }
    
    h->inChannels = new_inChannels;
    h->outChannels = new_outChannels;
    if (h->hybridMode){
        hyb_h->inChannels = new_inChannels;
        hyb_h->outChannels = new_outChannels;
    }
}

void afSTFTforward(void* handle, float** inTD, complexVector* outFD)
{
    afSTFT *h = (afSTFT*)(handle);
    int ch,k,hopIndex_this,hopIndex_this2;
    float *p1,*p2,*p3,*p4;
    int lr;
    
    for (ch=0;ch<h->inChannels;ch++) {
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
        vtClr(p1, h->hopSize*2);
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
            vtVma(p1, p2, p3, h->hopSize);  /* Vector multiply-add */
            hopIndex_this++;
            if (hopIndex_this >= h->totalHops)
            {
                hopIndex_this = 0;
            }
        }

        /* Apply FFT and copy the data to the output vector */
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

void afSTFTinverse(void* handle, complexVector* inFD, float** outTD)
{
    afSTFT *h = (afSTFT*)(handle);
    int ch,k,hopIndex_this,hopIndex_this2;
    float *p1,*p2,*p3,*p4;
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
        h->fftProcessFrameFD[0] = inFD[ch].re[0]; /* DC */
        h->fftProcessFrameFD[h->hopSize] = inFD[ch].re[h->hopSize]; /* Nyquist */
        p1 = inFD[ch].re + 1;
        p2 = inFD[ch].im + 1;
        p3 = h->fftProcessFrameFD + 1;
        p4 = h->fftProcessFrameFD + 1 + h->hopSize;
        memcpy((void*)p3,(void*)p1,sizeof(float)*(h->hopSize - 1));
        memcpy((void*)p4,(void*)p2,sizeof(float)*(h->hopSize - 1));
        
        /* The low delay mode requires this procedure corresponding to the circular shift of the data in the time domain */
        if (h->LDmode == 1)
        {
            for (k=1;k<h->hopSize;k+=2)
            {
                *p3 = -*p3;
                *p4 = -*p4;
                p3+=2;
                p4+=2;
            }
        }

        /* Inverse FFT */
        vtRunFFT(h->vtFFT, -1);

        /* Clear buffer at the pointer location and increment the pointer */
        p1 = &(h->outBuffer[ch][hopIndex_this2*h->hopSize]);
        vtClr(p1,h->hopSize);
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
            vtVma(p2, p3, p1, h->hopSize); /* Vector multiply-add */
            
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
    vtFreeFFT(h->vtFFT);
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
    h->analysisBuffer = (complexVector**)malloc(sizeof(complexVector*)* (h->inChannels));
    h->loopPointer=0;
    for (ch=0;ch< h->inChannels;ch++)
    {
        h->analysisBuffer[ch] = (complexVector*)malloc(sizeof(complexVector)*7);
        for (sample=0;sample<7;sample++)
        {
            h->analysisBuffer[ch][sample].re=(float*)calloc(sizeof(float),h->hopSize+1);
            h->analysisBuffer[ch][sample].im=(float*)calloc(sizeof(float),h->hopSize+1);
        }
    }
}

void afHybridForward(void* handle, complexVector* FD)
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
}

void afHybridInverse(void* handle, complexVector* FD)
{
    afHybrid *h = (afHybrid*)(handle);
    int ch,realImag;
    float *pr;

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
}
