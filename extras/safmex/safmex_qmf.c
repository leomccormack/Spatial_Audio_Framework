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
 * @file safmex_qmf.c
 * @brief MEX wrapper for qmf (see the .m file of the same name for
 *        documentation)
 * @author Leo McCormack
 * @date 28.08.2020
 */
 
#include "safmex.h"

/* ===================================================================== */
/*                                Config                                 */
/* ===================================================================== */

#define NUM_INPUT_ARGS_CREATE  ( 7 )
const MEX_DATA_TYPES inputDataTypes_create[NUM_INPUT_ARGS_CREATE] = {
    SM_INT32,                 
    SM_INT32,  
    SM_INT32,  
    SM_INT32,  
    SM_INT32,
    SM_INT32,
    SM_DOUBLE_REAL
}; 
#define NUM_INPUT_ARGS_FWD ( 1 )
#define NUM_OUTPUT_ARGS_FWD ( 1 )
#define NUM_INPUT_ARGS_BKWD ( 1 )
#define NUM_OUTPUT_ARGS_BKWD ( 1 )
const MEX_DATA_TYPES inputDataTypes_fwd[NUM_INPUT_ARGS_FWD] = {
    SM_DOUBLE_REAL_1D_OR_2D      
}; 
const MEX_DATA_TYPES inputDataTypes_bkwd[NUM_INPUT_ARGS_BKWD] = {
    SM_DOUBLE_COMPLEX_3D      
}; 
const MEX_DATA_TYPES outputDataTypes_fwd[NUM_OUTPUT_ARGS_FWD] = {
    SM_DOUBLE_COMPLEX_3D 
}; 
const MEX_DATA_TYPES outputDataTypes_bkwd[NUM_OUTPUT_ARGS_BKWD] = {
    SM_DOUBLE_REAL_1D_OR_2D      
}; 


/* ===================================================================== */
/*                                 Vars                                  */
/* ===================================================================== */

/* user arguments */
int nCHin;          /* Number of input channels */
int nCHout;         /* Number of output channels */
int hopsize;        /* Hop size, in samples */
int blocksize;      /* time domain samples to process at a time */
int hybridmode;     /* 0: disabled, 1: hybrid-filtering enabled */
int formatFlag;     /* 0: nBands x nCH x time, 1: time x nCH x nBands */
float fs;           /* sampling rate */

/* internal parameters */
void* hQMF = NULL;               /* qmf handle */
QMF_FDDATA_FORMAT format;        /* enum cast of "formatFlag" */
float* freqVector = NULL;
int nBands;
int procDelay;
int timeSlots;
float** dataTD_in = NULL;
float** dataTD_out = NULL;
float_complex*** dataFD_in = NULL;
float_complex*** dataFD_out = NULL;


/* ===================================================================== */
/*                              MEX Wrapper                              */
/* ===================================================================== */

void mexFunction
(
    int nlhs,             /* Number of output argments */
    mxArray *plhs[],      /* Pointers for output arguments */
    int nrhs,             /* Number of input argments */
    const mxArray *prhs[] /* Pointers for input arguments */
)
{  
    /* mex variables */
    int nDims;
    int *pDims = NULL;
     
    /* DESTROY */
    if(nrhs == 0){
        if(hQMF!=NULL){
            mexPrintf("Destroying QMF filterbank.\n");
            qmf_destroy(&hQMF);
            free(freqVector); freqVector = NULL;
            free(dataTD_in);  dataTD_in = NULL;
            free(dataFD_in);  dataFD_in = NULL;
            free(dataTD_out); dataTD_out = NULL;
            free(dataFD_out); dataFD_out = NULL;
            hQMF = NULL;
        } 
        else
            mexPrintf("QMF filterbank is already dead!\n"); 
    }
    
    /* CREATE */
    else if(nrhs==NUM_INPUT_ARGS_CREATE && (nlhs==0 || nlhs==1 || nlhs==2)){
        if(hQMF!=NULL)
            mexErrMsgIdAndTxt("MyToolbox:inputError","safmex_qmf is already initialised! First destroy it if you want to change its configuration.");
        
        /* Check input argument datatypes are as expected */ 
        checkArgDataTypes((mxArray**)prhs, (MEX_DATA_TYPES*)inputDataTypes_create, NUM_INPUT_ARGS_CREATE); 
        
        /* Copy user arguments */
        nCHin = (int)mxGetScalar(prhs[0]);
        nCHout = (int)mxGetScalar(prhs[1]);
        hopsize = (int)mxGetScalar(prhs[2]);
        blocksize = (int)mxGetScalar(prhs[3]); 
        hybridmode = (int)mxGetScalar(prhs[4]); 
        formatFlag = (int)mxGetScalar(prhs[5]); 
        fs = (float)mxGetScalar(prhs[6]); 
        switch(formatFlag){
            case 0: format = QMF_BANDS_CH_TIME; break;
            case 1: format = QMF_TIME_CH_BANDS; break;
            default:
                mexErrMsgIdAndTxt("MyToolbox:inputError","the value of the fifth argument should be 0 or 1");
        }
        
        /* Extra checks */
        if( !(hybridmode==0 || hybridmode==1) )
            mexErrMsgIdAndTxt("MyToolbox:inputError","'hybridmode' should be 0 (disabled) or 1 (enabled)");
        if( !(formatFlag==0 || formatFlag==1) )
            mexErrMsgIdAndTxt("MyToolbox:inputError","'formatFlag' should be 0 (bands x channels x time) or 1 (time x channels x bands)");
        if( !(hopsize==4 || hopsize==8 || hopsize==16 || hopsize==32 || hopsize==64 || hopsize==128) )
            mexErrMsgIdAndTxt("MyToolbox:inputError","the 'hopsize' should be 4, 8, 16, 32, 64, or 128");
        if( blocksize % hopsize != 0)
            mexErrMsgIdAndTxt("MyToolbox:inputError","'blocksize' must be a multiple of 'hopsize'");
                 
        /* Create an instance of the qmf filterbank */
        timeSlots = blocksize/hopsize;
        qmf_create(&hQMF, nCHin, nCHout, hopsize, hybridmode, format);
        nBands = qmf_getNBands(hQMF);
        procDelay = qmf_getProcDelay(hQMF);
        freqVector = malloc1d(nBands*sizeof(float));
        qmf_getCentreFreqs(hQMF, fs, nBands, freqVector);
        
        /* Allocate buffers */
        dataTD_in = (float**)malloc2d(nCHin, blocksize, sizeof(float)); 
        dataTD_out = (float**)malloc2d(nCHout, blocksize, sizeof(float));
        switch(formatFlag){
            case 0:
                dataFD_in = (float_complex***)malloc3d(nBands, nCHin, timeSlots, sizeof(float_complex));
                dataFD_out = (float_complex***)malloc3d(nBands, nCHout, timeSlots, sizeof(float_complex));
                break;
            case 1:
                dataFD_in = (float_complex***)malloc3d(timeSlots, nCHin, nBands, sizeof(float_complex));
                dataFD_out = (float_complex***)malloc3d(timeSlots, nCHout, nBands, sizeof(float_complex));
                break;
        }
        
        /* (optional) output frequency vector and processing delay */
        if(nlhs>0){
            nDims = 2;
            pDims = realloc1d(pDims, 2*sizeof(int));
            pDims[0] = nBands;
            pDims[1] = 1;
            SAFsingle2MEXdouble(freqVector, nDims, pDims, &plhs[0]); 
        } 
        if(nlhs>1){
            plhs[1] = mxCreateDoubleScalar(procDelay); 
        } 
         
        /* Mainly just for debugging... */
        mexPrintf("Creating QMF filterbank:");
        snprintf(message, MSG_STR_LENGTH, " %d input channels,", nCHin); mexPrintf(message);
        snprintf(message, MSG_STR_LENGTH, " %d output channels,", nCHout); mexPrintf(message);
        snprintf(message, MSG_STR_LENGTH, " %d hopsize,", hopsize); mexPrintf(message);
        snprintf(message, MSG_STR_LENGTH, " %d blocksize,", blocksize); mexPrintf(message);
        if(hybridmode)
            mexPrintf(" hybrid mode enabled,");
        else
            mexPrintf(" hybrid mode disabled,");
        if(formatFlag)
            mexPrintf(" format: time x channels x bands.\n");
        else
            mexPrintf(" format: bands x channels x time.\n");
    }
    
    /* TRANSFORM */
    else if(nrhs == 1 && nlhs == 1){
        if(hQMF==NULL)
            mexErrMsgIdAndTxt("MyToolbox:inputError","safmex_qmf is uninitialised!");
          
        /* Find dimensionality of input */
        mwSize nDims_mx;
        const mwSize *pDims_mx;
        nDims_mx = mxGetNumberOfDimensions(prhs[0]);
        pDims_mx = mxGetDimensions(prhs[0]); 
        
        /* FORWARD */
        if(!mxIsComplex(prhs[0])){ 
            /* Check input argument datatypes are as expected */ 
            checkArgDataTypes((mxArray**)prhs, (MEX_DATA_TYPES*)inputDataTypes_fwd, NUM_INPUT_ARGS_FWD); 
            
            /* extra checks */
            if( !(pDims_mx[0] == (mwSize)nCHin) ){
                snprintf(message, MSG_STR_LENGTH, "Was expecting %d input channels.", nCHin);
                mexErrMsgIdAndTxt("MyToolbox:inputError", message);
            }
            if( !(pDims_mx[1] == (mwSize)blocksize) ){
                snprintf(message, MSG_STR_LENGTH, "Was expecting a block size of %d samples.", blocksize);
                mexErrMsgIdAndTxt("MyToolbox:inputError", message);
            } 
            
            /* QMF analysis */
            MEXdouble2SAFsingle(prhs[0], &FLATTEN2D(dataTD_in), &nDims, &pDims);  
            qmf_analysis(hQMF, dataTD_in, blocksize, dataFD_in);
            
            /* output */
            nDims = 3;
            pDims = realloc1d(pDims, nDims*sizeof(int));
            switch(formatFlag){
                case 0: pDims[0] = nBands; pDims[1] = nCHin; pDims[2] = timeSlots; break;
                case 1: pDims[0] = timeSlots; pDims[1] = nCHin; pDims[2] = nBands; break;
            }
            SAFsingle2MEXdouble_complex(FLATTEN3D(dataFD_in), nDims, pDims, &plhs[0]);
            
            /* Check output argument datatypes are as expected */ 
            checkArgDataTypes((mxArray**)plhs, (MEX_DATA_TYPES*)outputDataTypes_fwd, NUM_OUTPUT_ARGS_FWD); 
        }
        
        /* BACKWARD */
        else if(mxIsComplex(prhs[0])){
            /* Check input argument datatypes are as expected */ 
            checkArgDataTypes((mxArray**)prhs, (MEX_DATA_TYPES*)inputDataTypes_bkwd, NUM_INPUT_ARGS_BKWD); 
            
            /* extra checks */
            if( !(pDims_mx[0] == (mwSize)nBands) && formatFlag==0 ){
                snprintf(message, MSG_STR_LENGTH, "Was expecting %d bands.", nBands);
                mexErrMsgIdAndTxt("MyToolbox:inputError", message);
            }
            if( !(pDims_mx[0] == (mwSize)timeSlots) && formatFlag==1 ){
                snprintf(message, MSG_STR_LENGTH, "Was expecting %d down-sampled time indices.", timeSlots);
                mexErrMsgIdAndTxt("MyToolbox:inputError", message);
            }
            if( !(pDims_mx[1] == (mwSize)nCHout) ){
                snprintf(message, MSG_STR_LENGTH, "Was expecting %d input channels.", nCHout);
                mexErrMsgIdAndTxt("MyToolbox:inputError", message);
            } 
            if( !(pDims_mx[2] == (mwSize)timeSlots) && formatFlag==0 ){
                snprintf(message, MSG_STR_LENGTH, "Was expecting %d down-sampled time indices.", timeSlots);
                mexErrMsgIdAndTxt("MyToolbox:inputError", message);
            } 
            if( !(pDims_mx[2] == (mwSize)nBands) && formatFlag==1 ){
                snprintf(message, MSG_STR_LENGTH, "Was expecting %d bands.", nBands);
                mexErrMsgIdAndTxt("MyToolbox:inputError", message);
            } 
            
            /* QMF synthesis */
            MEXdouble2SAFsingle_complex(prhs[0], &FLATTEN3D(dataFD_out), &nDims, &pDims); 
            qmf_synthesis(hQMF, dataFD_out, blocksize, dataTD_out);
             
            /* output */
            nDims = 2;
            pDims = realloc1d(pDims, nDims*sizeof(int));
            pDims[0] = nCHout;
            pDims[1] = blocksize; 
            SAFsingle2MEXdouble(FLATTEN2D(dataTD_out), nDims, pDims, &plhs[0]);
            
            /* Check output argument datatypes are as expected */ 
            checkArgDataTypes((mxArray**)plhs, (MEX_DATA_TYPES*)outputDataTypes_bkwd, NUM_OUTPUT_ARGS_BKWD);  
        }
        else
            mexErrMsgIdAndTxt("MyToolbox:inputError","Unrecognised input/output configuration, refer to help instructions.");
    }
    
    /* ERROR */
    else 
        mexErrMsgIdAndTxt("MyToolbox:inputError","Unrecognised input/output configuration, refer to help instructions.");
}
