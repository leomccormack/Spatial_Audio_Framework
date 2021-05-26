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
 * @file safmex_latticeDecorrelator.c
 * @brief MEX wrapper for latticeDecorrelator (see the .m file of the same name
 *        for documentation)
 * @author Leo McCormack
 * @date 06.08.2020
 */
 
#include "safmex.h"

/* ===================================================================== */
/*                                Config                                 */
/* ===================================================================== */

#define NUM_INPUT_ARGS_CREATE  ( 8 )
const MEX_DATA_TYPES inputDataTypes_create[NUM_INPUT_ARGS_CREATE] = { 
    SM_DOUBLE_REAL,
    SM_INT32,  
    SM_INT32,
    SM_DOUBLE_REAL_1D,  
    SM_DOUBLE_REAL_1D,  
    SM_INT32,
    SM_DOUBLE_REAL_1D, 
    SM_INT32 
}; 
#define NUM_INPUT_ARGS_APPLY ( 1 )
#define NUM_OUTPUT_ARGS_APPLY ( 1 ) 
const MEX_DATA_TYPES inputDataTypes_apply[NUM_INPUT_ARGS_APPLY] = {
    SM_DOUBLE_COMPLEX_3D      
};  
const MEX_DATA_TYPES outputDataTypes_apply[NUM_OUTPUT_ARGS_APPLY] = {
    SM_DOUBLE_COMPLEX_3D 
};  


/* ===================================================================== */
/*                                 Vars                                  */
/* ===================================================================== */

/* user arguments */
float fs;
int hopsize;
int nCH;                   /* Number of channels */
int* orders = NULL;        /* Lattice all-pass filter orders (2,3,4,6,8,10,12,14,
                              15,16 18, or 20) per band grouping; nCutoffs x 1 */
float* freqCutoffs = NULL; /* Frequency cut-offs defining the band groupings;
                              nCutoffs x 1 */
int maxDelay;              /* Maximum delay */
float* freqVector = NULL;  /* Frequency vector; nBands x 1 */ 
int nTimeSlots;            /* Number of TF frames to process at a time */

/* internal parameters */
void* hDecor = NULL;       /* decor handle */  
int nBands;
int nCutoffs; 
float_complex*** dataFD_in = NULL;
float_complex*** dataFD_out = NULL;

const int nSupportedOrders = 12;
const int supportedOrders[12] = {2,3,4,6,8,10,12,14,15,16,18,20};


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
    int i,j,nDims;
    int *pDims = NULL;
   
    /* DESTROY */
    if(nrhs == 0){
        if(hDecor!=NULL){
            mexPrintf("Destroying latticeDecorrelator.\n");
            latticeDecorrelator_destroy(&hDecor);
            free(orders);      orders = NULL;
            free(freqCutoffs); freqCutoffs = NULL;
            free(freqVector);  freqVector = NULL;
            free(dataFD_in);   dataFD_in = NULL;
            free(dataFD_out);  dataFD_out = NULL;
            hDecor = NULL;
        } 
        else
            mexPrintf("latticeDecorrelator is already dead!\n"); 
    }
    
    /* CREATE */
    else if(nrhs==NUM_INPUT_ARGS_CREATE){
        if(hDecor!=NULL)
            mexErrMsgIdAndTxt("MyToolbox:inputError","safmex_latticeDecorrelator is already initialised! First destroy it if you want to change its configuration.");
        
        /* Check input argument datatypes are as expected */ 
        checkArgDataTypes((mxArray**)prhs, (MEX_DATA_TYPES*)inputDataTypes_create, NUM_INPUT_ARGS_CREATE); 
        
        /* Copy user arguments */  
        fs = (float)mxGetScalar(prhs[0]); 
        hopsize = (int)mxGetScalar(prhs[1]); 
        nCH = (int)mxGetScalar(prhs[2]); 
        MEXdouble2SAFsingle_int(prhs[3], &orders, &nDims, &pDims);
         
        nCutoffs = pDims[0];
        if( nCutoffs<=1 )
            mexErrMsgIdAndTxt("MyToolbox:inputError","freqCutoffs vector must be longer than 1 element."); 
        for(i=0; i<nCutoffs; i++){
            int supportedOrder = 0;
            for(j=0; j<nSupportedOrders; j++)
                if(orders[i] == supportedOrders[j])
                    supportedOrder = 1;
            
            if(!supportedOrder)
                mexErrMsgIdAndTxt("MyToolbox:inputError","Supported 'orders' are: 2,3,4,6,8,10,12,14,15,16,18,20."); 
        }    
        
        MEXdouble2SAFsingle(prhs[4], &freqCutoffs, &nDims, &pDims);
        if(pDims[0] != nCutoffs )
            mexErrMsgIdAndTxt("MyToolbox:inputError","freqCutoffs vector must be the same length as orders vector.");
                
        maxDelay = (int)mxGetScalar(prhs[5]);

        MEXdouble2SAFsingle(prhs[6], &freqVector, &nDims, &pDims);
        nBands = pDims[0];  
        nTimeSlots = (int)mxGetScalar(prhs[7]);  
   
        /* Create an instance of latticeDecorrelator */ 
        latticeDecorrelator_create(&hDecor, fs, hopsize, freqVector, nBands, nCH, orders, freqCutoffs, nCutoffs, maxDelay, 0, 0.75f);

        /* Allocate buffers */ 
        dataFD_in = (float_complex***)malloc3d(nBands, nCH, nTimeSlots, sizeof(float_complex));
        dataFD_out = (float_complex***)malloc3d(nBands, nCH, nTimeSlots, sizeof(float_complex)); 
        
        /* Mainly just for debugging... */
        mexPrintf("Creating latticeDecorrelator:");
        snprintf(message, MSG_STR_LENGTH, " %d channels,", nCH); mexPrintf(message);
        mexPrintf(" filter orders = ["); 
        for(i=0; i<nCutoffs; i++){
            snprintf(message, MSG_STR_LENGTH, " %d ", orders[i]); 
            mexPrintf(message);
        }  
        mexPrintf("], ");
        mexPrintf(" cut-offs = ["); 
        for(i=0; i<nCutoffs; i++){
            snprintf(message, MSG_STR_LENGTH, " %.2f ", freqCutoffs[i]); 
            mexPrintf(message);
        }  
        mexPrintf("], "); 
        snprintf(message, MSG_STR_LENGTH, " %d nBands,", nBands); mexPrintf(message);
        snprintf(message, MSG_STR_LENGTH, " %d timeslots\n", nTimeSlots); mexPrintf(message);
    }
    
    /* APPLY */
    else if(nrhs == 1 && nlhs == 1){
        if(hDecor==NULL)
            mexErrMsgIdAndTxt("MyToolbox:inputError","safmex_latticeDecorrelator is uninitialised!");
          
        /* Find dimensionality of input */
        mwSize nDims_mx;
        const mwSize *pDims_mx;
        nDims_mx = mxGetNumberOfDimensions(prhs[0]);
        pDims_mx = mxGetDimensions(prhs[0]); 
 
        /* Check input argument datatypes are as expected */ 
        checkArgDataTypes((mxArray**)prhs, (MEX_DATA_TYPES*)inputDataTypes_apply, NUM_INPUT_ARGS_APPLY); 

        /* extra checks */ 
        if( !(pDims_mx[0] == (mwSize)nBands) ){
            snprintf(message, MSG_STR_LENGTH, "Was expecting %d bands.", nBands);
            mexErrMsgIdAndTxt("MyToolbox:inputError", message);
        } 
        if( !(pDims_mx[1] == (mwSize)nCH) ){
            snprintf(message, MSG_STR_LENGTH, "Was expecting %d channels.", nCH);
            mexErrMsgIdAndTxt("MyToolbox:inputError", message);
        } 
        if( !(pDims_mx[2] == (mwSize)nTimeSlots) ){
            snprintf(message, MSG_STR_LENGTH, "Was expecting %d down-sampled time indices.", nTimeSlots);
            mexErrMsgIdAndTxt("MyToolbox:inputError", message);
        }  

        /* Apply decorrelator */
        MEXdouble2SAFsingle_complex(prhs[0], &FLATTEN3D(dataFD_in), &nDims, &pDims);  
        latticeDecorrelator_apply(hDecor, dataFD_in, nTimeSlots, dataFD_out);

        /* output */
        nDims = 3;
        pDims = realloc1d(pDims, nDims*sizeof(int));
        pDims[0] = nBands; pDims[1] = nCH; pDims[2] = nTimeSlots; 
        SAFsingle2MEXdouble_complex(FLATTEN3D(dataFD_out), nDims, pDims, &plhs[0]);

        /* Check output argument datatypes are as expected */ 
        checkArgDataTypes((mxArray**)plhs, (MEX_DATA_TYPES*)outputDataTypes_apply, NUM_OUTPUT_ARGS_APPLY); 
    }
    
    /* ERROR */
    else 
        mexErrMsgIdAndTxt("MyToolbox:inputError","Unrecognised input/output configuration, refer to help instructions.");
}
