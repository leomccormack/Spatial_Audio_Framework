/*
 * This file is part of the saf_tracker module.
 * Copyright (c) 2020 - Leo McCormack
 *
 * The saf_tracker module is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at your
 * option) any later version.
 *
 * The saf_tracker module is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 * or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
 * more details.
 *
 * See <http://www.gnu.org/licenses/> for a copy of the GNU General Public
 * License.
 */
 
#include "safmex.h"

/* ===================================================================== */
/*                                Config                                 */
/* ===================================================================== */

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
tracker3d_config tpars;

/* internal parameters */
void* hT3d = NULL;                /* tracker3d handle */ 
int ifield, nfields;
 
/* ===================================================================== */
/*                              MEX Wrapper                              */
/* ===================================================================== */

void mexFunction
(
    int nlhs,             /* Number of input argments */
    mxArray *plhs[],      /* Pointers for input arguments */
    int nrhs,             /* Number of output argments */
    const mxArray *prhs[] /* Pointers for output arguments */
)
{  
    /* mex variables */
    int nDims;
    int *pDims = NULL;
     
    /* DESTROY */
    if(nrhs == 0){
        if(hT3d!=NULL){
            mexPrintf("Destroying tracker3d.\n");
            tracker3d_destroy(&hT3d); 
            hT3d = NULL;
        } 
        else
            mexPrintf("tracker3d is already dead!\n"); 
    }
    
    /* CREATE */
    else if(nrhs==1){
        if(hT3d!=NULL)
            mexErrMsgIdAndTxt("MyToolbox:inputError","tracker3d is already initialised! First destroy it if you want to change its configuration.");
         
          /* Check if prhs[1] is a struct */
          if(!mxIsStruct(prhs[0]))
              mexErrMsgIdAndTxt("MyToolbox:inputError","Input must be a struct"); 
          
          /* get the values from the struct */
          const mxArray  *mxTmp; 
          double         *tmp;

          /* tpars.Np */
          mxTmp = mxGetField(prhs[0],0,"Np"); 
          if(mxTmp==NULL)
              mexErrMsgIdAndTxt("MyToolbox:inputError","'Np' is not defined"); 
          if (mxIsComplex(mxTmp) || mxGetNumberOfElements(mxTmp)!=1 || (int)mxGetScalar(mxTmp)<1 || (int)mxGetScalar(mxTmp)>100)
              mexErrMsgIdAndTxt("MyToolbox:inputError","'Np' must be an integer between 1 and 100");
          tpars.Np = (int)mxGetScalar(mxTmp);
          
          /* tpars.maxNactiveTargets */
          mxTmp = mxGetField(prhs[0],0,"maxNactiveTargets"); 
          if(mxTmp==NULL)
              mexErrMsgIdAndTxt("MyToolbox:inputError","'maxNactiveTargets' is not defined"); 
          if (mxIsComplex(mxTmp) || mxGetNumberOfElements(mxTmp)!=1 || (int)mxGetScalar(mxTmp)<1 || (int)mxGetScalar(mxTmp)>100)
              mexErrMsgIdAndTxt("MyToolbox:inputError","'maxNactiveTargets' must be an integer between 1 and 100");
          tpars.maxNactiveTargets = (int)mxGetScalar(mxTmp);
          
          /* tpars.noiseLikelihood */
          mxTmp = mxGetField(prhs[0],0,"noiseLikelihood"); 
          if(mxTmp==NULL)
              mexErrMsgIdAndTxt("MyToolbox:inputError","'noiseLikelihood' is not defined"); 
          if (mxIsComplex(mxTmp) || mxGetNumberOfElements(mxTmp)!=1 || mxGetScalar(mxTmp)<0 || mxGetScalar(mxTmp)>1)
              mexErrMsgIdAndTxt("MyToolbox:inputError","'noiseLikelihood' must be a scalar between 0 and 1");
          tpars.noiseLikelihood = mxGetScalar(mxTmp);
          
          /* tpars.measNoiseSD */
          mxTmp = mxGetField(prhs[0],0,"measNoiseSD"); 
          if(mxTmp==NULL)
              mexErrMsgIdAndTxt("MyToolbox:inputError","'measNoiseSD' is not defined"); 
          if (mxIsComplex(mxTmp) || mxGetNumberOfElements(mxTmp)!=1 || mxGetScalar(mxTmp)<0)
              mexErrMsgIdAndTxt("MyToolbox:inputError","'measNoiseSD' must be a scalar");
          tpars.measNoiseSD = mxGetScalar(mxTmp);
          
          /* tpars.noiseSpecDen */
          mxTmp = mxGetField(prhs[0],0,"noiseSpecDen"); 
          if(mxTmp==NULL)
              mexErrMsgIdAndTxt("MyToolbox:inputError","'noiseSpecDen' is not defined"); 
          if (mxIsComplex(mxTmp) || mxGetNumberOfElements(mxTmp)!=1 || mxGetScalar(mxTmp)<0)
              mexErrMsgIdAndTxt("MyToolbox:inputError","'noiseSpecDen' must be a scalar");
          tpars.noiseSpecDen = mxGetScalar(mxTmp);
          
          /* tpars.ALLOW_MULTI_DEATH */
          mxTmp = mxGetField(prhs[0],0,"ALLOW_MULTI_DEATH"); 
          if(mxTmp==NULL)
              mexErrMsgIdAndTxt("MyToolbox:inputError","'ALLOW_MULTI_DEATH' is not defined"); 
          if (mxIsComplex(mxTmp) || mxGetNumberOfElements(mxTmp)!=1 || (int)mxGetScalar(mxTmp)<0 || (int)mxGetScalar(mxTmp)>1)
              mexErrMsgIdAndTxt("MyToolbox:inputError","'ALLOW_MULTI_DEATH' must be 0 or 1");
          tpars.ALLOW_MULTI_DEATH = (int)mxGetScalar(mxTmp);
                  
          /* tpars.init_birth */
          mxTmp = mxGetField(prhs[0],0,"init_birth"); 
          if(mxTmp==NULL)
              mexErrMsgIdAndTxt("MyToolbox:inputError","'init_birth' is not defined"); 
          if (mxIsComplex(mxTmp) || mxGetNumberOfElements(mxTmp)!=1 || mxGetScalar(mxTmp)<0 || mxGetScalar(mxTmp)>1)
              mexErrMsgIdAndTxt("MyToolbox:inputError","'init_birth' must be a scalar between 0 and 1");
          tpars.init_birth = mxGetScalar(mxTmp);
          
          /* tpars.alpha_death */
          mxTmp = mxGetField(prhs[0],0,"alpha_death"); 
          if(mxTmp==NULL)
              mexErrMsgIdAndTxt("MyToolbox:inputError","'alpha_death' is not defined"); 
          if (mxIsComplex(mxTmp) || mxGetNumberOfElements(mxTmp)!=1 || mxGetScalar(mxTmp)<0 )
              mexErrMsgIdAndTxt("MyToolbox:inputError","'alpha_death' must be a scalar");
          tpars.alpha_death = mxGetScalar(mxTmp);
          
          /* tpars.beta_death */
          mxTmp = mxGetField(prhs[0],0,"beta_death"); 
          if(mxTmp==NULL)
              mexErrMsgIdAndTxt("MyToolbox:inputError","'beta_death' is not defined"); 
          if (mxIsComplex(mxTmp) || mxGetNumberOfElements(mxTmp)!=1 || mxGetScalar(mxTmp)<0 )
              mexErrMsgIdAndTxt("MyToolbox:inputError","'beta_death' must be a scalar");
          tpars.beta_death = mxGetScalar(mxTmp);
          
          /* tpars.dt */
          mxTmp = mxGetField(prhs[0],0,"dt"); 
          if(mxTmp==NULL)
              mexErrMsgIdAndTxt("MyToolbox:inputError","'dt' is not defined"); 
          if (mxIsComplex(mxTmp) || mxGetNumberOfElements(mxTmp)!=1 || mxGetScalar(mxTmp)<0 )
              mexErrMsgIdAndTxt("MyToolbox:inputError","'dt' must be a scalar");
          tpars.dt = mxGetScalar(mxTmp);
      
          /* tpars.W_avg_coeff */
          mxTmp = mxGetField(prhs[0],0,"W_avg_coeff"); 
          if(mxTmp==NULL)
              mexErrMsgIdAndTxt("MyToolbox:inputError","'W_avg_coeff' is not defined"); 
          if (mxIsComplex(mxTmp) || mxGetNumberOfElements(mxTmp)!=1 || mxGetScalar(mxTmp)<0 || mxGetScalar(mxTmp)>1)
              mexErrMsgIdAndTxt("MyToolbox:inputError","'W_avg_coeff' must be a scalar between 0 and 1");
          tpars.W_avg_coeff = mxGetScalar(mxTmp);
          
          /* tpars.FORCE_KILL_TARGETS */
          mxTmp = mxGetField(prhs[0],0,"FORCE_KILL_TARGETS"); 
          if(mxTmp==NULL)
              mexErrMsgIdAndTxt("MyToolbox:inputError","'FORCE_KILL_TARGETS' is not defined"); 
          if (mxIsComplex(mxTmp) || mxGetNumberOfElements(mxTmp)!=1 || (int)mxGetScalar(mxTmp)<0 || (int)mxGetScalar(mxTmp)>1)
              mexErrMsgIdAndTxt("MyToolbox:inputError","'FORCE_KILL_TARGETS' must be 0 or 1");
          tpars.FORCE_KILL_TARGETS = (int)mxGetScalar(mxTmp);
          
          /* tpars.forceKillDistance */
          mxTmp = mxGetField(prhs[0],0,"forceKillDistance"); 
          if(mxTmp==NULL)
              mexErrMsgIdAndTxt("MyToolbox:inputError","'forceKillDistance' is not defined"); 
          if (mxIsComplex(mxTmp) || mxGetNumberOfElements(mxTmp)!=1 || mxGetScalar(mxTmp)<0 )
              mexErrMsgIdAndTxt("MyToolbox:inputError","'forceKillDistance' must be a scalar");
          tpars.forceKillDistance = mxGetScalar(mxTmp);
          
          /* tpars.M0 */
          mxTmp = mxGetField(prhs[0],0,"M0"); 
          if(mxTmp==NULL)
              mexErrMsgIdAndTxt("MyToolbox:inputError","'M0' is not defined"); 
          if (mxIsComplex(mxTmp) || mxGetNumberOfElements(mxTmp)!=6)
              mexErrMsgIdAndTxt("MyToolbox:inputError","'M0' must be a 6-element vector");
           
          
          /* Find dimensionality of input */
          //nDims_mx = mxGetNumberOfDimensions(tmp);
          //pDims_mx = mxGetDimensions(tmp);
          tmp = mxGetData(mxTmp);
          
//      
//     float M0[6];              /**< [0,1,2] Position of sound source PRIORs
//                                *   (x,y,z), [3,4,5] Mean velocity PRIORs (x,y,z)
//                                */
//     float P0[6][6];           /**< Diagonal matrix, [0,1,2] Variance PRIORs of
//                                *   estimates along the x,y,z axes; [3,4,5]
//                                *   Velocity PRIORs of estimates along the x,y,z
//                                *   axes */
 
          mxTmp = mxGetField(prhs[0],0,"cd"); 
          if(mxTmp==NULL)
              mexErrMsgIdAndTxt("MyToolbox:inputError","'cd' is not defined"); 
          if (mxIsComplex(mxTmp) || mxGetNumberOfElements(mxTmp)!=1 || mxGetScalar(mxTmp)<0 )
              mexErrMsgIdAndTxt("MyToolbox:inputError","'cd' must be a scalar");
          tpars.cd = mxGetScalar(mxTmp);
    
           
          
//           tmp   = mxGetPr(mxTmp);
 
//           int iPar1;
//           for (iPar1=0; iPar1<4; iPar1++)
//           {
//             mexPrintf("par1[%i] = %.2f \n",iPar1,tmp[iPar1]);
//           }
//         
        
//         /* Check input argument datatypes are as expected */ 
//         checkArgDataTypes((mxArray**)prhs, (MEX_DATA_TYPES*)inputDataTypes_create, NUM_INPUT_ARGS_CREATE); 
//         
//         /* Copy user arguments */
//         nCHin = (int)mxGetScalar(prhs[0]);
//         nCHout = (int)mxGetScalar(prhs[1]);
//         hopsize = (int)mxGetScalar(prhs[2]);
//         blocksize = (int)mxGetScalar(prhs[3]); 
//         hybridmode = (int)mxGetScalar(prhs[4]); 
//         formatFlag = (int)mxGetScalar(prhs[5]); 
//         fs = (float)mxGetScalar(prhs[6]); 
//         switch(formatFlag){
//             case 0: format = AFSTFT_BANDS_CH_TIME; break;
//             case 1: format = AFSTFT_TIME_CH_BANDS; break;
//             default:
//                 mexErrMsgIdAndTxt("MyToolbox:inputError","the value of the fifth argument should be 0 or 1");
//         }
//         
//         /* Extra checks */
//         if( !(hybridmode==0 || hybridmode==1) )
//             mexErrMsgIdAndTxt("MyToolbox:inputError","'hybridmode' should be 0 (disabled) or 1 (enabled)");
//         if( !(formatFlag==0 || formatFlag==1) )
//             mexErrMsgIdAndTxt("MyToolbox:inputError","'formatFlag' should be 0 (bands x channels x time) or 1 (time x channels x bands)");
//         if( !(hopsize==4 || hopsize==8 || hopsize==16 || hopsize==32 || hopsize==64 || hopsize==128) )
//             mexErrMsgIdAndTxt("MyToolbox:inputError","the 'hopsize' should be 4, 8, 16, 32, 64, or 128");
//         if( blocksize % hopsize != 0)
//             mexErrMsgIdAndTxt("MyToolbox:inputError","'blocksize' must be a multiple of 'hopsize'");
//                  
//         /* Create an instance of the afSTFT filterbank */
//         timeSlots = blocksize/hopsize;
//         afSTFT_create(&hSTFT, nCHin, nCHout, hopsize, 0, hybridmode, format);
//         nBands = afSTFT_getNBands(hSTFT);
//         procDelay = afSTFT_getProcDelay(hSTFT);
//         freqVector = malloc1d(nBands*sizeof(float));
//         afSTFT_getCentreFreqs(hSTFT, fs, nBands, freqVector);
//         
//         /* Allocate buffers */
//         dataTD_in = (float**)malloc2d(nCHin, blocksize, sizeof(float)); 
//         dataTD_out = (float**)malloc2d(nCHout, blocksize, sizeof(float));
//         switch(formatFlag){
//             case 0:
//                 dataFD_in = (float_complex***)malloc3d(nBands, nCHin, timeSlots, sizeof(float_complex));
//                 dataFD_out = (float_complex***)malloc3d(nBands, nCHout, timeSlots, sizeof(float_complex));
//                 break;
//             case 1:
//                 dataFD_in = (float_complex***)malloc3d(timeSlots, nCHin, nBands, sizeof(float_complex));
//                 dataFD_out = (float_complex***)malloc3d(timeSlots, nCHout, nBands, sizeof(float_complex));
//                 break;
//         }
//         
//         /* (optional) output frequency vector and processing delay */
//         if(nlhs>0){
//             nDims = 2;
//             pDims = realloc1d(pDims, 2*sizeof(int));
//             pDims[0] = nBands;
//             pDims[1] = 1;
//             SAFsingle2MEXdouble(freqVector, nDims, pDims, &plhs[0]); 
//         } 
//         if(nlhs>1){
//             plhs[1] = mxCreateDoubleScalar(procDelay); 
//         } 
         
//         /* Mainly just for debugging... */
//         mexPrintf("Creating afSTFT filterbank:");
//         snprintf(message, MSG_STR_LENGTH, " %d input channels,", nCHin); mexPrintf(message);
//         snprintf(message, MSG_STR_LENGTH, " %d output channels,", nCHout); mexPrintf(message);
//         snprintf(message, MSG_STR_LENGTH, " %d hopsize,", hopsize); mexPrintf(message);
//         snprintf(message, MSG_STR_LENGTH, " %d blocksize,", blocksize); mexPrintf(message);
//         if(hybridmode)
//             mexPrintf(" hybrid mode enabled,");
//         else
//             mexPrintf(" hybrid mode disabled,");
//         if(formatFlag)
//             mexPrintf(" format: time x channels x bands.\n");
//         else
//             mexPrintf(" format: bands x channels x time.\n");
    }
    
    /* TRANSFORM */
    else if(nrhs == 1 && nlhs == 1){
//         if(hSTFT==NULL)
//             mexErrMsgIdAndTxt("MyToolbox:inputError","safmex_afSTFT is uninitialised!");
//           
//         /* Find dimensionality of input */
//         mwSize nDims_mx;
//         const mwSize *pDims_mx;
//         nDims_mx = mxGetNumberOfDimensions(prhs[0]);
//         pDims_mx = mxGetDimensions(prhs[0]); 
//         
//         /* FORWARD */
//         if(!mxIsComplex(prhs[0])){ 
//             /* Check input argument datatypes are as expected */ 
//             checkArgDataTypes((mxArray**)prhs, (MEX_DATA_TYPES*)inputDataTypes_fwd, NUM_INPUT_ARGS_FWD); 
//             
//             /* extra checks */
//             if( !(pDims_mx[0] == (mwSize)nCHin) ){
//                 snprintf(message, MSG_STR_LENGTH, "Was expecting %d input channels.", nCHin);
//                 mexErrMsgIdAndTxt("MyToolbox:inputError", message);
//             }
//             if( !(pDims_mx[1] == (mwSize)blocksize) ){
//                 snprintf(message, MSG_STR_LENGTH, "Was expecting a block size of %d samples.", blocksize);
//                 mexErrMsgIdAndTxt("MyToolbox:inputError", message);
//             } 
//             
//             /* afSTFT forward */
//             MEXdouble2SAFsingle(prhs[0], &FLATTEN2D(dataTD_in), &nDims, &pDims);  
//             afSTFT_forward(hSTFT, dataTD_in, blocksize, dataFD_in);
//             
//             /* output */
//             nDims = 3;
//             pDims = realloc1d(pDims, nDims*sizeof(int));
//             switch(formatFlag){
//                 case 0: pDims[0] = nBands; pDims[1] = nCHin; pDims[2] = timeSlots; break;
//                 case 1: pDims[0] = timeSlots; pDims[1] = nCHin; pDims[2] = nBands; break;
//             }
//             SAFsingle2MEXdouble_complex(FLATTEN3D(dataFD_in), nDims, pDims, &plhs[0]);
//             
//             /* Check output argument datatypes are as expected */ 
//             checkArgDataTypes((mxArray**)plhs, (MEX_DATA_TYPES*)outputDataTypes_fwd, NUM_OUTPUT_ARGS_FWD); 
//        }
        
        /* BACKWARD */
//        else if(mxIsComplex(prhs[0])){
//             /* Check input argument datatypes are as expected */ 
//             checkArgDataTypes((mxArray**)prhs, (MEX_DATA_TYPES*)inputDataTypes_bkwd, NUM_INPUT_ARGS_BKWD); 
//             
//             /* extra checks */
//             if( !(pDims_mx[0] == (mwSize)nBands) && formatFlag==0 ){
//                 snprintf(message, MSG_STR_LENGTH, "Was expecting %d bands.", nBands);
//                 mexErrMsgIdAndTxt("MyToolbox:inputError", message);
//             }
//             if( !(pDims_mx[0] == (mwSize)timeSlots) && formatFlag==1 ){
//                 snprintf(message, MSG_STR_LENGTH, "Was expecting %d down-sampled time indices.", timeSlots);
//                 mexErrMsgIdAndTxt("MyToolbox:inputError", message);
//             }
//             if( !(pDims_mx[1] == (mwSize)nCHout) ){
//                 snprintf(message, MSG_STR_LENGTH, "Was expecting %d input channels.", nCHout);
//                 mexErrMsgIdAndTxt("MyToolbox:inputError", message);
//             } 
//             if( !(pDims_mx[2] == (mwSize)timeSlots) && formatFlag==0 ){
//                 snprintf(message, MSG_STR_LENGTH, "Was expecting %d down-sampled time indices.", timeSlots);
//                 mexErrMsgIdAndTxt("MyToolbox:inputError", message);
//             } 
//             if( !(pDims_mx[2] == (mwSize)nBands) && formatFlag==1 ){
//                 snprintf(message, MSG_STR_LENGTH, "Was expecting %d bands.", nBands);
//                 mexErrMsgIdAndTxt("MyToolbox:inputError", message);
//             } 
//             
//             /* afSTFT inverse */
//             MEXdouble2SAFsingle_complex(prhs[0], &FLATTEN3D(dataFD_out), &nDims, &pDims); 
//             afSTFT_backward(hSTFT, dataFD_out, blocksize, dataTD_out);
//              
//             /* output */
//             nDims = 2;
//             pDims = realloc1d(pDims, nDims*sizeof(int));
//             pDims[0] = nCHout;
//             pDims[1] = blocksize; 
//             SAFsingle2MEXdouble(FLATTEN2D(dataTD_out), nDims, pDims, &plhs[0]);
//             
//             /* Check output argument datatypes are as expected */ 
//             checkArgDataTypes((mxArray**)plhs, (MEX_DATA_TYPES*)outputDataTypes_bkwd, NUM_OUTPUT_ARGS_BKWD);  
//        }
//        else
//            mexErrMsgIdAndTxt("MyToolbox:inputError","Unrecognised input/output configuration, refer to help instructions.");
    }
    
    /* ERROR */
    else 
        mexErrMsgIdAndTxt("MyToolbox:inputError","Unrecognised input/output configuration, refer to help instructions.");
}