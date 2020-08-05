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
  
#include "mex.h"
#include "saf.h"

#define MSG_STR_LENGTH ( 256 )
char message[MSG_STR_LENGTH];
/* snprintf(message, MSG_STR_LENGTH, "mem required: %d", numElements);
   mexPrintf(message); */

typedef enum _MEX_DATA_TYPES{
    INT32 = 0, 
    
    /* Double-precision */
    DOUBLE_REAL, 
    DOUBLE_COMPLEX, 
    DOUBLE_REAL_1D, 
    DOUBLE_COMPLEX_1D, 
    DOUBLE_REAL_1D_OR_2D, 
    DOUBLE_COMPLEX_1D_OR_2D,
    DOUBLE_REAL_2D, 
    DOUBLE_COMPLEX_2D,
    DOUBLE_REAL_3D, 
    DOUBLE_COMPLEX_3D
            
}MEX_DATA_TYPES;

void checkNumInOutArgs(int nInputs, int nOutputs, int nInputs_expected, int nOutputs_expected)
{ 
    if(nInputs != nInputs_expected){
        snprintf(message, MSG_STR_LENGTH, "Number of inputs expected: %d", nInputs_expected);
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs", message); 
    }
    if(nOutputs != nOutputs_expected){
        snprintf(message, MSG_STR_LENGTH, "Number of outputs expected: %d", nOutputs_expected);
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nlhs", message);
    }
}

void checkArgDataTypes(mxArray** hData, MEX_DATA_TYPES* dataTypes, int nArgs)
{
    int i,j;
    mwSize nDims, true_nDims;
    const mwSize* pDims; 
    
    for(i=0; i<nArgs; i++){
        /* Find number of dimensions */
        nDims = mxGetNumberOfDimensions(hData[i]);
        pDims = mxGetDimensions(hData[i]); 
        true_nDims = 0; /* (note: true_nDims==0 infers argument is a scalar) */
        for(j=0; j<nDims; j++) 
            true_nDims = pDims[j]!=1 ? true_nDims+1 : true_nDims; 
         
        /* Check data types */
        switch(dataTypes[i]){
            case INT32:
                if (mxIsComplex(hData[i]) || mxGetNumberOfElements(hData[i])!=1){
                    snprintf(message, MSG_STR_LENGTH, "The following input argument must be an integer scalar: %d", i+1);
                    mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notScalar", message);
                }
                break;
            case DOUBLE_REAL: 
                if (!mxIsDouble(hData[i]) || mxIsComplex(hData[i]) || true_nDims!=0){
                    snprintf(message, MSG_STR_LENGTH, "The following input argument must be a real-valued double-precision scalar: %d", i+1);
                    mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notScalar", message);
                }
                break;
            case DOUBLE_COMPLEX:
                if (!mxIsDouble(hData[i]) || !mxIsComplex(hData[i]) || true_nDims!=0){
                    snprintf(message, MSG_STR_LENGTH, "The following input argument must be a complex-valued double-precision scalar: %d", i+1);
                    mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notScalar", message);
                }
                break;
            case DOUBLE_REAL_1D: 
                if (!mxIsDouble(hData[i]) || mxIsComplex(hData[i]) || true_nDims!=1){
                    snprintf(message, MSG_STR_LENGTH, "The following input argument must be a real-valued double-precision 1-D vector: %d", i+1);
                    mexErrMsgIdAndTxt("MyToolbox:inputError", message);
                }
                break;
            case DOUBLE_COMPLEX_1D:
                if (!mxIsDouble(hData[i]) || !mxIsComplex(hData[i]) || true_nDims!=1){
                    snprintf(message, MSG_STR_LENGTH, "The following input argument must be a complex-valued double-precision 1-D vector: %d", i+1);
                    mexErrMsgIdAndTxt("MyToolbox:inputError", message);
                }
                break;
            case DOUBLE_REAL_1D_OR_2D: 
                if (!mxIsDouble(hData[i]) || mxIsComplex(hData[i]) || true_nDims==0 || true_nDims>2){
                    snprintf(message, MSG_STR_LENGTH, "The following input argument must be a real-valued double-precision 1-D vector or 2-D matrix: %d", i+1);
                    mexErrMsgIdAndTxt("MyToolbox:inputError", message);
                }
                break;
            case DOUBLE_COMPLEX_1D_OR_2D:
                if (!mxIsDouble(hData[i]) || !mxIsComplex(hData[i]) || true_nDims==0 || true_nDims>2){
                    snprintf(message, MSG_STR_LENGTH, "The following input argument must be a complex-valued double-precision 1-D vector or 2-D matrix: %d", i+1);
                    mexErrMsgIdAndTxt("MyToolbox:inputError", message);
                }
                break;
            case DOUBLE_REAL_2D:
                if( !mxIsDouble(hData[i]) || mxIsComplex(hData[i]) || true_nDims!=2) {
                    snprintf(message, MSG_STR_LENGTH, "The following input argument must be a real-valued double-precision 2-D matrix: %d", i+1);
                    mexErrMsgIdAndTxt("MyToolbox:inputError", message); 
                }
                break;
            case DOUBLE_COMPLEX_2D:
                if( !mxIsDouble(hData[i]) || !mxIsComplex(hData[i]) || true_nDims!=2) {
                    snprintf(message, MSG_STR_LENGTH, "The following input argument must be a complex-valued double-precision 2-D matrix: %d", i+1);
                    mexErrMsgIdAndTxt("MyToolbox:inputError", message); 
                }
                break;
            case DOUBLE_REAL_3D:
                if( !mxIsDouble(hData[i]) || mxIsComplex(hData[i]) || nDims!=3) {
                    snprintf(message, MSG_STR_LENGTH, "The following input argument must be a real-valued double-precision 3-D matrix: %d", i+1);
                    mexErrMsgIdAndTxt("MyToolbox:inputError", message); 
                }
                break;
            case DOUBLE_COMPLEX_3D:
                if( !mxIsDouble(hData[i]) || !mxIsComplex(hData[i]) || nDims!=3) {
                    snprintf(message, MSG_STR_LENGTH, "The following input argument must be a complex-valued double-precision 3-D matrix: %d", i+1);
                    mexErrMsgIdAndTxt("MyToolbox:inputError", message); 
                }
                break;
        }
    } 
}

void MEXdouble2SAFsingle(const mxArray* in, float** out, int* nDims, int** pDims)
{ 
    int i, j, numElements;
    double *inMatrix;
    mwSize nDims_mx;
    const mwSize *pDims_mx;
    
    /* Find dimensionality of input */
    nDims_mx = mxGetNumberOfDimensions(in);
    pDims_mx = mxGetDimensions(in);
    
    /* convert mwSize->int */
    (*nDims) = (int)nDims_mx;
    (*pDims) = malloc1d(*nDims*sizeof(int));
    for(i=0; i<(*nDims); i++)
        (*pDims)[i] = (int)pDims_mx[i];
    
    /* Find number of elements */
    numElements = 1;
    for(i=0; i<(*nDims); i++)
        numElements *= (*pDims)[i];
     
    /* Convert input mex array to saf array */
#if MX_HAS_INTERLEAVED_COMPLEX 
    assert(0); /* need to implement switch case for number of dims, and interleave the reading of inMatrix somehow... */
#endif
    inMatrix = mxGetData(in);
    if((*out)==NULL)
        (*out) = malloc1d(numElements*sizeof(float)); 
    
    /* column-major -> row-major */
    switch(*nDims){
        case 0: /* scalar */ break;
        case 1: assert(0); break;
        case 2: 
            for(i=0; i<(*pDims)[0]; i++)
                for(j=0; j<(*pDims)[1]; j++)
                    (*out)[i*(*pDims)[1]+j] = (float)inMatrix[j*(*pDims)[0]+i];
            break;
    }
}

void MEXdouble2SAFsingle_complex(const mxArray* in, float_complex** out, int* nDims, int** pDims)
{ 
    int i, j, k, numElements;
 #if MX_HAS_INTERLEAVED_COMPLEX 
    mxComplexDouble* inMatrix;
#else
    double *inMatrix_r;
    double *inMatrix_i;
#endif
    mwSize nDims_mx;
    const mwSize *pDims_mx;
    
    /* Find dimensionality of input */
    nDims_mx = mxGetNumberOfDimensions(in);
    pDims_mx = mxGetDimensions(in);
    
    /* convert mwSize->int */
    (*nDims) = (int)nDims_mx;
    (*pDims) = malloc1d(*nDims*sizeof(int));
    for(i=0; i<(*nDims); i++)
        (*pDims)[i] = (int)pDims_mx[i];
    
    /* Find number of elements */
    numElements = 1;
    for(i=0; i<(*nDims); i++)
        numElements *= (*pDims)[i];
     
    /* Convert input mex array to saf array */
#if MX_HAS_INTERLEAVED_COMPLEX 
    inMatrix = mxGetData(in);
#else 
    inMatrix_r = mxGetPr(in);
    inMatrix_i = mxGetPi(in);
#endif 
    if((*out)==NULL)
        (*out) = malloc1d(numElements*sizeof(float_complex)); 
    
    /* column-major -> row-major */
    switch(*nDims){
        case 0: /* scalar */ break;
        case 1: assert(0); break;
        case 2: 
#if MX_HAS_INTERLEAVED_COMPLEX 
#else 
            for(i=0; i<(*pDims)[0]; i++)
                for(j=0; j<(*pDims)[1]; j++)
                    (*out)[i*(*pDims)[1]+j] = cmplxf((float)inMatrix_r[j*(*pDims)[0]+i], (float)inMatrix_i[j*(*pDims)[0]+i]);
#endif 
            break;
        case 3:
#if MX_HAS_INTERLEAVED_COMPLEX 
#else 
            for(i=0; i<(*pDims)[0]; i++)
                for(j=0; j<(*pDims)[1]; j++)
                    for(k=0; k<(*pDims)[2]; k++)
                        (*out)[i*(*pDims)[1]*(*pDims)[2]+j*(*pDims)[2]+k] = cmplxf((float)inMatrix_r[k*(*pDims)[1]*(*pDims)[0]+ j*(*pDims)[0]+i], (float)inMatrix_i[k*(*pDims)[1]*(*pDims)[0]+ j*(*pDims)[0]+i]);
#endif 
            break;
    }
}

void SAFsingle2MEXdouble(float* in, int nDims, int* dims, mxArray** out)
{
    int i,j;
    double* pData;
    mwSize nDims_mx;
    mwSize* pDims_mx;
    nDims_mx = (mwSize)nDims;
    
    /* Define dimensionality of output and convert to mwSize */
    pDims_mx = malloc1d(nDims*sizeof(mwSize));
    for(i=0; i<nDims; i++)
        pDims_mx[i] = (mwSize)dims[i];
     
    /* create and copy data to output */ 
    *out = mxCreateNumericArray(nDims_mx, pDims_mx, mxDOUBLE_CLASS, mxREAL);
    pData = mxGetData(*out); 
    
    /* row-major -> column-major */
    switch(nDims){
        case 0: /* scalar */ break;
        case 1: assert(0); break;
        case 2: 
            for(i=0; i<dims[0]; i++)
                for(j=0; j<dims[1]; j++)
                    pData[j*dims[0]+i] = (double)in[i*dims[1]+j];
            break;
    }
    
    /* clean-up */
    free(pDims_mx);
}

void SAFsingle2MEXdouble_complex(float_complex* in, int nDims, int* dims, mxArray** out)
{
    int i,j,k;
#if MX_HAS_INTERLEAVED_COMPLEX 
    mxComplexDouble* pData;
#else
    double *pData_r, *pData_i;
#endif
    mwSize nDims_mx;
    mwSize* pDims_mx;
    nDims_mx = (mwSize)nDims;
    
    /* Define dimensionality of output and convert to mwSize */
    pDims_mx = malloc1d(nDims*sizeof(mwSize));
    for(i=0; i<nDims; i++)
        pDims_mx[i] = (mwSize)dims[i];
      
    /* create and copy data to output */ 
    *out = mxCreateNumericArray(nDims_mx, pDims_mx, mxDOUBLE_CLASS, mxCOMPLEX);
#if MX_HAS_INTERLEAVED_COMPLEX 
    pData = mxGetData(*out); 
#else 
    pData_r = mxGetPr(*out);
    pData_i = mxGetPi(*out);
#endif
    
    /* row-major -> column-major */
    switch(nDims){
        case 0: /* scalar */ break;
        case 1: assert(0); break;
        case 2: 
            for(i=0; i<dims[0]; i++){
                for(j=0; j<dims[1]; j++){
#if MX_HAS_INTERLEAVED_COMPLEX 
                    pData[j*dims[0]+i].real = (double)crealf(in[i*dims[1]+j]);
                    pData[j*dims[0]+i].imag = (double)cimagf(in[i*dims[1]+j]);
#else 
                    pData_r[j*dims[0]+i] = (double)crealf(in[i*dims[1]+j]);
                    pData_i[j*dims[0]+i] = (double)cimagf(in[i*dims[1]+j]);
#endif
                }
            } 
        case 3: 
            for(i=0; i<dims[0]; i++){
                for(j=0; j<dims[1]; j++){
                    for(k=0; k<dims[2]; k++){
#if MX_HAS_INTERLEAVED_COMPLEX 
                        assert(0);
                        pData[k*dims[1]*dims[0]+j*dims[0]+i].real = (double)crealf(in[i*dims[1]*dims[2]+j*dims[2]+k]);
                        pData[k*dims[1]*dims[0]+j*dims[0]+i].imag = (double)cimagf(in[i*dims[1]*dims[2]+j*dims[2]+k]);
#else 
                        pData_r[k*dims[1]*dims[0]+j*dims[0]+i] = (double)crealf(in[i*dims[1]*dims[2]+j*dims[2]+k]);
                        pData_i[k*dims[1]*dims[0]+j*dims[0]+i] = (double)cimagf(in[i*dims[1]*dims[2]+j*dims[2]+k]); 
#endif
                    }
                }
            }
            break;
    }
     
    /* clean-up */
    free(pDims_mx);
}