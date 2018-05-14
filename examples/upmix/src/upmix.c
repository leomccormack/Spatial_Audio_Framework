/*
 Copyright 2018 Leo McCormack
 
 Permission to use, copy, modify, and/or distribute this software for any purpose with or
 without fee is hereby granted, provided that the above copyright notice and this permission
 notice appear in all copies.
 
 THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES WITH REGARD TO
 THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS. IN NO EVENT
 SHALL THE AUTHOR BE LIABLE FOR ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR
 ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF
 CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE
 OR PERFORMANCE OF THIS SOFTWARE.
*/
/*
 * Filename:
 *     upmix.c
 * Description:
 *     A (soon to be) collection of upmixing algorithms. However, currently, only stereo to
 *     5.x is supported, utilising a modified version of the direct-ambient decomposition
 *     approach described in: Faller, C. (2006). Multiple-loudspeaker playback of stereo
 *     signals. Journal of the Audio Engineering Society, 54(11), 1051-1064.
 * Dependencies:
 *     saf_utilities, afSTFTlib, saf_vbap, saf_sh
 * Author, date created:
 *     Leo McCormack, 04.04.2018
 */
 
#include "upmix_internal.h"

static inline float matlab_fmodf(float x, float y) {
    float tmp = fmodf(x, y);
    return tmp >= 0 ? tmp : tmp + y;
}

void upmix_create
(
    void ** const phUpmx
)
{
    upmix_data* pData = (upmix_data*)malloc(sizeof(upmix_data));
    if (pData == NULL) { return;/*error*/ }
    *phUpmx = (void*)pData;
    int t, ch;
    
    /* time-frequency transform + buffers */
    afSTFTinit(&(pData->hSTFT), HOP_SIZE, MAX_NUM_INPUT_CHANNELS, MAX_NUM_OUTPUT_CHANNELS, 0, 1);
    pData->STFTInputFrameTF = (complexVector**)malloc2d(TIME_SLOTS, MAX_NUM_INPUT_CHANNELS, sizeof(complexVector));
    pData->STFTOutputFrameTF = (complexVector**)malloc2d(TIME_SLOTS, MAX_NUM_OUTPUT_CHANNELS, sizeof(complexVector));
    for(t=0; t<TIME_SLOTS; t++) {
        for(ch=0; ch< MAX_NUM_INPUT_CHANNELS; ch++) {
            pData->STFTInputFrameTF[t][ch].re = (float*)calloc(HYBRID_BANDS, sizeof(float));
            pData->STFTInputFrameTF[t][ch].im = (float*)calloc(HYBRID_BANDS, sizeof(float));
        }
        for(ch=0; ch< MAX_NUM_OUTPUT_CHANNELS; ch++) {
            pData->STFTOutputFrameTF[t][ch].re = (float*)calloc(HYBRID_BANDS, sizeof(float));
            pData->STFTOutputFrameTF[t][ch].im = (float*)calloc(HYBRID_BANDS, sizeof(float));
        }
    }
    pData->tempHopFrameTD = (float**)malloc2d( MAX(MAX_NUM_OUTPUT_CHANNELS, MAX_NUM_INPUT_CHANNELS), HOP_SIZE, sizeof(float));
    
    /* internal parameters */
    pData->reInitCodec = 1;
    pData->pars = malloc(sizeof(codecPars));
    codecPars* pars = pData->pars;
    pars->grid_vbap_gtable = NULL;
    pars->grp_freqs = NULL;
    pars->grp_idx = NULL;
    pars->prev_est_dir = NULL; 
     
    /* user parameters */
    pData->pValueCoeff = 0.5f;
    pData->paramAvgCoeff = 0.0f;
    pData->scaleDoAwidth = 1.0f;
    pData->covAvg = 0.85f;
}


void upmix_destroy
(
    void ** const phUpmx
)
{
    upmix_data *pData = (upmix_data*)(*phUpmx);
    int t, ch;

    if (pData != NULL) {
        afSTFTfree(pData->hSTFT);
        for (t = 0; t<TIME_SLOTS; t++) {
            for(ch=0; ch< MAX_NUM_INPUT_CHANNELS; ch++) {
                free(pData->STFTInputFrameTF[t][ch].re);
                free(pData->STFTInputFrameTF[t][ch].im);
            }
            for (ch = 0; ch<MAX_NUM_OUTPUT_CHANNELS; ch++) {
                free(pData->STFTOutputFrameTF[t][ch].re);
                free(pData->STFTOutputFrameTF[t][ch].im);
            }
        }
        free2d((void**)pData->STFTInputFrameTF, TIME_SLOTS);
        free2d((void**)pData->STFTOutputFrameTF, TIME_SLOTS);
        free2d((void**)pData->tempHopFrameTD, MAX(MAX_NUM_INPUT_CHANNELS, MAX_NUM_OUTPUT_CHANNELS));
 
        free(pData);
        pData = NULL;
    }
}

void upmix_init
(
    void * const hUpmx,
    int          sampleRate
)
{
    upmix_data *pData = (upmix_data*)(hUpmx);
    int band;
    
    /* define frequency vectors */
    pData->fs = sampleRate;
    if(pData->fs==44100)
        for(band=0; band <HYBRID_BANDS; band++)
            pData->freqVector[band] = (float)__afCenterFreq44100[band];
    else
        for(band=0; band <HYBRID_BANDS; band++)
            pData->freqVector[band] = (float)__afCenterFreq48e3[band];
    
    /* define pValues */
    getPvalues(pData->pValueCoeff, (float*)pData->freqVector, HYBRID_BANDS, (float*)pData->pValues);
    
    /* default starting values */
    memset(pData->Cx, 0, HYBRID_BANDS*MAX_NUM_INPUT_CHANNELS*MAX_NUM_INPUT_CHANNELS*sizeof(float_complex));
    memset(pData->inputframeTF_buffer, 0, HYBRID_BANDS*MAX_NUM_INPUT_CHANNELS*DIFFUSE_DELAY_TIME_SLOTS*sizeof(float_complex));
    pData->buffer_rIdx = 1;
    pData->buffer_wIdx = 0;
}

void upmix_process
(
    void  *  const hUpmx,
    float ** const inputs,
    float ** const outputs,
    int            nInputs,
    int            nOutputs,
    int            nSamples,
    int            isPlaying
)
{
    upmix_data *pData = (upmix_data*)(hUpmx);
    codecPars* pars = pData->pars;
    int t, sample, ch, i, j, k, band, grpband, num_grpBands, idx2D, ls;
    int* grp_bands;
    float est_dir, dummy;
    double Cx_grp00, Cx_grp11, ICC_01, A1, A2, B, C, src_en, diff_en, src_diff_en, w_denom;
    double pv_f, gains2D_sum_pvf;
    float est_dir_xyz[3], prev_est_dir_xyz[3], est_dir_xyz_avg[3];
    double gains2D[MAX_NUM_OUTPUT_CHANNELS];
    double w_src[1][2], w_diff[2][2];
    double Ms_S[MAX_NUM_OUTPUT_CHANNELS][MAX_NUM_INPUT_CHANNELS];
    double Ms_N[MAX_NUM_OUTPUT_CHANNELS][MAX_NUM_INPUT_CHANNELS];
    double Md_N[MAX_NUM_OUTPUT_CHANNELS][MAX_NUM_INPUT_CHANNELS];
    float_complex new_Cx[MAX_NUM_INPUT_CHANNELS][MAX_NUM_INPUT_CHANNELS];
    float_complex Cx_grp[MAX_NUM_INPUT_CHANNELS][MAX_NUM_INPUT_CHANNELS];
    const float_complex calpha = cmplxf(1.0f, 0.0f), cbeta = cmplxf(0.0f, 0.0f);
    const double mix_LR[MAX_NUM_OUTPUT_CHANNELS][MAX_NUM_INPUT_CHANNELS] = { {sqrt(4.0), 0.0}, {0.0, sqrt(4.0)}, {0.0, 0.0f}, {0.0, 0.0}, {0.0, 0.0}};
    const double mix_LsRs[MAX_NUM_OUTPUT_CHANNELS][MAX_NUM_INPUT_CHANNELS] = { { 0.0, 0.0}, {0.0,  0.0}, {0.0, 0.0}, {sqrt(4.0), 0.0}, {0.0, sqrt(4.0)}};
    
    /* local copies of user parameters */
    int nLoudspeakers;
    float paramAvgCoeff, scaleDoAwidth, covAvg;
    
    /* reinitialise codec if needed */
    if(pData->reInitCodec){
        pData->reInitCodec = 2;
        upmix_initCodec(hUpmx);
        pData->reInitCodec = 0;
    }
    if ((nSamples == FRAME_SIZE) && (isPlaying == 1) && (pData->reInitCodec == 0) ) {
        nLoudspeakers = pData->nLoudspeakers;
        paramAvgCoeff = pData->paramAvgCoeff;
        scaleDoAwidth = pData->scaleDoAwidth;
        covAvg = pData->covAvg;
    
        /* Load time-domain data */
        for(i=0; i < MIN(MAX_NUM_INPUT_CHANNELS,nInputs); i++)
            memcpy(pData->inputFrameTD[i], inputs[i], FRAME_SIZE * sizeof(float));
        for(; i<MAX_NUM_INPUT_CHANNELS; i++)
            memset(pData->inputFrameTD[i], 0, FRAME_SIZE * sizeof(float));
        
        /* Apply time-frequency transform (TFT) */
        for ( t=0; t< TIME_SLOTS; t++) {
            for( ch=0; ch < MAX_NUM_INPUT_CHANNELS; ch++)
                for ( sample=0; sample < HOP_SIZE; sample++)
                    pData->tempHopFrameTD[ch][sample] = pData->inputFrameTD[ch][sample + t*HOP_SIZE];
            afSTFTforward(pData->hSTFT, (float**)pData->tempHopFrameTD, (complexVector*)pData->STFTInputFrameTF[t]);
        }
        for(band=0; band<HYBRID_BANDS; band++)
            for( ch=0; ch < MAX_NUM_INPUT_CHANNELS; ch++)
                for ( t=0; t<TIME_SLOTS; t++)
                    pData->inputframeTF[band][ch][t] = cmplxf(pData->STFTInputFrameTF[t][ch].re[band], pData->STFTInputFrameTF[t][ch].im[band]);
   
        /* update covarience matrix per band */
        for(band=0; band<HYBRID_BANDS; band++){
            cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasConjTrans, MAX_NUM_INPUT_CHANNELS, MAX_NUM_INPUT_CHANNELS, TIME_SLOTS, &calpha,
                        pData->inputframeTF[band], TIME_SLOTS,
                        pData->inputframeTF[band], TIME_SLOTS, &cbeta,
                        new_Cx, MAX_NUM_INPUT_CHANNELS);
            
            /* average over time */
            for(i=0; i<MAX_NUM_INPUT_CHANNELS; i++)
                for(j=0; j<MAX_NUM_INPUT_CHANNELS; j++)
                    pData->Cx[band][i][j] = new_Cx[i][j]*(1.0f-covAvg) + pData->Cx[band][i][j] * covAvg;
        }
        
        /* Calculate mixing matrices for upmixing */
        for (grpband = 0; grpband < pars->nGrpBands-1; grpband++) {
            /* define which bands make up the bark/erb scale grouping */
            num_grpBands = pars->grp_idx[grpband+1] - pars->grp_idx[grpband];
            grp_bands = (int*)malloc(num_grpBands*sizeof(int));
            for(i=pars->grp_idx[grpband], j=0; i<pars->grp_idx[grpband+1]; i++, j++)
                grp_bands[j] = i-1; /* -1 as the indices start from 1 */
            
            /* sum the covarience matrices for these grouped bands and calculate ICC between left and right channels */
            memset(Cx_grp, 0, MAX_NUM_INPUT_CHANNELS*MAX_NUM_INPUT_CHANNELS*sizeof(float_complex));
            for(i=0; i<num_grpBands; i++)
                for(j=0; j<MAX_NUM_INPUT_CHANNELS; j++)
                    for(k=0; k<MAX_NUM_INPUT_CHANNELS; k++)
                        Cx_grp[j][k] = ccaddf(Cx_grp[j][k], pData->Cx[grp_bands[i]][j][k]);
            Cx_grp00 = creal((double)Cx_grp[0][0]);
            Cx_grp11 = creal((double)Cx_grp[1][1]);
            ICC_01 = creal((double)Cx_grp[0][1])/(sqrt(Cx_grp00*Cx_grp11)+2.23e-9);
            
            /* estimate the short-time energy of the source and diffuse signals */
            /* Faller, C. (2006). Multiple-loudspeaker playback of stereo signals. Journal of the Audio Engineering Society, 54(11), 1051-1064. */
            C = crealf(Cx_grp[0][1]);
            B = Cx_grp11 - Cx_grp00 + sqrt( pow(Cx_grp00-Cx_grp11, 2.0) + 4.0 * Cx_grp00*Cx_grp11 * pow(ICC_01, 2.0));
            A1 = B/(2.0*C + 2.23e-9);
            A2 = (2.0*C)/(B+ 2.23e-9);
            src_en = (2.0f*C*C)/(B+2.23e-9);
            diff_en = Cx_grp00 - src_en;
            src_diff_en = src_en * diff_en;
            
            /* Determine source azimuth (-180..180) based on the real-valued amplitude ratio */
            if (A1<=1.0 && A1>=-1.0)
                est_dir = A1<0.0 ? 150.0*A1 - 30.0f : 30.0*A1 - 30.0;
            else if (A2<=1.0 && A2>=-1.0)
                est_dir = A2<0.0 ? -(150.0*A2 - 30.0f) : -(30.0*A2 - 30.0);
            else
                est_dir = 0.0;
            est_dir = est_dir*scaleDoAwidth; /* manipulate the width by scaling this estimate */
            
            /* Average source DoA over time */
            unitSph2Cart(est_dir*M_PI/180.0f, 0.0f, est_dir_xyz);
            unitSph2Cart(pars->prev_est_dir[grpband]*M_PI/180.0f, 0.0f, prev_est_dir_xyz);
            for(i=0; i<3; i++)
                est_dir_xyz_avg[i] = (1.0f-paramAvgCoeff)*est_dir_xyz[i] + paramAvgCoeff*prev_est_dir_xyz[i];
            unitCart2Sph_aziElev( est_dir_xyz_avg, &est_dir, &dummy);
            est_dir *= 180.0f/M_PI;
            pars->prev_est_dir[grpband] = est_dir;
            
            /* estimate the mixing weights reqiured to obtain source and diffuse components via a least-square approximation */
            w_denom = (A1*A1+1.0)*src_diff_en + diff_en*diff_en + 2.23e-9;
            w_src[0][0] = (src_diff_en)/w_denom;
            w_src[0][1] = w_src[0][0]*A1;
            w_diff[0][0] = ((A1*A1)*src_diff_en + diff_en*diff_en)/w_denom;
            w_diff[0][1] = (-A1*src_diff_en + diff_en*diff_en)/w_denom;
            w_diff[1][0] = w_diff[0][1];
            w_diff[1][1] = (src_diff_en + diff_en*diff_en)/w_denom;
            
            for(band=0; band<num_grpBands; band++){
                /* Pull loudspeaker gains from vbap table */
                idx2D = (int)((matlab_fmodf(est_dir+180.0f,360.0f)/pars->vbap_azi_res)+0.5f);
                for (ls = 0; ls < nLoudspeakers; ls++)
                    gains2D[ls] = (double)pars->grid_vbap_gtable[idx2D*nLoudspeakers+ls];
                 
                /* apply pValue normalisation (i.e. amplitude normalises the VBAP gains for low frequencies depending on room) */
                pv_f = pData->pValues[grp_bands[band]];
                if(pv_f != 2.0f){
                    gains2D_sum_pvf = 0.0f;
                    for (ls = 0; ls < nLoudspeakers; ls++)
                        gains2D_sum_pvf += pow(MAX(gains2D[ls], 0.0), pv_f);
                    gains2D_sum_pvf = pow(gains2D_sum_pvf, 1.0/(pv_f+2.23e-9));
                    for (ls = 0; ls < nLoudspeakers; ls++)
                        gains2D[ls] = gains2D[ls] / (gains2D_sum_pvf+2.23e-9);
                }

                /* formulate direct mixing matrix */
                cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, MAX_NUM_OUTPUT_CHANNELS, MAX_NUM_INPUT_CHANNELS, 1, 1.0,
                            (const double*)gains2D, 1,
                            (const double*)w_src, MAX_NUM_INPUT_CHANNELS, 0.0,
                            (double*)Ms_S, MAX_NUM_INPUT_CHANNELS);
                cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, MAX_NUM_OUTPUT_CHANNELS, MAX_NUM_INPUT_CHANNELS, MAX_NUM_INPUT_CHANNELS, 1.0,
                            (const double*)mix_LR, MAX_NUM_INPUT_CHANNELS,
                            (const double*)w_diff, MAX_NUM_INPUT_CHANNELS, 0.0,
                            (double*)Ms_N, MAX_NUM_INPUT_CHANNELS);
                for(i=0; i<MAX_NUM_OUTPUT_CHANNELS; i++)
                    for(j=0; j<MAX_NUM_INPUT_CHANNELS; j++)
                        pData->new_Ms[grp_bands[band]][i][j] = cmplxf((float)Ms_S[i][j] + (float)Ms_N[i][j], 0.0f);
                
                /* formulate diffuse mixing matrix */
                cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, MAX_NUM_OUTPUT_CHANNELS, MAX_NUM_INPUT_CHANNELS, MAX_NUM_INPUT_CHANNELS, 1.0,
                            (const double*)mix_LsRs, MAX_NUM_INPUT_CHANNELS,
                            (const double*)w_diff, MAX_NUM_INPUT_CHANNELS, 0.0f,
                            (double*)Md_N, MAX_NUM_INPUT_CHANNELS);
                for(i=0; i<MAX_NUM_OUTPUT_CHANNELS; i++)
                    for(j=0; j<MAX_NUM_INPUT_CHANNELS; j++)
                        pData->new_Md[grp_bands[band]][i][j] = cmplxf(pars->diff_lpf[grp_bands[band]] * (float)Md_N[i][j], 0.0f);
            }
            free(grp_bands);
        }
        
        /* obtain delayed inputframe */
        for(t=0; t<TIME_SLOTS; t++){
            for(band=0; band<HYBRID_BANDS; band++){
                for(ch=0; ch<MAX_NUM_INPUT_CHANNELS; ch++){
                    pData->inputframeTF_buffer[band][ch][pData->buffer_wIdx] = pData->inputframeTF[band][ch][t];
                    pData->inputframeTF_del[band][ch][t] = pData->inputframeTF_buffer[band][ch][pData->buffer_rIdx];
                }
            }
            /* increment circular buffer indices */
            pData->buffer_wIdx++;
            pData->buffer_rIdx++;
            if(pData->buffer_wIdx==DIFFUSE_DELAY_TIME_SLOTS)
                pData->buffer_wIdx = 0;
            if(pData->buffer_rIdx==DIFFUSE_DELAY_TIME_SLOTS)
                pData->buffer_rIdx = 0;
        }
        
        /* Apply mixing matrices to current and delayed intputframe */
        for(band=0; band<HYBRID_BANDS; band++){
            cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, MAX_NUM_OUTPUT_CHANNELS, TIME_SLOTS, MAX_NUM_INPUT_CHANNELS, &calpha,
                        (const float*)pData->new_Ms[band], MAX_NUM_INPUT_CHANNELS,
                        (const float*)pData->inputframeTF[band], TIME_SLOTS, &cbeta,
                        (float*)pData->directframeTF[band], TIME_SLOTS);
            cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, MAX_NUM_OUTPUT_CHANNELS, TIME_SLOTS, MAX_NUM_INPUT_CHANNELS, &calpha,
                        (const float*)pData->new_Md[band], MAX_NUM_INPUT_CHANNELS,
                        (const float*)pData->inputframeTF_del[band], TIME_SLOTS, &cbeta,
                        (float*)pData->diffuseframeTF[band], TIME_SLOTS);
            /* combine */
            for(ch=0; ch<MAX_NUM_OUTPUT_CHANNELS; ch++)
                for(t=0; t<TIME_SLOTS; t++)
                    pData->outputframeTF[band][ch][t] = ccaddf(pData->directframeTF[band][ch][t], pData->diffuseframeTF[band][ch][t]);
        }
        
        /* inverse-TFT */
        for (band = 0; band < HYBRID_BANDS; band++) {
            for (ch = 0; ch < MAX_NUM_OUTPUT_CHANNELS; ch++) {
                for (t = 0; t < TIME_SLOTS; t++) {
                    pData->STFTOutputFrameTF[t][ch].re[band] = crealf(pData->outputframeTF[band][ch][t]);
                    pData->STFTOutputFrameTF[t][ch].im[band] = cimagf(pData->outputframeTF[band][ch][t]);
                }
            }
        }
        for (t = 0; t < TIME_SLOTS; t++) {
            afSTFTinverse(pData->hSTFT, pData->STFTOutputFrameTF[t], pData->tempHopFrameTD);
            for (ch = 0; ch < MIN(MAX_NUM_OUTPUT_CHANNELS, nOutputs); ch++)
                for (sample = 0; sample < HOP_SIZE; sample++)
                    outputs[ch][sample + t* HOP_SIZE] = pData->tempHopFrameTD[ch][sample];
            for (; ch < nOutputs; ch++)
                for (sample = 0; sample < HOP_SIZE; sample++)
                    outputs[ch][sample + t* HOP_SIZE] = 0.0f;
        }
    }
    else{
        for (ch=0; ch < nOutputs; ch++)
            memset(outputs[ch],0, FRAME_SIZE*sizeof(float));
    } 
}


/* Set Functions */

void upmix_setPValueCoeff(void* const hUpmx, float newValue)
{
    upmix_data *pData = (upmix_data*)(hUpmx);
    pData->pValueCoeff = newValue;
    getPvalues(pData->pValueCoeff, pData->freqVector, HYBRID_BANDS, pData->pValues);
}

void upmix_setParamAvgCoeff(void* const hUpmx, float newValue)
{
    upmix_data *pData = (upmix_data*)(hUpmx);
    pData->paramAvgCoeff = newValue;
}

void upmix_setScaleDoAwidth(void* const hUpmx, float newValue)
{
    upmix_data *pData = (upmix_data*)(hUpmx);
    pData->scaleDoAwidth = newValue;
}

void upmix_setCovAvg(void* const hUpmx, float newValue)
{
    upmix_data *pData = (upmix_data*)(hUpmx);
    pData->covAvg = newValue;
}


/* Get Functions */

float upmix_getPValueCoeff(void* const hUpmx)
{
    upmix_data *pData = (upmix_data*)(hUpmx);
    return pData->pValueCoeff;
}

float upmix_getParamAvgCoeff(void* const hUpmx)
{
    upmix_data *pData = (upmix_data*)(hUpmx);
    return pData->paramAvgCoeff;
}

float upmix_getScaleDoAwidth(void* const hUpmx)
{
    upmix_data *pData = (upmix_data*)(hUpmx);
    return pData->scaleDoAwidth;
}

float upmix_getCovAvg(void* const hUpmx)
{
    upmix_data *pData = (upmix_data*)(hUpmx);
    return pData->covAvg;
}







    
    
