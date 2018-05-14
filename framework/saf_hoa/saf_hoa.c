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
 *     saf_hoa.c
 * Description:
 *     A Higher-order Ambisonics C library; largely derived from the MatLab library by
 *     Archontis Politis: https://github.com/polarch/Higher-Order-Ambisonics
 * Dependencies:
 *     saf_utilities, saf_sh, saf_vbap
 * Author, date created:
 *     Leo McCormack, 19.03.2018
 */

#include "saf_hoa.h"
#include "saf_hoa_internal.h"

/* Zotter, F., Frank, M. (2012). All-Round Ambisonic Panning and Decoding. Journal of the Audio Engineering Society, 60(10), 807?820. */
void getMaxREweights
(
    int order,
    float* a_n /* (order+1)^2 x (order+1)^2 */
)
{
    int n, i, idx, nSH;
    float* ppm;
    
    nSH = (order+1)*(order+1);
    memset(a_n, 0, nSH*nSH*sizeof(float));
    ppm = calloc((order+1),sizeof(float));
    idx = 0;
    for(n=0; n<=order; n++){
        legendreP(n, cosf(137.9f*(M_PI/180.0f)/((float)order+1.51f)), ppm);
        
        /* store the first Legendre polynomial value for each order along the diagonal of a_n */
        for(i = 0; i<2*n+1; i++)
            a_n[(idx+i)*nSH + (idx+i)] = ppm[0];
        idx += 2*n+1;
    }
    free(ppm);
}

void getAmbiDecoder
(
    float* ls_dirs_deg,
    int nLS,
    AMBI_DECODER_METHODS method,
    int order,
    float** decMtx
)
{
    int i, j, nSH;
    float* Y_ls;
    
    nSH = (order+1) * (order+1);
    Y_ls = NULL;
    free(*decMtx);
    (*decMtx) = malloc(nLS*nSH*sizeof(float));
    switch(method){
        default:
        case DECODER_DEFAULT:
        case DECODER_SAD:
            /* Sampling Ambisonic Decoder (SAD) is simply the loudspeaker spherical harmonic matrix scaled by the
             * number of loudspeakers. */
            getRSH(order, ls_dirs_deg, nLS, &Y_ls);
            for(i=0; i<nLS; i++)
                for(j=0; j<nSH; j++)
                    (*decMtx)[i*nSH+j] = Y_ls[j*nLS + i]/(float)nLS;
            free(Y_ls);
            break;
           
        case DECODER_MMD:
            /* Mode-Matching Decoder (MMD) is simply the psuedo inverse of the loudspeaker spherical harmonic
             * matrix. */
            getRSH(order, ls_dirs_deg, nLS, &Y_ls);
            utility_spinv(Y_ls, nSH, nLS, (*decMtx));
            free(Y_ls);
            break;
            
        case DECODER_EPAD:
            getEPAD(order, ls_dirs_deg, nLS, decMtx);
            break;
            
        case DECODER_ALLRAD:
            getAllRAD(order, ls_dirs_deg, nLS, decMtx);
            break;
    }
}






















