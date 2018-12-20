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
 *     saf_hoa_internal.c
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

/* Zotter, F., Pomberger, H., Noisternig, M. (2012). Energy-Preserving Ambisonic Decoding. Acta Acustica United with Acustica, 98(1), 37:47.
* The function has been written to also work when the number of spherical harmonic components exceeds the number of loudspeakers.
* In which case, the 'U' matrix from the SVD is truncated instead. However, ideally... nLS > nSH, like in the paper */
void getEPAD
(
    int order,
    float* ls_dirs_deg,
    int nLS,
    float **decMtx
)
{
    int i, j, nSH;
    float* Y_ls, *U, *S, *V, *U_tr, *V_tr;
    
    nSH = (order+1)*(order+1);
    Y_ls = U = S = V = NULL;
    getRSH(order, ls_dirs_deg, nLS, &Y_ls);
    utility_ssvd(Y_ls, nSH, nLS, &U, &S, &V);
    if(nSH>nLS){
        /* truncate the U matrix */
        U_tr = malloc(nSH*nLS*sizeof(float));
        for(i=0; i<nSH; i++)
            for(j=0; j<nLS; j++)
                U_tr[i*nLS+j] = U[i*nSH+j];
        cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasTrans, nLS, nSH, nLS, 1.0f,
                    V, nLS,
                    U_tr, nLS, 0.0f,
                    (*decMtx), nSH);
        free(U_tr);
    }
    else{
        /* truncate the V matrix (NOT V^T!) */
        V_tr = malloc(nLS*nSH*sizeof(float));
        for(i=0; i<nLS; i++)
            for(j=0; j<nSH; j++)
                V_tr[i*nSH+j] = V[i*nLS+j];
        cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasTrans, nLS, nSH, nSH, 1.0f,
                    V_tr, nSH,
                    U, nSH, 0.0f,
                    (*decMtx), nSH);
        free(V_tr);
    }
    for(i=0; i<nLS*nSH; i++)
        (*decMtx)[i] /= (float)nLS;
    
    free(U);
    free(S);
    free(V);
    free(Y_ls);
}

/* Zotter, F., Frank, M. (2012). All-Round Ambisonic Panning and Decoding. Journal of the Audio Engineering Society, 60(10), 807:820. */
void getAllRAD(int order, float* ls_dirs_deg, int nLS, float **decMtx)
{
    int i, t, nDirs_td, N_gtable, nGroups, nSH;
    float* Y_td, *G_td, *t_dirs;
    
    nSH = (order+1)*(order+1);
    
    /* define a sufficiently dense t-design for this decoding order, as to conserve omni energy */
    t = 4*order;
    Y_td = G_td = NULL;
    if(t<=21){
        /* suitable for up to 5th order */
        nDirs_td = __Tdesign_nPoints_per_degree[t-1];
        t_dirs = (float*)__HANDLES_Tdesign_dirs_deg[t-1];
    }
    else if (order > 7){
        /* suitable for >7th order */
        nDirs_td = 5100; /* Minimum t-design of degree 100 has 5100 points */
        t_dirs = (float*)__Tdesign_degree_100_dirs_deg;
    }
    else{
        /* suitable for 6th & 7th order */
        nDirs_td = 480; /* Minimum t-design of degree 30 has 480 points (sufficient for up to 7th order) */
        t_dirs = (float*)__Tdesign_degree_30_dirs_deg;
    }
    
    /* calculate vbap gains and SH matrix for this t-design */
    generateVBAPgainTable3D_srcs(t_dirs, nDirs_td, ls_dirs_deg, nLS, 0, 0, 0.0f, &G_td, &N_gtable, &nGroups);
    getRSH(order, t_dirs, nDirs_td, &Y_td);
    
    /* AllRAD decoder is simply (G_td * T_td * 1/nDirs_td) */
    cblas_sgemm(CblasRowMajor, CblasTrans, CblasTrans, nLS, nSH, nDirs_td, 1.0f,
                G_td, nLS,
                Y_td, nDirs_td, 0.0f,
                (*decMtx), nSH);
    for(i=0; i<nLS*nSH; i++)
        (*decMtx)[i] /= (float)nDirs_td;
 
    free(Y_td);
    free(G_td);
}


