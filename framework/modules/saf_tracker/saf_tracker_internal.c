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

/**
 * @file saf_tracker_internal.c
 * @brief Particle filtering based tracker
 *
 * Based on the RBMCDA Matlab toolbox (GPLv2 license) by Simo Särkkä and Jouni
 * Hartikainen (Copyright (C) 2003-2008):
 *     https://users.aalto.fi/~ssarkka/#softaudio
 *
 * And also inspired by the work of Sharath Adavanne, Archontis Politis, Joonas
 * Nikunen, and Tuomas Virtanen (GPLv2 license):
 *     https://github.com/sharathadavanne/multiple-target-tracking
 *
 * @author Leo McCormack
 * @date 12.08.2020
 */

#include "saf_tracker.h"
#include "saf_tracker_internal.h"

/* ========================================================================== */
/*                             Static Prototypes                              */
/* ========================================================================== */

/**
 * Natural logarithm of gamma function
 *
 * @param[in]     x      input
 * @param[in,out] sgngam (&) sign(Gamma(X))
 *
 * @returns logarithm of the absolute value of the Gamma(X).
 *
 * Adapted from: https://www.alglib.net/download.php#cpp
 * Original Copyright 1984, 1987, 1989, 1992, 2000 by Stephen L. Moshier (GPLv2)
 */
static double lngamma(double x, double* sgngam);

/**
 * Incomplete gamma integral
 *
 * The function is defined by:
 * \verbatim
 *                          x
 *                           -
 *                  1       | |  -t  a-1
 * igam(a,x)  =   -----     |   e   t   dt.
 *                 -      | |
 *                | (a)    -
 *                          0
 * \endverbatim
 *
 * In this implementation both arguments must be positive. The integral is
 * evaluated by either a power series or continued fraction expansion, depending
 * on the relative values of a and x.
 *
 * Adapted from: https://www.alglib.net/download.php#cpp
 * Original Copyright 1985, 1987, 2000 by Stephen L. Moshier (GPLv2)
 */
static double incompletegamma(double a, double x);

/**
 * Complemented incomplete gamma integral
 *
 * The function is defined by
 * \verbatim
 * igamc(a,x)   =   1 - igam(a,x)
 *
 *                           inf.
 *                             -
 *                    1       | |  -t  a-1
 *              =   -----     |   e   t   dt.
 *                   -      | |
 *                  | (a)    -
 *                            x
 * \endverbatim
 *
 * In this implementation both arguments must be positive. The integral is
 * evaluated by either a power series or continued fraction expansion, depending
 * on the relative values of a and x.
 *
 * Adapted from: https://www.alglib.net/download.php#cpp
 * Original Copyright 1985, 1987, 2000 by Stephen L. Moshier  (GPLv2)
 */
static double incompletegammac(double a, double x);


/* ========================================================================== */
/*                             Internal Functions                             */
/* ========================================================================== */

void tracker3d_predict
(
    void* const hT3d,
    int Tinc
)
{
    tracker3d_data *pData = (tracker3d_data*)(hT3d);
    int i, j, k, n, nDead, isDead, ind;
    int* dead;
    float dt0, dt1, p_death, angle_diff, norm, dot, rand01;
    float cross_Mjk[3];
    MCS_data* S;
    tracker3d_config* tpars = &(pData->tpars);
#ifdef TRACKER_VERBOSE
    char event[256];
#endif

    /* prep */
    nDead = 0;
    dead = NULL;

    /* Loop over particles */
    for (i=0; i<tpars->Np; i++){
        S = (MCS_data*)pData->SS[i];
#ifdef TRACKER_VERBOSE
        sprintf(S->evstr, "Particle: %d ", i);
#endif
        /* Loop over targets */
        for (j=0; j<S->nTargets; j++){

            /* No target has died yet or multiple targets are allowed to die in
             * one prediction step */
            if (nDead==0 || tpars->ALLOW_MULTI_DEATH){
                /* Probability of death */
                dt0 = (float)S->Tcount[j] * S->dt;
                dt1 = dt0 + S->dt * (float)Tinc;
                if (dt0 == 0)
                    p_death = gamma_cdf(dt1, tpars->alpha_death, tpars->beta_death, 0.0f);
                else
                    p_death = 1.0f - (1.0f-gamma_cdf(dt1, tpars->alpha_death,tpars->beta_death, 0.0f)) /
                    (1.0f-gamma_cdf(dt0,tpars->alpha_death,tpars->beta_death, 0.0f));

                /* Force probability of death to 1, if this target is too close
                   another target that has been alive longer. */
                if (tpars->FORCE_KILL_TARGETS){
                    for(k=0; k<S->nTargets; k++){
                        if (k!=j){
                            crossProduct3(S->M[j].M, S->M[k].M, cross_Mjk);
                            norm = L2_norm3(cross_Mjk);
                            utility_svvdot(S->M[j].M, S->M[k].M, 3, &dot);
                            angle_diff = atan2f(norm, dot);

                            if (angle_diff < tpars->forceKillAngle_rad && S->Tcount[j] <= S->Tcount[k])
                                p_death = 1;
                        }
                    }
                }

                /* Decide whether target should die */
                rand_0_1(&rand01, 1);
                if (rand01 < p_death){
                    nDead++; /* Target dies */
                    dead = realloc1d(dead, nDead*sizeof(float));
                    dead[nDead-1] = j;
                }
            }

            /* Prediction step if target is alive */
            if( tpars->ALLOW_MULTI_DEATH ){
                isDead = 0;
                for(n=0; n<nDead; n++)
                    if(j==dead[n])
                        isDead = 1;

                /* Kalman Filter prediction for the target if alive */
                if(!isDead)
                    kf_predict6(S->M[j].M, S->P[j].P, pData->A, pData->Q);
            }
            else {
                /* Kalman Filter prediction for the target if alive */
                if ( (nDead==0) || (j != dead[0]) )
                    kf_predict6(S->M[j].M, S->P[j].P, pData->A, pData->Q);
            }
        }

        /* Remove the dead target */
        if (tpars->ALLOW_MULTI_DEATH){
            for(j=0; j<nDead; j++){
                /* Find index of the target to remove */
                ind = -1;
                for(k=0; k<S->nTargets; k++)
                    if(dead[j]==k)
                        ind = k;
                assert(ind != -1);

                /* Shimy target data down by 1... overriding the dead target */
                S->nTargets--;
                if(ind!=S->nTargets){
                    memmove(&S->M[ind], &S->M[ind+1], (S->nTargets-ind)*sizeof(M6));
                    memmove(&S->P[ind], &S->P[ind+1], (S->nTargets-ind)*sizeof(P66));
                    memmove(&S->Tcount[ind], &S->Tcount[ind+1], (S->nTargets-ind)*sizeof(int));
                    memmove(&S->targetIDs[ind], &S->targetIDs[ind+1], (S->nTargets-ind)*sizeof(int));
                }

                /* resize */
                S->M = realloc1d(S->M, S->nTargets*sizeof(M6));
                S->P = realloc1d(S->P, S->nTargets*sizeof(P66));
                S->Tcount = realloc1d(S->Tcount, S->nTargets*sizeof(int));
                S->targetIDs = realloc1d(S->targetIDs, S->nTargets*sizeof(int));

#ifdef TRACKER_VERBOSE
                sprintf(event, ", Target %d died ", ind);
                strcpy(S->evstr, event);
#endif
            }
        }
        else {
            if (nDead==1){
                /* Find index of the target to remove */
                ind = -1;
                for(k=0; k<S->nTargets; k++)
                    if(dead[0]==k)
                        ind = k;
                assert(ind != -1);

                S->nTargets--;
                if(ind!=S->nTargets){
                    memmove(&S->M[ind], &S->M[ind+1], (S->nTargets-ind)*sizeof(M6));
                    memmove(&S->P[ind], &S->P[ind+1], (S->nTargets-ind)*sizeof(P66));
                    memmove(&S->Tcount[ind], &S->Tcount[ind+1], (S->nTargets-ind)*sizeof(int));
                    memmove(&S->targetIDs[ind], &S->targetIDs[ind+1], (S->nTargets-ind)*sizeof(int));
                }

                /* resize */
                S->M = realloc1d(S->M, S->nTargets*sizeof(M6));
                S->P = realloc1d(S->P, S->nTargets*sizeof(P66));
                S->Tcount = realloc1d(S->Tcount, S->nTargets*sizeof(int));
                S->targetIDs = realloc1d(S->targetIDs, S->nTargets*sizeof(int));

#ifdef TRACKER_VERBOSE
                sprintf(event, ", Target %d died ", ind);
                strcpy(S->evstr, event);
#endif
            }
        }

#ifdef TRACKER_VERBOSE
        printf("%s\n", S->evstr);
#endif
    }
}

#define TRACKER3D_MAX_NUM_EVENTS ( 24 )

void tracker3d_update
(
    void* const hT3d,
    float* Y,
    int Tinc
)
{
    tracker3d_data *pData = (tracker3d_data*)(hT3d);
    int i, n_events, count;
    float TP0;
    MCS_data* S;
    tracker3d_config* tpars = &(pData->tpars);

    /* Loop over particles */
//    ev_strs = cell(size(S));
//    ev_IDs = cell(size(S));
    for (i=0; i<tpars->Np; i++){
        S = (MCS_data*)pData->SS[i];

        /* Association priors to targets */
        TP0 = (1.0f-tpars->noiseLikelihood)/(S->nTargets+2.23e-10f);

        /* Possible events: */
        if (tpars->MULTI_ACTIVE){
//            n_events = 1; % clutter
//            nTargets = min(S{i}.nTargets,tpars.maxNactiveTargets);
//            if nTargets > 0
//                n_events = n_events + tpars.ATinds{nTargets}.nCombinations;   % combinations of active targets
//            end
        }
        else{
            /* clutter or 1 of the targets is active */
            n_events = S->nTargets + 2;
        }
        if( S->nTargets < tpars->maxNactiveTargets)
            n_events++; /* chance of a new target */

        assert(n_events<=TRACKER3D_MAX_NUM_EVENTS);

        evt = cell(1,n_events);   % Event descriptions
        evta = cell(1,n_events);  % Event descriptions
        evp = zeros(1,n_events);  % Event priors
        evl = zeros(1,n_events);  % Event likelhoods
        str = cell(1,n_events);   % Structure after each event
        count = 0; /* Event counter */

        /* Association to clutter */
        count = count+1;
        evt{count} = 'Clutter';
        evta{count} = [];
        evp(count) = (1-tpars.init_birth)*tpars.noiseLikelihood;
        evl(count) = tpars.cd;
        str{count} = S{i};
        str{count}.B = 0;
        str{count}.D = 0;

        /* Loop over associations to targets */
        for (j=0; j<S->nTargets; j++){
            /* Compute update result and likelihood for association to signal j */
            [M,P,~,~,~,lhood] = kf_update(S{i}.M{j},S{i}.P{j},Y,tpars.H,tpars.R);

            /* Assocation to target j */
            count = count+1;
            evt{count} = sprintf('Target %d',S{i}.targetIDs(j));
            evta{count} = S{i}.targetIDs(j);
            evp(count) = (1-tpars.init_birth)*TP0;
            evl(count) = lhood;
            str{count} = S{i};
            str{count}.M{j} = M;
            str{count}.P{j} = P;
            %str{count}.Tcount(j) = S{i}.Tcount(j) + Tinc;
            str{count}.Tcount = S{i}.Tcount + Tinc;
            str{count}.B = 0;
            str{count}.D = 0;

            /*  Handles for MULTI_ACTIVE */
//            if tpars.MULTI_ACTIVE
//                strTARGETS{j} = str{count}; %#ok
//                evpTARGETS(j) = evp(count); %#ok
//                evlTARGETS(j) = evl(count); %#ok
//            end
        }

        /* Loop over associations to multiple active target combinations */
        if (tpars->MULTI_ACTIVE == 1) && (nTargets > 0){
//            for nA=1:tpars.ATinds{nTargets}.nCombinationsAbove2
//                count = count+1;
//                inds = tpars.ATinds{nTargets}.indsAbove2{nA};
//                evt{count} = sprintf(['MultiTargets' strcat(int2str(inds))]);
//                evp(count) = prod(evpTARGETS(inds)); % assuming targets are independent of one another
//                evl(count) = prod(evlTARGETS(inds)); % assuming targets are independent of one another
//                evta{count} = [];
//
//                % copy already computed track info
//                str{count} = S{i};
//                for j=inds
//                    evta{count} = [evta{count} S{i}.targetIDs(j)];
//
//                    str{count}.M{j} = strTARGETS{j}.M{j};
//                    str{count}.P{j} = strTARGETS{j}.P{j};
//                    str{count}.Tcount(j) = strTARGETS{j}.Tcount(j);
//                end
//            end
        }

        /* Association to new target */
        if (S->nTargets < tpars->maxNactiveTargets){
            /* Initialization of new target */
            [M,P,~,~,~,lhood] =  kf_update(S{i}.M0,S{i}.P0,Y,tpars.H,tpars.R);

            /* find a free ID */
            for ss = 1:tpars.maxNactiveTargets
                if all(ss~=S{i}.targetIDs)% && all(ss~=S{i}.D)
                    j_new = ss;
                    break;
                end
            end
            count = count+1;
            j = S{i}.nTargets+1;
            evt{count} = sprintf('New target %d',j);
            evta{count} = j;
            evp(count) = tpars.init_birth;
            evl(count) = lhood;
            str{count} = S{i};
            str{count}.nTargets = j;
            str{count}.M{j} = M;
            str{count}.P{j} = P;
            str{count}.Tcount(j) = 0;
            str{count}.targetIDs(j,1) = j_new;
            str{count}.B = j_new;
            str{count}.D = 0;
        }

        /* Draw sample from importance distribution */
        evp = evp./(sum(evp));
        imp = evp.*evl;
        imp = imp./(sum(imp));
        ev = categ_rnd(imp);  /* Event index */

        /* Update particle */
        S{i} = str{ev};       % Copy the updated structure
        ev_comment = evt{ev}; % Event description string

        S{i}.W = S{i}.W*evl(ev)*evp(ev)/imp(ev);
        ev_strs{i} = ev_comment;
        ev_IDs{i} = evta{ev};
    }
}

void tracker3d_particleCreate
(
     void** phPart,
     float W0,
     float M0[6],
     float P0[6][6],
     float dt
)
{
    *phPart = malloc1d(sizeof(MCS_data));
    MCS_data *p = (MCS_data*)(*phPart);

    p->W = W0;
    p->W_prev = W0;
    p->W0 = W0;
    p->nTargets = 0;
    p->nActive = 0;
    p->dt = dt;
    p->M = NULL;
    p->P = NULL;
    memcpy(p->M0.M, M0, 6*sizeof(float));
    memcpy(p->P0.P, P0, 6*6*sizeof(float));
    p->targetIDs = NULL;
    p->activeIDs = NULL;
    p->Tcount = NULL;
}

void tracker3d_particleDestroy
(
     void** phPart
)
{
    MCS_data *p = (MCS_data*)(*phPart); 

    if(p!=NULL){
        free(p->M);
        free(p->P);
        free(p->targetIDs);
        free(p->activeIDs);
        free(p->Tcount);
        p=NULL;
        *phPart = NULL;
    }
}


/* ========================================================================== */
/*                              RBMCDA Functions                              */
/* ========================================================================== */

/* hard-coded for length(M)=6 ... */
void kf_predict6(float M[6], float P[6][6], float A[6][6], float Q[6][6])
{
    float AM[6], AP[6][6], APAT[6][6];

    /* Perform prediction */
    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 6, 1, 6, 1.0f,
                (float*)A, 6,
                (float*)M, 6, 0.0f,
                (float*)AM, 1);
    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 6, 6, 6, 1.0f,
                (float*)A, 6,
                (float*)P, 6, 0.0f,
                (float*)AP, 6);
    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasTrans, 6, 6, 6, 1.0f,
                (float*)AP, 6,
                (float*)A, 6, 0.0f,
                (float*)APAT, 6);

    /* Override M and P, with new M and P */
    memcpy(M, AM, 6*sizeof(float));
    utility_svvadd((float*)APAT, (float*)Q, 36, (float*)P);
}

float gamma_cdf(float x, float gam, float beta, float mu)
{
    /* Convert to standard form */
    x = (x-mu)/beta;

    /* Compute the probability using the imcomplete gamma function */
    return (float)incompletegamma((double)gam, (double)x)/tgamma(x);
}

void lti_disc
(
    float* F,
    int len_N,
    int len_Q,
    float* opt_L,
    float* opt_Qc,
    float dt,
    float* A,
    float* Q
)
{
    int i, j;
    float* L, *Qc;
    float* Fdt;
    float** L_Qc, **L_Qc_LT, **Phi, **ZE, **B, **AB, **AB1_T, **AB2_T, **Q_T;

    /* Defaults: */
    if(opt_L==NULL){ /* Identity */
        L = calloc1d(len_N*len_Q, sizeof(float));
        for(i=0; i<MIN(len_N, len_Q); i++)
            L[i*len_Q+i] = 1.0f;
    }
    else
        L = opt_L;
    if(opt_Qc==NULL) /* zeros */
        Qc = calloc1d(len_Q*len_Q, sizeof(float));
    else
        Qc = opt_Qc;
    Fdt = malloc1d(len_N*len_N*sizeof(float));

    /* Closed form integration of transition matrix */
    for(i=0; i<len_N*len_N; i++)
        Fdt[i] = F[i]*dt;
    gexpm(Fdt, len_N, 0, A);

    /* Closed form integration of covariance by matrix fraction decomposition */
    L_Qc = (float**)malloc2d(len_N, len_Q, sizeof(float));
    L_Qc_LT = (float**)malloc2d(len_N, len_N, sizeof(float));
    Phi = (float**)calloc2d(len_N*2, len_N*2, sizeof(float)); // ca
    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, len_N, len_Q, len_Q, 1.0f,
                L, len_Q,
                Qc, len_Q, 0.0f,
                FLATTEN2D(L_Qc), len_Q);
    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasTrans, len_N, len_N, len_Q, 1.0f,
                FLATTEN2D(L_Qc), len_Q,
                L, len_Q, 0.0f,
                FLATTEN2D(L_Qc_LT), len_N);
    for(i=0; i<len_N; i++){
        for(j=0; j<len_N; j++){
            Phi[i][j] = F[i*len_N+j];
            Phi[i][j+len_N] = L_Qc_LT[i][j];
            Phi[i+len_N][j+len_N] = -F[j*len_N+i];
        }
    }
    utility_svsmul(FLATTEN2D(Phi), &dt, (len_N*2)*(len_N*2), NULL);
    ZE = (float**)calloc2d(len_N*2, len_N, sizeof(float));
    for(i=0; i<len_N; i++)
        ZE[i+len_N][i] = 1.0f;
    B = (float**)malloc2d(len_N*2, len_N*2, sizeof(float));
    AB = (float**)malloc2d(len_N*2, len_N, sizeof(float));
    gexpm(FLATTEN2D(Phi), len_N*2, 0, FLATTEN2D(B)); 
    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, len_N*2, len_N, len_N*2, 1.0f,
                FLATTEN2D(B), len_N*2,
                FLATTEN2D(ZE), len_N, 0.0f,
                FLATTEN2D(AB), len_N);
    AB1_T = (float**)malloc2d(len_N, len_N, sizeof(float));
    AB2_T = (float**)malloc2d(len_N, len_N, sizeof(float));
    Q_T = (float**)malloc2d(len_N, len_N, sizeof(float));
    for(i=0; i<len_N; i++){
        for(j=0; j<len_N; j++){
            AB1_T[j][i] = AB[i][j];
            AB2_T[j][i] = AB[i+len_N][j];
        }
    }
    utility_sglslv(FLATTEN2D(AB2_T), len_N, FLATTEN2D(AB1_T), len_N, FLATTEN2D(Q_T));

    /* transpose back */
    for(i=0; i<len_N; i++)
        for(j=0; j<len_N; j++)
            Q[i*len_N+j] = Q_T[j][i];
 
    /* clean-up */
    if(opt_L==NULL)
        free(L);
    if(opt_Qc==NULL)
        free(Qc);
    free(Fdt);
    free(L_Qc);
    free(L_Qc_LT);
    free(Phi);
    free(ZE);
    free(B);
    free(AB);
    free(AB1_T);
    free(AB2_T);
    free(Q_T);
}

/** GAUSS_PDF  Multivariate Gaussian PDF
%
% Syntax:
%   P = GAUSS_PDF(X,M,S)
%
% Author:
%   Simo Särkkä, 2002
%
% In:
%   X - Dx1 value or N values as DxN matrix
%   M - Dx1 mean of distibution or N values as DxN matrix.
%   S - DxD covariance matrix
%
% Out:
%   P - Probability of X.
%
% Description:
%   Calculate values of PDF (Probability Density
%   Function) of multivariate Gaussian distribution
%
%    N(X | M, S)
%
%   Function returns probability of X in PDF. If multiple
%   X's or M's are given (as multiple columns), function
%   returns probabilities for each of them. X's and M's are
%   repeated to match each other, S must be the same for all.
%
% See also:
%   GAUSS_RND

% History:
%   14.05.2003  Returns also the energy
%   20.11.2002  The first official version.
%
% Copyright (C) 2002 Simo Särkkä
%
% $Id: gauss_pdf.m,v 1.1.1.1 2003/09/15 10:54:34 ssarkka Exp $
%
% This software is distributed under the GNU General Public
% Licence (version 2 or later); please refer to the file
% Licence.txt, included with the software, for details.
 */
/* hard-coded for length(M)=3 ... */
float gauss_pdf3(float X[3], float M[3], float S[3][3])
{
    float E;
    float DX


    DX = X-repmat(M,1,size(X,2));
    E = 0.5*sum(DX.*(S\DX),1);
    d = size(M,1);
    E = E + 0.5 * d * log(2*pi) + 0.5 * log(det(S));
    P = exp(-E);
}


/* ========================================================================== */
/*                              Static Functions                              */
/* ========================================================================== */

static double lngamma(double x, double* sgngam)
{
    double a, b, c, p, q, u, w, z, logpi, ls2pi, tmp, result;
    int i;

    *sgngam = 0;
    *sgngam = (double)(1);
    logpi = 1.14472988584940017414;
    ls2pi = 0.91893853320467274178;
    if( x<-34.0 ) {
        q = -x;
        w = lngamma(q, &tmp);
        p = floor(q);
        i = floor(p+0.5);
        if( i%2==0 )
            *sgngam = (double)(-1);
        else
            *sgngam = (double)(1);
        z = q-p;
        if( z>0.5 ) {
            p = p+1.0;
            z = p-q;
        }
        z = q*sin(SAF_PId*z);
        result = logpi-log(z)-w;
        return result;
    }
    if( x<13.0 ) {
        z = (double)(1);
        p = (double)(0);
        u = x;
        while(u>=(double)3) {
            p = p-1.0;
            u = x+p;
            z = z*u;
        }
        while(u<(double)2) {
            z = z/u;
            p = p+1;
            u = x+p;
        }
        if( z<(double)0 ) {
            *sgngam = (double)(-1);
            z = -z;
        }
        else
            *sgngam = (double)(1);
        if( u==(double)(2) ) {
            result = log(z);
            return result;
        }
        p = p-2.0;
        x = x+p;
        b = -1378.25152569120859100;
        b = -38801.6315134637840924+x*b;
        b = -331612.992738871184744+x*b;
        b = -1162370.97492762307383+x*b;
        b = -1721737.00820839662146+x*b;
        b = -853555.664245765465627+x*b;
        c = (double)(1);
        c = -351.815701436523470549+x*c;
        c = -17064.2106651881159223+x*c;
        c = -220528.590553854454839+x*c;
        c = -1139334.44367982507207+x*c;
        c = -2532523.07177582951285+x*c;
        c = -2018891.41433532773231+x*c;
        p = x*b/c;
        result = log(z)+p;
        return result;
    }
    q = (x-0.5)*log(x)-x+ls2pi;
    if( x>(double)(100000000) ) {
        result = q;
        return result;
    }
    p = 1.0/(x*x);
    if( x>=1000.0 )
        q = q+((7.9365079365079365079365*0.0001*p-2.7777777777777777777778*0.001)*p+0.0833333333333333333333)/x;
    else {
        a = 8.11614167470508450300*0.0001;
        a = -5.95061904284301438324*0.0001+p*a;
        a = 7.93650340457716943945*0.0001+p*a;
        a = -2.77777777730099687205*0.001+p*a;
        a = 8.33333333333331927722*0.01+p*a;
        q = q+a/x;
    }
    result = q;
    return result;
}

static double incompletegammac(double a, double x)
{
    double igammaepsilon, igammabignumber, igammabignumberinv, result;
    double ans, ax, c, yc, r, t, y, z, pk, pkm1, pkm2, qk, qkm1, qkm2, tmp;

    igammaepsilon = 0.000000000000001;
    igammabignumber = 4503599627370496.0;
    igammabignumberinv = 2.22044604925031308085*0.0000000000000001;
    if((x<=(double)(0))|| (a<=(double)(0)) ) {
        result = (double)(1);
        return result;
    }
    if( (x<(double)(1))||(x<a) ) {
        result = 1.0-incompletegamma(a, x);
        return result;
    }
    ax = a*log(x)-x-lngamma(a, &tmp);
    if( (ax<-709.78271289338399) ) {
        result = (double)(0);
        return result;
    }
    ax = exp(ax);
    y = 1.0-a;
    z = x+y+1.0;
    c = (double)(0);
    pkm2 = (double)(1);
    qkm2 = x;
    pkm1 = x+1;
    qkm1 = z*x;
    ans = pkm1/qkm1;
    do {
        c = c+1.0;
        y = y+1.0;
        z = z+2.0;
        yc = y*c;
        pk = pkm1*z-pkm2*yc;
        qk = qkm1*z-qkm2*yc;
        if( !(qk==(double)(0)) ) {
            r = pk/qk;
            t = fabs((ans-r)/r);
            ans = r;
        }
        else
            t = (double)(1);
        pkm2 = pkm1;
        pkm1 = pk;
        qkm2 = qkm1;
        qkm1 = qk;
        if( (fabs(pk)>igammabignumber) ) {
            pkm2 = pkm2*igammabignumberinv;
            pkm1 = pkm1*igammabignumberinv;
            qkm2 = qkm2*igammabignumberinv;
            qkm1 = qkm1*igammabignumberinv;
        }
    }
    while((t>igammaepsilon));
    result = ans*ax;
    return result;
}

static double incompletegamma(double a, double x)
{
    double igammaepsilon, ans, ax, c, r, tmp, result;

    igammaepsilon = 0.000000000000001;
    if( (x<=(double)(0))||(a<=(double)(0)) ) {
        result = (double)(0);
        return result;
    }
    if( (x>(double)(1))&&(x>a) ) {
        result = 1.0-incompletegammac(a, x);
        return result;
    }
    ax = a*log(x)-x-lngamma(a, &tmp);
    if( (ax<-709.78271289338399) ) {
        result = (double)(0);
        return result;
    }
    ax = exp(ax);
    r = a;
    c = (double)(1);
    ans = (double)(1);
    do {
        r = r+1.0;
        c = c*x/r;
        ans = ans+c;
    }
    while((c/ans>igammaepsilon));
    result = ans*ax/a;
    return result;
}
