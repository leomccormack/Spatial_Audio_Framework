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
 * @file saf_utility_geometry.c
 * @ingroup Utilities
 * @brief A collection of computational geometry related functions
 *
 * @author Leo McCormack
 * @date 03.07.2020
 * @license ISC
 */

#include "saf_utilities.h"
#include "saf_externals.h" 

/** Helper function for euler2rotationMatrix() */
static void getRx
(
    float theta_rad,
    float Rx[3][3]
)
{
    Rx[0][0] = 1.0f;
    Rx[0][1] = 0.0f;
    Rx[0][2] = 0.0f;
    Rx[1][0] = 0.0f;
    Rx[1][1] = cosf(theta_rad);
    Rx[1][2] = sinf(theta_rad);
    Rx[2][0] = 0.0f;
    Rx[2][1] = -sinf(theta_rad);
    Rx[2][2] = cosf(theta_rad);
}

/** Helper function for euler2rotationMatrix() */
static void getRy
(
    float theta_rad,
    float Ry[3][3]
)
{
    Ry[0][0] = cosf(theta_rad);
    Ry[0][1] = 0.0f;
    Ry[0][2] = -sinf(theta_rad);
    Ry[1][0] = 0.0f;
    Ry[1][1] = 1.0f;
    Ry[1][2] = 0.0f;
    Ry[2][0] = sinf(theta_rad);
    Ry[2][1] = 0.0f;
    Ry[2][2] = cosf(theta_rad);
}

/** Helper function for euler2rotationMatrix() */
static void getRz
(
    float theta_rad,
    float Rz[3][3]
)
{
    Rz[0][0] = cosf(theta_rad);
    Rz[0][1] = sinf(theta_rad);
    Rz[0][2] = 0.0f;
    Rz[1][0] = -sinf(theta_rad);
    Rz[1][1] = cosf(theta_rad);
    Rz[1][2] = 0.0f;
    Rz[2][0] = 0.0f;
    Rz[2][1] = 0.0f;
    Rz[2][2] = 1.0f;
}


/* ========================================================================== */
/*                        Basic Geometrical Functions                         */
/* ========================================================================== */

void quaternion2rotationMatrix
(
    quaternion_data* Q,
    float R[3][3]
)
{
    R[0][0] = 2.0f * (Q->w * Q->w + Q->z * Q->z) - 1.0f;
    R[0][1] = 2.0f * (Q->z * Q->y - Q->w * Q->x);
    R[0][2] = 2.0f * (Q->z * Q->x + Q->w * Q->y);
    R[1][0] = 2.0f * (Q->z * Q->y + Q->w * Q->x);
    R[1][1] = 2.0f * (Q->w * Q->w + Q->y * Q->y) - 1.0f;
    R[1][2] = 2.0f * (Q->y * Q->x - Q->w * Q->z);
    R[2][0] = 2.0f * (Q->z * Q->x - Q->w * Q->y);
    R[2][1] = 2.0f * (Q->y * Q->x + Q->w * Q->z);
    R[2][2] = 2.0f * (Q->w * Q->w + Q->x * Q->x) - 1.0f; 
}

/* Adapted from: https://www.euclideanspace.com/maths/geometry/rotations/conversions/matrixToQuaternion/forum.htm */
void rotationMatrix2quaternion
(
    float R[3][3],
    quaternion_data* Q
)
{
    Q->w = sqrtf( SAF_MAX( 0.0f, 1.0f + R[0][0] + R[1][1] + R[2][2] ) ) / 2.0f;
    Q->z = sqrtf( SAF_MAX( 0.0f, 1.0f + R[0][0] - R[1][1] - R[2][2] ) ) / 2.0f;
    Q->y = sqrtf( SAF_MAX( 0.0f, 1.0f - R[0][0] + R[1][1] - R[2][2] ) ) / 2.0f;
    Q->x = sqrtf( SAF_MAX( 0.0f, 1.0f - R[0][0] - R[1][1] + R[2][2] ) ) / 2.0f;
    Q->z = copysignf( Q->z, R[2][1] - R[1][2] );
    Q->y = copysignf( Q->y, R[0][2] - R[2][0] );
    Q->x = copysignf( Q->x, R[1][0] - R[0][1] );
}

/* Adapted from (ISC License): https://github.com/MartinWeigel/Quaternion  */
void euler2Quaternion
(
    float alpha,
    float beta,
    float gamma,
    int degreesFlag,
    EULER_ROTATION_CONVENTIONS convention,
    quaternion_data* Q
)
{
    float cy, sy, cr, sr, cp, sp;

    cy = sy = cr = sr = cp = sp = 0.0f; 
    switch(convention){
        case EULER_ROTATION_Y_CONVENTION: /* fall through*/
        case EULER_ROTATION_X_CONVENTION: saf_print_error("This convention is not supported"); return;
        case EULER_ROTATION_YAW_PITCH_ROLL:
            cy = cosf((degreesFlag ? alpha*SAF_PI/180.0f : alpha)  * 0.5f); /* x */
            sy = sinf((degreesFlag ? alpha*SAF_PI/180.0f : alpha)  * 0.5f); /* x */
            cp = cosf((degreesFlag ? beta*SAF_PI/180.0f : beta) * 0.5f);    /* y */
            sp = sinf((degreesFlag ? beta*SAF_PI/180.0f : beta) * 0.5f);    /* y */
            cr = cosf((degreesFlag ? gamma*SAF_PI/180.0f : gamma)  * 0.5f); /* z */
            sr = sinf((degreesFlag ? gamma*SAF_PI/180.0f : gamma)  * 0.5f); /* z */
            break;
        case EULER_ROTATION_ROLL_PITCH_YAW:
            cy = cosf((degreesFlag ? gamma*SAF_PI/180.0f : gamma)  * 0.5f); /* x */
            sy = sinf((degreesFlag ? gamma*SAF_PI/180.0f : gamma)  * 0.5f); /* x */
            cp = cosf((degreesFlag ? beta*SAF_PI/180.0f : beta) * 0.5f);    /* y */
            sp = sinf((degreesFlag ? beta*SAF_PI/180.0f : beta) * 0.5f);    /* y */
            cr = cosf((degreesFlag ? alpha*SAF_PI/180.0f : alpha)  * 0.5f); /* z */
            sr = sinf((degreesFlag ? alpha*SAF_PI/180.0f : alpha)  * 0.5f); /* z */
            break;
    }
    Q->w = cy * cr * cp + sy * sr * sp;
    Q->x = cy * sr * cp - sy * cr * sp;
    Q->y = cy * cr * sp + sy * sr * cp;
    Q->z = sy * cr * cp - cy * sr * sp;
}

/* Adapted from (ISC License): https://github.com/MartinWeigel/Quaternion  */
void quaternion2euler
(
    quaternion_data* Q,
    int degreesFlag,
    EULER_ROTATION_CONVENTIONS convention,
    float* alpha,
    float* beta,
    float* gamma
)
{
    float sinr_cosp, cosr_cosp, sinp, siny_cosp, cosy_cosp;

    sinr_cosp = 2.0f * (Q->w * Q->x + Q->y * Q->z);
    cosr_cosp = 1.0f - 2.0f * (Q->x * Q->x + Q->y * Q->y);
    sinp = 2.0f * (Q->w * Q->y - Q->z * Q->x);
    siny_cosp = 2.0f * (Q->w * Q->z + Q->x * Q->y);
    cosy_cosp = 1.0f - 2.0f * (Q->y * Q->y + Q->z * Q->z);
    switch(convention){
        case EULER_ROTATION_Y_CONVENTION: /* fall through */
        case EULER_ROTATION_X_CONVENTION: saf_print_error("This convention is not supported"); break;
        case EULER_ROTATION_YAW_PITCH_ROLL:
            /* Yaw (z-axis rotation) */
            (*gamma) = atan2f(sinr_cosp, cosr_cosp);
            /* Pitch (y-axis rotation) */
            if (fabsf(sinp) >= 1.0f)
                (*beta) = copysignf(SAF_PI / 2.0f, sinp); /* use 90 degrees if out of range */
            else
                (*beta) = asinf(sinp);
            /* Roll (x-axis rotation) */
           (*alpha) = atan2f(siny_cosp, cosy_cosp);
            break;
        case EULER_ROTATION_ROLL_PITCH_YAW:
            /* Roll (x-axis rotation) */
            (*alpha) = atan2f(sinr_cosp, cosr_cosp);
            /* Pitch (y-axis rotation) */
            if (fabs(sinp) >= 1.0f)
                (*beta) = copysignf(SAF_PI / 2.0f, sinp); /* use 90 degrees if out of range */
            else
                (*beta) = asinf(sinp);
            /* Yaw (z-axis rotation) */
            (*gamma) = atan2f(siny_cosp, cosy_cosp); 
            break;
    }
    if(degreesFlag){
        (*alpha) *= 180.0f/SAF_PI;
        (*beta)  *= 180.0f/SAF_PI;
        (*gamma) *= 180.0f/SAF_PI;
    }
}

void euler2rotationMatrix
(
    float alpha,
    float beta,
    float gamma,
    int degreesFlag,
    EULER_ROTATION_CONVENTIONS convention,
    float R[3][3]
)
{
    float R1[3][3], R2[3][3], R3[3][3], Rtmp[3][3];

    switch(convention){
        case EULER_ROTATION_Y_CONVENTION:
            getRz(degreesFlag ? alpha*SAF_PI/180.0f : alpha, R1);
            getRy(degreesFlag ? beta*SAF_PI/180.0f  : beta,  R2);
            getRz(degreesFlag ? gamma*SAF_PI/180.0f : gamma, R3);
            break;
        case EULER_ROTATION_X_CONVENTION:
            getRz(degreesFlag ? alpha*SAF_PI/180.0f : alpha, R1);
            getRx(degreesFlag ? beta*SAF_PI/180.0f  : beta,  R2);
            getRz(degreesFlag ? gamma*SAF_PI/180.0f : gamma, R3);
            break;
        case EULER_ROTATION_YAW_PITCH_ROLL:
            getRz(degreesFlag ? alpha*SAF_PI/180.0f : alpha, R1);
            getRy(degreesFlag ? beta*SAF_PI/180.0f  : beta,  R2);
            getRx(degreesFlag ? gamma*SAF_PI/180.0f : gamma, R3);
            break;
        case EULER_ROTATION_ROLL_PITCH_YAW:
            getRx(degreesFlag ? alpha*SAF_PI/180.0f : alpha, R1);
            getRy(degreesFlag ? beta*SAF_PI/180.0f  : beta,  R2);
            getRz(degreesFlag ? gamma*SAF_PI/180.0f : gamma, R3);
            break;
    }
    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 3, 3, 3, 1.0f,
                (float*)R2, 3,
                (float*)R1, 3, 0.0f,
                (float*)Rtmp, 3);
    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 3, 3, 3, 1.0f,
                (float*)R3, 3,
                (float*)Rtmp, 3, 0.0f,
                (float*)R, 3);
}

void yawPitchRoll2Rzyx
(
    float yaw,
    float pitch,
    float roll,
    int rollPitchYawFLAG,
    float R[3][3]
)
{
    if(rollPitchYawFLAG)
        euler2rotationMatrix(yaw, pitch, roll, 0, EULER_ROTATION_ROLL_PITCH_YAW, R);
    else
        euler2rotationMatrix(yaw, pitch, roll, 0, EULER_ROTATION_YAW_PITCH_ROLL, R);
}

void sph2cart(float* sph,
              int nDirs,
              int anglesInDegreesFLAG,
              float* cart)
{
    int i;
    float tmp_rad[2];

    if(anglesInDegreesFLAG){
        for(i=0; i<nDirs; i++){
            tmp_rad[0] = sph[i*3]   * SAF_PI/180.0f;
            tmp_rad[1] = sph[i*3+1] * SAF_PI/180.0f;
            cart[i*3]   = sph[i*3+2] * cosf(tmp_rad[1]) * cosf(tmp_rad[0]);
            cart[i*3+1] = sph[i*3+2] * cosf(tmp_rad[1]) * sinf(tmp_rad[0]);
            cart[i*3+2] = sph[i*3+2] * sinf(tmp_rad[1]);
        }
    }
    else { /* Angles given in radians */
        for(i=0; i<nDirs; i++){
            cart[i*3]   = sph[i*3+2] * cosf(sph[i*3+1]) * cosf(sph[i*3]);
            cart[i*3+1] = sph[i*3+2] * cosf(sph[i*3+1]) * sinf(sph[i*3]);
            cart[i*3+2] = sph[i*3+2] * sinf(sph[i*3+1]);
        }
    } 
}

void cart2sph(float* cart,
              int nDirs,
              int anglesInDegreesFLAG,
              float* sph)
{
    int i;
    float hypotxy;

    for(i=0; i<nDirs; i++){
        hypotxy = sqrtf(cart[i*3]*cart[i*3] + cart[i*3+1]*cart[i*3+1]);
        sph[i*3]   = atan2f(cart[i*3+1], cart[i*3]);
        sph[i*3+1] = atan2f(cart[i*3+2], hypotxy);
        sph[i*3+2] = L2_norm3(&cart[i*3]);
    }

    /* Return in degrees instead... */
    if(anglesInDegreesFLAG){
        for(i=0; i<nDirs; i++){
            sph[i*3] *= (180.0f/SAF_PI);
            sph[i*3+1] *= (180.0f/SAF_PI);
        }
    }
}

void unitSph2cart
(
    float* dirs,
    int nDirs,
    int anglesInDegreesFLAG,
    float* dirs_xyz
)
{
    int i;
    float tmp_rad[2];

    if(anglesInDegreesFLAG){
        for(i=0; i<nDirs; i++){
            tmp_rad[0] = dirs[i*2]   * SAF_PI/180.0f;
            tmp_rad[1] = dirs[i*2+1] * SAF_PI/180.0f;
            dirs_xyz[i*3]   = cosf(tmp_rad[1]) * cosf(tmp_rad[0]);
            dirs_xyz[i*3+1] = cosf(tmp_rad[1]) * sinf(tmp_rad[0]);
            dirs_xyz[i*3+2] = sinf(tmp_rad[1]);
        }
    }
    else { /* Angles given in radians */
        for(i=0; i<nDirs; i++){
            dirs_xyz[i*3]   = cosf(dirs[i*2+1]) * cosf(dirs[i*2]);
            dirs_xyz[i*3+1] = cosf(dirs[i*2+1]) * sinf(dirs[i*2]);
            dirs_xyz[i*3+2] = sinf(dirs[i*2+1]);
        }
    } 
}

void unitCart2sph
(
    float* dirs_xyz,
    int nDirs,
    int anglesInDegreesFLAG,
    float* dirs
)
{
    int i;

    for(i=0; i<nDirs; i++){
        dirs[i*2]   = atan2f(dirs_xyz[i*3+1], dirs_xyz[i*3]);
        dirs[i*2+1] = atan2f(dirs_xyz[i*3+2], sqrtf(dirs_xyz[i*3]*dirs_xyz[i*3] + dirs_xyz[i*3+1]*dirs_xyz[i*3+1]));
    }

    /* Return in degrees instead... */
    if(anglesInDegreesFLAG)
        for(i=0; i<nDirs*2; i++)
            dirs[i] *= (180.0f/SAF_PI);
}

float L2_norm3
(
    float v[3]
)
{
    return sqrtf(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}

float L2_norm
(
    float* v,
    int lenV
)
{
    int i;
    float res;
    res = 0.0f;
    for(i=0; i<lenV; i++)
        res += (v[i]*v[i]);
    return sqrtf(res);
}

float Frob_norm
(
    float* M,
    int lenX,
    int lenY
)
{
    int i;
    float res;
    float* MMT;
    MMT = malloc1d(lenX*lenX*sizeof(float));
    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasTrans, lenX, lenX, lenY, 1.0f,
                M, lenY,
                M, lenY, 0.0f,
                MMT, lenX);
    res = 0.0f;
    for(i=0; i<lenX; i++)
        res += MMT[i*lenX+i];
    free(MMT); 
    return sqrtf(res);
}

void crossProduct3
(
    float a[3],
    float b[3],
    float c[3]
)
{
    c[0] = a[1]*b[2]-a[2]*b[1];
    c[1] = a[2]*b[0]-a[0]*b[2];
    c[2] = a[0]*b[1]-a[1]*b[0];
}

float getDistBetweenPointAndLine
(
    float point[3],
    float v1[3],
    float v2[3]
)
{
    float a[3], b[3], cross_a_ab[3];
    a[0] = v1[0] - v2[0];
    a[1] = v1[1] - v2[1];
    a[2] = v1[2] - v2[2];
    b[0] = point[0] - v2[0];
    b[1] = point[1] - v2[1];
    b[2] = point[2] - v2[2];
    crossProduct3(a, b, cross_a_ab);
    return L2_norm3(cross_a_ab)/(L2_norm3(a)+2.3e-9f);
}

float getDistBetween2Points
(
    float point_a[3],
    float point_b[3]
)
{
#if defined(SAF_USE_APPLE_ACCELERATE)
    float dist;
    vDSP_distancesq((const float*)point_a, 1, (const float*)point_b, 1, &dist, 3);
    return dist;
#else
    float a_b[3];
    a_b[0] = point_a[0] - point_b[0];
    a_b[1] = point_a[1] - point_b[1];
    a_b[2] = point_a[2] - point_b[2];
    return L2_norm3(a_b);
#endif
}


/* ========================================================================== */
/*                     Computational Geometry Functions                       */
/* ========================================================================== */

void convhull3d
(
    const float* vertices,
    const int nVert,
    int** faces,
    int* nFaces
)
{
    int i;
    ch_vertex* ch_vertices;

    /* convert vertices to use "ch_vertex" format used by convhull_3d_build() */
    ch_vertices = malloc1d(nVert*sizeof(ch_vertex));
    for(i = 0; i < nVert; i++) {
        ch_vertices[i].z = (CH_FLOAT)vertices[i*3+2];
        ch_vertices[i].x = (CH_FLOAT)vertices[i*3];
        ch_vertices[i].y = (CH_FLOAT)vertices[i*3+1];
    }

    /* build convex hull */
    saf_assert(*faces == NULL, "nFaces not known yet, and so shouldn't be pre-allocated...");
    convhull_3d_build(ch_vertices, nVert, faces, NULL, NULL, nFaces);

    /* clean-up */
    free(ch_vertices);
}

void convhullnd
(
    const float* points,
    const int nPoints,
    const int nd,
    int** faces,
    int* nFaces
)
{
    int i, j;
    CH_FLOAT* ch_points;

    /* convert vertices to use CH_FLOAT used by convhull_nd_build() */
    ch_points = malloc1d(nPoints*nd*sizeof(CH_FLOAT));
    for(i = 0; i < nPoints; i++) {
        for(j=0; j<nd; j++)
            ch_points[i*nd+j] = (CH_FLOAT)points[i*nd+j];
    }

    /* build convex hull */
    saf_assert(*faces == NULL, "nFaces not known yet, and so shouldn't be pre-allocated...");
    convhull_nd_build(ch_points, nPoints, nd, faces, NULL, NULL, nFaces);

    /* clean-up */
    free(ch_points);
}

void delaunaynd
(
    const float* points,
    const int nPoints,
    const int nd,
    int** DT,
    int* nDT
)
{
    int i, j, k, nHullFaces, maxW_idx, nVisible;
    int* hullfaces;
    CH_FLOAT w0, w_optimal, w_optimal2;
    CH_FLOAT* projpoints, *cf, *df, *p0, *p, *visible;

    /* Project the N-dimensional points onto a N+1-dimensional paraboloid */
    projpoints = malloc1d(nPoints*(nd+1)*sizeof(CH_FLOAT));
    for(i = 0; i < nPoints; i++) {
        projpoints[i*(nd+1)+nd] = 0.0;
        for(j=0; j<nd; j++){
            projpoints[i*(nd+1)+j] = (CH_FLOAT)points[i*nd+j] + 0.0000001*(CH_FLOAT)rand()/(CH_FLOAT)RAND_MAX;  
            projpoints[i*(nd+1)+nd] += (projpoints[i*(nd+1)+j]*projpoints[i*(nd+1)+j]); /* w vector */
        }
    }

    /* The N-dimensional delaunay triangulation requires first computing the convex hull of this N+1-dimensional paraboloid */
    hullfaces = NULL;
    cf = df = NULL;
    convhull_nd_build(projpoints, nPoints, nd+1, &hullfaces, &cf, &df, &nHullFaces);

    /* Find the coordinates of the point with the maximum (N+1 dimension) coordinate (i.e. the w vector) */
    if(sizeof(CH_FLOAT)==sizeof(double))
        maxW_idx = (int)cblas_idamax(nPoints, (double*)&projpoints[nd], nd+1);
    else
        maxW_idx = (int)cblas_isamax(nPoints, (float*)&projpoints[nd], nd+1);
    w0 = projpoints[maxW_idx*(nd+1)+nd];
    p0 = malloc1d(nd*sizeof(CH_FLOAT));
    for(j=0; j<nd; j++)
        p0[j] = projpoints[maxW_idx*(nd+1)+j];

    /* Find the point where the plane tangent to the point (p0,w0) on the paraboloid crosses the w axis.
     * This is the point that can see the entire lower hull. */
    w_optimal = 0.0;
    for(j=0; j<nd; j++)
       w_optimal += (2.0*pow(p0[j], 2.0));
    w_optimal = w0-w_optimal;

    /* Subtract 1000 times the absolute value of w_optimal to ensure that the point where the tangent plane
     * crosses the w axis will see all points on the lower hull. This avoids numerical roundoff errors. */
    w_optimal2=w_optimal-1000.0*fabs(w_optimal);

    /* Set the point where the tangent plane crosses the w axis */
    p = calloc1d((nd+1),sizeof(CH_FLOAT));
    p[nd] = w_optimal2;

    /* Find all faces that are visible from this point */
    visible = malloc1d(nHullFaces*sizeof(CH_FLOAT));
    if(sizeof(CH_FLOAT)==sizeof(double)){
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nHullFaces, 1, nd+1, 1.0,
                    (double*)cf, nd+1,
                    (double*)p, 1, 0.0,
                    (double*)visible, 1);
    }
    else{
        cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nHullFaces, 1, nd+1, 1.0f,
                    (float*)cf, nd+1,
                    (float*)p, 1, 0.0f,
                    (float*)visible, 1);
    }
    nVisible = 0;
    for(j=0; j<nHullFaces; j++){
        visible[j] += df[j];
        if(visible[j]>0.0)
            nVisible++;
    }
 
    /* Output */
    (*nDT) = nVisible;
    if(nVisible>0){
        (*DT) = malloc1d(nVisible*(nd+1)*sizeof(int));
        for(i=0, j=0; i<nHullFaces; i++){
            if(visible[i]>0.0){
                for(k=0; k<nd+1; k++)
                    (*DT)[j*(nd+1)+k] = hullfaces[i*(nd+1)+k];
                j++;
            }
        }
        saf_assert(j==nVisible, "Ugly error");
    }

    /* clean up */
    free(projpoints);
    free(hullfaces);
    free(cf);
    free(df);
    free(p0);
    free(p);
    free(visible);
}

void sphDelaunay
(
    const float* dirs_deg,
    const int nDirs,
    int** faces,
    int* nFaces,
    float* vertices
)
{
    int i;
    float* vertices_tmp;
    float rcoselev;

    /* sph to unit Cart: */
    vertices_tmp = malloc1d(nDirs*3*sizeof(float));
    for(i = 0; i < nDirs; i++) {
        vertices_tmp[i*3+2] = sinf(dirs_deg[i*2+1]* SAF_PI/180.0f);
        rcoselev = cosf(dirs_deg[i*2+1]*SAF_PI/180.0f);
        vertices_tmp[i*3] = rcoselev * cosf(dirs_deg[i*2+0]*SAF_PI/180.0f);
        vertices_tmp[i*3+1] = rcoselev * sinf(dirs_deg[i*2+0]*SAF_PI/180.0f);
    }

    /* Delaunay triangulation of a spherical grid is equivalent to the computing
     * the 3d convex hull */
    convhull3d(vertices_tmp, nDirs, faces, nFaces);

    /* optionally, also output the vertices */
    if(vertices!=NULL)
        memcpy(vertices, vertices_tmp, nDirs*3*sizeof(float));

    /* clean up */
    free(vertices_tmp);
}

void sphVoronoi
(
    int* faces,
    int nFaces,
    float* vertices,
    int nDirs,
    voronoi_data* voronoi
)
{
    int n, m;
    int* duplicates;
    float r_12[3], r_13[3], r_normal[3];
    float norm;

    /* prep */
    voronoi->nVert = nFaces;
    voronoi->nFaces = nDirs;
    voronoi->vert = (float**)malloc2d(voronoi->nVert, 3, sizeof(float));
    voronoi->faces = (int**)malloc1d(voronoi->nFaces*sizeof(int*));
    voronoi->nPointsPerFace = malloc1d(voronoi->nFaces*sizeof(int));

    /* Calculate the voronoi vertices for each triangle, which for the unit
     * sphere are given by the unit normal vector of the triangle */
    for(n = 0; n<nFaces; n++){
        r_12[0] = vertices[faces[n*3+1]*3]   - vertices[faces[n*3]*3];
        r_12[1] = vertices[faces[n*3+1]*3+1] - vertices[faces[n*3]*3+1];
        r_12[2] = vertices[faces[n*3+1]*3+2] - vertices[faces[n*3]*3+2];
        r_13[0] = vertices[faces[n*3+2]*3]   - vertices[faces[n*3]*3];
        r_13[1] = vertices[faces[n*3+2]*3+1] - vertices[faces[n*3]*3+1];
        r_13[2] = vertices[faces[n*3+2]*3+2] - vertices[faces[n*3]*3+2];
        crossProduct3(r_12, r_13, r_normal);

        norm = 1.0f/L2_norm3(r_normal);
        utility_svsmul(r_normal, &norm, 3, r_normal);
        memcpy(voronoi->vert[n], r_normal, 3*sizeof(float));
    }

    /* Find duplicate vertices if any, due to two triangles sharing the same
     * circumscribed circle */
    duplicates = calloc1d(voronoi->nVert, sizeof(int));
    for(n = 0; n<voronoi->nVert; n++){
        if (duplicates[n] == 0 ){
            for (m = 0; m<voronoi->nVert; m++){
                if (n != m){
                    if (fabsf(voronoi->vert[n][0] - voronoi->vert[m][0]) < 1.0e-5f &&
                        fabsf(voronoi->vert[n][1] - voronoi->vert[m][1]) < 1.0e-5f &&
                        fabsf(voronoi->vert[n][2] - voronoi->vert[m][2]) < 1.0e-5f ){
                        duplicates[m] = n;
                    }
                }
            }
        }
    }

    int NOT_SORTED;
    int i, j, k, l, nFaceIdx, currentfaceIdx, currentvertIdx, currentvert, nSorted;
    int currentface[3];
    int* faceIdx, *sorted, *tempfacelist;
    faceIdx = NULL;
    sorted = NULL;
    tempfacelist = NULL;

    /* Calculate the voronoi polygons
     *
     * find the an ordered sequence of the triangles around each vertex and get
     * the proper sequence of the voronoi vertices that constitute a polygon */
    for (n = 0; n<voronoi->nFaces; n++){
        nFaceIdx=0;
        for(m=0; m<voronoi->nVert; m++)
            if(faces[m*3+0]==n || faces[m*3+1]==n || faces[m*3+2]==n)
                nFaceIdx++;
        faceIdx = realloc1d(faceIdx, nFaceIdx*sizeof(int));

        /* list of triangles that contain this specific vertex */
        i=0;
        for(m=0; m<voronoi->nVert; m++){
            if(faces[m*3+0]==n || faces[m*3+1]==n || faces[m*3+2]==n){
                faceIdx[i] = m;
                i++;
            }
        }

        /* Each triangle from the list contain the common vertex and two other
         * vertices - each triangle has two common vertices with each other. One
         * (brute) way of sorting the sequence is to pick one triangle, find the
         * neighbour triangle by finding their common
         * vertex, move to that triangle and iterate till all the number of
         * triangles have been checked. */
        k = 0;
        currentfaceIdx = faceIdx[k]; /* pick-up one of the triangles in the list */
        memcpy(currentface, &faces[currentfaceIdx*3], 3*sizeof(int)); /* the triangle's vertex indices */

        /* pick-up one of the vertices that is not the central one */
        currentvertIdx = -1;
        for(j=0; j<3; j++){
            if(currentface[j] != n){
                currentvertIdx = j;
                break;
            }
        }
        saf_assert(currentvertIdx!=-1, "Ugly error");
        currentvert = currentface[currentvertIdx];

        /* Prep */
        nSorted = 1;
        sorted = realloc1d(sorted, nSorted*sizeof(int));
        sorted[0] = faceIdx[k]; /* this is the list that keeps the the ordered triangles */
        NOT_SORTED = 1;
        while(NOT_SORTED){
            /* exclude the current triangle from the temporary list */
            tempfacelist = realloc1d(tempfacelist, (nFaceIdx-1)*sizeof(int));
            l=0;
            for(i=0; i<nFaceIdx; i++){
                if(faceIdx[i]!=currentfaceIdx){
                    tempfacelist[l] = faceIdx[i];
                    l++;
                }
            }
            saf_assert(l==nFaceIdx-1, "Ugly error");

            for (l = 0; l<nFaceIdx-1; l++){
                currentfaceIdx = tempfacelist[l];
                memcpy(currentface, &faces[currentfaceIdx*3], 3*sizeof(int));

                /* if the currentvert exists in the current triangles vertices... */
                if (currentface[0] == currentvert || currentface[1] == currentvert || currentface[2] == currentvert){
                    /* then it's the neighbour triangle, so store its index */
                    nSorted++;
                    sorted = realloc1d(sorted, nSorted*sizeof(int));
                    sorted[nSorted-1] = currentfaceIdx;

                    /* if the sorted list has the length of faceIdx, then we're done */
                    if (nSorted == nFaceIdx){
                        NOT_SORTED = 0;
                        break;
                    }

                    /* find the next vertex from current triangle that excludes the central one and the previous one */
                    for(j=0; j<3; j++)
                        if(currentface[j] != n && currentface[j] != currentvert)
                            currentvertIdx = j;
                    currentvert = currentface[currentvertIdx];
                    break;
                }
            }
        }

        /* remove the duplicate vertices from the list */
        for (i=0; i<nSorted; i++)
            if (duplicates[sorted[i]] != 0)
                sorted[i] = duplicates[sorted[i]];


        /* Identify unique IDs, and sort them in accending order */
        int nUnique;
        int* uniqueIDs;
        uniqueIDs = NULL;
        unique_i(sorted, nSorted, NULL, &uniqueIDs, &nUnique);
        sorti(uniqueIDs, uniqueIDs, NULL, nUnique, 0);

        /* Save */
        voronoi->faces[n] = malloc1d(nUnique*sizeof(int));
        for(i=0; i<nUnique; i++)
            voronoi->faces[n][i] = sorted[uniqueIDs[i]];
        voronoi->nPointsPerFace[n] = nUnique;

        /* clean-up */
        free(uniqueIDs);
    }

    /* clean-up */
    free(duplicates);
    free(faceIdx);
    free(sorted);
    free(tempfacelist);
}

void sphVoronoiAreas
(
    voronoi_data* voronoi,
    float* areas
)
{
    int i, m, n, N_poly, tmp_i;
    int* face;
    float tmp, r_21_norm, r_23_norm;
    float *theta;
    float r_01[3], r_02[3], r_2x1[3], r_21[3], r_03[3], r_2x3[3], r_23[3];

    face = NULL;
    theta = NULL;
    for(m=0; m<voronoi->nFaces; m++){
        N_poly = voronoi->nPointsPerFace[m]; /* number of vertices in the polygon */
        face = realloc1d(face, N_poly*sizeof(int));
        theta = realloc1d(theta, N_poly*sizeof(float));
        memcpy(face, voronoi->faces[m], N_poly*sizeof(int)); /* current face */
 
        for(n=0; n<N_poly; n++){
            /* find vector between vertex origin and vertices 1 & 2 */
            memcpy(r_01, voronoi->vert[face[0]], 3*sizeof(float));
            memcpy(r_02, voronoi->vert[face[1]], 3*sizeof(float));

            /* find normal vector to the great circle of 1 & 2 */
            crossProduct3(r_02, r_01, r_2x1);

            /* find tangent vector to great circle at 2 */
            crossProduct3(r_2x1, r_02, r_21);

            /* find vector between vertex origin and vertex 3 */
            memcpy(r_03, voronoi->vert[face[2]], 3*sizeof(float));

            /* find normal vector to the great circle of 2 & 3 */
            crossProduct3(r_02, r_03, r_2x3);

            /* find tangent vector to great circle at 2 */
            crossProduct3(r_2x3, r_02, r_23);

            /* normalise tangent vectors */
            r_21_norm = 1.0f/L2_norm3(r_21);
            utility_svsmul(r_21, &r_21_norm, 3, r_21);
            r_23_norm = 1.0f/L2_norm3(r_23);
            utility_svsmul(r_23, &r_23_norm, 3, r_23);

            /* get angle between the normals */
            utility_svvdot(r_21, r_23, 3, &tmp);
            theta[n] = acosf(tmp);

            /* shift the vertex list by one position and repeat */
            tmp_i = face[0];
            for(i=1; i<N_poly; i++)
                face[i-1] = face[i];
            face[N_poly-1] = tmp_i;
        }

        /* compute area of spherical polygon by spherical excess */
        tmp = 0.0f;
        for(i=0; i<N_poly; i++)
            tmp += theta[i];
        areas[m] = tmp - ((float)N_poly-2.0f)*SAF_PI; 
    }
    free(face);
    free(theta);
}

void getVoronoiWeights
(
    float* dirs_deg,
    int nDirs,
    int diagFLAG,
    float* weights
)
{
    int i, nFaces;
    int* faces;
    float* vertices, *areas;
    voronoi_data voronoi;

    /* Perform delaunay triangulation */
    faces = NULL;
    vertices = malloc1d(nDirs*3*sizeof(float));
    sphDelaunay(dirs_deg, nDirs, &faces, &nFaces, vertices);

    /* Get voronoi diagrams */
    sphVoronoi(faces, nFaces, vertices, nDirs, &voronoi);

    /* Compute areas of spherical voronoi polygons */
    areas = malloc1d(voronoi.nFaces * sizeof(float));
    sphVoronoiAreas(&voronoi, areas);

    /* Store weights */
    if(diagFLAG){
        /* along the diagonal... */
        memset(weights, 0, nDirs*nDirs*sizeof(float));
        for(i=0; i<nDirs; i++)
            weights[i*nDirs+i] = areas[i];
    }
    else
        memcpy(weights, areas, nDirs*sizeof(float));

    /* clean-up */
    free(faces);
    free(vertices);
    free(areas);
    for(i=0; i<voronoi.nFaces; i++)
        free(voronoi.faces[i]);
    free(voronoi.faces);
    free(voronoi.vert);
    free(voronoi.nPointsPerFace);
}
