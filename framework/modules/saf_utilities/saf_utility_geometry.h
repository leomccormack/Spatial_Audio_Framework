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
 *@addtogroup Utilities
 *@{
 * @file saf_utility_geometry.h
 * @brief A collection of computational geometry related functions
 *
 * @author Leo McCormack
 * @date 03.07.2020
 * @license ISC
 */

#ifndef SAF_GEOMETRY_H_INCLUDED
#define SAF_GEOMETRY_H_INCLUDED

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/** Quaternion data structure */
typedef struct _quaternion_data {
    union {
        struct { float w; /**< W value of the quaternion [-1..1] */
                 float x; /**< X value of the quaternion [-1..1] */
                 float y; /**< Y value of the quaternion [-1..1] */
                 float z; /**< Z value of the quaternion [-1..1] */
        };
        float Q[4]; /**< WXYZ values of the quaternion [-1..1] */
    };
} quaternion_data;

/** Available euler2rotationMatrix() conventions */
typedef enum {
    EULER_ROTATION_Y_CONVENTION,   /**< y-convention, 'zyz' */
    EULER_ROTATION_X_CONVENTION,   /**< x-convention, 'zxz' */
    EULER_ROTATION_YAW_PITCH_ROLL, /**< yaw-pitch-roll, 'zyx' */
    EULER_ROTATION_ROLL_PITCH_YAW  /**< roll-pitch-yaw, 'xyz' */

} EULER_ROTATION_CONVENTIONS;

/**
 * Data structure for Voronoi diagrams
 *
 * @warning 'faces' is NOT contiguously allocated! [i.e., you cannot do
 *          free(faces), you must loop over the i=0:nFaces-1: free(faces[i]) ]
 */
typedef struct _voronoi_data{
    int nVert;           /**< Number of vertices */
    int nFaces;          /**< Number of faces/polygons */
    float** vert;        /**< Vertices; nVert x 3 */
    int** faces;         /**< faces; nFaces x nPointsPerFace[i]  */
    int* nPointsPerFace; /**< Number of points for each face; nFaces x 1 */

}voronoi_data;


/* ========================================================================== */
/*                        Basic Geometrical Functions                         */
/* ========================================================================== */

/** Constructs a 3x3 rotation matrix based on a quaternion */
void quaternion2rotationMatrix(/* Input Arguments */
                               quaternion_data* Q,
                               /* Output Arguments */
                               float R[3][3]);

/** Calculates the quaternion corresponding to a 3x3 rotation matrix */
void rotationMatrix2quaternion(/* Input Arguments */
                               float R[3][3],
                               /* Output Arguments */
                               quaternion_data* Q);

/** Converts Euler angles to a quaternion */
void euler2Quaternion(/* Input Arguments */
                      float alpha,
                      float beta,
                      float gamma,
                      int degreesFlag,
                      EULER_ROTATION_CONVENTIONS convention,
                      /* Output Arguments */
                      quaternion_data* Q);

/** Converts a quaternion to Euler angles */
void quaternion2euler(/* Input Arguments */
                      quaternion_data* Q,
                      int degreesFlag,
                      EULER_ROTATION_CONVENTIONS convention,
                      /* Output Arguments */
                      float* alpha,
                      float* beta,
                      float* gamma);

/**
 * Constructs a 3x3 rotation matrix from the Euler angles
 *
 * @param[in]  alpha       first rotation angle
 * @param[in]  beta        first rotation angle
 * @param[in]  gamma       first rotation angle
 * @param[in]  degreesFlag 1: angles are in degrees, 0: angles are in radians
 * @param[in]  convention  see #EULER_ROTATION_CONVENTIONS enum
 * @param[out] R           resulting 3x3 rotation matrix
 */
void euler2rotationMatrix(/* Input Arguments */
                          float alpha,
                          float beta,
                          float gamma,
                          int degreesFlag,
                          EULER_ROTATION_CONVENTIONS convention,
                          /* Output Arguments */
                          float R[3][3]);

/**
 * Constructs a 3x3 rotation matrix from the Euler angles, using the
 * yaw-pitch-roll (zyx) convention
 *
 * @note DEPRICATED. This function now just calls: euler2rotationMatrix()
 *
 * @param[in]  yaw              Yaw angle in radians
 * @param[in]  pitch            Pitch angle in radians
 * @param[in]  roll             Roll angle in radians
 * @param[in]  rollPitchYawFLAG '1' to use Rxyz, i.e. apply roll, pitch and then
 *                              yaw, '0' Rzyx / y-p-r
 * @param[out] R                zyx rotation matrix; 3 x 3
 */
void yawPitchRoll2Rzyx (/* Input Arguments */
                        float yaw,
                        float pitch,
                        float roll,
                        int rollPitchYawFLAG,
                        /* Output Arguments */
                        float R[3][3]);

/**
 * Converts spherical coordinates to Cartesian coordinates
 *
 * @param[in]  sph                 Spherical coordinates; FLAT: nDirs x 3
 * @param[in]  nDirs               Number of directions/coordinates
 * @param[in]  anglesInDegreesFLAG 0: dirs given in radians, 1: degrees instead
 * @param[out] cart                Cartesian coordinates; FLAT: nDirs x 3
 */
void sph2cart(/* Input Arguments */
              float* sph,
              int nDirs,
              int anglesInDegreesFLAG,
              /* Output Arguments */
              float* cart);

/**
 * Converts Cartesian coordinates to spherical coordinates
 *
 * @param[in]  cart                Cartesian coordinates; FLAT: nDirs x 3
 * @param[in]  nDirs               Number of directions/coordinates
 * @param[in]  anglesInDegreesFLAG 0: dirs wanted in radians, 1: degrees instead
 * @param[out] sph                 Spherical coordinates; FLAT: nDirs x 3
 */
void cart2sph(/* Input Arguments */
              float* cart,
              int nDirs,
              int anglesInDegreesFLAG,
              /* Output Arguments */
              float* sph);

/**
 * Converts spherical coordinates to Cartesian coordinates of unit length
 *
 * @param[in]  dirs                Spherical coordinates; FLAT: nDirs x 2
 * @param[in]  nDirs               Number of directions/coordinates
 * @param[in]  anglesInDegreesFLAG 0: dirs given in radians, 1: degrees instead
 * @param[out] dirs_xyz            Cartesian coordinates; FLAT: nDirs x 3
 */
void unitSph2cart(/* Input Arguments */
                  float* dirs,
                  int nDirs,
                  int anglesInDegreesFLAG,
                  /* Output Arguments */
                  float* dirs_xyz);

/**
 * Converts Cartesian coordinates of unit length to spherical coordinates
 *
 * @param[in]  dirs_xyz            Cartesian coordinates; FLAT: nDirs x 3
 * @param[in]  nDirs               Number of directions/coordinates
 * @param[in]  anglesInDegreesFLAG 0: dirs wanted in radians, 1: degrees instead
 * @param[out] dirs                Spherical coordinates; FLAT: nDirs x 2
 */
void unitCart2sph(/* Input Arguments */
                  float* dirs_xyz,
                  int nDirs,
                  int anglesInDegreesFLAG,
                  /* Output Arguments */
                  float* dirs);

/** Returns the L2 (Euclidean) norm of a 3-element vector */
float L2_norm3(float v[3]);

/** Returns the L2 (Euclidean) norm of an arbitrary length vector */
float L2_norm(float* v, int lenV);

/** Returns the Frobenius Norm of a matrix M, of dimensions: lenX x lenY */
float Frob_norm(float* M, int lenX, int lenY);

/** Cross product between two 3-element floating point vectors (c = a x b) */
void crossProduct3(/* Input Arguments */
                   float a[3],
                   float b[3],
                   /* Output Arguments */
                   float c[3]);

/**
 * Returns the distance between a "point" and an infinite line described by the
 * two points "v1" and "v2"
 */
float getDistBetweenPointAndLine(/* Input Arguments */
                                 float point[3],
                                 float v1[3],
                                 float v2[3]);

/** Returns the distance between "point_a" and "point_b" */
float getDistBetween2Points(/* Input Arguments */
                            float point_a[3],
                            float point_b[3]);


/* ========================================================================== */
/*                     Computational Geometry Functions                       */
/* ========================================================================== */

/**
 * Builds the convex hull of an arrangement of vertices in 3-dimensional space
 *
 * This function employs algorithms originally implemented in MATLAB by George
 * Papazafeiropoulos [1] (BSD 2-clause license), which are based on [2].
 *
 * @warning Currently, this does not check if there are duplicate vertices or
 *          whether any of them are co-linear!
 *
 * @param[in]  vertices The vertices; FLAT: nDirs x 3
 * @param[in]  nVert    Number of vertices
 * @param[out] faces    (&) The face indices; FLAT: nFaces x 3
 * @param[out] nFaces   (&) Number of faces found
 *
 * @see [1] https://www.mathworks.com/matlabcentral/fileexchange/48509-computational-geometry-toolbox?focused=3851550&tab=example
 * @see [2] The Quickhull Algorithm for Convex Hull, C. Bradford Barber, David
 *          P. Dobkin and Hannu Huhdanpaa, Geometry Center Technical Report
 *          GCG53, July 30, 1993
 */
void convhull3d(/* Input Arguments */
                const float* vertices,
                const int nVert,
                /* Output Arguments */
                int** faces,
                int* nFaces);

/**
 * Builds the convex hull of an arrangement of points in N-dimensional space
 *
 * This function employs algorithms originally implemented in MATLAB by George
 * Papazafeiropoulos [1] (BSD 2-clause license), which are based on [2].
 *
 * @param[in]  points  The input points; FLAT: nDirs x nd
 * @param[in]  nPoints Number of points
 * @param[in]  nd      The number of dimensions (max=5)
 * @param[out] faces   (&) The face indices; FLAT: nFaces x nd
 * @param[out] nFaces  (&) Number of faces found
 *
 * @see [1] https://www.mathworks.com/matlabcentral/fileexchange/48509-computational-geometry-toolbox?focused=3851550&tab=example
 * @see [2] The Quickhull Algorithm for Convex Hull, C. Bradford Barber, David
 *          P. Dobkin and Hannu Huhdanpaa, Geometry Center Technical Report
 *          GCG53, July 30, 1993
 */
void convhullnd(/* Input Arguments */
                const float* points,
                const int nPoints,
                const int nd,
                /* Output Arguments */
                int** faces,
                int* nFaces);

/**
 * Computes the Delaunay triangulation of an arrangement of points in
 * N-dimensional space
 *
 * This function employs algorithms originally implemented in MATLAB by George
 * Papazafeiropoulos [1] (BSD 2-clause license), which are based on [2].
 *
 * @note If you know that your points all reside on a sphere, then you should
 *       use sphDelaunay() instead; as it is faster and more accurate.
 *
 * @param[in]  points  The intput points; FLAT: nDirs x nd
 * @param[in]  nPoints Number of points
 * @param[in]  nd      The number of dimensions (max=5)
 * @param[out] DT      (&) The indices defining the Delaunay triangulation of
 *                     the points; FLAT: nDT x (nd+1)
 * @param[out] nDT     (&) Number of triangulations
 *
 * @see [1] https://www.mathworks.com/matlabcentral/fileexchange/48509-computational-geometry-toolbox?focused=3851550&tab=example
 * @see [2] The Quickhull Algorithm for Convex Hull, C. Bradford Barber, David
 *          P. Dobkin and Hannu Huhdanpaa, Geometry Center Technical Report
 *          GCG53, July 30, 1993
 */
void delaunaynd(/* Input Arguments */
                const float* points,
                const int nPoints,
                const int nd,
                /* Output Arguments */
                int** DT,
                int* nDT);

/**
 * Delaunay triangulation of a spherical arrangement of points
 *
 * @param[in]  dirs_deg Coordinates for spherically arranged points, in degrees;
 *                      FLAT: nDirs x 2
 * @param[in]  nDirs    Number of points/directions
 * @param[out] faces    (&) The face indices; FLAT: nFaces x 3
 * @param[out] nFaces   (&) Number of faces found
 * @param[out] vertices (Optional) the vertices (x,y,z) of the points (set to
 *                      NULL if not wanted); FLAT: nDirs x 3
 */
void sphDelaunay(/* Input Arguments */
                 const float* dirs_deg,
                 const int nDirs,
                 /* Output Arguments */
                 int** faces,
                 int* nFaces,
                 float* vertices);

/**
 * Computes the Voronoi diagram for a spherical arrangement of points
 *
 * @param[in]  faces    The face indices; FLAT: nFaces x 3
 * @param[in]  nFaces   Number of faces
 * @param[in]  vertices The vertices (x,y,z) of the points; FLAT: nDirs x 3
 * @param[in]  nDirs    Number of points/directions
 * @param[out] voronoi  (&) The Voronoi diagram
 */
void sphVoronoi(/* Input Arguments */
                int* faces,
                int nFaces,
                float* vertices,
                int nDirs,
                /* Output Arguments */
                voronoi_data* voronoi);

/**
 * Computes the areas of a Voronoi diagram on the unit sphere [sum(areas)=4pi]
 *
 * @param[in]  voronoi The Voronoi diagram
 * @param[out] areas   The areas; voronoi.nFaces x 1
 */
void sphVoronoiAreas(/* Input Arguments */
                     voronoi_data* voronoi,
                     /* Output Arguments */
                     float* areas);

/**
 * Computes the integration weights, based on the areas of each face of the
 * corresponding Voronoi diagram [sum(weights)=4pi]
 *
 * @test test__getVoronoiWeights()
 *
 * @param[in]  dirs_deg Coordinates for spherically arranged points, in degrees;
 *                      FLAT: nDirs x 2
 * @param[in]  nDirs    Number of points/directions
 * @param[in]  diagFLAG 0: weights returned as vector, 1: weights given along
 *                      the diagonal of a square matrix
 * @param[out] weights  The weights; nDirs x 1, or FLAT: nDirs x nDirs
 */
void getVoronoiWeights(/* Input Arguments */
                       float* dirs_deg,
                       int nDirs,
                       int diagFLAG,
                       /* Output Arguments */
                       float* weights);


#ifdef __cplusplus
}/* extern "C" */
#endif /* __cplusplus */

#endif /* SAF_GEOMETRY_H_INCLUDED */

/**@} */ /* doxygen addtogroup Utilities */
