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
 */

#ifndef SAF_GEOMETRY_H_INCLUDED
#define SAF_GEOMETRY_H_INCLUDED

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/**
 * Data structure for Voronoi diagrams
 *
 * @warning 'faces' is NOT contiguously allocated! [i.e., you cannot do
 *          free(faces), you must loop over nFaces: free(faces[i]) ]
 */
typedef struct _voronoi_data{
    int nVert;           /**< Number of vertices */
    int nFaces;          /**< Number of faces/polygons */
    float** vert;        /**< Vertices; nVert x 3 */
    int** faces;         /**< faces; nFaces x nPointsPerFace[i]  */
    int* nPointsPerFace; /**< Number of points for each face; nFaces x 1 */

}voronoi_data;

/**
 * Constructs a 3x3 rotation matrix from the Euler angles, using the
 * yaw-pitch-roll (zyx) convention
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

/**
 * L2 (Euclidean) norm of a 3-element vector
 */
float L2_norm(float v[3]);

/**
 * Cross product between two 3-element floating point vectors (c = a x b)
 */
void crossProduct(/* Input Arguments */
                  float a[3],
                  float b[3],
                  /* Output Arguments */
                  float c[3]);

/**
 * Builds the 3-D convex hull given a list of vertices
 *
 * @warning Currently, this does not check if there are duplicate vertices or
 *          whether any of them are co-linear.
 *
 * @param[in]  vertices The vertices; FLAT: nDirs x 3
 * @param[in]  nDirs    Number of vertices
 * @param[out] faces    (&) The face indices; FLAT: nFaces x 3
 * @param[out] nFaces   (&) Number of faces found
 */
void convhull3d(/* Input Arguments */
                const float* vertices,
                const int nDirs,
                /* Output Arguments */
                int** faces,
                int* nFaces);

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
