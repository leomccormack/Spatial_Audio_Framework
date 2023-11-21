/*
 Copyright (c) 2016-2017, Symonics GmbH, Christian Hoene
 All rights reserved.

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are
 met:

     (1) Redistributions of source code must retain the above copyright
     notice, this list of conditions and the following disclaimer.

     (2) Redistributions in binary form must reproduce the above copyright
     notice, this list of conditions and the following disclaimer in
     the documentation and/or other materials provided with the
     distribution.

     (3)The name of the author may not be used to
     endorse or promote products derived from this software without
     specific prior written permission.

 THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR
 IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 DISCLAIMED. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT,
 INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
 STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
 IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 POSSIBILITY OF SUCH DAMAGE.
*/

/*
 * tools.c
 *
 *  Created on: 13.01.2017
 *      Author: hoene
 */

#if defined(SAF_ENABLE_SOFA_READER_MODULE)

#ifdef _MSC_VER
#pragma warning( disable : 4244)
#endif

#ifndef _USE_MATH_DEFINES
# define _USE_MATH_DEFINES
#endif // !_USE_MATH_DEFINES

#include "mysofa_internal.h"
#include "hdf_reader.h"
#include "../mysofa.h"
#include "../../../saf_utilities/saf_utilities.h"

/* ========================================================================== */
/*                                   CHECK                                    */
/* ========================================================================== */

static int compareValues(struct MYSOFA_ARRAY *array, const float *compare,
                         int elements, int size) {
  int i, j;
  if (array->values == NULL || (int)array->elements != elements * size)
    return 0;
  for (j = 0; j < (int)array->elements;)
    for (i = 0; i < elements; i++, j++)
      if (!fequals(array->values[j], compare[i]))
        return 0;
  return 1;
}

static const float array000[] = {0, 0, 0};
static const float array001[] = {0, 0, 1};
static const float array100[] = {1, 0, 0};

MYSOFA_EXPORT int mysofa_check(struct MYSOFA_HRTF *hrtf) {

  /* check for valid parameter ranges */
  /*
   Attributes":{
   "APIName":"ARI SOFA API for Matlab\/Octave",
   "APIVersion":"0.4.0",
   "ApplicationName":"Demo of the SOFA API",
   "ApplicationVersion":"0.4.0",
   "AuthorContact":"piotr@majdak.com",
   "Comment":"",
   "Conventions":"SOFA",
   "DataType":"FIR",
   "DatabaseName":"ARI",
   "DateCreated":"2014-03-20 17:35:22",
   "DateModified":"2014-03-20 17:35:22",
   "History":"Converted from the ARI format",
   "License":"No license provided, ask the author for permission",
   "ListenerShortName":"",
   "Organization":"Acoustics Research Institute",
   "Origin":"",
   "References":"",
   "RoomType":"free field",
   "SOFAConventions":"SimpleFreeFieldHRIR",
   "SOFAConventionsVersion":"0.4",
   "Title":"",
   "Version":"0.6"
   },
   */
  if (!verifyAttribute(hrtf->attributes, "Conventions", "SOFA") ||
      !verifyAttribute(hrtf->attributes, "SOFAConventions",
                       "SimpleFreeFieldHRIR") ||
      /* TODO: Support FT too */
      !verifyAttribute(hrtf->attributes, "DataType", "FIR"))
    return MYSOFA_INVALID_ATTRIBUTES; // LCOV_EXCL_LINE

  if (!verifyAttribute(hrtf->attributes, "RoomType", "free field") &&
      !verifyAttribute(hrtf->attributes, "RoomType", "reverberant") &&
      !verifyAttribute(hrtf->attributes, "RoomType", "shoebox"))
    return MYSOFA_INVALID_ATTRIBUTES; // LCOV_EXCL_LINE

  /*==============================================================================
   dimensions
   ==============================================================================
 */

  if (hrtf->C != 3 || hrtf->_I != 1 || hrtf->E != 1 || hrtf->R != 2 ||
      hrtf->M == 0)
    return MYSOFA_INVALID_DIMENSIONS; // LCOV_EXCL_LINE

  /* verify format */

  if (hrtf->ListenerView.values) {
    int m = 1;
    if (!verifyAttribute(hrtf->ListenerView.attributes, "DIMENSION_LIST",
                         "I,C")) {
      if (!verifyAttribute(hrtf->ListenerView.attributes, "DIMENSION_LIST",
                           "M,C")) {
        return MYSOFA_INVALID_DIMENSION_LIST; // LCOV_EXCL_LINE
      }
      m = hrtf->M;
    }
    if (verifyAttribute(hrtf->ListenerView.attributes, "Type", "cartesian")) {
      if (!compareValues(&hrtf->ListenerView, array100, 3, m))
        return MYSOFA_INVALID_FORMAT; // LCOV_EXCL_LINE
    } else if (verifyAttribute(hrtf->ListenerView.attributes, "Type",
                               "spherical")) {
      if (!compareValues(&hrtf->ListenerView, array001, 3, m))
        return MYSOFA_INVALID_FORMAT; // LCOV_EXCL_LINE
    } else
      return MYSOFA_INVALID_COORDINATE_TYPE; // LCOV_EXCL_LINE
  }

#if 0
    if(hrtf->ListenerUp.values) {
        if(!verifyAttribute(hrtf->ListenerUp.attributes,"DIMENSION_LIST","I,C"))
        return MYSOFA_INVALID_FORMAT;
        if(verifyAttribute(hrtf->ListenerUp.attributes,"Type","cartesian")) {
            if(!compareValues(&hrtf->ListenerUp,array001,3))
            return MYSOFA_INVALID_FORMAT;
        }
        else if(verifyAttribute(hrtf->ListenerUp.attributes,"Type","spherical")) {
            if(!compareValues(&hrtf->ListenerUp,array0901,3))
            return MYSOFA_INVALID_FORMAT;
        }
    }

    /* TODO. support M,C too */
    if(!verifyAttribute(hrtf->ListenerPosition.attributes,"DIMENSION_LIST","I,C"))
    return MYSOFA_INVALID_FORMAT;
    if(!compareValues(&hrtf->ListenerPosition,array000,3))
    return MYSOFA_INVALID_FORMAT;
#endif

  int m = 1;
  if (!verifyAttribute(hrtf->EmitterPosition.attributes, "DIMENSION_LIST",
                       "E,C,I")) {
    if (!verifyAttribute(hrtf->EmitterPosition.attributes, "DIMENSION_LIST",
                         "E,C,M")) {
      return MYSOFA_ONLY_EMITTER_WITH_ECI_SUPPORTED; // LCOV_EXCL_LINE
    }
    m = hrtf->M;
  }

  if (!compareValues(&hrtf->EmitterPosition, array000, 3, m))
    return MYSOFA_ONLY_EMITTER_WITH_ECI_SUPPORTED; // LCOV_EXCL_LINE

  if (hrtf->DataDelay.values) {
    if (!verifyAttribute(hrtf->DataDelay.attributes, "DIMENSION_LIST", "I,R") &&
        !verifyAttribute(hrtf->DataDelay.attributes, "DIMENSION_LIST", "M,R"))
      return MYSOFA_ONLY_DELAYS_WITH_IR_OR_MR_SUPPORTED; // LCOV_EXCL_LINE
  }
  /* TODO: Support different sampling rate per measurement, support default
   sampling rate of 48000 However, so far, I have not seen any sofa files with
   an format other and I */
  if (!verifyAttribute(hrtf->DataSamplingRate.attributes, "DIMENSION_LIST",
                       "I"))
    return MYSOFA_ONLY_THE_SAME_SAMPLING_RATE_SUPPORTED; // LCOV_EXCL_LINE

  if (verifyAttribute(hrtf->ReceiverPosition.attributes, "DIMENSION_LIST",
                      "R,C,I")) {
    // do nothing
  } else if (verifyAttribute(hrtf->ReceiverPosition.attributes,
                             "DIMENSION_LIST", "R,C,M")) {
    for (int i = 0; i < 6; i++) {
      int offset = i * hrtf->M;
      double receiverPosition = hrtf->ReceiverPosition.values[offset];
      for (int j = 1; j < (int)hrtf->M; j++)
        if (!fequals(receiverPosition,
                     hrtf->ReceiverPosition.values[offset + j]))
          return MYSOFA_RECEIVERS_WITH_RCI_SUPPORTED; // LCOV_EXCL_LINE
    }
  } else {
    return MYSOFA_RECEIVERS_WITH_RCI_SUPPORTED; // LCOV_EXCL_LINE
  }

  if (!verifyAttribute(hrtf->ReceiverPosition.attributes, "Type", "cartesian"))
    return MYSOFA_RECEIVERS_WITH_CARTESIAN_SUPPORTED; // LCOV_EXCL_LINE

  if (hrtf->ReceiverPosition.elements < 6 ||
      !fequals(hrtf->ReceiverPosition.values[0], 0.f) ||
      !fequals(hrtf->ReceiverPosition.values[2], 0.f) ||
      !fequals(hrtf->ReceiverPosition.values[3], 0.f) ||
      !fequals(hrtf->ReceiverPosition.values[5], 0.f)) {
    return MYSOFA_INVALID_RECEIVER_POSITIONS; // LCOV_EXCL_LINE
  }
  if (!fequals(hrtf->ReceiverPosition.values[4],
               -hrtf->ReceiverPosition.values[1]))
    return MYSOFA_INVALID_RECEIVER_POSITIONS; // LCOV_EXCL_LINE
  if (hrtf->ReceiverPosition.values[1] < 0) {
    if (!verifyAttribute(hrtf->attributes, "APIName",
                         "ARI SOFA API for Matlab/Octave"))
      return MYSOFA_INVALID_RECEIVER_POSITIONS; // LCOV_EXCL_LINE

    const char *version = mysofa_getAttribute(hrtf->attributes, "APIVersion");
    if (version == NULL)
      return MYSOFA_INVALID_RECEIVER_POSITIONS; // LCOV_EXCL_LINE

    int a, b, c;
    int res = sscanf(version, "%d.%d.%d", &a, &b, &c);
    if (res != 3)
      return MYSOFA_INVALID_RECEIVER_POSITIONS; // LCOV_EXCL_LINE
    if (a > 1)
      return MYSOFA_INVALID_RECEIVER_POSITIONS; // LCOV_EXCL_LINE
    if (a == 1 && b > 1)
      return MYSOFA_INVALID_RECEIVER_POSITIONS; // LCOV_EXCL_LINE
    if (a == 1 && b == 1 && c > 0)
      return MYSOFA_INVALID_RECEIVER_POSITIONS; // LCOV_EXCL_LINE

    if (hrtf->ReceiverPosition.values[1] >= 0)
      return MYSOFA_INVALID_RECEIVER_POSITIONS; // LCOV_EXCL_LINE

    // old versions of sofaapi sometimes has been used wrongly. Thus, they wrote
    // left and right ears to different possitions
    mylog("WARNING: SOFA file is written with wrong receiver positions. %d "
          "%d.%d.%d %f<>%f\n",
          res, a, b, c, hrtf->ReceiverPosition.values[1],
          hrtf->ReceiverPosition.values[4]);
  }

  /* read source positions */
  if (!verifyAttribute(hrtf->SourcePosition.attributes, "DIMENSION_LIST",
                       "M,C"))
    return MYSOFA_ONLY_SOURCES_WITH_MC_SUPPORTED; // LCOV_EXCL_LINE

  return MYSOFA_OK;
}


/* ========================================================================== */
/*                                INTERPOLATE                                 */
/* ========================================================================== */

MYSOFA_EXPORT float *mysofa_interpolate(struct MYSOFA_HRTF *hrtf,
                                        float *cordinate, int nearest,
                                        int *neighborhood, float *fir,
                                        float *delays) {
  int i, use[6];
  float d, d6[6];
  float weight;
  int size = hrtf->N * hrtf->R;

  d = distance(cordinate, hrtf->SourcePosition.values + nearest * hrtf->C);
  if (fequals(d, 0)) {
    if (hrtf->DataDelay.elements > hrtf->R) {
      delays[0] = hrtf->DataDelay.values[nearest * hrtf->R];
      delays[1] = hrtf->DataDelay.values[nearest * hrtf->R + 1];
    } else {
      delays[0] = hrtf->DataDelay.values[0];
      delays[1] = hrtf->DataDelay.values[1];
    }
    float *ret = hrtf->DataIR.values + nearest * size;
    copyFromFloat(fir, ret, size);
    return ret;
  }

  for (i = 0; i < 6; i++) {
    use[i] = 0;
    d6[i] = 1;
  }

  if (neighborhood[0] >= 0 && neighborhood[1] >= 0) {
    d6[0] = distance(cordinate,
                     hrtf->SourcePosition.values + neighborhood[0] * hrtf->C);
    d6[1] = distance(cordinate,
                     hrtf->SourcePosition.values + neighborhood[1] * hrtf->C);

    if (!fequals(d6[0], d6[1])) {
      if (d6[0] < d6[1])
        use[0] = 1;
      else
        use[1] = 1;
    }
  } else if (neighborhood[0] >= 0) {
    d6[0] = distance(cordinate,
                     hrtf->SourcePosition.values + neighborhood[0] * hrtf->C);
    use[0] = 1;
  } else if (neighborhood[1] >= 0) {
    d6[1] = distance(cordinate,
                     hrtf->SourcePosition.values + neighborhood[1] * hrtf->C);
    use[1] = 1;
  }

  if (neighborhood[2] >= 0 && neighborhood[3] >= 0) {
    d6[2] = distance(cordinate,
                     hrtf->SourcePosition.values + neighborhood[2] * hrtf->C);
    d6[3] = distance(cordinate,
                     hrtf->SourcePosition.values + neighborhood[3] * hrtf->C);
    if (!fequals(d6[2], d6[3])) {
      if (d6[2] < d6[3])
        use[2] = 1;
      else
        use[3] = 1;
    }
  } else if (neighborhood[2] >= 0) {
    d6[2] = distance(cordinate,
                     hrtf->SourcePosition.values + neighborhood[2] * hrtf->C);
    use[2] = 1;
  } else if (neighborhood[3] >= 0) {
    d6[3] = distance(cordinate,
                     hrtf->SourcePosition.values + neighborhood[3] * hrtf->C);
    use[3] = 1;
  }

  if (neighborhood[4] >= 0 && neighborhood[5] >= 0) {
    d6[4] = distance(cordinate,
                     hrtf->SourcePosition.values + neighborhood[4] * hrtf->C);
    d6[5] = distance(cordinate,
                     hrtf->SourcePosition.values + neighborhood[5] * hrtf->C);
    if (!fequals(d6[4], d6[5])) {
      if (d6[4] < d6[5])
        use[4] = 1;
      else
        use[5] = 1;
    }
  } else if (neighborhood[4] >= 0) {
    d6[4] = distance(cordinate,
                     hrtf->SourcePosition.values + neighborhood[4] * hrtf->C);
    use[4] = 1;
  } else if (neighborhood[5] >= 0) {
    d6[5] = distance(cordinate,
                     hrtf->SourcePosition.values + neighborhood[5] * hrtf->C);
    use[5] = 1;
  }

  weight = 1 / d;
  copyArrayWeighted(fir, hrtf->DataIR.values + nearest * size, size, weight);
  if (hrtf->DataDelay.elements > hrtf->R) {
    delays[0] = hrtf->DataDelay.values[nearest * hrtf->R] * weight;
    delays[1] = hrtf->DataDelay.values[nearest * hrtf->R + 1] * weight;
  } else {
    delays[0] = hrtf->DataDelay.values[0] * weight;
    delays[1] = hrtf->DataDelay.values[1] * weight;
  }
#ifdef VDEBUG
  printf("%d ! %f ", nearest, d);
#endif
  for (i = 0; i < 6; i++) {
    if (use[i]) {
      float w = 1 / d6[i];
#ifdef VDEBUG
      printf("%d - %f ", neighborhood[i], d6[i]);
#endif
      addArrayWeighted(fir, hrtf->DataIR.values + neighborhood[i] * size, size,
                       w);
      weight += w;
      if (hrtf->DataDelay.elements > hrtf->R) {
        delays[0] += hrtf->DataDelay.values[neighborhood[i] * hrtf->R] * w;
        delays[1] += hrtf->DataDelay.values[neighborhood[i] * hrtf->R + 1] * w;
      }
    }
  }
#ifdef VDEBUG
  printf("\n");
#endif
  weight = 1 / weight;
  scaleArray(fir, size, weight);
  delays[0] *= weight;
  delays[1] *= weight;
  return fir;
}


/* ========================================================================== */
/*                                  LOUDNESS                                  */
/* ========================================================================== */

MYSOFA_EXPORT float mysofa_loudness(struct MYSOFA_HRTF *hrtf) {
  float c[3], factor;
  float min = FLT_MAX;
  int radius = 0;
  unsigned int i, index = 0;
  int cartesian =
      verifyAttribute(hrtf->SourcePosition.attributes, "Type", "cartesian");

  /*
   * find frontal source position
   */
  for (i = 0; i < hrtf->SourcePosition.elements; i += hrtf->C) {
    c[0] = hrtf->SourcePosition.values[i];
    c[1] = hrtf->SourcePosition.values[i + 1];
    c[2] = hrtf->SourcePosition.values[i + 2];

    if (cartesian)
      mysofa_c2s(c);

    if (min > c[0] + c[1]) {
      min = c[0] + c[1];
      radius = c[2];
      index = i;
    } else if (min == c[0] + c[1] && radius < c[2]) {
      radius = c[2];
      index = i;
    }
  }

  /* get loudness of frontal fir filter, for both channels */
  factor = loudness(hrtf->DataIR.values + (index / hrtf->C) * hrtf->N * hrtf->R,
                    hrtf->N * hrtf->R);
  factor = sqrtf(2 / factor);
  if (fequals(factor, 1.f))
    return 1.f;

  scaleArray(hrtf->DataIR.values, hrtf->DataIR.elements, factor);

  return factor;
}

/* ========================================================================== */
/*                                 SPHERICAL                                  */
/* ========================================================================== */

static void convertArray(struct MYSOFA_ARRAY *array) {
  if (!changeAttribute(array->attributes, "Type", "cartesian", "spherical"))
    return;

  changeAttribute(array->attributes, "Units", NULL, "degree, degree, meter");

  convertCartesianToSpherical(array->values, array->elements);
}

MYSOFA_EXPORT void mysofa_tospherical(struct MYSOFA_HRTF *hrtf) {
  convertArray(&hrtf->ListenerView);
  convertArray(&hrtf->ListenerUp);
  convertArray(&hrtf->ListenerPosition);
  convertArray(&hrtf->EmitterPosition);
  convertArray(&hrtf->ReceiverPosition);
  convertArray(&hrtf->SourcePosition);
}

static void convertArray2(struct MYSOFA_ARRAY *array) {
  if (!changeAttribute(array->attributes, "Type", "spherical", "cartesian"))
    return;

  changeAttribute(array->attributes, "Units", NULL, "meter");

  convertSphericalToCartesian(array->values, array->elements);
}

MYSOFA_EXPORT void mysofa_tocartesian(struct MYSOFA_HRTF *hrtf) {
  convertArray2(&hrtf->ListenerView);
  convertArray2(&hrtf->ListenerUp);
  convertArray2(&hrtf->ListenerPosition);
  convertArray2(&hrtf->EmitterPosition);
  convertArray2(&hrtf->ReceiverPosition);
  convertArray2(&hrtf->SourcePosition);
}

/* ========================================================================== */
/*                                NEIGHBORS                                   */
/* ========================================================================== */

MYSOFA_EXPORT struct MYSOFA_NEIGHBORHOOD *
mysofa_neighborhood_init(struct MYSOFA_HRTF *hrtf,
                         struct MYSOFA_LOOKUP *lookup) {
  return mysofa_neighborhood_init_withstepdefine(
      hrtf, lookup, MYSOFA_DEFAULT_NEIGH_STEP_ANGLE,
      MYSOFA_DEFAULT_NEIGH_STEP_RADIUS);
}

MYSOFA_EXPORT struct MYSOFA_NEIGHBORHOOD *
mysofa_neighborhood_init_withstepdefine(struct MYSOFA_HRTF *hrtf,
                                        struct MYSOFA_LOOKUP *lookup,
                                        float angleStep, float radiusStep) {
  int i, index;
  float *origin, *test;
  float radius, radius2;
  float theta;
  float phi;

  // distance (degree) beyond which neighbor search is abandonned
  float maxNeighborSearchAngle = 45;

  struct MYSOFA_NEIGHBORHOOD *neighbor =
      malloc(sizeof(struct MYSOFA_NEIGHBORHOOD));
  if (!neighbor)
    return NULL;

  neighbor->elements = hrtf->M;
  neighbor->index = malloc(sizeof(int) * neighbor->elements * 6);
  if (!neighbor->index) {
    free(neighbor);
    return NULL;
  }
  for (i = 0; i < neighbor->elements * 6; i++)
    neighbor->index[i] = -1;

  origin = malloc(sizeof(float) * hrtf->C);
  test = malloc(sizeof(float) * hrtf->C);

  for (i = 0; i < (int)hrtf->M; i++) {
    memcpy(origin, hrtf->SourcePosition.values + i * hrtf->C,
           sizeof(float) * hrtf->C);
    convertCartesianToSpherical(origin, hrtf->C);

    if ((lookup->phi_max - lookup->phi_min) > FLT_MIN) {
      phi = angleStep;
      do {
        test[0] = origin[0] + phi;
        test[1] = origin[1];
        test[2] = origin[2];
        convertSphericalToCartesian(test, 3);
        index = mysofa_lookup(lookup, test);
        if (index != i) {
          neighbor->index[i * 6 + 0] = index;
          break;
        }
        phi += angleStep;
      } while (phi <= maxNeighborSearchAngle);

      phi = -angleStep;
      do {
        test[0] = origin[0] + phi;
        test[1] = origin[1];
        test[2] = origin[2];
        convertSphericalToCartesian(test, 3);
        index = mysofa_lookup(lookup, test);
        if (index != i) {
          neighbor->index[i * 6 + 1] = index;
          break;
        }
        phi -= angleStep;
      } while (phi >= -maxNeighborSearchAngle);
    }

    if ((lookup->theta_max - lookup->theta_min) > FLT_MIN) {
      theta = angleStep;
      do {
        test[0] = origin[0];
        test[1] = origin[1] + theta;
        test[2] = origin[2];
        convertSphericalToCartesian(test, 3);
        index = mysofa_lookup(lookup, test);
        if (index != i) {
          neighbor->index[i * 6 + 2] = index;
          break;
        }
        theta += angleStep;
      } while (theta <= maxNeighborSearchAngle);

      theta = -angleStep;
      do {
        test[0] = origin[0];
        test[1] = origin[1] + theta;
        test[2] = origin[2];
        convertSphericalToCartesian(test, 3);
        index = mysofa_lookup(lookup, test);
        if (index != i) {
          neighbor->index[i * 6 + 3] = index;
          break;
        }
        theta -= angleStep;
      } while (theta >= -maxNeighborSearchAngle);
    }

    if ((lookup->radius_max - lookup->radius_min) > FLT_MIN) {
      radius = radiusStep;
      do {
        test[0] = origin[0];
        test[1] = origin[1];
        radius2 = test[2] = origin[2] + radius;
        convertSphericalToCartesian(test, 3);
        index = mysofa_lookup(lookup, test);
        if (index != i) {
          neighbor->index[i * 6 + 4] = index;
          break;
        }
        radius += radiusStep;
      } while (radius2 <= lookup->radius_max + radiusStep);

      radius = -radiusStep;
      do {
        test[0] = origin[0];
        test[1] = origin[1];
        radius2 = test[2] = origin[2] + radius;
        convertSphericalToCartesian(test, 3);
        index = mysofa_lookup(lookup, test);
        if (index != i) {
          neighbor->index[i * 6 + 5] = index;
          break;
        }
        radius -= radiusStep;
      } while (radius2 >= lookup->radius_min - radiusStep);
    }
  }
  free(test);
  free(origin);
  return neighbor;
}

MYSOFA_EXPORT int *mysofa_neighborhood(struct MYSOFA_NEIGHBORHOOD *neighborhood,
                                       int index) {
  if (index < 0 || index >= neighborhood->elements)
    return NULL;
  return neighborhood->index + index * 6;
}

MYSOFA_EXPORT void
mysofa_neighborhood_free(struct MYSOFA_NEIGHBORHOOD *neighborhood) {
  if (neighborhood) {
    free(neighborhood->index);
    free(neighborhood);
  }
}


/* ========================================================================== */
/*                                 MINPHASE                                   */
/* ========================================================================== */

static void trunk(float *in, int size, int *start, int *end, float threshold) {
  float energy = 0;
  int s = 0;
  int e = size - 1;
  float ss, ee;

  float l = loudness(in, size);
  threshold = threshold * l;

  ss = in[s] * in[s];
  ee = in[e] * in[e];
  while (s < e) {
    if (ss <= ee) {
      if (energy + ss > threshold)
        break;
      energy += ss;
      s++;
      ss = in[s] * in[s];
    } else {
      if (energy + ee > threshold)
        break;
      energy += ee;
      e--;
      ee = in[e] * in[e];
    }
  }
  *start = s;
  *end = e + 1;
}

MYSOFA_EXPORT int mysofa_minphase(struct MYSOFA_HRTF *hrtf, float threshold) {
  int i;
  int max = 0;
  int filters;
  int *start;
  int *end;
  float samplerate;
  float d[2];

  if (hrtf->DataDelay.elements != 2)
    return -1;

  filters = hrtf->M * hrtf->R;
  start = malloc(filters * sizeof(int));
  end = malloc(filters * sizeof(int));

  /*
   * find maximal length of a filter
   */
  for (i = 0; i < filters; i++) {
    trunk(hrtf->DataIR.values + i * hrtf->N, hrtf->N, start + i, end + i,
          threshold);
    if (end[i] - start[i] > max)
      max = end[i] - start[i];
  }

  if (max == (int)hrtf->N) {
    free(start);
    free(end);
    return max;
  }

  /*
   * update delay and filters
   */
  samplerate = hrtf->DataSamplingRate.values[0];
  d[0] = hrtf->DataDelay.values[0];
  d[1] = hrtf->DataDelay.values[1];
  hrtf->DataDelay.elements = filters;
  hrtf->DataDelay.values =
      realloc(hrtf->DataDelay.values, sizeof(float) * filters);
  for (i = 0; i < filters; i++) {
    if (start[i] + max > (int)hrtf->N)
      start[i] = hrtf->N - max;
    hrtf->DataDelay.values[i] = d[i % 1] + (start[i] / samplerate);
    memmove(hrtf->DataIR.values + i * max,
            hrtf->DataIR.values + i * hrtf->N + start[i], max * sizeof(float));
  }

  /*
   * update hrtf structure
   */
  hrtf->N = max;
  hrtf->DataIR.elements = max * filters;
  hrtf->DataIR.values =
      realloc(hrtf->DataIR.values, sizeof(float) * hrtf->DataIR.elements);

  free(start);
  free(end);
  return max;
}

/* ========================================================================== */
/*                                 RESAMPLE                                   */
/* ========================================================================== */

MYSOFA_EXPORT int mysofa_resample(struct MYSOFA_HRTF *hrtf, float samplerate) {
  unsigned int i;
  int err;
  float factor;
  unsigned newN;
  float *values;
  SpeexResamplerState *resampler;
  float *out;
  float zero[10] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

  if (hrtf->DataSamplingRate.elements != 1 || samplerate < 8000.)
    return MYSOFA_INVALID_FORMAT;

  if (samplerate == hrtf->DataSamplingRate.values[0])
    return MYSOFA_OK;

  factor = samplerate / hrtf->DataSamplingRate.values[0];
  newN = ceil(hrtf->N * factor);

  /*
   * resample FIR filter
   */
  values = malloc(newN * hrtf->R * hrtf->M * sizeof(float));
  if (values == NULL)
    return MYSOFA_NO_MEMORY;

  resampler = speex__resampler_init(1, hrtf->DataSamplingRate.values[0],
                                   samplerate, 10, &err);
  if (resampler == NULL) {
    free(values);
    return err;
  }

  out = malloc(sizeof(float) *
               (newN + speex__resampler_get_output_latency(resampler)));
  for (i = 0; i < hrtf->R * hrtf->M; i++) {
    unsigned inlen = hrtf->N;
    unsigned outlen = newN;
    speex__resampler_reset_mem(resampler);
    speex__resampler_skip_zeros(resampler);
    speex__resampler_process_float(resampler, 0,
                                  hrtf->DataIR.values + i * hrtf->N, &inlen,
                                  values + i * newN, &outlen);
    assert(inlen == hrtf->N);
    while (outlen < newN) {
      unsigned difflen = newN - outlen;
      inlen = 10;
      speex__resampler_process_float(resampler, 0, zero, &inlen,
                                    values + i * newN + outlen, &difflen);
      outlen += difflen;
    }
  }
  free(out);
  speex__resampler_destroy(resampler);

  free(hrtf->DataIR.values);
  hrtf->DataIR.values = values;
  hrtf->DataIR.elements = newN * hrtf->R * hrtf->M;

  /*
   * update delay values
   */
  for (i = 0; i < hrtf->DataDelay.elements; i++)
    hrtf->DataDelay.values[i] *= factor;

  /*
   * update sample rate
   */
  hrtf->DataSamplingRate.values[0] = samplerate;
  hrtf->N = newN;

  return MYSOFA_OK;
}


/* ========================================================================== */
/*                                  TOOLS                                     */
/* ========================================================================== */

char *mysofa_strdup(const char *str) {
  size_t size = strlen(str) + 1;
  char *copy = malloc(size);
  if (copy)
    memcpy(copy, str, size);
  return copy;
}

int verifyAttribute(struct MYSOFA_ATTRIBUTE *attr, char *name, char *value) {
  while (attr) {
    if (attr->name && !strcmp(name, attr->name) && attr->value &&
        !strcmp(value, attr->value))
      return 1;
    attr = attr->next;
  }
  return 0;
}

int changeAttribute(struct MYSOFA_ATTRIBUTE *attr, char *name, char *value,
                    char *newvalue) {
  while (attr) {
    if (!strcmp(name, attr->name) &&
        (value == NULL || attr->value == NULL || !strcmp(value, attr->value))) {
      free(attr->value);
      attr->value = mysofa_strdup(newvalue);
      return 1;
    }
    attr = attr->next;
  }
  return 0;
}

MYSOFA_EXPORT
char *mysofa_getAttribute(struct MYSOFA_ATTRIBUTE *attr, char *name) {
  while (attr) {
    if (attr->name && !strcmp(name, attr->name)) {
      return attr->value;
    }
    attr = attr->next;
  }
  return NULL;
}

MYSOFA_EXPORT void mysofa_c2s(float values[3]) {
  float x, y, z, r, theta, phi;
  x = values[0];
  y = values[1];
  z = values[2];
  r = radius(values);

  theta = atan2f(z, sqrtf(x * x + y * y));
  phi = atan2f(y, x);

  values[0] = fmodf(phi * (180.0 / SAF_PId) + 360, 360);
  values[1] = theta * (180.0 / SAF_PId);
  values[2] = r;
}

MYSOFA_EXPORT void mysofa_s2c(float values[3]) {
  float x, r, theta, phi;
  phi = values[0] * (SAF_PId / 180.0);
  theta = values[1] * (SAF_PId / 180.0);
  r = values[2];
  x = cosf(theta) * r;
  values[2] = sinf(theta) * r;
  values[0] = cosf(phi) * x;
  values[1] = sinf(phi) * x;
}

void convertCartesianToSpherical(float *values, int elements) {
  int i;

  for (i = 0; i < elements - 2; i += 3) {
    mysofa_c2s(values + i);
  }
}

void convertSphericalToCartesian(float *values, int elements) {

  int i;

  for (i = 0; i < elements - 2; i += 3) {
    mysofa_s2c(values + i);
  }
}

float radius(float *cartesian) {
  return sqrtf(powf(cartesian[0], 2.f) + powf(cartesian[1], 2.f) +
               powf(cartesian[2], 2.f));
}

/*
 * search of the nearest
 */

void nsearch(const void *key, const char *base, size_t num, size_t size,
             int (*cmp)(const void *key, const void *elt), int *lower,
             int *higher) {
  size_t start = 0, end = num;
  int result;

  while (start < end) {
    size_t mid = start + (end - start) / 2;

    result = cmp(key, base + mid * size);
    if (result < 0)
      end = mid;
    else if (result > 0)
      start = mid + 1;
    else {
      *lower = (int)mid;
      *higher = (int)mid;
      return;
    }
  }

  if (start == num) {
    *lower = (int)start - 1;
    *higher = -1;
  } else if (start == 0) {
    *lower = -1;
    *higher = 0;
  } else {
    *lower = (int)start - 1;
    *higher = (int)start;
  }
}

void copyToFloat(float *out, float *in, int size) {
  while (size > 0) {
    *out++ = *in++;
    size--;
  }
}

void copyFromFloat(float *out, float *in, int size) {
  while (size > 0) {
    *out++ = *in++;
    size--;
  }
}

void copyArrayWeighted(float *dst, float *src, int size, float w) {
  while (size > 0) {
    *dst++ = *src++ * w;
    size--;
  }
}

void addArrayWeighted(float *dst, float *src, int size, float w) {
  while (size > 0) {
    *dst++ += *src++ * w;
    size--;
  }
}

void scaleArray(float *dst, int size, float w) {
  while (size > 0) {
    *dst++ *= w;
    size--;
  }
}

float loudness(float *in, int size) {
  float res = 0;
  while (size > 0) {
    res += *in * *in;
    in++;
    size--;
  }
  return res;
}

#else
extern int to_avoid_iso_compiler_warning_when_there_are_no_symbols;
#endif /* SAF_ENABLE_SOFA_READER_MODULE */
