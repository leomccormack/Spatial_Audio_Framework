/*
 Copyright (c) 2016, Symonics GmbH, Christian Hoene
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

#if defined(SAF_ENABLE_SOFA_READER_MODULE)

#ifdef _MSC_VER
#pragma warning( disable : 4244)
#endif

#include "internal/hdf_reader.h"
#include "internal/kdtree.h"
#include "internal/mysofa_internal.h"
#include "mysofa.h"

#include <ctype.h>
#include <errno.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <float.h>

#define CMAKE_INSTALL_PREFIX ""
#define CPACK_PACKAGE_VERSION_MAJOR 0
#define CPACK_PACKAGE_VERSION_MINOR 0
#define CPACK_PACKAGE_VERSION_PATCH 0

/* ========================================================================== */
/*                                   MYSOFA                                   */
/* ========================================================================== */

static int mystrcmp(char *s1, char *s2) {
  if (s1 == NULL && s2 == NULL)
    return 0;
  if (s1 == NULL)
    return -1;
  if (s2 == NULL)
    return 1;
  return strcmp(s1, s2);
}

static int checkAttribute(struct MYSOFA_ATTRIBUTE *attribute, char *name,
                          char *value) {
  while (attribute) {
    if (!mystrcmp(attribute->name, name) && !mystrcmp(attribute->value, value))
      return MYSOFA_OK;
    attribute = attribute->next;
  }

  return MYSOFA_INVALID_FORMAT;
}

static int getDimension(unsigned *dim, struct DATAOBJECT *dataobject) {
  int err;
  struct MYSOFA_ATTRIBUTE *attr = dataobject->attributes;

  if (!!(err = checkAttribute(dataobject->attributes, "CLASS",
                              "DIMENSION_SCALE")))
    return err;

  while (attr) {
    mylog(" %s=%s\n", attr->name, attr->value);

    if (!strcmp(attr->name, "NAME") && attr->value &&
        !strncmp(attr->value,
                 "This is a netCDF dimension but not a netCDF variable.", 53)) {
      char *p = attr->value + strlen(attr->value) - 1;
      while (isdigit(*p)) {
        p--;
      }
      p++;
      *dim = atoi(p);
      mylog("NETCDF DIM %u\n", *dim);
      return MYSOFA_OK;
    }
    attr = attr->next;
  }
  return MYSOFA_INVALID_FORMAT;
}

static int getArray(struct MYSOFA_ARRAY *array, struct DATAOBJECT *dataobject) {
  float *p1;
  double *p2;
  unsigned int i;

  struct MYSOFA_ATTRIBUTE *attr = dataobject->attributes;
  while (attr) {
    mylog(" %s=%s\n", attr->name ? attr->name : "(null)",
          attr->value ? attr->value : "(null)");
    attr = attr->next;
  }

  if (dataobject->dt.u.f.bit_precision != 64)
    return MYSOFA_UNSUPPORTED_FORMAT;

  array->attributes = dataobject->attributes;
  dataobject->attributes = NULL;
  array->elements = dataobject->data_len / 8;

  p1 = dataobject->data;
  p2 = dataobject->data;
  for (i = 0; i < array->elements; i++)
    *p1++ = (float)*p2++;
  array->values = realloc(dataobject->data, array->elements * sizeof(float));

  dataobject->data = NULL;

  return MYSOFA_OK;
}

static void arrayFree(struct MYSOFA_ARRAY *array) {
  while (array->attributes) {
    struct MYSOFA_ATTRIBUTE *next = array->attributes->next;
    free(array->attributes->name);
    free(array->attributes->value);
    free(array->attributes);
    array->attributes = next;
  }
  free(array->values);
}

static int addUserDefinedVariable(struct MYSOFA_HRTF *hrtf,
                                  struct DATAOBJECT *dataobject) {
  int err;

  // init variable
  struct MYSOFA_VARIABLE *var = malloc(sizeof(struct MYSOFA_VARIABLE));

  if (!var) {
    return errno;
  }

  memset(var, 0, sizeof(struct MYSOFA_VARIABLE));

  // init values array
  var->value = malloc(sizeof(struct MYSOFA_ARRAY));
  if (!var->value) {
    free(var);
    return errno;
  }
  memset(var->value, 0, sizeof(struct MYSOFA_ARRAY));

  var->next = NULL;
  // copy name
  var->name = malloc(strlen(dataobject->name) + 1);
  if (!var->name) {
    free(var->value);
    free(var);
    return errno;
  }
  strcpy(var->name, dataobject->name);

  err = getArray(var->value, dataobject);
  if (err != MYSOFA_OK) {
    arrayFree(var->value);
    free(var->value);
    free(var->name);
    free(var);
    return err;
  }

  var->next = hrtf->variables;
  hrtf->variables = var;

  return MYSOFA_OK;
}

static struct MYSOFA_HRTF *getHrtf(struct READER *reader, int *err) {
  int dimensionflags = 0;
  struct DIR *dir = reader->superblock.dataobject.directory;

  struct MYSOFA_HRTF *hrtf = malloc(sizeof(struct MYSOFA_HRTF));
  if (!hrtf) {
    *err = errno;
    return NULL;
  }
  memset(hrtf, 0, sizeof(struct MYSOFA_HRTF));

  /* copy SOFA file attributes */
  hrtf->attributes = reader->superblock.dataobject.attributes;
  reader->superblock.dataobject.attributes = NULL;

  /* check SOFA file attributes */
  if (!!(*err = checkAttribute(hrtf->attributes, "Conventions", "SOFA"))) {
    mylog("no Conventions=SOFA attribute\n");
    goto error;
  }

  /* read dimensions */
  while (dir) {
    if (dir->dataobject.name && dir->dataobject.name[0] &&
        dir->dataobject.name[1] == 0) {
      switch (dir->dataobject.name[0]) {
      case 'I':
        *err = getDimension(&hrtf->_I, &dir->dataobject);
        dimensionflags |= 1;
        break;
      case 'C':
        *err = getDimension(&hrtf->C, &dir->dataobject);
        dimensionflags |= 2;
        break;
      case 'R':
        *err = getDimension(&hrtf->R, &dir->dataobject);
        dimensionflags |= 4;
        break;
      case 'E':
        *err = getDimension(&hrtf->E, &dir->dataobject);
        dimensionflags |= 8;
        break;
      case 'N':
        *err = getDimension(&hrtf->N, &dir->dataobject);
        dimensionflags |= 0x10;
        break;
      case 'M':
        *err = getDimension(&hrtf->M, &dir->dataobject);
        dimensionflags |= 0x20;
        break;
      case 'S':
        break; /* be graceful, some issues with API version 0.4.4 */
      default:
        mylog("UNKNOWN SOFA VARIABLE %s", dir->dataobject.name);
        goto error;
      }
      if (*err)
        goto error;
    }
    dir = dir->next;
  }

  if (dimensionflags != 0x3f || hrtf->_I != 1 || hrtf->C != 3) {
    mylog("dimensions are missing or wrong\n");
    goto error;
  }

  dir = reader->superblock.dataobject.directory;
  while (dir) {

    if (!dir->dataobject.name) {
      mylog("SOFA VARIABLE IS NULL.\n");
    } else if (!strcmp(dir->dataobject.name, "ListenerPosition")) {
      *err = getArray(&hrtf->ListenerPosition, &dir->dataobject);
    } else if (!strcmp(dir->dataobject.name, "ReceiverPosition")) {
      *err = getArray(&hrtf->ReceiverPosition, &dir->dataobject);
    } else if (!strcmp(dir->dataobject.name, "SourcePosition")) {
      *err = getArray(&hrtf->SourcePosition, &dir->dataobject);
    } else if (!strcmp(dir->dataobject.name, "EmitterPosition")) {
      *err = getArray(&hrtf->EmitterPosition, &dir->dataobject);
    } else if (!strcmp(dir->dataobject.name, "ListenerUp")) {
      *err = getArray(&hrtf->ListenerUp, &dir->dataobject);
    } else if (!strcmp(dir->dataobject.name, "ListenerView")) {
      *err = getArray(&hrtf->ListenerView, &dir->dataobject);
    } else if (!strcmp(dir->dataobject.name, "Data.IR")) {
      *err = getArray(&hrtf->DataIR, &dir->dataobject);
    } else if (!strcmp(dir->dataobject.name, "Data.SamplingRate")) {
      *err = getArray(&hrtf->DataSamplingRate, &dir->dataobject);
    } else if (!strcmp(dir->dataobject.name, "Data.Delay")) {
      *err = getArray(&hrtf->DataDelay, &dir->dataobject);
    } else if (!(dir->dataobject.name[0] && !dir->dataobject.name[1])) {
      *err = addUserDefinedVariable(hrtf, &dir->dataobject);
    }
    dir = dir->next;
  }

  return hrtf;

error:
  free(hrtf);
  if (!*err)
    *err = MYSOFA_INVALID_FORMAT;
  return NULL;
}

MYSOFA_EXPORT struct MYSOFA_HRTF *mysofa_load(const char *filename, int *err) {
  struct READER reader;
  struct MYSOFA_HRTF *hrtf = NULL;

  if (filename == NULL)
    filename = CMAKE_INSTALL_PREFIX "/share/libmysofa/default.sofa";

  if (strcmp(filename, "-"))
    reader.fhd = fopen(filename, "rb");
  else
    reader.fhd = stdin;

  if (!reader.fhd) {
    mylog("cannot open file %s\n", filename);
    *err = errno;
    return NULL;
  }
  reader.gcol = NULL;
  reader.all = NULL;
  reader.recursive_counter = 0;

  *err = superblockRead(&reader, &reader.superblock);

  if (!*err) {
    hrtf = getHrtf(&reader, err);
  }

  superblockFree(&reader, &reader.superblock);
  gcolFree(reader.gcol);
  if (strcmp(filename, "-"))
    fclose(reader.fhd);

  return hrtf;
}

MYSOFA_EXPORT void mysofa_free(struct MYSOFA_HRTF *hrtf) {
  if (!hrtf)
    return;

  while (hrtf->attributes) {
    struct MYSOFA_ATTRIBUTE *next = hrtf->attributes->next;
    free(hrtf->attributes->name);
    free(hrtf->attributes->value);
    free(hrtf->attributes);
    hrtf->attributes = next;
  }

  while (hrtf->variables) {
    struct MYSOFA_VARIABLE *next = hrtf->variables->next;
    free(hrtf->variables->name);
    arrayFree(hrtf->variables->value);
    free(hrtf->variables->value);
    free(hrtf->variables);
    hrtf->variables = next;
  }

  arrayFree(&hrtf->ListenerPosition);
  arrayFree(&hrtf->ReceiverPosition);
  arrayFree(&hrtf->SourcePosition);
  arrayFree(&hrtf->EmitterPosition);
  arrayFree(&hrtf->ListenerUp);
  arrayFree(&hrtf->ListenerView);
  arrayFree(&hrtf->DataIR);
  arrayFree(&hrtf->DataSamplingRate);
  arrayFree(&hrtf->DataDelay);
  free(hrtf);
}

MYSOFA_EXPORT void mysofa_getversion(int *major, int *minor, int *patch) {
  *major = CPACK_PACKAGE_VERSION_MAJOR;
  *minor = CPACK_PACKAGE_VERSION_MINOR;
  *patch = CPACK_PACKAGE_VERSION_PATCH;
}

/* ========================================================================== */
/*                                    EASY                                    */
/* ========================================================================== */

static struct MYSOFA_EASY *
mysofa_open_default(const char *filename, float samplerate, int *filterlength,
                    int *err, bool applyNorm, float neighbor_angle_step,
                    float neighbor_radius_step) {

  struct MYSOFA_EASY *easy = malloc(sizeof(struct MYSOFA_EASY));

  if (!easy) {
    *err = MYSOFA_NO_MEMORY;
    return NULL;
  }

  // set all values of struct to their default "0" (to avoid freeing unallocated
  // values in mysofa_free)
  *easy = (struct MYSOFA_EASY){0};

  easy->hrtf = mysofa_load(filename, err);
  if (!easy->hrtf) {
    mysofa_close(easy);
    return NULL;
  }

  *err = mysofa_check(easy->hrtf);
  if (*err != MYSOFA_OK) {
    mysofa_close(easy);
    return NULL;
  }

  *err = mysofa_resample(easy->hrtf, samplerate);
  if (*err != MYSOFA_OK) {
    mysofa_close(easy);
    return NULL;
  }

  if (applyNorm) {
    mysofa_loudness(easy->hrtf);
  }

  /* does not sound well:
   mysofa_minphase(easy->hrtf,0.01);
   */

  mysofa_tocartesian(easy->hrtf);

  easy->lookup = mysofa_lookup_init(easy->hrtf);
  if (easy->lookup == NULL) {
    *err = MYSOFA_INTERNAL_ERROR;
    mysofa_close(easy);
    return NULL;
  }

  easy->neighborhood = mysofa_neighborhood_init_withstepdefine(
      easy->hrtf, easy->lookup, neighbor_angle_step, neighbor_radius_step);

  *filterlength = easy->hrtf->N;

  easy->fir = malloc(easy->hrtf->N * easy->hrtf->R * sizeof(float));
  assert(easy->fir);

  return easy;
}

MYSOFA_EXPORT struct MYSOFA_EASY *mysofa_open(const char *filename,
                                              float samplerate,
                                              int *filterlength, int *err) {
  return mysofa_open_default(filename, samplerate, filterlength, err, true,
                             MYSOFA_DEFAULT_NEIGH_STEP_ANGLE,
                             MYSOFA_DEFAULT_NEIGH_STEP_RADIUS);
}

MYSOFA_EXPORT struct MYSOFA_EASY *mysofa_open_no_norm(const char *filename,
                                                      float samplerate,
                                                      int *filterlength,
                                                      int *err) {
  return mysofa_open_default(filename, samplerate, filterlength, err, false,
                             MYSOFA_DEFAULT_NEIGH_STEP_ANGLE,
                             MYSOFA_DEFAULT_NEIGH_STEP_RADIUS);
}

MYSOFA_EXPORT struct MYSOFA_EASY *
mysofa_open_advanced(const char *filename, float samplerate, int *filterlength,
                     int *err, bool norm, float neighbor_angle_step,
                     float neighbor_radius_step) {
  return mysofa_open_default(filename, samplerate, filterlength, err, norm,
                             neighbor_angle_step, neighbor_radius_step);
}

MYSOFA_EXPORT struct MYSOFA_EASY *mysofa_open_cached(const char *filename,
                                                     float samplerate,
                                                     int *filterlength,
                                                     int *err) {
  struct MYSOFA_EASY *res = mysofa_cache_lookup(filename, samplerate);
  if (res) {
    *filterlength = res->hrtf->N;
    return res;
  }
  res = mysofa_open(filename, samplerate, filterlength, err);
  if (res) {
    res = mysofa_cache_store(res, filename, samplerate);
  }
  return res;
}

MYSOFA_EXPORT void mysofa_getfilter_short(struct MYSOFA_EASY *easy, float x,
                                          float y, float z, short *IRleft,
                                          short *IRright, int *delayLeft,
                                          int *delayRight) {
  float c[3];
  float delays[2];
  float *fl;
  float *fr;
  int nearest;
  int *neighbors;
  unsigned int i;

  c[0] = x;
  c[1] = y;
  c[2] = z;
  nearest = mysofa_lookup(easy->lookup, c);
  assert(nearest >= 0);
  neighbors = mysofa_neighborhood(easy->neighborhood, nearest);

  mysofa_interpolate(easy->hrtf, c, nearest, neighbors, easy->fir, delays);

  *delayLeft = delays[0] * easy->hrtf->DataSamplingRate.values[0];
  *delayRight = delays[1] * easy->hrtf->DataSamplingRate.values[0];

  fl = easy->fir;
  fr = easy->fir + easy->hrtf->N;
  for (i = easy->hrtf->N; i > 0; i--) {
    *IRleft++ = (short)(*fl++ * 32767.);
    *IRright++ = (short)(*fr++ * 32767.);
  }
}

MYSOFA_EXPORT void mysofa_getfilter_float_advanced(
    struct MYSOFA_EASY *easy, float x, float y, float z, float *IRleft,
    float *IRright, float *delayLeft, float *delayRight, bool interpolate) {
  float c[3];
  float delays[2];
  float *fl;
  float *fr;
  int nearest;
  int *neighbors;
  int i;

  c[0] = x;
  c[1] = y;
  c[2] = z;
  nearest = mysofa_lookup(easy->lookup, c);
  assert(nearest >= 0);
  neighbors = mysofa_neighborhood(easy->neighborhood, nearest);

  // bypass interpolate by forcing current cooordinates to nearest's
  if (!interpolate) {
    memcpy(c, easy->hrtf->SourcePosition.values + nearest * easy->hrtf->C,
           sizeof(float) * easy->hrtf->C);
  }

  float *res =
      mysofa_interpolate(easy->hrtf, c, nearest, neighbors, easy->fir, delays);

  *delayLeft = delays[0];
  *delayRight = delays[1];

  fl = res;
  fr = res + easy->hrtf->N;
  for (i = easy->hrtf->N; i > 0; i--) {
    *IRleft++ = *fl++;
    *IRright++ = *fr++;
  }
}

MYSOFA_EXPORT void mysofa_getfilter_float(struct MYSOFA_EASY *easy, float x,
                                          float y, float z, float *IRleft,
                                          float *IRright, float *delayLeft,
                                          float *delayRight) {
  mysofa_getfilter_float_advanced(easy, x, y, z, IRleft, IRright, delayLeft,
                                  delayRight, true);
}

MYSOFA_EXPORT void
mysofa_getfilter_float_nointerp(struct MYSOFA_EASY *easy, float x, float y,
                                float z, float *IRleft, float *IRright,
                                float *delayLeft, float *delayRight) {
  mysofa_getfilter_float_advanced(easy, x, y, z, IRleft, IRright, delayLeft,
                                  delayRight, false);
}

MYSOFA_EXPORT void mysofa_close(struct MYSOFA_EASY *easy) {
  if (easy) {
    if (easy->fir)
      free(easy->fir);
    if (easy->neighborhood)
      mysofa_neighborhood_free(easy->neighborhood);
    if (easy->lookup)
      mysofa_lookup_free(easy->lookup);
    if (easy->hrtf)
      mysofa_free(easy->hrtf);
    free(easy);
  }
}

MYSOFA_EXPORT void mysofa_close_cached(struct MYSOFA_EASY *easy) {
  mysofa_cache_release(easy);
}

/* ========================================================================== */
/*                                   LOOKUP                                   */
/* ========================================================================== */

MYSOFA_EXPORT struct MYSOFA_LOOKUP *
mysofa_lookup_init(struct MYSOFA_HRTF *hrtf) {
  int i;
  struct MYSOFA_LOOKUP *lookup;

  /*
   * alloc memory structure
   */
  if (!verifyAttribute(hrtf->SourcePosition.attributes, "Type", "cartesian"))
    return NULL;

  lookup = malloc(sizeof(struct MYSOFA_LOOKUP));
  if (!lookup)
    return NULL;

  /*
   * find smallest and largest phi, theta, and radius (to reduce neighbors table
   * init)
   */
  float *origin;
  origin = malloc(sizeof(float) * hrtf->C);
  lookup->phi_min = FLT_MAX;
  lookup->phi_max = FLT_MIN;
  lookup->theta_min = FLT_MAX;
  lookup->theta_max = FLT_MIN;
  lookup->radius_min = FLT_MAX;
  lookup->radius_max = FLT_MIN;
  for (i = 0; i < (int)hrtf->M; i++) {
    memcpy(origin, hrtf->SourcePosition.values + i * hrtf->C,
           sizeof(float) * hrtf->C);
    convertCartesianToSpherical(origin, hrtf->C);
    if (origin[0] < lookup->phi_min) {
      lookup->phi_min = origin[0];
    }
    if (origin[0] > lookup->phi_max) {
      lookup->phi_max = origin[0];
    }
    if (origin[1] < lookup->theta_min) {
      lookup->theta_min = origin[1];
    }
    if (origin[1] > lookup->theta_max) {
      lookup->theta_max = origin[1];
    }
    if (origin[2] < lookup->radius_min) {
      lookup->radius_min = origin[2];
    }
    if (origin[2] > lookup->radius_max) {
      lookup->radius_max = origin[2];
    }
  }
  free(origin);

  /*
   * Allocate kd tree
   */
  lookup->kdtree = kd_create();
  if (!lookup->kdtree) {
    free(lookup);
    return NULL;
  }

  /*
   * add coordinates to the tree
   */
  for (i = 0; i < (int)hrtf->M; i++) {
    float *f = hrtf->SourcePosition.values + i * hrtf->C;
    kd_insert((struct kdtree *)lookup->kdtree, f, (void *)(intptr_t)i);
  }

  return lookup;
}

/*
 * looks for a filter that is similar to the given Cartesian coordinate
 * BE AWARE: The coordinate vector will be normalized if required
 * A return value of -1 = MYSOFA_INTERNAL_ERROR indicates an error
 */
MYSOFA_EXPORT int mysofa_lookup(struct MYSOFA_LOOKUP *lookup,
                                float *coordinate) {

  int index;
  void *res;
  int success;
  float r = radius(coordinate);
  if (r > lookup->radius_max) {
    r = lookup->radius_max / r;
    coordinate[0] *= r;
    coordinate[1] *= r;
    coordinate[2] *= r;
  } else if (r < lookup->radius_min) {
    r = lookup->radius_min / r;
    coordinate[0] *= r;
    coordinate[1] *= r;
    coordinate[2] *= r;
  }

  success = kd_nearest((struct kdtree *)lookup->kdtree, coordinate, &res);
  if (success != 0) {
    return MYSOFA_INTERNAL_ERROR;
  }
  index = (int)((uintptr_t)res);
  return index;
}

MYSOFA_EXPORT void mysofa_lookup_free(struct MYSOFA_LOOKUP *lookup) {
  if (lookup) {
    kd_free((struct kdtree *)lookup->kdtree);
    free(lookup);
  }
}

/* ========================================================================== */
/*                                   CACHE                                    */
/* ========================================================================== */

static struct MYSOFA_CACHE_ENTRY {
  struct MYSOFA_CACHE_ENTRY *next;
  struct MYSOFA_EASY *easy;
  char *filename;
  float samplerate;
  int count;
} * cache;

static int compare_filenames(const char *a, const char *b) {
  if (a == NULL && b == NULL)
    return 0;
  if (a == NULL)
    return -1;
  else if (b == NULL)
    return 1;
  return strcmp(a, b);
}

MYSOFA_EXPORT struct MYSOFA_EASY *mysofa_cache_lookup(const char *filename,
                                                      float samplerate) {
  struct MYSOFA_CACHE_ENTRY *p;
  struct MYSOFA_EASY *res = NULL;

  p = cache;

  while (p) {
    if (samplerate == p->samplerate &&
        !compare_filenames(filename, p->filename)) {
      res = p->easy;
      p->count++;
      break;
    }
    p = p->next;
  }

  return res;
}

MYSOFA_EXPORT struct MYSOFA_EASY *mysofa_cache_store(struct MYSOFA_EASY *easy,
                                                     const char *filename,
                                                     float samplerate) {
  struct MYSOFA_CACHE_ENTRY *p;

  assert(easy);

  p = cache;

  while (p) {
    if (samplerate == p->samplerate &&
        !compare_filenames(filename, p->filename)) {
      mysofa_close(easy);
      return p->easy;
    }
    p = p->next;
  }

  p = malloc(sizeof(struct MYSOFA_CACHE_ENTRY));
  if (p == NULL) {
    return NULL;
  }
  p->next = cache;
  p->samplerate = samplerate;
  p->filename = NULL;
  if (filename != NULL) {
    p->filename = mysofa_strdup(filename);
    if (p->filename == NULL) {
      free(p);
      return NULL;
    }
  }
  p->easy = easy;
  p->count = 1;
  cache = p;
  return easy;
}

MYSOFA_EXPORT void mysofa_cache_release(struct MYSOFA_EASY *easy) {
  struct MYSOFA_CACHE_ENTRY **p;
  int count;

  assert(easy);
  assert(cache);

  p = &cache;

  for (count = 0;; count++) {
    if ((*p)->easy == easy)
      break;
    p = &((*p)->next);
    assert(*p);
  }

  if ((*p)->count == 1 && (count > 0 || (*p)->next != NULL)) {
    struct MYSOFA_CACHE_ENTRY *gone = *p;
    free(gone->filename);
    mysofa_close(easy);
    *p = (*p)->next;
    free(gone);
  } else {
    (*p)->count--;
  }
}

MYSOFA_EXPORT void mysofa_cache_release_all() {
  struct MYSOFA_CACHE_ENTRY *p;

  p = cache;
  while (p) {
    struct MYSOFA_CACHE_ENTRY *gone = p;
    p = p->next;
    free(gone->filename);
    free(gone->easy);
    free(gone);
  }
  cache = NULL;
}

#else
extern int to_avoid_iso_compiler_warning_when_there_are_no_symbols;
#endif /* SAF_ENABLE_SOFA_READER_MODULE */
