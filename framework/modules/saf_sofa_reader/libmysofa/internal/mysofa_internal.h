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
 * tools.h
 *
 *  Created on: 13.01.2017
 *      Author: hoene
 */

#ifndef MYSOFA_SRC_TOOLS_H_
#define MYSOFA_SRC_TOOLS_H_

#ifdef __cplusplus
extern "C" {
#endif

#if defined(SAF_ENABLE_SOFA_READER_MODULE)

#include "../mysofa.h"
#include <math.h>
#include <stdlib.h>

int changeAttribute(struct MYSOFA_ATTRIBUTE *attr, char *name, char *value,
                    char *newvalue);
int verifyAttribute(struct MYSOFA_ATTRIBUTE *attr, char *name, char *value);
char *getAttribute(struct MYSOFA_ATTRIBUTE *attr, char *name);

void convertCartesianToSpherical(float *values, int elements);
void convertSphericalToCartesian(float *values, int elements);

#define fequals(a, b) (fabs(a - b) < 0.00001)

float radius(float *cartesian);

#define distance(cartesian1, cartesian2)                                       \
  (sqrtf(powf((cartesian1)[0] - (cartesian2)[0], 2.f) +                        \
         powf((cartesian1)[1] - (cartesian2)[1], 2.f) +                        \
         powf((cartesian1)[2] - (cartesian2)[2], 2.f)))

void copyToFloat(float *out, float *in, int size);
void copyFromFloat(float *out, float *in, int size);

void copyArrayWeighted(float *dst, float *src, int size, float w);
void addArrayWeighted(float *dst, float *src, int size, float w);
void scaleArray(float *dst, int size, float w);
float loudness(float *in, int size);

void nsearch(const void *key, const char *base, size_t num, size_t size,
             int (*cmp)(const void *key, const void *elt), int *lower,
             int *higher);


#endif /* SAF_ENABLE_SOFA_READER_MODULE */

#ifdef __cplusplus
}
#endif

#endif /* MYSOFA_SRC_TOOLS_H_ */
