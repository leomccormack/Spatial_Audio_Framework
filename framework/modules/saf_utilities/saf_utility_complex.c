/*
 * Copyright 2016-2018 Leo McCormack
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
 * @file saf_utility_complex.c
 * @ingroup Utilities
 * @brief Contains wrappers for handling complex numbers across both
 *        C99-compliant compilers and Microsoft Visual Compiler (MSVC)
 
 * @author Leo McCormack
 * @date 11.07.2016
 * @license ISC
 */

#include "saf_utility_complex.h"

#if __STDC_VERSION__ >= 199901L
  /* for C99 (and above) compliant compilers */

  /*
  Single-Precision Complex Operations
  */
  inline float_complex cmplxf(float re, float im) { return re + im * I; }
  inline float_complex ccaddf(float_complex x, float_complex y) { return x + y; }
  inline float_complex craddf(float_complex x, float y) { return x + y; }
  inline float_complex ccsubf(float_complex x, float_complex y) { return x - y; }
  inline float_complex crsubf(float_complex x, float y) { return x - y; }
  inline float_complex ccmulf(float_complex x, float_complex y) { return x * y; }
  inline float_complex cccmulf(float_complex x, float_complex y, float_complex z) { return x * y * z; }
  inline float_complex crmulf(float_complex x, float y) { return x * y; }
  inline float_complex ccdivf(float_complex x, float_complex y) { return x / y; }
  inline float_complex crdivf(float_complex x, float y) { return x / y; }

  /*
  Double-Precision Complex Operations
  */
  inline double_complex cmplx(double re, double im) { return re + im * I; }
  inline double_complex ccadd(double_complex x, double_complex y) { return x + y; }
  inline double_complex cradd(double_complex x, double y) { return x + y; }
  inline double_complex ccsub(double_complex x, double_complex y) { return x - y; }
  inline double_complex crsub(double_complex x, double y) { return x - y; }
  inline double_complex ccmul(double_complex x, double_complex y) { return x * y; }
  inline double_complex cccmul(double_complex x, double_complex y, double_complex z) { return x * y * z; }
  inline double_complex crmul(double_complex x, double y) { return x * y; }
  inline double_complex ccdiv(double_complex x, double_complex y) { return x / y; }
  inline double_complex crdiv(double_complex x, double y) { return x / y; }

#endif
