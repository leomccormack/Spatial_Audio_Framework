/*
 Copyright 2016-2018 Leo McCormack
 
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
 *     saf_complex.h
 * Description:
 *     Contains a collection of useful memory allocation functions and cross-platform
 *     complex number wrappers. Optimised linear algebra routines utilising BLAS and LAPACK
 *     are also included.
 * Dependencies:
 *     Windows users only: Intel's MKL must be installed, which can be freely aquired via:
 *     https://software.intel.com/en-us/articles/free-ipsxe-tools-and-libraries
 *     Mac users only: saf_utilities will utilise Apple's Accelerate library.
 * Author, date created:
 *     Leo McCormack, 11.07.2016
 */

#ifndef SAF_COMPLEX_H_INCLUDED
#define SAF_COMPLEX_H_INCLUDED

#ifdef __cplusplus
extern "C" {
#endif
    
#if __STDC_VERSION__ >= 199901L
  /* for fully C99+ compliant compilers */

    #include <complex.h>
    typedef float _Complex float_complex;
    typedef double _Complex double_complex;

  /*
  Single-Precision Complex Operations
  */
  float_complex cmplxf(float re, float im);
  float_complex ccaddf(float_complex x, float_complex y);
  float_complex craddf(float_complex x, float y);
  float_complex ccsubf(float_complex x, float_complex y);
  float_complex crsubf(float_complex x, float y);
  float_complex ccmulf(float_complex x, float_complex y);
  float_complex cccmulf(float_complex x, float_complex y, float_complex z);
  float_complex crmulf(float_complex x, float y);
  float_complex ccdivf(float_complex x, float_complex y);
  float_complex crdivf(float_complex x, float y);

  /*
  Double-Precision Complex Operations
  */
  double_complex cmplx(double re, double im);
  double_complex ccadd(double_complex x, double_complex y);
  double_complex cradd(double_complex x, double y);
  double_complex ccsub(double_complex x, double_complex y);
  double_complex crsub(double_complex x, double y);
  double_complex ccmul(double_complex x, double_complex y);
  double_complex cccmul(double_complex x, double_complex y, double_complex z);
  double_complex crmul(double_complex x, double y);
  double_complex ccdiv(double_complex x, double_complex y);
  double_complex crdiv(double_complex x, double y);


#elif _MSC_VER >= 1900
  /* for Microsoft's ancient C compiler */
    
    #if _MSC_VER >= 1900
      #include <complex.h>
      typedef _Fcomplex float_complex;
      typedef _Dcomplex double_complex;
    #else
    
    #endif

  /*
  Single-Precision Complex Operations
  */
  inline float_complex cmplxf(float re, float im) {
    float_complex z;
    z._Val[0] = re; z._Val[1] = im;
    return z;
  }

  inline float_complex ccaddf(float_complex x, float_complex y) {
    float_complex z;
    z._Val[0] = x._Val[0] + y._Val[0]; z._Val[1] = x._Val[1] + y._Val[1];
    return z;
  }

  inline float_complex craddf(float_complex x, float y) {
    float_complex z;
    z._Val[0] = x._Val[0] + y; z._Val[1] = x._Val[1];
    return z;
  }

  inline float_complex ccsubf(float_complex x, float_complex y) {
    float_complex z;
    z._Val[0] = x._Val[0] - y._Val[0]; z._Val[1] = x._Val[1] - y._Val[1];
    return z;
  }

  inline float_complex crsubf(float_complex x, float y) {
    float_complex z;
    z._Val[0] = x._Val[0] - y; z._Val[1] = x._Val[1];
    return z;
  }

  inline float_complex ccmulf(float_complex x, float_complex y) {
    return _FCmulcc(x, y);
  }

  inline float_complex cccmulf(float_complex x, float_complex y, float_complex z) {
    return _FCmulcc(_FCmulcc(x, y), z);
  }

  inline float_complex crmulf(float_complex x, float y) {
    return _FCmulcr(x, y);
  }

  inline float_complex ccdivf(float_complex x, float_complex y) {
    float_complex z;
    z._Val[0] = (x._Val[0] * y._Val[0] + x._Val[1] * y._Val[1]) / (y._Val[0] * y._Val[0] + y._Val[1] * y._Val[1] );
    z._Val[1] = (x._Val[1] * y._Val[0] - x._Val[0] * y._Val[1]) / (y._Val[0] * y._Val[0] + y._Val[1] * y._Val[1] );
    return z;
  }

  inline float_complex crdivf(float_complex x, float y) {
    float_complex z;
    z._Val[0] = x._Val[0] / y;
    z._Val[1] = x._Val[1] / y;
    return z;
  }

  /* 
  Double-Precision Complex Operations 
  */
  inline double_complex cmplx(double re, double im) {
    double_complex z;
    z._Val[0] = re; z._Val[1] = im;
    return z;
  }

  inline double_complex ccadd(double_complex x, double_complex y) {
    double_complex z;
    z._Val[0] = x._Val[0] + y._Val[0]; z._Val[1] = x._Val[1] + y._Val[1];
    return z;
  }

  inline double_complex cradd(double_complex x, double y) {
    double_complex z;
    z._Val[0] = x._Val[0] + y; z._Val[1] = x._Val[1];
    return z;
  }

  inline double_complex ccsub(double_complex x, double_complex y) {
    double_complex z;
    z._Val[0] = x._Val[0] - y._Val[0]; z._Val[1] = x._Val[1] - y._Val[1];
    return z;
  }

  inline double_complex crsub(double_complex x, double y) {
    double_complex z;
    z._Val[0] = x._Val[0] - y; z._Val[1] = x._Val[1];
    return z;
  }

  inline double_complex ccmul(double_complex x, double_complex y) {
    return _Cmulcc(x, y);
  }

  inline double_complex cccmul(double_complex x, double_complex y, double_complex z) {
    return _Cmulcc(_Cmulcc(x, y), z);
  }

  inline double_complex crmul(double_complex x, double y) {
    return _Cmulcr(x, y);
  }

  inline double_complex ccdiv(double_complex x, double_complex y) {
    double_complex z;
    z._Val[0] = (x._Val[0] * y._Val[0] + x._Val[1] * y._Val[1]) / (y._Val[0] * y._Val[0] + y._Val[1] * y._Val[1]);
    z._Val[1] = (x._Val[1] * y._Val[0] - x._Val[0] * y._Val[1]) / (y._Val[0] * y._Val[0] + y._Val[1] * y._Val[1]);
    return z;
  }
    
  inline double_complex crdiv(double_complex x, double y) {
      double_complex z;
      z._Val[0] = x._Val[0] / y;
      z._Val[1] = x._Val[1] / y;
      return z;
  }

#endif

    
#ifdef __cplusplus
}/* extern "C" */
#endif


#endif /* SAF_COMPLEX_H_INCLUDED */
