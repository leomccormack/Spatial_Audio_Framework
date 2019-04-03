/* FFT library adopted from http://www.kurims.kyoto-u.ac.jp/~ooura/fft.html 
 
Original copyright:
 Copyright(C) 1996-2001 Takuya OOURA
 email: ooura@mmm.t.u-tokyo.ac.jp
 download: http://momonga.t.u-tokyo.ac.jp/~ooura/fft.html
 You may use, copy, modify this code for any purpose and
 without fee. You may distribute this ORIGINAL package.
 */

#include <math.h>

void cdft(int n, int isgn, float *a, int *ip, float *w);
void rdft(int n, int isgn, float *a, int *ip, float *w);
void ddct(int n, int isgn, float *a, int *ip, float *w);
void ddst(int n, int isgn, float *a, int *ip, float *w);
void dfct(int n, float *a, float *t, int *ip, float *w);
void dfst(int n, float *a, float *t, int *ip, float *w);
void makewt(int nw, int *ip, float *w);
void makect(int nc, int *ip, float *c);
void bitrv2(int n, int *ip, float *a);
void bitrv2conj(int n, int *ip, float *a);
void cftfsub(int n, float *a, float *w);
void cftbsub(int n, float *a, float *w);
void cft1st(int n, float *a, float *w);
void cftmdl(int n, int l, float *a, float *w);
void rftfsub(int n, float *a, int nc, float *c);
void rftbsub(int n, float *a, int nc, float *c);
void dctsub(int n, float *a, int nc, float *c);
void dstsub(int n, float *a, int nc, float *c);

