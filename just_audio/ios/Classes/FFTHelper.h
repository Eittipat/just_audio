//
//  FFTHelper.h
//  just_audio
//
//  Created by Eittipat on 10/10/2564 BE.
//

#ifndef FFTHelper_h
#define FFTHelper_h

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define MIN(a, b) (((a) < (b)) ? (a) : (b))
#define MAX(a, b) (((a) > (b)) ? (a) : (b))
#define CLAMP(a, b, c) MIN(MAX((a), (b)), (c))
#define SCALE(x, xa, xb, ya, yb) ya + ((x - xa) * (yb - ya) / (xb - xa));

void irfft(int n, int *a);
void rfft(int n, int *a);
void ifft(int n, int *a, int *b);
void fft(int n, int *a, int *b);
void rearrange(int n, int *a, int *b, int *c);

#endif /* FFTHelper_h */
