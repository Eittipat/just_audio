/*

MIT License

Copyright (c) 2018 fukuroda

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

*/

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
