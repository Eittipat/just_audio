//
//  FFTHelper.c
//  just_audio
//
//  Created by Eittipat on 10/10/2564 BE.
//
//  Original source codes
//  - https://github.com/fukuroder/intfft/blob/master/src/intfft.cpp
//  - https://developer.android.com/reference/android/media/audiofx/Visualizer#getFft(byte[])

#include "FFTHelper.h"

int lift(int swap, int real, int image, float c, float s)
{
    if (s == 0.0f)
    {
        return swap == 0 ? real : image;
    }

    if (s > c)
    {
        if (s > -c)
        { // (0.25pi, 0.75pi)
            // swap
            int t = real;
            real = image;
            image = t;
            real += (int)(image * (s - 1) / c);
            image += (int)(real * c);
            real += (int)(image * (s - 1) / c);
            real = -real;
        }
        else
        { // (0.75pi, 1.25pi)
            image = -image;
            real += (int)(image * (-c - 1) / s);
            image += (int)(real * s);
            real += (int)(image * (-c - 1) / s);
            real = -real;
        }
    }
    else
    {
        if (s < -c)
        { // (-0.75pi, -0.25pi)
            real += (int)(image * (-s - 1) / c);
            image += (int)(real * c);
            real += (int)(image * (-s - 1) / c);
            int t = real;
            real = image;
            image = -t;
        }
        else
        { // (-0.25pi, 0.25pi)
            real += (int)(image * (c - 1) / s);
            image += (int)(real * s);
            real += (int)(image * (c - 1) / s);
        }
    }
    return swap == 0 ? real : image;
}

int ilift(int swap, int real, int image, float c, float s)
{

    if (s == 0.0f)
    {
        return swap == 0 ? real : image;
    }

    if (s > c)
    {
        if (s > -c)
        { // (0.25pi, 0.75pi)
            real = -real;
            real -= (int)(image * (s - 1) / c);
            image -= (int)(real * c);
            real -= (int)(image * (s - 1) / c);
            int t = real;
            real = image;
            image = t;
        }
        else
        { // (0.75pi, 1.25pi)
            real = -real;
            real -= (int)(image * (-c - 1) / s);
            image -= (int)(real * s);
            real -= (int)(image * (-c - 1) / s);
            image = -image;
        }
    }
    else
    {
        if (s < -c)
        { // (-0.75pi, -0.25pi)
            int t = real;
            real = -image;
            image = t;
            real -= (int)(image * (-s - 1) / c);
            image -= (int)(real * c);
            real -= (int)(image * (-s - 1) / c);
        }
        else
        { // (-0.25pi, 0.25pi)
            real -= (int)(image * (c - 1) / s);
            image -= (int)(real * s);
            real -= (int)(image * (c - 1) / s);
        }
    }
    return swap == 0 ? real : image;
}

void rfft(int n, int *a)
{
    // scrambler
    for (int i = 0, j = 1; j < n - 1; j++)
    {
        for (int k = n >> 1; k > (i ^= k); k >>= 1)
            ;
        if (j < i)
        {
            int xr = a[j];
            a[j] = a[i];
            a[i] = xr;
        }
    }

    for (int mh = 1, m; (m = mh << 1) <= n; mh = m)
    {
        double theta = -2 * M_PI / m;
        int mq = mh >> 1;

        // real to real butterflies (W == 1)
        for (int jr = 0; jr < n; jr += m)
        {
            int kr = jr + mh;
            int xr = a[kr];
            a[kr] = a[jr] - xr;
            a[jr] += xr;
        }

        // complex to complex butterflies (W != 1)
        for (int i = 1; i < mq; i++)
        {
            double wr = cos(theta * i);
            double wi = sin(theta * i);
            for (int j = 0; j < n; j += m)
            {
                int jr = j + i;
                int ji = j + mh - i;
                int kr = j + mh + i;
                int ki = j + m - i;
                int xr = lift(0, a[kr], -a[ki], wr, wi);
                int xi = lift(1, a[kr], -a[ki], wr, wi);
                a[kr] = -a[ji] - xi;
                a[ki] = a[ji] - xi;
                a[ji] = a[jr] - xr;
                a[jr] = a[jr] + xr;
            }
        }
    }
}

void irfft(int n, int *a)
{
    for (int m = n, mh; (mh = m >> 1) >= 1; m = mh)
    {
        float theta = -2 * M_PI / m;
        int mq = mh >> 1;

        // complex to complex butterflies (W != 1)
        for (int i = 1; i < mq; i++)
        {
            float wr = cos(theta * i);
            float wi = sin(theta * i);
            for (int j = 0; j < n; j += m)
            {
                int jr = j + i;
                int ji = j + mh - i;
                int kr = j + mh + i;
                int ki = j + m - i;
                int xr = -(a[ji] - a[jr]) / 2;
                int xi = -(a[kr] + a[ki]) / 2;
                a[jr] = (a[ji] + a[jr]) / 2;
                a[ji] = -(a[kr] - a[ki]) / 2;
                a[kr] = ilift(0, xr, xi, wr, wi);
                a[ki] = ilift(1, xr, xi, wr, wi);
                a[ki] = -a[ki];
            }
        }

        // real to real butterflies (W == 1)
        for (int jr = 0; jr < n; jr += m)
        {
            int kr = jr + mh;
            int xr = a[jr];
            a[jr] = (xr + a[kr]) / 2;
            a[kr] = (xr - a[kr]) / 2;
        }
    }

    // unscrambler
    for (int i = 0, j = 1; j < n - 1; j++)
    {
        for (int k = n >> 1; k > (i ^= k); k >>= 1)
            ;
        if (j < i)
        {
            int xr = a[j];
            a[j] = a[i];
            a[i] = xr;
        }
    }
}

void fft(int n, int *ar, int *ai)
{
    // L shaped butterflies
    for (int m = n; m > 2; m >>= 1)
    {
        float theta = -2 * M_PI / m;
        int mq = m >> 2;
        for (int i = 0; i < mq; i++)
        {
            float s1 = sin(theta * i);
            float c1 = cos(theta * i);
            float s3 = sin(theta * 3 * i);
            float c3 = cos(theta * 3 * i);
            for (int k = m; k <= n; k <<= 2)
            {
                for (int j0 = k - m + i; j0 < n; j0 += 2 * k)
                {
                    int j1 = j0 + mq;
                    int j2 = j1 + mq;
                    int j3 = j2 + mq;
                    int x1r = ar[j0] - ar[j2];
                    int x1i = ai[j0] - ai[j2];
                    ar[j0] += ar[j2];
                    ai[j0] += ai[j2];
                    int x3r = ar[j1] - ar[j3];
                    int x3i = ai[j1] - ai[j3];
                    ar[j1] += ar[j3];
                    ai[j1] += ai[j3];
                    //std::tie(ar[j2], ai[j2]) = lift(x1r + x3i, x1i - x3r, c1, s1);
                    //std::tie(ar[j3], ai[j3]) = lift(x1r - x3i, x1i + x3r, c3, s3);
                    ar[j2] = lift(0, x1r + x3i, x1i - x3r, c1, s1);
                    ai[j2] = lift(1, x1r + x3i, x1i - x3r, c1, s1);
                    ar[j3] = lift(0, x1r - x3i, x1i + x3r, c3, s3);
                    ai[j3] = lift(1, x1r - x3i, x1i + x3r, c3, s3);
                }
            }
        }
    }

    // radix 2 butterflies
    for (int k = 2; k <= n; k <<= 2)
    {
        for (int j = k - 2; j < n; j += 2 * k)
        {
            int x0r = ar[j] - ar[j + 1];
            int x0i = ai[j] - ai[j + 1];
            ar[j] += ar[j + 1];
            ai[j] += ai[j + 1];
            ar[j + 1] = x0r;
            ai[j + 1] = x0i;
        }
    }
    // unscrambler
    for (int i = 0, j = 1; j < n - 1; j++)
    {
        for (int k = n >> 1; k > (i ^= k); k >>= 1)
            ;
        if (j < i)
        {
            int x0r = ar[j];
            int x0i = ai[j];
            ar[j] = ar[i];
            ai[j] = ai[i];
            ar[i] = x0r;
            ai[i] = x0i;
        }
    }
}

void ifft(int n, int *ar, int *ai)
{
    // scrambler
    for (int i = 0, j = 1; j < n - 1; j++)
    {
        for (int k = n >> 1; k > (i ^= k); k >>= 1)
            ;
        if (j < i)
        {
            int x0r = ar[j];
            int x0i = ai[j];
            ar[j] = ar[i];
            ai[j] = ai[i];
            ar[i] = x0r;
            ai[i] = x0i;
        }
    }

    // radix 2 butterflies
    for (int k = 2; k <= n; k <<= 2)
    {
        for (int j = k - 2; j < n; j += 2 * k)
        {
            int x0r = ar[j];
            int x0i = ai[j];
            ar[j] = (x0r + ar[j + 1]) / 2;
            ai[j] = (x0i + ai[j + 1]) / 2;
            ar[j + 1] = x0r - ar[j];
            ai[j + 1] = x0i - ai[j];
        }
    }

    // L shaped butterflies
    for (int m = 4; m <= n; m <<= 1)
    {
        double theta = -2 * M_PI / m;
        int mq = m >> 2;
        for (int i = 0; i < mq; i++)
        {
            float s1 = sin(theta * i);
            float c1 = cos(theta * i);
            float s3 = sin(theta * 3 * i);
            float c3 = cos(theta * 3 * i);
            for (int k = m; k <= n; k <<= 2)
            {
                for (int j0 = k - m + i; j0 < n; j0 += 2 * k)
                {
                    int j1 = j0 + mq;
                    int j2 = j1 + mq;
                    int j3 = j2 + mq;
                    int x0r = ar[j0];
                    int x0i = ai[j0];
                    int x1r = ar[j1];
                    int x1i = ai[j1];
                    //auto [x2r, x2i] = ilift(ar[j2], ai[j2], c1, s1);
                    //auto [x3r, x3i] = ilift(ar[j3], ai[j3], c3, s3);
                    int x2r = ilift(0, ar[j2], ai[j2], c1, s1);
                    int x2i = ilift(1, ar[j2], ai[j2], c1, s1);
                    int x3r = ilift(0, ar[j3], ai[j3], c3, s3);
                    int x3i = ilift(1, ar[j3], ai[j3], c3, s3);
                    int x2r_ = (x2r + x3r) / 2;
                    int x2i_ = (x2i + x3i) / 2;
                    int x3r_ = -(x2i - x2i_);
                    int x3i_ = (x2r - x2r_);
                    ar[j0] = (x0r + x2r_) / 2;
                    ai[j0] = (x0i + x2i_) / 2;
                    ar[j1] = (x1r + x3r_) / 2;
                    ai[j1] = (x1i + x3i_) / 2;
                    ar[j2] = (x0r - ar[j0]);
                    ai[j2] = (x0i - ai[j0]);
                    ar[j3] = (x1r - ar[j1]);
                    ai[j3] = (x1i - ai[j1]);
                }
            }
        }
    }
}

// rearrange array to match with android
void rearrange(int n, int *out, int *real, int *image)
{
    int nOver2 = n / 2;
    out[0] = real[0];
    out[1] = real[nOver2];
    int rIndex = 1;
    int iIndex = 1;
    for (int i = 2; i < n - 1; i += 2)
    {
        out[i] = real[rIndex++];
        out[i + 1] = image[iIndex++];
    }
}

/*
//
//  SignalProcessing.m
//
//  Created by Eittipat on 7/10/2564 BE.
//
//  Original source code
//  - http://www.myuiviews.com/2016/03/04/visualizing-audio-frequency-spectrum-on-ios-via-accelerate-vdsp-fast-fourier-transform.html
//


#import "SignalProcessing.h"
#include <Accelerate/Accelerate.h>
#include "FFTHelper.h"

@interface SignalProcessing()
{
    FFTSetup _fftSetup;
    int _fftBuffSize;
    int _fftBuffSizeLog2;
    float _fftNormFactor;
}
@end

@implementation SignalProcessing


- (void) integerFFT:(int*) buffer :(int) length {
    rfft(length,buffer);
    magnitudes(length,buffer);
}

- (void) initFFT:(int)length {
    _fftBuffSize = length;
    _fftNormFactor = 1.0/( 2 * length);
    _fftBuffSizeLog2 = round(log2(length));
    _fftSetup = vDSP_create_fftsetup(_fftBuffSizeLog2, kFFTRadix2);
}

-(void) performFFT:(float*)audioInput :(float*)fftOutput {
    int nOver2 = _fftBuffSize/2;
    float outReal[nOver2];
    float outImaginary[nOver2];
    COMPLEX_SPLIT output = { .realp = outReal, .imagp = outImaginary };
    vDSP_ctoz((COMPLEX *)audioInput, 2, &output, 1, nOver2);
    vDSP_fft_zrip(_fftSetup, &output, 1, _fftBuffSizeLog2, FFT_FORWARD);
    vDSP_vsmul(output.realp, 1, &_fftNormFactor, output.realp, 1, nOver2);
    vDSP_vsmul(output.imagp, 1, &_fftNormFactor, output.imagp, 1, nOver2);
    vDSP_zvabs(&output, 1,fftOutput, 1, nOver2);
}

-(void) destoryFFT {
    vDSP_destroy_fftsetup(_fftSetup);
}

@end


*/