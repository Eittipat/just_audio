//
//  FFTHelper.c
//  just_audio
//
//  Created by Eittipat on 10/10/2564 BE.
//
//  Original source codes
//  - https://developer.android.com/reference/android/media/audiofx/Visualizer#getFft(byte[])
//  - http://www.myuiviews.com/2016/03/04/visualizing-audio-frequency-spectrum-on-ios-via-accelerate-vdsp-fast-fourier-transform.html

#include "FFTHelper.h"

@interface FFTHelper()
{
    FFTSetup _fftSetup;
    int _fftBuffSize;
    int _fftBuffSizeLog2;
    float _fftNormFactor;
}
@end

@implementation FFTHelper


-(void) initialFFT:(int)length {
    _fftBuffSize = length;
    _fftNormFactor = 0.5;
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
    // At this point, result is equal to MATLAB
    
    // rearrange to match with android
    fftOutput[0] = outReal[0];
    fftOutput[1] = outReal[nOver2];
    int rIndex = 1;
    int iIndex = 1;
    for (int i = 2; i < _fftBuffSize - 1; i += 2)
    {
        fftOutput[i] = outReal[rIndex++];
        fftOutput[i + 1] = outImaginary[iIndex++];
    }
}

-(void) destoryFFT {
    vDSP_destroy_fftsetup(_fftSetup);
}

    
@end



