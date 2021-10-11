//
//  FFTHelper.h
//  just_audio
//
//  Created by Eittipat on 10/10/2564 BE.
//

#import <Foundation/Foundation.h>
#include <Accelerate/Accelerate.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>


#define CLAMP(a, b, c) MIN(MAX((a), (b)), (c))
#define SCALE(x, xa, xb, ya, yb) ya + ((x - xa) * (yb - ya) / (xb - xa));


@interface FFTHelper : NSObject

-(void) initialFFT:(int)length;
-(void) performFFT:(float*)audioInput :(float*)fftOutput;
-(void) destoryFFT;

@end
