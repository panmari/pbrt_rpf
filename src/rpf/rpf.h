
/*
    Copyright(c) 2012-2013 Tzu-Mao Li
    All rights reserved.

    The code is based on PBRT: http://www.pbrt.org

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are
    met:

    - Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
    IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
    PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
    HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 */


#if defined(_MSC_VER)
#pragma once
#endif

#ifndef PBRT_RPF_H
#define PBRT_RPF_H

#include "pbrt.h"
#include "memory.h"
#include "rng.h"
#include "filter_utils/VectorNf.h"
#include "filter_utils/TwoDArray.h"
#include "SampleData.h"
#include "RandomParameterFilter.h"

class RPF {
public:       

    RPF(int xs, int ys, int w, int h, float jouni, string qual);

    void AddSample(const CameraSample &sample, const Spectrum &L, 
            const Intersection &isect);
    void GetAdaptPixels(int spp, vector<vector<int> > &pixels);
    void WriteImage(const string &filename, int xres, int yres, bool dump);
    void AssembleImages(bool dump);

    void SetSPP(int spp) {
    	this->spp = spp;
    	printf("Set spp to %d", spp);
    	allSamples.resize(xPixelCount * yPixelCount * spp);
    }
	static void dumpAsBinary(const string &filenameBase, const int w, const int h,
				const int spp, const vector<SampleData> &allSamples);

private:
    void WriteImage(const string &filename, const TwoDArray<Color> &image, int xres, int yres) const;
    TwoDArray<Color> FloatImageToColor(const TwoDArray<float> &image) const;
    int spp;
    long volatile sampleCount;

    vector<SampleData> allSamples;
    int xPixelStart, yPixelStart;
    int xPixelCount, yPixelCount;
    const float jouni;
    const string quality;

    // Storing the image, features and their variance 
    // reconstructed by default filter
    TwoDArray<Color> colImg;

    // Filtered image
    TwoDArray<Color> fltImg;

    // These images are stored for debug and visualization
    TwoDArray<Color> rhoImg;
    TwoDArray<Color> normalImg;
    RNG rng;

    //new debug img
    TwoDArray<Color> secNormalImg;
	TwoDArray<Color> secOrigImg;
	TwoDArray<Color> thirdOrigImg;
	TwoDArray<Color> lensImg;
	TwoDArray<float> timeImg;
};

#endif //PBRT_RPF_H

