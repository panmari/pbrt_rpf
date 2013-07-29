
/*
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

#ifndef SBF_RANDOM_PARAMETER_FILTER_H__
#define SBF_RANDOM_PARAMETER_FILTER_H__

#include "filter_utils/TwoDArray.h"
#include "pbrt.h"
#include "SampleData.h"
#include "rng.h"
#include "MutualInformation.h"

class RandomParameterFilter {
public:

	enum Quality {
	        HIGH,
	        MEDIUM
	    };

    RandomParameterFilter(const int width, const int height,
    		const int spp, const float jouni, const Quality quality,
    		vector<SampleData> &allSamples);

    void Apply();

private:
    int h, w, spp;

	FILE *debugLog;
	vector<SampleData> &allSamples;
	const float jouni;

	void preprocessSamples();
	void dumpIntermediateResults(int iterStep);
    vector<SampleData> determineNeighbourhood(const int boxsize, const int maxSamples, const int pixelIdx);
    void computeWeights(vector<float> &alpha, vector<float> &beta, float &W_r_c, vector<SampleData> &neighbourhood,int iterStep);
    void filterColorSamples(vector<float> &alpha, vector<float> &beta, float W_r_c, vector<SampleData> &neighbourhood, int currentPixelIdx);
    //some helpers
    inline float sqr(float a) const {return a*a;};
    inline float rcp(const float a) const { return (a) ? 1.f/ a : 0.f; };
    void getPixelMeanAndStd(int pixelIdx, SampleData &sampleMean, SampleData &sampleStd) const;
    void getGaussian(const float stddev, float &x, float &y, RNG &rng) const;
    SampleData& getRandomSampleAt(const int x, const int y, int &idx, RNG &rng) const; //TODO: long for veeery big images?
};

#endif
