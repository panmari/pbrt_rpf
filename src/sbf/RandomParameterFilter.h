
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

#include "SBFCommon.h"
#include "TwoDArray.h"
#include "VectorNf.h"
#include "pbrt.h"
#include "SampleData.h"
#include "rng.h"

class RandomParameterFilter {
public:
    RandomParameterFilter(const int width, const int height,
    		const int spp, vector<SampleData> &allSamples);

    void Apply();

private:
    int h, w, spp;
	FILE *debugLog;
	vector<SampleData> &allSamples;
	const RNG rng; //random generator

	void preprocessSamples();
	void dumpIntermediateResults(int iterStep);
    vector<SampleData> determineNeighbourhood(const int boxsize, const int maxSamples, const int pixelIdx);
    void computeWeights(vector<float> &alpha, vector<float> &beta, float &W_r_c, vector<SampleData> &neighbourhood,int iterStep);
    void filterColorSamples(vector<float> &alpha, vector<float> &beta, float W_r_c, vector<SampleData> &neighbourhood, int currentPixelIdx);
    //some helpers
    inline float sqr(float a) {return a*a;};
    inline float rcp(const float a) { return (a) ? 1.f/ a : 0.f; };
    void getPixelMeanAndStd(int pixelIdx, SampleData &sampleMean, SampleData &sampleStd);
    void getGaussian(float stddev, int &x, int &y) const;
    SampleData& getRandomSampleAt(const int x, const int y, int &idx); //TODO: long for veeery big images?
};

/**
 * Again heavily inspired by jklethinens code.
 * You must make an instance of this (instead of static), so memory for histograms only needs to be assigned once.
 */
#define NR_BUCKETS 5
class MutualInformation {
public:
	float mutualinfo(const vector<SampleData> &neighbourhood, const int firstChannel, const int secondChannel) {
		clearHistograms();
		for (const SampleData& s: neighbourhood) {
			int a = quantize(s[firstChannel]);
			int b = quantize(s[secondChannel]);
			hist_a[a]++;
			hist_b[b]++;
			hist_ab[a*NR_BUCKETS+b]++;
		}
		//compute entropies
		float ent_a = 0.f;
		float ent_b = 0.f;
		for (int i = 0; i < NR_BUCKETS; i++) {
			if(hist_a[i]) {
				float prob_a = hist_a[i]/neighbourhood.size();
				ent_a += -prob_a*log2f(prob_a);
			}
			if(hist_b[i]) {
				float prob_b = hist_b[i]/neighbourhood.size();
				ent_b += -prob_b*log2f(prob_b);
			}
		}
		float ent_ab = 0.f;
		for (int i = 0; i < NR_BUCKETS*NR_BUCKETS; i++) {
			if(hist_ab[i]) {
				float prob_ab = hist_ab[i]/neighbourhood.size();
				ent_ab += -prob_ab*log2f(prob_ab);
			}
		}
		return ent_a + ent_b - ent_ab;
	}

private:
	void clearHistograms() {
		for (int i = 0; i < NR_BUCKETS; i++) {
			hist_a[i] = hist_b[i] = 0.f;
		}
		for (int i = 0; i < NR_BUCKETS*NR_BUCKETS; i++) {
			hist_ab[i] = 0.f;
		}
	}
	float hist_a[NR_BUCKETS];
	float hist_b[NR_BUCKETS];
	float hist_ab[NR_BUCKETS*NR_BUCKETS];

	inline int quantize(float v) {
		v = (v+2)/4;
		v *= NR_BUCKETS-1;
		int bucket = (int)(v + 0.5f);
		return min(max(bucket, 0), NR_BUCKETS - 1);
	}
};

#endif
