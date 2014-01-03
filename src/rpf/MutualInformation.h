/*
 * MutualInfo.h
 *
 *  Created on: Jul 20, 2013
 *      Author: moser
 */

#ifndef MUTUALINFO_H_
#define MUTUALINFO_H_
#include "SampleData.h"
#include <vector>
#include "filter_utils/fmath.hpp"
/**
 * Again heavily inspired by jklethinens code.
 * You must make an instance of this (instead of static), so memory for histograms only needs to be assigned once.
 */
#define NR_BUCKETS 5
#define NORMED false
class MutualInformation {
public:


	float mutualinfo(const vector<SampleData> &neighbourhood, const int firstChannel, const int secondChannel) {
		clearHistograms();
		const uint stride = SampleData::getSize();
		const float* A = (float*)(&neighbourhood[0][firstChannel]);
		const float* B = (float*)(&neighbourhood[0][secondChannel]);
		for (uint i=0,j=0; j<neighbourhood.size(); i+=stride,j++) {
			int a = quantize(A[i]);
			int b = quantize(B[i]);
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
				ent_a += -prob_a * fmath::log2(prob_a);
			}
			if(hist_b[i]) {
				float prob_b = hist_b[i]/neighbourhood.size();
				ent_b += -prob_b * fmath::log2(prob_b);
			}
		}
		float ent_ab = 0.f;
		for (int i = 0; i < NR_BUCKETS*NR_BUCKETS; i++) {
			if(hist_ab[i]) {
				float prob_ab = hist_ab[i]/neighbourhood.size();
				ent_ab += -prob_ab * fmath::log2(prob_ab);
			}
		}
		float mi = (ent_a+ent_b-ent_ab);
		if (NORMED)
			return mi*rcp(ent_ab);
		else
			return mi;
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

    inline float rcp(const float a) const { return (a) ? 1.f/ a : 0.f; };
	inline int quantize(float v) const {
		v = (v+2)/4;
		v *= NR_BUCKETS-1;
		int bucket = (int)(v + 0.5f);
		return min(max(bucket, 0), NR_BUCKETS - 1);
	}
};

#endif /* MUTUALINFO_H_ */
