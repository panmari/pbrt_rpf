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
//debugging stuff
#define DEBUG false
#define DEBUG_PIXEL_NR 173+100*w
//910344/spp
//200 + 200*w
#define DUMP_INTERMEDIATE_RESULTS true

//some parameters that should stay true for most cases
#define CROP_BOX true                 				  	// jlehtinen => true, sen => false?
#define HDR_CLAMP true									// both true
// For some scenes this is very problematic, because spikes are not properly removed if activated...
// But does also make tone mapping necessary!
#define REINSERT_ENERGY_HDR_CLAMP false					// jlehtinen => true, sen => false
#define PER_CHANNEL_ALPHA false							// jlehtinen => false, sen => true
#define PREAPPLY_GAMMA 0								// Set to 0 if should not be preapplied, usually a bad idea
// You'll most likely want to change this:
#define RANDOM_PARAMS_SIZE 5 // TODO: make this a user-defined parameter.

#include "RandomParameterFilter.h"

#include "filter_utils/fmath.hpp"
#include "filter_utils/VectorNf.h"

#include "parallel.h"
#include "progressreporter.h"
#include "imageio.h"
#include "MutualInformation.h"
#include "SampleData.h"
#include <boost/range/numeric.hpp>
#include <sys/time.h>

const int BOX_SIZE[] = { 55, 35, 17, 7 };
const float MAX_SAMPLES_FACTOR_MEDIUM[] = { 0.02f, 0.04f, 0.3f, 0.5f }; // for fast prototyping, by jklethinen
const float MAX_SAMPLES_FACTOR_HIGH[] = { 0.5f, 0.5f, 0.5f, 0.5f }; // by sen
int MAX_SAMPLES[4];

RandomParameterFilter::RandomParameterFilter(const int width, const int height,
		const int spp, const float _jouni, vector<SampleData> &_allSamples) :
	allSamples(_allSamples), jouni(_jouni) {
	this->w = width;
	this->h = height;
	this->spp = spp;
	if (DEBUG) {
		this->debugLog = fopen("rpf.log", "w");
		fprintf(debugLog, "Number of samples: %lu, size: %dx%d, spp: %d \n",
				allSamples.size(), w, h, spp);
	}
}

void RandomParameterFilter::Apply() {
	timeval startTime, endTime;
	gettimeofday(&startTime, NULL);
	preprocessSamples();
	TwoDArray<Color> fltImg = TwoDArray<Color>(w, h);
	for (int iterStep = 0; iterStep < 1; iterStep++) {
		ProgressReporter reporter(w*h, "Applying RPF filter, pass " + std::to_string(iterStep + 1) + " of 4");
#pragma omp parallel for num_threads(PbrtOptions.nCores)
		for (int pixel_nr = 0; pixel_nr < w * h; pixel_nr++) {
			float D_r_c = 0.f, D_p_c = 0.f, D_f_c = 0.f;
			vector<SampleData> neighbourhood;
			neighbourhood.resize(spp);
			for (int i=0; i < spp; i++) {
				neighbourhood.push_back(allSamples[pixel_nr*spp + i]);
			}
			MutualInformation mi;
			vector<float> m_D_fk_c = vector<float>(FEATURES_SIZE);
			std::fill(m_D_fk_c.begin(), m_D_fk_c.end(), 0.f);
			for(int l=0; l < COLOR_SIZE; l++) {
				float m_D_r_cl = 0.f;
				float m_D_p_cl = 0.f;
				float m_D_f_cl = 0.f;
				for(int k=0; k < RANDOM_PARAMS_SIZE; k++) {
					m_D_r_cl += mi.mutualinfo(neighbourhood,
							l + COLOR_OFFSET, k + RANDOM_PARAMS_OFFSET);
				}
				for(int k=0; k < IMG_POS_SIZE; k++) {
					m_D_p_cl += mi.mutualinfo(neighbourhood,
							l + COLOR_OFFSET, k + IMG_POS_OFFSET);
				}
				for(int k=0; k < FEATURES_SIZE; k++) {
					// needs to be saved per feature and per color
					const float m_D_fk_cl = mi.mutualinfo(neighbourhood,
							l + COLOR_OFFSET, k + FEATURES_OFFSET);
					m_D_fk_c[k] += m_D_fk_cl;
					m_D_f_cl += m_D_fk_cl;
				}
				D_r_c += m_D_r_cl;
				D_p_c += m_D_p_cl;
				D_f_c += m_D_f_cl;
			}
			SampleData &s = allSamples[pixel_nr*8];
			//Todo: use D_a_c?
			const float D_a_c = D_r_c + D_p_c + D_f_c;
			fltImg(s.x, s.y) = Color(D_r_c*rcp(D_r_c + D_p_c));
		}
	}
	WriteImage("lens_dependancy_color.exr", (float*)fltImg.GetRawPtr(), NULL, w, h,
						 w, h, 0, 0);
	gettimeofday(&endTime, NULL);
	int duration(endTime.tv_sec - startTime.tv_sec);
	printf("The whole rendering process took %d minutes and %d seconds \n", duration/60, duration%60);
}

/**
 * Or something like this... probably I'd have to overwrite all features of pixelMean
 */
void RandomParameterFilter::preprocessSamples() {
	printf("Preprocessing... \n");
	vector<int> pixelWithInvalidSamplesCount(spp);
	RNG rng(42);
	for (uint pixelOffset = 0; pixelOffset < allSamples.size(); pixelOffset+= spp) {
		SampleData pixelValidSamplesMean;
		pixelValidSamplesMean.reset();
		vector<uint> validSamplesIdx, invalidSamplesIdx;
		for (int sampleOffset = 0; sampleOffset < spp; sampleOffset++) {
			uint idx = pixelOffset + sampleOffset;
			SampleData &s = allSamples[idx];
			if (PREAPPLY_GAMMA) {
				for(int i=0; i<3; i++) {
					s.inputColors[i] = pow(s.rgb[i], 1/PREAPPLY_GAMMA);
				}
			}
			bool valid = true;
			// invalid, if second or third origin have one component very far off
			for (int f = FEATURES_OFFSET; f < 6; f++) {
				if (s[f] > 1e10f) {
					valid = false;
					break;
				}
			}
			if (valid) {
				pixelValidSamplesMean += s;
				validSamplesIdx.push_back(idx);
			} else {
				invalidSamplesIdx.push_back(idx);
			}
		}
		if (validSamplesIdx.size() > 0) {
			pixelValidSamplesMean.divide(validSamplesIdx.size());
		}
		// need to set x/y separately bc they were interpreted as float before
		pixelValidSamplesMean.x = allSamples[pixelOffset].x;
		pixelValidSamplesMean.y = allSamples[pixelOffset].y;
		if (invalidSamplesIdx.size() > 0)
			pixelWithInvalidSamplesCount[invalidSamplesIdx.size() - 1]++;
		for (uint invalidSampleIdx: invalidSamplesIdx) {
			SampleData &s = allSamples[invalidSampleIdx];
			if (validSamplesIdx.size() > 0) {
				//replace invalid sample with random valid sample from same pixel
				int replaceIdx = (int) (rng.RandomFloat()*validSamplesIdx.size());
				s = allSamples[validSamplesIdx[replaceIdx]];
				//put radiance to black of overwritten sample
				for (int i = 0; i < 3; i++) {
					s.rgb[i] = 0;
					s.inputColors[i] = 0;
					s.outputColors[i] = 0;
				}
			}
			else {
				s = pixelValidSamplesMean; //set everything to 0 basically
			}
		}
	}
	bool fixedInvalidSamples = false;
	for (int i=0; i < spp; i++) {
		if (pixelWithInvalidSamplesCount[i]) {
			printf("%d pixels with %d invalid samples \n", pixelWithInvalidSamplesCount[i], i + 1);
			fixedInvalidSamples = true;
		}
	}
	if (!fixedInvalidSamples) {
		printf("No invalid samples found. \n");
	}
	printf("Done! \n");
}

void RandomParameterFilter::getGaussian(const float stddev, float &x, float &y, RNG &rng) const {
	// Box-Muller method, adapted from @ jtlehtin's code.
	float S, V1, V2;
	do {
		V1 = 2 * rng.RandomFloat() - 1;
		V2 = 2 * rng.RandomFloat() - 1;
		S = V1 * V1 + V2 * V2;
	} while (S >= 1);

	x = sqrt(-2 * log(S) / S) * V1 * stddev;
	y = sqrt(-2 * log(S) / S) * V2 * stddev;
}

SampleData& RandomParameterFilter::getRandomSampleAt(const int x, int y, int &idx, RNG &rng) const {
	idx = (x + y*w)*spp + (int)(spp*rng.RandomFloat());
	return allSamples[idx];
}

void RandomParameterFilter::setQuality(string quality_string){
	Quality quality;
	if (quality_string == "high" || quality_string == "sen") {
		quality = Quality::HIGH;
		printf("Filter quality set to high\n");
	}
	else {
		quality = RandomParameterFilter::Quality::MEDIUM;
		printf("Filter quality set to medium\n");
	}
	for (int i = 0; i < 4; i++) {
		MAX_SAMPLES[i] = sqr(BOX_SIZE[i]) * spp;
		if (quality == Quality::MEDIUM)
			MAX_SAMPLES[i] *= MAX_SAMPLES_FACTOR_MEDIUM[i];
		else
			MAX_SAMPLES[i] *= MAX_SAMPLES_FACTOR_HIGH[i];
	}
}
