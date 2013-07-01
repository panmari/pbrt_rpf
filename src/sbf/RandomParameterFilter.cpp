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
#define DEBUG true

#include "RandomParameterFilter.h"

#include "fmath.hpp"
#include "parallel.h"
#include "progressreporter.h"

const int BOX_SIZE[] = { 55, 35, 17, 7 };
const float MAX_SAMPLES_FACTOR[] = { 0.02f, 0.04f, 0.3f, 0.5f };
const vector<SampleData> allSamples;
int MAX_SAMPLES[4];

RandomParameterFilter::RandomParameterFilter(const int width, const int height,
		const int spp, const vector<SampleData> allSamples) {
	this->w = width;
	this->h = height;
	this->spp = spp;
	this->rng = RNG(42);
	this->allSamples = allSamples;
	this->debugLog = fopen("rpf.log", "w");
	fprintf(debugLog, "Number of samples: %lu", allSamples.size());
	for (int i = 0; i < 4; i++) {
		MAX_SAMPLES[i] = pow(BOX_SIZE[i], 2) * spp * MAX_SAMPLES_FACTOR[0];
	}
}

void RandomParameterFilter::Apply() {
	ProgressReporter reporter(4, "Applying RPF filter");
	for (int iterStep = 0; iterStep < 4; iterStep++) {
		reporter.Update(iterStep);
		for (int pixel_nr = 0; pixel_nr < 1; pixel_nr++) {
		//for (int pixel_nr = 0; pixel_nr < w * h; pixel_nr++) {
			const int pixel_idx = pixel_nr * spp;
			vector<SampleData> neighbourhood = determineNeighbourhood(
					BOX_SIZE[iterStep], MAX_SAMPLES[iterStep], pixel_idx);

		}
	}

	reporter.Done();
}

vector<SampleData> RandomParameterFilter::determineNeighbourhood(
		const int boxsize, const int maxSamples, const int pixelIdx) {
	vector<SampleData> neighbourhood;
	neighbourhood.reserve(maxSamples);
	// add all samples of current pixel
	for (int i = 0; i < spp; i++) {
		neighbourhood.push_back(allSamples[pixelIdx + i]);
	}

	// add more samples from neighbourhood
	const float stdv = boxsize / 4.f;

	SampleData pixelMean, pixelStd;
	getPixelMeanAndStd(pixelIdx, pixelMean, pixelStd);
	for (int i = 0; i < maxSamples - spp; i++) {
		int x, y;
		//retry, as long as its not in picture or original pixel
		do {
			getGaussian(stdv, pixelMean.x, pixelMean.y, x, y);
		} while(x == pixelMean.x || y == pixelMean.y || x < 0 || y < 0 || x >= w || y >= h);
		SampleData &sample = getRandomSampleAt(x, y);
		// to check if sample from right location was retrieved
		//if (DEBUG) { fprintf(debugLog, "[%d,%d vs %d,%d]", x, y, sample.x, sample.y); }
		printf("\n checking samle");
		bool flag = true;
		for (int f = 0; f < SampleData::getFeaturesSize() && flag; f++) {
			//printf("\n %f vs %f", sample[f], pixelMean[f]);
			const float lim = (f < 6) ? 30.f : 3.f;
			if( fabs(sample[f] - pixelMean[f]) > lim*pixelStd[f] &&
					(fabs(sample[f] - pixelMean[f]) > 0.1f || pixelStd[f] > 0.1f)) {
				//printf(" rejected!");
				flag = false;
			}
		}

		if (flag) {
			//need to add a copy to neighbourhood, since it will be changed!
			neighbourhood.push_back(sample);
		}
	}

	if (DEBUG) {
		fprintf(debugLog, "\n Samples in Neighbourhood: \n");
		for (unsigned int i=0;i<neighbourhood.size();i++) {fprintf(debugLog, "[%d,%d]",neighbourhood[i].x, neighbourhood[i].y); }
	}

	return neighbourhood;
}

void RandomParameterFilter::getPixelMeanAndStd(int pixelIdx,
		SampleData &pixelMean, SampleData &pixelStd) {
	SampleData pixelMeanSquare;
	for (int sampleOffset = 0; sampleOffset < spp; sampleOffset++) {
		const SampleData &currentSample = allSamples[pixelIdx + sampleOffset];
		for(int f=0;f<SampleData::getSize();f++)
		{
			pixelMean[f] += currentSample[f];
			pixelMeanSquare[f] += currentSample[f]*currentSample[f];
		}
	}
		for(int f=0;f<SampleData::getSize();f++)
		{
			pixelMean[f] /= spp;
			pixelMeanSquare[f] /= spp;

			pixelStd[f] = sqrt(max(0.f,pixelMeanSquare[f] - pixelMean[f]*pixelMean[f]));	// max() avoids accidental NaNs
		}
	}

void RandomParameterFilter::getGaussian(float stddev, int meanX, int meanY,
		int &x, int &y) {
	// Box-Muller method, adapted from @ jtlehtin's code.
	float S, V1, V2;
	do {
		V1 = 2 * rng.RandomFloat() - 1;
		V2 = 2 * rng.RandomFloat() - 1;
		S = V1 * V1 + V2 * V2;
	} while (S >= 1);

	x = sqrt(-2 * log(S) / S) * V1 * stddev + meanX;
	y = sqrt(-2 * log(S) / S) * V2 * stddev + meanY;
}
