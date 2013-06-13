
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


#include "RandomParameterFilter.h"

#include "fmath.hpp"
#include "parallel.h"
#include "progressreporter.h"

const int BOX_SIZE[] = {55, 35, 17, 7};
const float MAX_SAMPLES_FACTOR[] = {0.02f, 0.04f, 0.3f, 0.5f};
int MAX_SAMPLES[4];

RandomParameterFilter::RandomParameterFilter(const int width, const int height, const int spp) {
	this->w = width;
	this->h = height;
	this->spp = spp;
	this->rng = RNG(42);
	for (int i = 0; i < 4; i++) {
		MAX_SAMPLES[i] = pow(BOX_SIZE[i],2)*spp*MAX_SAMPLES_FACTOR[0];
	}
}

void RandomParameterFilter::Apply(const vector<SampleData> &allSamples) {
	this->allSamples = allSamples;
	ProgressReporter reporter(4, "Applying RPF filter");
	for (int iterStep = 0; iterStep < 4; iterStep++) {
		reporter.Update(iterStep);
		for (int pixel_nr = 0; pixel_nr < w*h; pixel_nr++) {
			const int pixel_idx = pixel_nr*spp;
			vector<SampleData> neighbourhood = determineNeighbourhood(BOX_SIZE[iterStep], MAX_SAMPLES[iterStep], pixel_idx, allSamples);

		}
	}

	reporter.Done();
}

vector<SampleData> RandomParameterFilter::determineNeighbourhood(const int boxsize, int maxSamples,
		 const int pixelIdx, const vector<SampleData> &allSamples) {
	 	vector<SampleData> neighbourhood;
	 	neighbourhood.reserve(maxSamples);
	 	// add all samples of current pixel
	 	for (int i = 0; i < spp; i++) {
	 		neighbourhood.push_back(allSamples[pixelIdx + i]);
	 	}

	 	// add more samples from neighbourhood
	 	const float stdv = boxsize/4.f;

	 	SampleData pixelMean;
	 	SampleData pixelStd;
	 	getPixelMeanAndStd(pixelIdx, pixelMean, pixelStd);

	 	for (int i = 0; i < maxSamples - spp; i++) {
	 		int x, y;
	 		getGaussian(stdv, pixelMean.x, pixelMean.y, x, y);
	 		getRandomSampleAt(x, y)
	 	}
	return neighbourhood;
}

void RandomParameterFilter::getGaussian(float stddev, int meanX, int meanY, int &x, int &y) {
	// Box-Muller method, adapted from @ jtlehtin's code.
	float S, V1, V2;
	do {
		V1 = 2*rng.RandomFloat() - 1;
		V2 = 2*rng.RandomFloat() - 1;
		S = V1*V1 + V2*V2;
	} while(S >= 1);

	x = sqrt(-2 * log(S)/S) * V1 *stddev + meanX;
	y = sqrt(-2 * log(S)/S) * V2 * stddev + meanY;
}
