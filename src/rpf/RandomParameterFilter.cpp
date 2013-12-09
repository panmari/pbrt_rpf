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
// But it does also darken these scenes very heavily!
#define REINSERT_ENERGY_HDR_CLAMP true					// jlehtinen => true, sen => false
#define PER_CHANNEL_ALPHA false							// jlehtinen => false, sen => true
#define PREAPPLY_GAMMA 0								// Set to 0 if should not be preapplied, usually a bad idea

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
	printf("Number of samples: %lu, size: %dx%d, spp: %d \n",
			allSamples.size(), w, h, spp);
}

void RandomParameterFilter::Apply() {
	timeval startTime, endTime;
	gettimeofday(&startTime, NULL);
	preprocessSamples();

	for (int iterStep = 0; iterStep < 4; iterStep++) {
		ProgressReporter reporter(w*h, "Applying RPF filter, pass " + std::to_string(iterStep + 1) + " of 4");
		if (DEBUG) fprintf(debugLog, "\n*** Starting pass number %d ***\n", iterStep);
#if DEBUG
		for (int pixel_nr = DEBUG_PIXEL_NR; pixel_nr <= DEBUG_PIXEL_NR; pixel_nr++) {
			fprintf(debugLog, "Debugging pixel nr %d, at %d, %d \n", pixel_nr, pixel_nr%w, (int)pixel_nr/w);
#else
#pragma omp parallel for num_threads(PbrtOptions.nCores)
		for (int pixel_nr = 0; pixel_nr < w * h; pixel_nr++) {
#endif
			const int pixel_idx = pixel_nr * spp;
			vector<SampleData> neighbourhood = determineNeighbourhood(BOX_SIZE[iterStep], MAX_SAMPLES[iterStep], pixel_idx);

			if (DEBUG) {
				fprintf(debugLog, "\nNormalized feature vectors in neighbourhood: \n");
				for (SampleData& s: neighbourhood) {
					//verified with matlab, has mean 0 and std 1
					for (int f=FEATURES_OFFSET; f < FEATURES_SIZE; f++) {fprintf(debugLog, "%-.3f ", s[f]); }
					fprintf(debugLog, "\n");
				}
				fflush(debugLog);
			}

			vector<float> alpha = vector<float>(COLOR_SIZE);
			vector<float> beta = vector<float>(FEATURES_SIZE);
			float W_r_c;
			computeWeights(alpha, beta, W_r_c, neighbourhood, iterStep);

			if (DEBUG) {
				fprintf(debugLog, "\nalpha: ");
				for(uint i=0; i<alpha.size(); i++) { fprintf(debugLog, "%-.3f, ", alpha[i]); }
				fprintf(debugLog, "\nbeta: ");
				for(uint i=0; i<beta.size(); i++) { fprintf(debugLog, "%-.3f, ", beta[i]); }
				fflush(debugLog);
			}
			filterColorSamples(alpha, beta, W_r_c, neighbourhood, pixel_idx);

			if (pixel_nr % (20*w) == 0) {
				reporter.Update(20*w);
			}
		}

		//write output to input
		for (SampleData &s: allSamples) {
			for (int k=0; k<3; k++) {
				s.inputColors[k] = s.outputColors[k];
			}
		}

		reporter.Done();
		if (DUMP_INTERMEDIATE_RESULTS)
			dumpIntermediateResults(iterStep);
	}
	gettimeofday(&endTime, NULL);
	int duration(endTime.tv_sec - startTime.tv_sec);
	printf("The whole rendering process took %d minutes and %d seconds \n", duration/60, duration%60);
}

void RandomParameterFilter::dumpIntermediateResults(int iterStep) {
	TwoDArray<Color> fltImg = TwoDArray<Color>(w, h);
	for (uint i=0; i < allSamples.size(); i+=spp) {
		Color c;
		for (int j=0; j<spp; j++)
			for(int k=0; k<3;k++){
				//TODO: if jkl_dump is used, this should also be multiplied with rho
				c[k] += allSamples[i+j].outputColors[k]; //*allSamples[i+j].rho[k];
			}
		c /= spp;
		fltImg(allSamples[i].x, allSamples[i].y) = c;
	}
	// passing null as alpha makes it 1.f for pixel
	WriteImage("pass" + to_string(iterStep+1) + ".exr", (float*)fltImg.GetRawPtr(), NULL, w, h,
					 w, h, 0, 0);
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
	RNG rng(pixelIdx);
	for (int i = 0; i < maxSamples - spp; i++) {
		int x = 0, y = 0, idx; // x, y are only set to prevent warning
		//retry, as long as its not in picture or original pixel
		do {
			float offsetX, offsetY;
			getGaussian(stdv, offsetX, offsetY, rng);
			if (CROP_BOX && (fabs(offsetX) >= boxsize/2.f || fabs(offsetY) >= boxsize/2.f))		// get only pixels inside of 'box'
				continue;
			x = pixelMean.x + int(floor(offsetX+0.5f));
			y = pixelMean.y + int(floor(offsetY+0.5f));
		} while((x == pixelMean.x && y == pixelMean.y) || 			// can not be same pixel
				x < 0 || y < 0 || x >= w || y >= h);				// or outside of image
		SampleData &sample = getRandomSampleAt(x, y, idx, rng);
		bool flag = true;
		for (int f = FEATURES_OFFSET; f < FEATURES_SIZE && flag; f++) {
			const float lim = (f < 6) ? 30.f : 3.f;
			if( fabs(sample[f] - pixelMean[f]) > lim*pixelStd[f] &&
					(fabs(sample[f] - pixelMean[f]) > 0.1f || pixelStd[f] > 0.1f)) {
					flag = false;
			}
		}
		if (flag) {
			//by default, this pushes a copy there
			neighbourhood.push_back(sample);
		}
	}

	if (DEBUG) {
		fprintf(debugLog, "\nSamples in Neighbourhood (%ld): \n", neighbourhood.size());
		for (unsigned int i=0;i<neighbourhood.size();i++) {
			fprintf(debugLog, "[%d,%d]", neighbourhood[i].x, neighbourhood[i].y);
		}
	}
	
	// Normalization of neighbourhood
	SampleData nMean, nMeanSquare, nStd;
	nMean.reset(); nMeanSquare.reset();
	for (int f = 0; f < LAST_NORMALIZED_OFFSET; f++) {
		for (SampleData& s: neighbourhood) {
			nMean[f] += s[f];
			nMeanSquare[f] += sqr(s[f]);
		}
		nMean[f] /= neighbourhood.size();
		nMeanSquare[f] /= neighbourhood.size();
		nStd[f] = sqrt(max(0.f, nMeanSquare[f] - sqr(nMean[f])));
	}
	for (int f = 0; f < LAST_NORMALIZED_OFFSET; f++) {
		float overStd = rcp(nStd[f]);
		for (SampleData& s: neighbourhood) {
			s[f] = (s[f] - nMean[f])*overStd;
		}
	}
	return neighbourhood;
}

void RandomParameterFilter::computeWeights(vector<float> &alpha, vector<float> &beta,
		float &W_r_c, vector<SampleData> &neighbourhood,int iterStep) {
	MutualInformation mi;
	// dependency for colors

	W_r_c = 0.f;
	float D_r_c = 0.f, D_p_c = 0.f, D_f_c = 0.f;
	vector<float> m_D_fk_c = vector<float>(FEATURES_SIZE);
	std::fill(m_D_fk_c.begin(), m_D_fk_c.end(), 0.f);

	for(int l=0; l < COLOR_SIZE; l++) {
		float m_D_r_cl = 0.f;
		float m_D_p_cl = 0.f;
		float m_D_f_cl = 0.f;
		for(int k=0; k < randomParamsSize; k++) {
			m_D_r_cl += mi.mutualinfo(neighbourhood,
					l + COLOR_OFFSET, k + randomParamsOffset);
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

		// sets alpha channel dependent
		if (PER_CHANNEL_ALPHA) {
			const float W_r_cl = m_D_r_cl*rcp(m_D_r_cl + m_D_p_cl);
			// yields better results with (dubious) factor 2
			alpha[l] = max(1 - 2*(1 + 0.1f*iterStep)*W_r_cl, 0.f);
			W_r_c += W_r_cl;
		}
	}

	// sets alpha for every channel to the same value
	if (!PER_CHANNEL_ALPHA) {
		W_r_c = D_r_c*rcp(D_r_c + D_p_c);
		std::fill(alpha.begin(), alpha.end(), max(1 - (1 + 0.1f*iterStep)*W_r_c, 0.f));
	}

	const float D_a_c = D_r_c + D_p_c + D_f_c;

	for(int k = 0; k < FEATURES_SIZE; k++) {
		float m_D_fk_r = 0.f, m_D_fk_p = 0.f;
		for(int l = 0; l < randomParamsSize; l++) {
			m_D_fk_r += mi.mutualinfo(neighbourhood,
					l + randomParamsOffset, k + FEATURES_OFFSET);
		}
		for(int l=0; l < IMG_POS_SIZE; l++) {
			m_D_fk_p += mi.mutualinfo(neighbourhood,
					l + IMG_POS_OFFSET, k + FEATURES_OFFSET);
		}
		const float W_fk_r = m_D_fk_r * rcp(m_D_fk_r + m_D_fk_p);
		const float W_fk_c = m_D_fk_c[k] * rcp(D_a_c);
		beta[k] = W_fk_c * max(1-(1+0.1f*iterStep)*W_fk_r, 0.f);
	}
}

void RandomParameterFilter::filterColorSamples(vector<float> &alpha, vector<float> &beta, float W_r_c,
		vector<SampleData> &neighbourhood, int pixelIdx) {
	const float var = 8*jouni/spp;

	const float scale_f = -sqr(1 - W_r_c) / (2*var);
	const float scale_c = scale_f;
	if (DEBUG) fprintf(debugLog, "\nInput colors vs Output colors (before HDR Clamp):\n");
	for (int i=0; i<spp; i++) {
		float color[3];
		for (int j=0; j<3;j++) {color[j] = 0.f; }
		float sum_relative_weights = 0.f;
		for (uint j=0; j<neighbourhood.size(); j++) {
			float dist_c = 0.f;
			for (int k=0; k<COLOR_SIZE; k++) {
				const int offset = k + COLOR_OFFSET;
				dist_c += alpha[k] * sqr(neighbourhood[i][offset] - neighbourhood[j][offset]);
			}

			float dist_f = 0.f;
			for (int k=0; k<FEATURES_SIZE; k++) {
				const int offset = k + FEATURES_OFFSET;
				dist_f += beta[k] * sqr(neighbourhood[i][offset] - neighbourhood[j][offset]);
			}

			const float w_ij = fmath::exp(scale_c*dist_c + scale_f*dist_f);
			sum_relative_weights += w_ij;
			for (int k=0; k < 3; k++) {
				color[k] += neighbourhood[j].inputColors[k]*w_ij; //should not be normalized, check?
			}
		}
		SampleData &s = allSamples[pixelIdx + i];
		for (int k = 0; k <3; k++) { //can I assign the whole array at once?
			s.outputColors[k] = color[k]/sum_relative_weights;
			if (DEBUG) fprintf(debugLog, "%-.4f, %-.4f\n", s.inputColors[k], s.outputColors[k]);
		}
	}

	if (HDR_CLAMP) {
		float colorMean[3], colorMeanSquare[3], colorStd[3], colorMeanAfter[3];;
		for (int i=0; i<3; i++) { colorMean[i] = colorMeanSquare[i] = colorMeanAfter[i] = 0.f; }
		for (int i=0; i<spp; i++) {
			SampleData &s = allSamples[pixelIdx + i];
			for(int j=0; j<3; j++) {
				colorMean[j] += s.outputColors[j];
				colorMeanSquare[j] += sqr(s.outputColors[j]);
			}
		}
		for (int i=0; i<3; i++) {
			colorMean[i] /= spp;
			colorMeanSquare[i] /= spp;
			colorStd[i] = sqrt(max(0.f, colorMeanSquare[i] - sqr(colorMean[i])));
		}
	#define STD_FACTOR 1
		for (int i=0; i<spp; i++) {
			SampleData &s = allSamples[pixelIdx + i];
			if( fabs(s.outputColors[0] - colorMean[0]) > STD_FACTOR*colorStd[0] ||
				fabs(s.outputColors[1] - colorMean[1]) > STD_FACTOR*colorStd[1] ||
				fabs(s.outputColors[2] - colorMean[2]) > STD_FACTOR*colorStd[2]) {
				for (int j=0; j<3;j++) {
					s.outputColors[j] = colorMean[j];
				}
			}
			for (int j=0; j<3; j++) {
				colorMeanAfter[j] += s.outputColors[j];
			}
		}

		if (REINSERT_ENERGY_HDR_CLAMP) {
			for (int j=0; j<3; j++) {
				colorMeanAfter[j] /= spp;
			}

			// reinsert energy from HDR clamp
			float lostEnergyPerSample[3];
			for (int i=0; i<3; i++) { lostEnergyPerSample[i] = colorMean[i] - colorMeanAfter[i]; }
			for (int i=0; i<spp; i++) {
				SampleData &s = allSamples[pixelIdx + i];
				for (int j=0; j<3; j++) {
					s.outputColors[j] += lostEnergyPerSample[j];
				}
			}
		}
	}

	if (DEBUG) {
		fprintf(debugLog, "After HDR-clamp and reinsertion of energy: \n");
		for (int i = 0; i < spp; i++) { //can I assign the whole array at once?
			SampleData &s = allSamples[pixelIdx + i];
			for (int k=0; k<3; k++) fprintf(debugLog, "%-.4f, %-.4f\n", s.inputColors[k], s.outputColors[k]);
		}
	}
}

/**
 * Only the features of the returned SampleData contain the desired values.
 * Everything else is set to 0 (rgb, random params etc.)
 * Only x, y are taken from the first sample of the pixel and assigned to pixelMean
 */
void RandomParameterFilter::getPixelMeanAndStd(int pixelIdx,
		SampleData &pixelMean, SampleData &pixelStd) const {
	SampleData pixelMeanSquare;
	pixelMean.reset(); pixelMeanSquare.reset();
	//set x and y separately
	pixelMean.x = allSamples[pixelIdx].x;
	pixelMean.y = allSamples[pixelIdx].y;
	// Online algorithm for computing variance from
	// http://en.wikipedia.org/wiki/Algorithms_for_calculating_variance
	vector<float> M2 = vector<float>(FEATURES_SIZE);
	std::fill(M2.begin(), M2.end(), 0.f);
	int n = 0;
	for (int sampleOffset = 0; sampleOffset < spp; sampleOffset++) {
		const SampleData &currentSample = allSamples[pixelIdx + sampleOffset];
		n += 1;
		for(int f=0;f<FEATURES_SIZE;f++)
		{
			float delta = currentSample[f] - pixelMean[f];
			pixelMean[f] += delta/n;
			M2[f] += delta*(currentSample[f] - pixelMean[f]);
		}
	}
	for(int f=0;f<FEATURES_SIZE;f++) {
		pixelStd[f] = max(0.f, M2[f]/(spp - 1));	// max() avoids accidental NaNs
	}
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
void RandomParameterFilter::setRandomParams(string randomParamsString) {
	if (randomParamsString == "frd") {
		randomParamsOffset = 20;
		randomParamsSize = 3;
	} else if (randomParamsString == "lens") {
		randomParamsOffset = 23;
		randomParamsSize = 2;
	} else if (randomParamsString == "time") {
		randomParamsOffset = 25;
		randomParamsSize = 1;
	} else if (std::string::npos != randomParamsString.find("frd") &&
			std::string::npos != randomParamsString.find("lens")) {
		randomParamsOffset = 20;
		randomParamsSize = 5;
	} else if (std::string::npos != randomParamsString.find("lens") &&
			std::string::npos != randomParamsString.find("time")) {
		randomParamsOffset = 23;
		randomParamsSize = 3;
	} else {
		randomParamsOffset = 20;
		randomParamsSize = 6;
	}
	printf("Using features from %d until %d as random parameters \n",
			randomParamsOffset, randomParamsOffset + randomParamsSize);
}
