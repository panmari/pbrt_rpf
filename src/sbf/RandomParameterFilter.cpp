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
#define DEBUG_PIXEL_NR 163 + 280*w

//TODO: make this configurable in scene file
#define JOUNI 0.1f

//some parameters that should stay true for most things
#define CROP_BOX true
#define HDR_CLAMP true
#define REINSERT_ENERGY_HDR_CLAMP true

//used to prevent div by zero
#define EPSILON 1e-10

#include "RandomParameterFilter.h"

#include "fmath.hpp"
#include "parallel.h"
#include "progressreporter.h"
#include "imageio.h"
#include<boost/range/numeric.hpp>
const int BOX_SIZE[] = { 55, 35, 17, 7 };
const float MAX_SAMPLES_FACTOR[] = { 0.02f, 0.04f, 0.3f, 0.5f }; // for fast prototyping, by jklethinen
//const float MAX_SAMPLES_FACTOR[] = { 0.1f, 0.2f, 0.3f, 0.5f }; // by me
//const float MAX_SAMPLES_FACTOR[] = { 0.5f, 0.5f, 0.5f, 0.5f }; // by sen
int MAX_SAMPLES[4];

RandomParameterFilter::RandomParameterFilter(const int width, const int height,
		const int spp, vector<SampleData> &allSamples) :
	allSamples(allSamples), rng(RNG(42)) {
	this->w = width;
	this->h = height;
	this->spp = spp;
	this->debugLog = fopen("rpf.log", "w");
	fprintf(debugLog, "Number of samples: %lu", allSamples.size());
	for (int i = 0; i < 4; i++) {
		MAX_SAMPLES[i] = pow(BOX_SIZE[i], 2) * spp * MAX_SAMPLES_FACTOR[0];
	}
}

void RandomParameterFilter::Apply() {
	preprocessSamples();

	for (int iterStep = 0; iterStep < 4; iterStep++) {
		ProgressReporter reporter(w*h, "Applying RPF filter, pass " + std::to_string(iterStep + 1) + " of 4");
		if (DEBUG) fprintf(debugLog, "\n*** Starting pass number %d ***\n", iterStep);
#pragma omp parallel for num_threads(PbrtOptions.nCores)
		for (int pixel_nr = 0; pixel_nr < w * h; pixel_nr++) {
			const int pixel_idx = pixel_nr * spp;
			vector<SampleData> neighbourhood = determineNeighbourhood(BOX_SIZE[iterStep], MAX_SAMPLES[iterStep], pixel_idx);

			if (DEBUG) {
				fprintf(debugLog, "\nNormalized feature vectors in neighbourhood: \n");
				for (SampleData& s: neighbourhood) {
					//verified with matlab, has mean 0 and std 1
					for (int f=SampleData::getFeaturesOffset(); f < SampleData::getFeaturesSize(); f++) {fprintf(debugLog, "%-.3f ", s[f]); }
					fprintf(debugLog, "\n");
				}
				fflush(debugLog);
			}

			vector<float> alpha = vector<float>(SampleData::getColorSize());
			vector<float> beta = vector<float>(SampleData::getFeaturesSize());
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

		//write output to input (might be faster non-parallelized)
//#pragma omp parallel for num_threads(PbrtOptions.nCores)
		for (SampleData &s: allSamples) {
			for (int k=0; k<3; k++) {
				s.inputColors[k] = s.outputColors[k];
			}
		}
		reporter.Done();
	}


}

/**
 * Or something like this... probably I'd have to overwrite all features of pixelMean
 */
void RandomParameterFilter::preprocessSamples() {
	printf("Preprocessing... \n");
	for (uint pixelOffset = 0; pixelOffset < allSamples.size(); pixelOffset+= spp) {
		SampleData pixelValidSamplesMean;
		pixelValidSamplesMean.reset();
		vector<uint> validSamplesIdx, invalidSamplesIdx;
		for (int sampleOffset = 0; sampleOffset < spp; sampleOffset++) {
			uint idx = pixelOffset + sampleOffset;
			SampleData &s = allSamples[idx];
			bool valid = true;
			for (int f=0; f < 6; f++) {
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
		pixelValidSamplesMean.divide(validSamplesIdx.size());
		if (validSamplesIdx.size() < spp/2)
			printf("Pixel has only %lu valid samples \n", validSamplesIdx.size());
		for (uint invalidSampleIdx: invalidSamplesIdx) {
			//replace invalid sample with random valid sample from same pixel
			SampleData &s = allSamples[invalidSampleIdx];
			int replaceIdx = (int) (rng.RandomFloat()*validSamplesIdx.size());
			SampleData &s2 = allSamples[validSamplesIdx[replaceIdx]];
			s = s2;
			//put radiance to black of overwritten sample
			for (int i = 0; i < 3; i++) {
				s.rgb[i] = 0;
				s.inputColors[i] = 0;
				s.outputColors[i] = 0;
			}
		}
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
	for (int i = 0; i < maxSamples - spp; i++) {
		int x, y, idx;
		//retry, as long as its not in picture or original pixel
		do {
			getGaussian(stdv, x, y);
			if (CROP_BOX && (abs(x) >= boxsize/2.f || abs(y) >= boxsize/2.f))		// get only pixels inside of 'box'
				continue;
			x += pixelMean.x;
			y += pixelMean.y;
		} while(x == pixelMean.x || y == pixelMean.y || 			// can not be sampe pixel
				x < 0 || y < 0 || x >= w || y >= h);				// or outside of image
		SampleData &sample = getRandomSampleAt(x, y, idx);
		// to check if sample from right location was retrieved
		//if (DEBUG) { fprintf(debugLog, "[%d,%d vs %d,%d]", x, y, sample.x, sample.y); }
		bool flag = true;
		for (int f = SampleData::getFeaturesOffset(); f < SampleData::getFeaturesSize() && flag; f++) {
			//printf("\n %f vs %f", sample[f], pixelMean[f]);
			const float lim = (f < 6) ? 30.f : 3.f;
			if( fabs(sample[f] - pixelMean[f]) > lim*pixelStd[f] &&
					(fabs(sample[f] - pixelMean[f]) > 0.1f || pixelStd[f] > 0.1f)) {
				flag = false;
			}
		}

		if (flag) {
			//by standard, this pushes a copy there
			neighbourhood.push_back(sample);
		}
	}

	if (DEBUG) {
		fprintf(debugLog, "\nSamples in Neighbourhood (%ld): \n", neighbourhood.size());
		for (unsigned int i=0;i<neighbourhood.size();i++) {
			fprintf(debugLog, "[%d,%d]", neighbourhood[i].x, neighbourhood[i].y);
		}
		//that's very verbose...
		/*
		fprintf(debugLog, "\nRgb vs input_colors: \n");
		for (unsigned int i=0;i<neighbourhood.size();i++) {
			for (int j=0; j <3;j++) {
				fprintf(debugLog, "%-.3f\t%-.3f\n",neighbourhood[i].rgb[j], neighbourhood[i].inputColors[j]);
			}
		}
		*/
	}

	// Normalization of neighbourhood
	SampleData nMean, nMeanSquare, nStd;
	nMean.reset(); nMeanSquare.reset();
	for (int f = 0; f < SampleData::getLastNormalizedOffset(); f++) {
		for (SampleData& s: neighbourhood) {
			nMean[f] += s[f];
			nMeanSquare[f] += s[f]*s[f];
		}
		nMean[f] /= neighbourhood.size();
		nMeanSquare[f] /= neighbourhood.size();
		nStd[f] = sqrt(max(0.f,nMeanSquare[f] - nMean[f]*nMean[f]));
	}
	for (SampleData& s: neighbourhood) {
		for (int f = 0; f < SampleData::getLastNormalizedOffset(); f++) {
			//todo: could optimize this to not divide if std is 0
			s[f] = (s[f] - nMean[f])/(EPSILON + nStd[f]);
		}
	}
	return neighbourhood;
}

void RandomParameterFilter::computeWeights(vector<float> &alpha, vector<float> &beta,
		float &W_r_c, vector<SampleData> &neighbourhood,int iterStep) {
	MutualInformation mi;
	// dependency for colors

	vector<float> m_D_rk_c = vector<float>(SampleData::getRandomParametersSize());
	vector<float> m_D_pk_c = vector<float>(SampleData::getImgPosSize());
	vector<float> m_D_fk_c = vector<float>(SampleData::getFeaturesSize());
	std::fill(m_D_rk_c.begin(), m_D_rk_c.end(), 0.f );
	std::fill(m_D_pk_c.begin(), m_D_pk_c.end(), 0.f );
	std::fill(m_D_fk_c.begin(), m_D_fk_c.end(), 0.f );

	for(int l=0; l < SampleData::getColorSize(); l++) {
		for(int k=0; k < SampleData::getRandomParametersSize(); k++) {
			m_D_rk_c[k] += mi.mutualinfo(neighbourhood,
					l + SampleData::getColorOffset(), k + SampleData::getRandomParamsOffset());
		}
		for(int k=0; k < SampleData::getImgPosSize(); k++) {
			m_D_pk_c[k] += mi.mutualinfo(neighbourhood,
					l + SampleData::getColorOffset(), k + SampleData::getImgPosOffset());
		}
		for(int k=0; k < SampleData::getFeaturesSize(); k++) {
			m_D_fk_c[k] += mi.mutualinfo(neighbourhood,
					l + SampleData::getColorOffset(), k + SampleData::getFeaturesOffset());
		}
	}

	// dependency for scene features
	vector<vector<float>> m_D_fk_rl = vector<vector<float>>(SampleData::getFeaturesSize());
	vector<vector<float>> m_D_fk_pl = vector<vector<float>>(SampleData::getFeaturesSize());
	vector<vector<float>> m_D_fk_cl = vector<vector<float>>(SampleData::getFeaturesSize());
	std::fill(m_D_fk_rl.begin(), m_D_fk_rl.end(), vector<float>(m_D_rk_c.size()));
	std::fill(m_D_fk_pl.begin(), m_D_fk_pl.end(), vector<float>(m_D_pk_c.size()));
	std::fill(m_D_fk_cl.begin(), m_D_fk_cl.end(), vector<float>(3)); //three color channels

	for(int k = 0; k < SampleData::getFeaturesSize(); k++) {
		for(int l = 0; l < SampleData::getRandomParametersSize(); l++) {
			m_D_fk_rl[k][l] = mi.mutualinfo(neighbourhood,
					l + SampleData::getRandomParamsOffset(), k + SampleData::getFeaturesOffset());
		}
		for(int l=0; l < SampleData::getImgPosSize(); l++) {
			m_D_fk_pl[k][l] = mi.mutualinfo(neighbourhood,
					l + SampleData::getImgPosOffset(), k + SampleData::getFeaturesOffset());
		}
		for(int l=0; l < SampleData::getColorSize(); l++) {
			m_D_fk_cl[k][l] = m_D_fk_c[k]/3; // average of color channels
		}
	}
	const float D_r_c = boost::accumulate(m_D_rk_c, 0.f);
	const float D_p_c = boost::accumulate(m_D_pk_c, 0.f);
	const float D_f_c = boost::accumulate(m_D_fk_c, 0.f);
	const float D_a_c = D_r_c + D_p_c + D_f_c;

	W_r_c = D_r_c /(D_r_c + D_p_c + EPSILON);

	// set alpha for every channel to the same value
	std::fill(alpha.begin(), alpha.end(), max(1 - (1 + 0.1f*iterStep)*W_r_c, 0.f));

	for(int k=0;k<SampleData::getFeaturesSize();k++) {
		const float D_fk_r = boost::accumulate(m_D_fk_rl[k], 0.f);
		const float D_fk_p = boost::accumulate(m_D_fk_pl[k], 0.f);
		const float D_fk_c = boost::accumulate(m_D_fk_cl[k], 0.f);

		const float W_fk_r = D_fk_r / (D_fk_r + D_fk_p + EPSILON);
		const float W_fk_c = D_fk_c / (D_a_c + EPSILON);

		beta[k] = W_fk_c * max(1-(1+0.1f*iterStep)*W_fk_r, 0.f);
	}
}

void RandomParameterFilter::filterColorSamples(vector<float> &alpha, vector<float> &beta, float W_r_c,
		vector<SampleData> &neighbourhood, int pixelIdx) {
	const float var_8 = JOUNI;
	const float var = 8*var_8/spp;

	const float scale_f = -sqr(1 - W_r_c) / (2*var);
	const float scale_c = scale_f;
	if (DEBUG) fprintf(debugLog, "\nInput colors vs Output colors (before HDR Clamp):\n");
	for (int i=0; i<spp; i++) {
		float color[3];
		for (int j=0; j<3;j++) {color[j] = 0.f; }
		float sum_relative_weights = 0.f;
		for (uint j=0; j<neighbourhood.size(); j++) {
			float dist_c = 0.f;
			for (int k=0; k<SampleData::getColorSize(); k++) {
				const int offset = k + SampleData::getColorOffset();
				dist_c += alpha[k] * sqr(neighbourhood[i][offset] - neighbourhood[j][offset]);
			}

			float dist_f = 0.f;
			for (int k=0; k<SampleData::getFeaturesSize(); k++) {
				const int offset = k + SampleData::getFeaturesOffset();
				dist_f += beta[k] * sqr(neighbourhood[i][offset] - neighbourhood[j][offset]);
			}

			const float w_ij = exp(scale_c*dist_c + scale_f*dist_f);
			sum_relative_weights += w_ij;
			for (int k=0; k < 3; k++) {
				color[k] += neighbourhood[j].inputColors[k]*w_ij; //should not be normalized, check?
			}
		}
		SampleData &s = allSamples[pixelIdx + i];
		for (int k = 0; k <3; k++) { //can I assign the whole array at once?
			s.outputColors[k] = color[k]/sum_relative_weights;
			if (DEBUG) fprintf(debugLog, "%-.3f, %-.3f\n", s.inputColors[k], s.outputColors[k]);
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
}

/**
 * Only the features of the returned SampleData contain the desired values.
 * Everything else is set to 0 (rgb, random params etc.)
 * Only x, y are taken from the first sample of the pixel and assigned to pixelMean
 */
void RandomParameterFilter::getPixelMeanAndStd(int pixelIdx,
		SampleData &pixelMean, SampleData &pixelStd) {
	SampleData pixelMeanSquare;
	pixelMean.reset(); pixelStd.reset(); pixelMeanSquare.reset();
	//set x and y separately
	pixelMean.x = allSamples[pixelIdx].x;
	pixelMean.y = allSamples[pixelIdx].y;
	for (int sampleOffset = 0; sampleOffset < spp; sampleOffset++) {
		const SampleData &currentSample = allSamples[pixelIdx + sampleOffset];
		for(int f=0;f<SampleData::getFeaturesSize();f++)
		{
			pixelMean[f] += currentSample[f];
			pixelMeanSquare[f] += sqr(currentSample[f]);
		}
	}
	for(int f=0;f<SampleData::getFeaturesSize();f++)
		{
			pixelMean[f] /= spp;
			pixelMeanSquare[f] /= spp;

			pixelStd[f] = sqrt(max(0.f, pixelMeanSquare[f] - sqr(pixelMean[f]) ));	// max() avoids accidental NaNs
		}
	}

void RandomParameterFilter::getGaussian(float stddev, int &x, int &y) const {
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

SampleData& RandomParameterFilter::getRandomSampleAt(const int x, int y, int &idx) {
	idx = (x + y*w)*spp + (int)(spp*rng.RandomFloat());
	return allSamples[idx];
}
