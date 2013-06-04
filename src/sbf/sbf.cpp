
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

#include "sbf.h"

#include "sampler.h"
#include "spectrum.h"
#include "intersection.h"
#include "imageio.h"
#include "montecarlo.h"
#include "progressreporter.h"
#include "CrossBilateralFilter.h"
#include "CrossNLMFilter.h"
#include "RandomParameterFilter.h"

#include "fmath.hpp"

#include <limits>
#include <algorithm>
#include <omp.h>
#include <cmath>

// Range sigma for bilateral filter, we found that with range term the result will be noisy,
// so we set the sigma to infinite to drop the range term(0 indicates infinite in our implementation)
const float c_SigmaC = 0.f;

SBF::SBF(int xs, int ys, int w, int h, 
         const Filter *filt, FilterType type,
         const vector<float> &_interParams,
         const vector<float> &_finalParams,
         float _sigmaN, float _sigmaR, float _sigmaD,
         float _interMseSigma, float _finalMseSigma) :
    fType(type), rFilter(filt),  
    interParams(_interParams), finalParams(_finalParams),
    sigmaN(_sigmaN), sigmaR(_sigmaR), sigmaD(_sigmaD),
    interMseSigma(_interMseSigma), finalMseSigma(_finalMseSigma) {
    xPixelStart = xs;
    yPixelStart = ys;
    xPixelCount = w;
    yPixelCount = h;
    pixelInfos = new BlockedArray<PixelInfo>(xPixelCount, yPixelCount);
    // TODO replace 8 with spp
    allSamples = vector<SampleData> (xPixelCount*yPixelCount*8);
    colImg = TwoDArray<Color>(xPixelCount, yPixelCount);
    varImg = TwoDArray<Color>(xPixelCount, yPixelCount);
    featureImg = TwoDArray<Feature>(xPixelCount, yPixelCount);
    featureVarImg = TwoDArray<Feature>(xPixelCount, yPixelCount);

    norImg = TwoDArray<Color>(xPixelCount, yPixelCount);
    rhoImg = TwoDArray<Color>(xPixelCount, yPixelCount);
    depthImg = TwoDArray<float>(xPixelCount, yPixelCount);
    rhoVarImg = TwoDArray<Color>(xPixelCount, yPixelCount);
    norVarImg = TwoDArray<Color>(xPixelCount, yPixelCount);
    depthVarImg = TwoDArray<float>(xPixelCount, yPixelCount);
    
    fltImg = TwoDArray<Color>(xPixelCount, yPixelCount);
    minMseImg = TwoDArray<float>(xPixelCount, yPixelCount);
    adaptImg = TwoDArray<float>(xPixelCount, yPixelCount);
    sigmaImg = TwoDArray<Color>(xPixelCount, yPixelCount);

    //new:
    secNormalImg = TwoDArray<Color>(xPixelCount, yPixelCount);
    secOrigImg = TwoDArray<Color>(xPixelCount, yPixelCount);
    thirdOrigImg = TwoDArray<Color>(xPixelCount, yPixelCount);
    lensImg = TwoDArray<Color>(xPixelCount, yPixelCount);
    timeImg = TwoDArray<float>(xPixelCount, yPixelCount);
    sampleCount = -1;
}

void SBF::AddSample(const CameraSample &sample, const Spectrum &L, 
                    const Intersection &isect) {    
    int x = Floor2Int(sample.imageX)-xPixelStart;
    int y = Floor2Int(sample.imageY)-yPixelStart;
    // Check if the sample is in the image
    if (x < 0 || y < 0 || x >= xPixelCount || y >= yPixelCount) 
        return;    

    // Convert to 3d color space from Spectrum
    float xyz[3];
    L.ToRGB(xyz);
    float rhoXYZ[3];
    isect.rho.ToRGB(rhoXYZ);

    int idx = AtomicAdd(&sampleCount, (long)1);
    SampleData& sd = allSamples[idx];
    sd.x = x;
    sd.y = y;
    for(int i = 0; i < 3; i++) {
    	sd.inputColors[i] = sd.rgb[i] = xyz[i];
    	sd.rho[i] = rhoXYZ[i];
    	sd.normal[i] = isect.shadingN[i];
    	sd.secondNormal[i] = isect.secondNormal[i];
    	sd.secondOrigin[i] = isect.secondOrigin[i];
    	sd.thirdOrigin[i] = isect.thirdOrigin[i];
    }
    sd.imgPos[0] = sample.imageX;
    sd.imgPos[1] = sample.imageY;
    sd.lensPos[0] = sample.lensU;
    sd.lensPos[1] = sample.lensV;
	sd.time = sample.time;
}

void SBF::GetAdaptPixels(int spp, vector<vector<int> > &pixels) {
    Update(false);
    
    // We use long long here since int will overflow for very extreme case
    // (e.g. for a very big image with size 2560x1920, unsigned int can only afford 873 spp)
    long long totalSamples = (long long)xPixelCount*
                             (long long)yPixelCount*
                             (long long)spp;

    long double probSum = 0.0L;
    for(int y = 0; y < yPixelCount; y++)
        for(int x = 0; x < xPixelCount; x++) {
            probSum += (long double)adaptImg(x, y);
        }
    long double invProbSum = 1.0L/probSum;

    // Clear pixels
    vector<vector<int> >().swap(pixels);

    pixels.resize(yPixelCount);
    for(int y = 0; y < yPixelCount; y++) {
        pixels[y].resize(xPixelCount);
        for(int x = 0; x < xPixelCount; x++) {
            pixels[y][x] = 
                max(Ceil2Int((long double)totalSamples * 
                             (long double)adaptImg(x, y) * invProbSum), 1);
        }
    }
}

float SBF::CalculateAvgSpp() const {
    unsigned long long totalSamples = 0;
    for(int y = 0; y < yPixelCount; y++)
        for(int x = 0; x < xPixelCount; x++) {
            totalSamples += (unsigned long long)(*pixelInfos)(x, y).sampleCount;
        }
    long double avgSpp = (long double)totalSamples/(long double)(xPixelCount*yPixelCount);
    return (float)avgSpp;
}

void SBF::WriteImage(const string &filename, int xres, int yres, bool dump) {
    Update(true);

    string filenameBase = filename.substr(0, filename.rfind("."));
    string filenameExt  = filename.substr(filename.rfind("."));

    printf("Avg spp: %.2f\n", CalculateAvgSpp());

    WriteImage(filenameBase+"_sbf_img"+filenameExt, colImg, xres, yres);
    WriteImage(filenameBase+"_sbf_flt"+filenameExt, fltImg, xres, yres);
    TwoDArray<Color> sImg = TwoDArray<Color>(xPixelCount, yPixelCount);
    for(int y = 0; y < yPixelCount; y++)
        for(int x = 0; x < xPixelCount; x++) {
            float sc = (float)(*pixelInfos)(x, y).sampleCount;
            sImg(x, y) = Color(sc, sc, sc);
        }        
    WriteImage(filenameBase+"_sbf_smp"+filenameExt, sImg, xres, yres);
    WriteImage(filenameBase+"_sbf_param"+filenameExt, sigmaImg, xres, yres);

    if(dump) { // Write debug images
        WriteImage(filenameBase+"_sbf_var"+filenameExt, varImg, xres, yres);
        
        // Normals contain negative values, normalize them here
        for(int y = 0; y < norImg.GetRowNum(); y++)
            for(int x = 0; x < norImg.GetColNum(); x++) {
                norImg(x, y) += Color(1.f, 1.f, 1.f);
                norImg(x, y) /= 2.f;
            }
        WriteImage(filenameBase+"_sbf_nor"+filenameExt, norImg, xres, yres);

        WriteImage(filenameBase+"_sbf_rho"+filenameExt, rhoImg, xres, yres);

        //new:
        WriteImage(filenameBase+"_sbf_second_normal"+filenameExt, secNormalImg, xres, yres);
        WriteImage(filenameBase+"_sbf_second_orig"+filenameExt, secOrigImg, xres, yres);
        WriteImage(filenameBase+"_sbf_third_orig"+filenameExt, thirdOrigImg, xres, yres);
        WriteImage(filenameBase+"_sbf_lens"+filenameExt, lensImg, xres, yres);
        TwoDArray<Color> timeColImg = FloatImageToColor(timeImg);
        WriteImage(filenameBase+"_sbf_time"+filenameExt, timeColImg, xres, yres);
    }
}

void SBF::WriteImage(const string &filename, const TwoDArray<Color> &image, int xres, int yres) const {
    ::WriteImage(filename, (float*)image.GetRawPtr(), NULL, xPixelCount, yPixelCount,
                 xres, yres, xPixelStart, yPixelStart);
}

TwoDArray<Color> SBF::FloatImageToColor(const TwoDArray<float> &image) const {
    TwoDArray<Color> colorImg(image.GetColNum(), image.GetRowNum());
    for(int y = 0; y < yPixelCount; y++)
        for(int x = 0; x < xPixelCount; x++) {
            float val = image(x, y);
            colorImg(x, y) = Color(val, val, val);
        } 
    return colorImg;
}

bool SBF::comparator(SampleData sd1, SampleData sd2) {
	if (sd1.y == sd2.y)
			return sd1.x < sd2.x;
	else return sd1.y < sd2.y;
}

void SBF::Update(bool final) {
    ProgressReporter reporter(1, "Dumping debug images");
    std::sort(allSamples.begin(), allSamples.end(), SBF::comparator);
#pragma omp parallel for num_threads(PbrtOptions.nCores)
    for(int i = 0; i < yPixelCount*xPixelCount*8; i++) {
		SampleData sd = allSamples[i];

		int x = Floor2Int(sd.imgPos[0])-xPixelStart;
		int y = Floor2Int(sd.imgPos[1])-yPixelStart;

		Color rgbC = Color(sd.rgb); //color in RGB
		Color normalC = Color(sd.normal);

		Color rhoC = Color(sd.rho);

		//new
		Color secNormalC = Color(sd.secondNormal);
		Color secOriginC = Color(sd.secondOrigin);
		Color thirdOriginC = Color(sd.thirdOrigin);
		Color lensC = Color(sd.lensPos[0], sd.lensPos[1], 0.f);

		colImg(x, y) = rgbC;
		norImg(x, y) = normalC;
		rhoImg(x, y) = rhoC;
		//new
		secNormalImg(x, y) = secNormalC;
		secOrigImg(x, y) = secOriginC;
		thirdOrigImg(x, y) = thirdOriginC;
		lensImg(x, y) = lensC;
	}
    reporter.Update();
    reporter.Done();

    TwoDArray<Color> rColImg = colImg;
    /**
     *  We use the image filtered by the reconstruction filter for MSE estimation,
     *  but apply filtering on the 1x1 box filtered image.
     *  We found that this gives sharper result and smoother filter selection
     */
    rFilter.Apply(rColImg);
    /**
     *  Theoratically, we should use squared kernel to filter variance,
     *  however we found that it will produce undersmoothed image(this is
     *  because we assumed Gaussian white noise when performing MSE 
     *  estimation, so we did not consider the covariances between pixels)
     *  Therefore we reconstruct the variance with the original filter.      
     */
    rFilter.Apply(varImg);

    // We reconstruct feature buffers with 1x1 box filter as it gives us sharper result
    // In the case that the feature buffers are very noisy like heavy DOF or very fast
    // motion, it might be a good idea to filter the feature buffer. But as the variance
    // will be very local, we will have to apply some adaptive filters.

    vector<float> sigma = final ? finalParams : interParams;
    Feature sigmaF;
    sigmaF[0] = sigmaF[1] = sigmaF[2] = sigmaN;
    sigmaF[3] = sigmaF[4] = sigmaF[5] = sigmaR;
    sigmaF[6] = sigmaD;

    vector<TwoDArray<Color> > fltArray;
    vector<TwoDArray<float> > mseArray;
    vector<TwoDArray<float> > fltMseArray;
    for(size_t i = 0; i < sigma.size(); i++) {
        fltArray.push_back(TwoDArray<Color>(xPixelCount, yPixelCount));
        mseArray.push_back(TwoDArray<float>(xPixelCount, yPixelCount));
        fltMseArray.push_back(TwoDArray<float>(xPixelCount, yPixelCount));
    }

    if(fType == CROSS_BILATERAL_FILTER) {
        for(size_t i = 0; i < sigma.size(); i++) {
            CrossBilateralFilter cbFilter(sigma[i], c_SigmaC, sigmaF, xPixelCount, yPixelCount); 
            TwoDArray<Color> flt(xPixelCount, yPixelCount);
            TwoDArray<float> mse(xPixelCount, yPixelCount);
            cbFilter.Apply(colImg, featureImg, featureVarImg, rColImg, varImg, flt, mse);
            mseArray[i] = mse;
            fltArray[i] = flt;
        }    

        CrossBilateralFilter mseFilter(final ? finalMseSigma : interMseSigma, 0.f, 
                                       sigmaF, xPixelCount, yPixelCount); 
        mseFilter.ApplyMSE(mseArray, featureImg, featureVarImg, fltMseArray);
    } else if (fType == CROSS_NLM_FILTER) {
        CrossNLMFilter nlmFilter(final ? 20 : 10, 2, sigma, sigmaF, 
                xPixelCount, yPixelCount);
        nlmFilter.Apply(colImg, featureImg, featureVarImg, 
                rColImg, varImg, fltArray, mseArray);
        // We use cross bilateral filter to filter MSE estimation even for NLM filters.
        CrossBilateralFilter mseFilter(final ? finalMseSigma : interMseSigma, 0.f, 
                                       sigmaF, xPixelCount, yPixelCount);   
        mseFilter.ApplyMSE(mseArray, featureImg, featureVarImg, fltMseArray);
        //filter.ApplyMSE(0.04f, mseArray, rColImg, featureImg, featureVarImg, fltMseArray);
    } else { //fType == RANDOM_PARAMETER_FILTER
    	RandomParameterFilter rpf;
    	rpf.Apply(allSamples);

    }

    minMseImg = numeric_limits<float>::infinity();   
    for(size_t i = 0; i < sigma.size(); i++) {
#pragma omp parallel for num_threads(PbrtOptions.nCores)
        for(int y = 0; y < yPixelCount; y++)
            for(int x = 0; x < xPixelCount; x++) {
                float error = fltMseArray[i](x, y);
                if(error < minMseImg(x, y)) {
                    Color c = fltArray[i](x, y);
                    adaptImg(x, y) = error/(c.Y()*c.Y()+1e-3f);
                    minMseImg(x, y) = error;
                    fltImg(x, y) = c;
                    sigmaImg(x, y) = Color((float)i/(float)sigma.size());
                }
            }
    }
}

