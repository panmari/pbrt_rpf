
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

#include "rpf.h"

#include "sampler.h"
#include "spectrum.h"
#include "intersection.h"
#include "imageio.h"
#include "progressreporter.h"
#include <fstream>
#include "filter_utils/fmath.hpp"

RPF::RPF(int xs, int ys, int w, int h,
          float _jouni, string _qual, string _randomParams) :
          jouni(_jouni), quality(_qual), randomParams(_randomParams) {
    xPixelStart = xs;
    yPixelStart = ys;
    xPixelCount = w;
    yPixelCount = h;
    spp = 8;
    allSamples = vector<SampleData> (spp*xPixelCount*yPixelCount);
    colImg = TwoDArray<Color>(xPixelCount, yPixelCount);

    //for debugging
    rhoImg = TwoDArray<Color>(xPixelCount, yPixelCount);
    normalImg = TwoDArray<Color>(xPixelCount, yPixelCount);
    secNormalImg = TwoDArray<Color>(xPixelCount, yPixelCount);
    secOrigImg = TwoDArray<Color>(xPixelCount, yPixelCount);
    thirdOrigImg = TwoDArray<Color>(xPixelCount, yPixelCount);
    lensImg = TwoDArray<Color>(xPixelCount, yPixelCount);
    timeImg = TwoDArray<float>(xPixelCount, yPixelCount);

    fltImg = TwoDArray<Color>(xPixelCount, yPixelCount);

    sampleCount = -1;
}

void RPF::AddSample(const CameraSample &sample, const Spectrum &L,
                    const Intersection &isect) {    
    int x = sample.x-xPixelStart;
    int y = sample.y-yPixelStart;
    // Check if the sample is in the image
    if (x < 0 || y < 0 || x >= xPixelCount || y >= yPixelCount)  {
        return;    
    }
    // Convert to 3d color space from Spectrum
    float xyz[3];
    L.ToRGB(xyz);
    float rhoXYZ[3];
    isect.rho.ToRGB(rhoXYZ);

    int idx = AtomicAdd(&sampleCount, 1);
    SampleData& sd = allSamples[idx];
    sd.x = x;
    sd.y = y;
    float frdLength = 0.f;
    for(int i = 0; i < 3; i++) {
    	//bit dangerous to also apply this to output color, but useful for debugging
    	sd.outputColors[i] = sd.inputColors[i] = sd.rgb[i] = xyz[i];
    	sd.rho[i] = rhoXYZ[i];
    	sd.normal[i] = isect.shadingN[i];
    	sd.secondNormal[i] = isect.secondNormal[i];
    	sd.secondOrigin[i] = isect.secondOrigin[i];
    	sd.thirdOrigin[i] = isect.thirdOrigin[i];
    	sd.firstReflectionDir[i] = sd.thirdOrigin[i] - sd.secondOrigin[i];
    	frdLength += sd.firstReflectionDir[i]*sd.firstReflectionDir[i];
    }
    // normalize firstReflectionDir
    frdLength = sqrt(frdLength);
    for (int i = 0; i < 3; i++)
    	sd.firstReflectionDir[i] /= frdLength;
    sd.imgPos[0] = sample.imageX;
    sd.imgPos[1] = sample.imageY;
    sd.lensPos[0] = sample.lensU;
    sd.lensPos[1] = sample.lensV;
	sd.time = sample.time;
}

void RPF::GetAdaptPixels(int spp, vector<vector<int> > &pixels) {
    //THis doesn't happen
}

void RPF::WriteImage(const string &filename, int xres, int yres, bool dump) {

	ProgressReporter reporter(2, "Sorting samples...");
	std::sort(allSamples.begin(), allSamples.end());
	reporter.Done();

    string filenameBase = filename.substr(0, filename.rfind(".")) + "_jouni_" + std::to_string(jouni).substr(2,3);
    string filenameExt  = filename.substr(filename.rfind("."));

    if (dump) {
		dumpAsBinary(filenameBase, xPixelCount, yPixelCount, spp, allSamples);
    }

	RandomParameterFilter rpf(xPixelCount, yPixelCount, spp, jouni, allSamples);
	rpf.setQuality(quality);
    rpf.setRandomParams(randomParams);
	rpf.Apply();

    AssembleImages(dump);

    //printf("Avg spp: %.2f\n", CalculateAvgSpp());

    WriteImage(filenameBase+"_rpf_img"+filenameExt, colImg, xres, yres);
    WriteImage(filenameBase+"_rpf_flt"+filenameExt, fltImg, xres, yres);

    if(dump) { // Write debug images
        // Normals contain negative values, normalize them here
        for(int y = 0; y < normalImg.GetRowNum(); y++)
            for(int x = 0; x < normalImg.GetColNum(); x++) {
                normalImg(x, y) += Color(1.f, 1.f, 1.f);
                normalImg(x, y) /= 2.f;
                secNormalImg(x, y) += Color(1.f, 1.f, 1.f);
                secNormalImg(x, y) /= 2.f;
            }
        WriteImage(filenameBase+"_rpf_nor"+filenameExt, normalImg, xres, yres);

        WriteImage(filenameBase+"_rpf_rho"+filenameExt, rhoImg, xres, yres);

        //new:
        WriteImage(filenameBase+"_rpf_second_normal"+filenameExt, secNormalImg, xres, yres);
        WriteImage(filenameBase+"_rpf_second_orig"+filenameExt, secOrigImg, xres, yres);
        WriteImage(filenameBase+"_rpf_third_orig"+filenameExt, thirdOrigImg, xres, yres);
        WriteImage(filenameBase+"_rpf_lens"+filenameExt, lensImg, xres, yres);
        TwoDArray<Color> timeColImg = FloatImageToColor(timeImg);
        WriteImage(filenameBase+"_rpf_time"+filenameExt, timeColImg, xres, yres);
    }
}

void RPF::WriteImage(const string &filename, const TwoDArray<Color> &image, int xres, int yres) const {
    ::WriteImage(filename, (float*)image.GetRawPtr(), NULL, xPixelCount, yPixelCount,
                 xres, yres, xPixelStart, yPixelStart);
}

TwoDArray<Color> RPF::FloatImageToColor(const TwoDArray<float> &image) const {
    TwoDArray<Color> colorImg(image.GetColNum(), image.GetRowNum());
    for(int y = 0; y < yPixelCount; y++)
        for(int x = 0; x < xPixelCount; x++) {
            float val = image(x, y);
            colorImg(x, y) = Color(val, val, val);
        } 
    return colorImg;
}

void RPF::AssembleImages(bool dump) {
#pragma omp parallel for num_threads(PbrtOptions.nCores)
    for(uint i = 0; i < allSamples.size(); i++) {
		SampleData sd = allSamples[i];

		int x = sd.x;
		int y = sd.y;

		Color rgbC = Color(sd.rgb); //color in RGB
		Color normalC = Color(sd.normal);

		Color rhoC = Color(sd.rho);

		//new
		Color secNormalC = Color(sd.secondNormal);
		Color secOriginC = Color(sd.secondOrigin);
		Color thirdOriginC = Color(sd.thirdOrigin);
		Color lensC = Color(sd.lensPos[0], sd.lensPos[1], 0.f);

		colImg(x, y) += rgbC;
		normalImg(x, y) += normalC;
		rhoImg(x, y) += rhoC;
		//new
		secNormalImg(x, y) += secNormalC;
		secOrigImg(x, y) += secOriginC;
		thirdOrigImg(x, y) += thirdOriginC;
		lensImg(x, y) += lensC;
	}

    for (int y=0; y < yPixelCount; y++) {
    	for (int x = 0; x < xPixelCount; x++) {
    		colImg(x, y) /= spp;
			normalImg(x, y) /= spp;
			rhoImg(x, y) /= spp;
			//new
			secNormalImg(x, y) /= spp;
			secOrigImg(x, y) /= spp;
			thirdOrigImg(x, y) /= spp;
			lensImg(x, y) /= spp;
    	}
    }

    for (uint i=0; i < allSamples.size(); i+=spp) {
    	Color c;
    	for (int j=0; j<spp; j++)
    		for(int k=0; k<3;k++){
    			c[k] += allSamples[i+j].outputColors[k];
    		}
    	c /= spp;
    	fltImg(allSamples[i].x, allSamples[i].y) = c;
    }

}

void RPF::dumpAsBinary(const string &filenameBase, const int w, const int h,
		const int spp, const vector<SampleData> &allSamples) {
	std::ofstream dump(filenameBase + ".bin",
			std::ifstream::out | std::ifstream::binary);
	// DUMP NUMBER allSamples.size()
	dump.write((char*) (&w), sizeof(int));
	dump.write((char*) (&h), sizeof(int));
	dump.write((char*) (&spp), sizeof(int));
	dump.write((char*) (&(allSamples[0])),
			allSamples.size() * sizeof(SampleData));
	dump.close();
}

