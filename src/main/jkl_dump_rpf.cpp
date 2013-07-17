/*
 * File:   rec.cpp
 * Author: fabrice
 * Modified: sm
 * Created on 25. avril 2011, 07:19
 */
#pragma GCC diagnostic ignored "-Wunused-result"

#include <cstdlib>
#include <iostream>
#include <fstream>

#include "core/pbrt.h"
#include "sbf/RandomParameterFilter.h"
#include "sbf/SampleData.h"
#include "core/imageio.h"

using namespace std;

void readDump(char** argv, vector<SampleData> &allSamples, int &w, int &h, int &spp) {
	string filename(argv[1]);
	int	m_width, m_height, m_numSamplesPerPixel;
	float m_version;
	/**
	 * largely taken from jkl, since it's his dump format
	 */
	FILE* fp  = fopen(filename.c_str(), "rb");

	FILE* fph = fopen((filename + ".header").c_str(),"rt");
	bool separateHeader = (fph!=NULL);
	if(!separateHeader)
		fph = fp;

	printf("Importing sample buffer... ");

	// Parse version and image size.

	fscanf(fph, "Version %f\n", &m_version);
	fscanf(fph, "Width %d\n", &m_width);
	fscanf(fph, "Height %d\n", &m_height);
	fscanf(fph, "Samples per pixel %d\n", &m_numSamplesPerPixel);
	printf(" Version: %.1f, size: %dx%d, spp: %d \n", m_version, m_width, m_height, m_numSamplesPerPixel);
	allSamples.resize(m_width*m_height*m_numSamplesPerPixel);
	w = m_width;
	h = m_height;
	spp = m_numSamplesPerPixel;

	//only handle this one version
	if(m_version == 2.2f) {
		// Parse the rest of the header.

		for (uint i=0; i< allSamples.size(); i++)
		{
			SampleData &s = allSamples[i];
			fread(&(s.imgPos),sizeof(float),2,fp); //x and y
			//fread(0 , sizeof(int),1,fp); //w, don't need this value
			fseek(fp, sizeof(int), SEEK_CUR);
			fread(&(s.lensPos),sizeof(float),2,fp); //u, v
			fread(&(s.time), sizeof(float), 1, fp); // t
			fread(&(s.rgb), sizeof(float), 3, fp); // rgb
			//fread(0, sizeof(float), 3, fp); // pri_mv
			fseek(fp, sizeof(float)*3, SEEK_CUR);
			fread(&(s.normal), sizeof(float), 3, fp); // pri_normal
			fread(&(s.rho), sizeof(float), 3, fp); // albedo
			fread(&(s.secondOrigin), sizeof(float), 3, fp); // sec_origin
			fread(&(s.thirdOrigin), sizeof(float), 3, fp); // sec_hitpoint
			//fread(0, sizeof(float), 3, fp); // sec_mv
			fseek(fp, sizeof(float)*3, SEEK_CUR);
			fread(&(s.secondNormal), sizeof(float), 3, fp); // sec_normal
			fseek(fp, sizeof(float)*9, SEEK_CUR);
			/*
			fread(0 , sizeof(float), 3, fp); // direct
			fread(0, sizeof(float), 3, fp); // sec_albedo
			fread(0, sizeof(float), 3, fp); // sec_direct
			*/
			s.x = (int)s.imgPos[0];
			s.y = (int)s.imgPos[1];
			float frdLength = 0.f;
			for (int i = 0; i < 3; i++) {
				s.firstReflectionDir[i] = s.thirdOrigin[i] - s.secondOrigin[i];
				frdLength += s.firstReflectionDir[i]*s.firstReflectionDir[i];
				s.inputColors[i] = s.outputColors[i] = s.rgb[i];
			}
			//normalize frd
			frdLength = sqrt(frdLength);
			for (int i=0; i<3; i++)
				s.firstReflectionDir[i] /= frdLength;
		}
	} else {
		printf("unsupported header format");
	}
	fclose(fp);
	if(separateHeader)
		fclose(fph);
	printf("done\n");
	SampleData &s = allSamples[(500*m_width + 500)*m_numSamplesPerPixel + 2];
	for (int f=0; f<SampleData::getLastNormalizedOffset(); f++)
		printf("%-.3f \n", s[f]);
	printf("pos: %d, %d \n", s.x, s.y);
}

int main(int argc, char** argv)
{
    if (argc == 1)
        Severe("No base name provided!");

    vector<SampleData> allSamples;
    int w, h, spp;
    readDump(argv, allSamples, w, h, spp);
    RandomParameterFilter rpf(w, h, spp, allSamples);
    rpf.Apply();

    TwoDArray<Color> fltImg = TwoDArray<Color>(w, h);
    // Dumping img (and multiply with rho/albedo
    for (uint i=0; i < allSamples.size(); i+=spp) {
    	Color c;
    	for (int j=0; j<spp; j++)
    		for(int k=0; k<3;k++){
    			c[k] += allSamples[i+j].outputColors[k]*allSamples[i+j].rho[k];
    		}
    	c /= spp;
    	fltImg(allSamples[i].x, allSamples[i].y) = c;
    }
    ::WriteImage("test.exr", (float*)fltImg.GetRawPtr(), NULL, w, h,
                     w, h, 0, 0);
   // SBF::WriteImage(string(argv[1]) + "_sbf_img.exr", fltImg, w, h);

    return 0;
}
