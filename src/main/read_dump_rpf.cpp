/*
 * File:   rec.cpp
 * Author: fabrice
 * Modified: sm
 */

#include <fstream>

#include "core/pbrt.h"
#include "rpf/RandomParameterFilter.h"
#include "rpf/SampleData.h"
#include "core/imageio.h"
#include "filter_utils/VectorNf.h"
using namespace std;

int main(int argc, char** argv)
{
    if (argc == 1)
        Severe("No base name provided!");
    string filename(argv[1]);

    vector<SampleData> allSamples;
    int w, h, spp;
    std::ifstream dump(filename, std::ifstream::in | std::ifstream::binary);

    dump.read((char*)&w, sizeof(int));
	dump.read((char*)&h, sizeof(int));
	dump.read((char*)&spp, sizeof(int));
	allSamples.resize(w*h*spp);
	dump.read((char*)&(allSamples[0]), allSamples.size() * sizeof(SampleData));
	dump.close();

    RandomParameterFilter rpf(w, h, spp, 0.02f, RandomParameterFilter::Quality::MEDIUM, allSamples);
    rpf.Apply();

    TwoDArray<Color> fltImg = TwoDArray<Color>(w, h);
    // Dumping img (and multiply with rho/albedo
    for (uint i=0; i < allSamples.size(); i+=spp) {
    	Color c;
    	for (int j=0; j<spp; j++)
    		for(int k=0; k<3;k++){
    			c[k] += allSamples[i+j].outputColors[k];
    		}
    	c /= spp;
    	fltImg(allSamples[i].x, allSamples[i].y) = c;
    }
    string filenameBase = filename.substr(0, filename.rfind("."));
    WriteImage(filenameBase + "_flt.exr", (float*)fltImg.GetRawPtr(), NULL, w, h,
                     w, h, 0, 0);

    return 0;
}
