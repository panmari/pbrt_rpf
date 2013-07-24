
/*
    Copyright(c) 1998-2012 Matt Pharr and Greg Humphreys.

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


// film/sbfimage.cpp*
#include "stdafx.h"
#include "spectrum.h"
#include "parallel.h"
#include "imageio.h"
#include "intersection.h"
#include "progressreporter.h"
#include "rpfimage.h"
#include "rpf/rpf.h"

// SBFImageFilm Method Definitions
RPFImageFilm::RPFImageFilm(int xres, int yres, Filter *filt, const float crop[4],
                     const string &fn, bool dp, float jouni)
    : Film(xres, yres) {
	filter = filt;
    memcpy(cropWindow, crop, 4 * sizeof(float));
    filename = fn;
    // Compute film image extent
    xPixelStart = Ceil2Int(xResolution * cropWindow[0]);
    xPixelCount = max(1, Ceil2Int(xResolution * cropWindow[1]) - xPixelStart);
    yPixelStart = Ceil2Int(yResolution * cropWindow[2]);
    yPixelCount = max(1, Ceil2Int(yResolution * cropWindow[3]) - yPixelStart);

    dump = dp;
    rpf = new RPF(xPixelStart, yPixelStart, xPixelCount, yPixelCount, jouni);
}


void RPFImageFilm::AddSample(const CameraSample &sample,
                          const Spectrum &L,
                          const Intersection &isect) {
    rpf->AddSample(sample, L, isect);
}


void RPFImageFilm::Splat(const CameraSample &sample, const Spectrum &L) {
    // TODO: Implement splatting
    Warning("[SBFImageFilm] Splatting is currently not supported");
}


void RPFImageFilm::GetSampleExtent(int *xstart, int *xend,
                                int *ystart, int *yend) const {
    *xstart = Floor2Int(xPixelStart + 0.5f - filter->xWidth);
    *xend   = Ceil2Int(xPixelStart + 0.5f + xPixelCount  +
                        filter->xWidth);

    *ystart = Floor2Int(yPixelStart + 0.5f - filter->yWidth);
    *yend   = Ceil2Int(yPixelStart + 0.5f + yPixelCount +
                        filter->yWidth);
}


void RPFImageFilm::GetPixelExtent(int *xstart, int *xend,
                               int *ystart, int *yend) const {
    *xstart = xPixelStart;
    *xend   = xPixelStart + xPixelCount;
    *ystart = yPixelStart;
    *yend   = yPixelStart + yPixelCount;
}


void RPFImageFilm::WriteImage(float splatScale) {
    rpf->WriteImage(filename, xResolution, yResolution, dump);
}

void RPFImageFilm::GetAdaptPixels(int spp, vector<vector<int> > &pixels) {
    rpf->GetAdaptPixels(spp, pixels);
}

RPFImageFilm *CreateRPFImageFilm(const ParamSet &params, Filter *filt) {
    string filename = params.FindOneString("filename", PbrtOptions.imageFile);
    if (filename == "")
#ifdef PBRT_HAS_OPENEXR
        filename = "pbrt.exr";
#else
        filename = "pbrt.tga";
#endif

    int xres = params.FindOneInt("xresolution", 640);
    int yres = params.FindOneInt("yresolution", 480);
    if (PbrtOptions.quickRender) xres = max(1, xres / 4);
    if (PbrtOptions.quickRender) yres = max(1, yres / 4);
    float crop[4] = { 0, 1, 0, 1 };
    int cwi;
    const float *cr = params.FindFloat("cropwindow", &cwi);
    if (cr && cwi == 4) {
        crop[0] = Clamp(min(cr[0], cr[1]), 0., 1.);
        crop[1] = Clamp(max(cr[0], cr[1]), 0., 1.);
        crop[2] = Clamp(min(cr[2], cr[3]), 0., 1.);
        crop[3] = Clamp(max(cr[2], cr[3]), 0., 1.);
    }
    bool debug = params.FindOneBool("dumpfeaturebuffer", false);

    float jouni = params.FindOneFloat("jouni", 0.02f);
    //TODO: add more parameters?

    return new RPFImageFilm(xres, yres, filt, crop, filename,
                            debug, jouni);
}

