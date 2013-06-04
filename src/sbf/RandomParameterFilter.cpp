
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


RandomParameterFilter::RandomParameterFilter() {

}

const int BOX_SIZE[] = {55, 35, 17, 7};
const float MAX_SAMPLES_FACTOR = {0.02f, 0.04f, 0.3f, 0.5f};

void RandomParameterFilter::Apply(const vector<SampleData> &allSamples, int w, int h, int spp) const {
	ProgressReporter reporter(4, "Applying RPF filter");
	for (int iterStep = 0; iterStep < 4; iterStep++) {
		reporter.Update(iterStep);
		for (int pixel = 0; pixel < w*h; pixel++) {
			SampleData &firstSamplePixel = allSamples[pixel*spp];
			determineNeighbourhood(BOX_SIZE[iterStep], w, h, firstSamplePixel);

		}
	}

	reporter.Done();
}
