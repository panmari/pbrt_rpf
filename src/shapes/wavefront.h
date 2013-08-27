
/*
    pbrt source code Copyright(c) 1998-2010 Matt Pharr and Greg Humphreys.

    This file is part of pbrt.

    pbrt is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.  Note that the text contents of
    the book "Physically Based Rendering" are *not* licensed under the
    GNU GPL.

    pbrt is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

 */

#ifndef PBRT_SHAPES_WAVEFRONT_H
#define PBRT_SHAPES_WAVEFRONT_H

// shapes/sphere.h*
#include "trianglemesh.h"
#include "shape.h"
#include "paramset.h"

// Sphere Declarations
class Wavefront : public Shape {
	public:
		Wavefront(const char* filename, const Transform *o2w, const Transform *w2o, bool reverseOrientation) ;
		~Wavefront() {
			if (vertexIndex) delete[] vertexIndex;
			if (p) delete[] p;
			if (n) delete[] n;
			if (uvs) delete[] uvs;
		}
		
		BBox ObjectBound() const;
		BBox WorldBound() const;
		bool CanIntersect() const { return false; }
		void Refine(vector<Reference<Shape> > &refined) const;

	private:
		void MergeIndicies(vector<Point> &points, vector<Normal> &normals, vector<float> &uvVec,
							vector<int> &vIndex, vector<int> &normalIndex, vector<int> &uvIndex);
		int ntris, nverts;
		int nVerticesTotal;
		int *vertexIndex;
		Point *p;
		Normal *n;
		float *uvs;
};


Wavefront *CreateWaveFrontShape(const Transform *o2w, const Transform *w2o,
        bool reverseOrientation, const ParamSet &params);

#endif 
