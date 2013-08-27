
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


// textures/scale.cpp*
#include "textures/divide.h"

// ScaleTexture Method Definitions
DivideTexture<float, float> *CreateDivideFloatTexture(const Transform &tex2world,
        const TextureParams &tp) {
    return new DivideTexture<float, float>(tp.GetFloatTexture("tex1", 1.f),
        tp.GetFloatTexture("tex2", 1.f));
}



DivideTexture<Spectrum, Spectrum> *CreateDivideSpectrumTexture(const Transform &tex2world,
        const TextureParams &tp) {
    return new DivideTexture<Spectrum, Spectrum>(
        tp.GetSpectrumTexture("tex1", Spectrum(1.f)),
        tp.GetSpectrumTexture("tex2", Spectrum(1.f)));
}


