#ifndef SAMPLE_DATA_H
#define SAMPLE_DATA_H

struct SampleData {
	SampleData() {
		for (int i = 0; i < 3; i++) {
			rgb[i] = normal[i] = secondNormal[i] = rho[i] = secondOrigin[i] =
					thirdOrigin[i] = inputColors[i] = outputColors[i] = 0.f;
		}
		lensPos[0] = lensPos[1] = time = 0.f;
		x = y = 0;
	}
	// position features (the first 6 values)
	float secondOrigin[3];
	float thirdOrigin[3];

	int x, y;

	float rgb[3];
	//features
	float normal[3];
	float secondNormal[3];
	float rho[3];

	float imgPos[2];
	//random parameters
	float lensPos[2];
	float time;

	// input/output colors:
	float inputColors[3];
	float outputColors[3];

	// Some handy accessor methods, thx @jklethinen
	// asserts that float and int are the same length
	static int getSize()					{ return sizeof(SampleData)/sizeof(float); }\
	void   operator+=(const SampleData& s)	{ for(int i=0;i<getSize();i++) (*this)[i] += s[i]; } \
	void   divide(int s)					{ for(int i=0;i<getSize();i++) (*this)[i] /= float(s); } \
	float& operator[](int i)				{ return ((float*)this)[i]; } \
	const float& operator[](int i) const	{ return ((float*)this)[i]; } \
	float sum() const						{ float s=0; for(int i=0;i<getSize();i++) s+=(*this)[i]; return s; } \
	float avg() const						{ return sum()/getSize(); } \
};

#endif //SAMPLE_DATA_H
