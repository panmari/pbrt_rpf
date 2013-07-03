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
	float secondOrigin[3];  //0, 1, 2
	float thirdOrigin[3];	//3, 4, 5
	//features
	float normal[3];		//6, 7, 8
	float secondNormal[3];	//9, 10, 11
	float rho[3];			//12, 13,14

	float rgb[3];			//15, 16, 17
	float imgPos[2];		//18, 19

	//random parameters
	float lensPos[2];		//20, 21
	float time;				//22

	// input/output colors:
	float inputColors[3];	//whatev
	float outputColors[3];
	int x, y;
	// Some handy accessor methods, thx @jklethinen
	// asserts that float and int are the same length
	static int getSize()					{ return sizeof(SampleData)/sizeof(float); }
	static int getFeaturesStart() 			{ return 0; }
	static int getFeaturesEnd()				{ return 15; }
	static int getRandomParametersStart()	{ return 20; }
	static int getRandomParametersEnd()		{ return 22; }
	//color can be accessed through .rgb
	void   operator+=(const SampleData& s)	{ for(int i=0;i<getSize();i++) (*this)[i] += s[i]; }
	void   divide(int s)					{ for(int i=0;i<getSize();i++) (*this)[i] /= float(s); }
	float& operator[](int i)				{ return ((float*)this)[i]; }
	const float& operator[](int i) const	{ return ((float*)this)[i]; }
	float sum() const						{ float s=0; for(int i=0;i<getSize();i++) s+=(*this)[i]; return s; }
	float avg() const						{ return sum()/getSize(); }
};

#endif //SAMPLE_DATA_H
