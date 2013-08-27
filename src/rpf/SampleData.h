#ifndef SAMPLE_DATA_H
#define SAMPLE_DATA_H

struct SampleData {

	void reset() {
		for (int i = 0; i < 3; i++) {
			rgb[i] = normal[i] = secondNormal[i] = rho[i] = secondOrigin[i] =
					thirdOrigin[i] = inputColors[i] = outputColors[i] =
					firstReflectionDir[i] = 0.f;
		}
		for (int i = 0; i < 2; i++) {
			lensPos[i] = imgPos[i]= 0.f;
		}
		time = 0.f;
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
	float firstReflectionDir[3]; //	20, 21, 22
	float lensPos[2];		//		23, 24
	float time;				//		25

	// input/output colors:
	float inputColors[3];	//whatev
	float outputColors[3];
	int x, y;
#define FEATURES_OFFSET 			0
#define FEATURES_SIZE 				15
#define COLOR_OFFSET				15
#define COLOR_SIZE					3
#define IMG_POS_OFFSET				18
#define IMG_POS_SIZE				2
#define LAST_NORMALIZED_OFFSET		26

	// Some handy accessor methods, thx @jlehtinen
	// asserts that float and int are the same length
	static uint getSize()					{ return sizeof(SampleData)/sizeof(float); }
	void   operator+=(const SampleData& s)	{ for(uint i=0;i<getSize();i++) (*this)[i] += s[i]; }
	void   divide(int s)					{ for(uint i=0;i<getSize();i++) (*this)[i] /= float(s); }
	float& operator[](int i)				{ return ((float*)this)[i]; }
	const float& operator[](int i) const	{ return ((float*)this)[i]; }
	float sum() const						{ float s=0; for(uint i=0;i<getSize();i++) s+=(*this)[i]; return s; }
	float avg() const						{ return sum()/getSize(); }

	/**
	 * For ordering in sbf.h. Hopefully, this produces the same ordering every time.
	 */
	bool operator <(const SampleData& other) const {
		if (y == other.y)
			return imgPos[0] < other.imgPos[0];
		else return imgPos[1] < other.imgPos[1];
	}
};

#endif //SAMPLE_DATA_H
