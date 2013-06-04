#ifndef SAMPLE_DATA_H
#define SAMPLE_DATA_H

struct SampleData {
            SampleData() {
                for(int i = 0; i < 3; i++) {
                    rgb[i] =
					normal[i] = secondNormal[i] =
					rho[i] =
					secondOrigin[i] = thirdOrigin[i] =
					inputColors[i] = outputColors[i] = 0.f;
                }
                lensPos[0] = lensPos[1] = time = 0.f;
                x = y = 0;
            }
            int x, y;

            float rgb[3];
            //features
            float normal[3];
            float secondNormal[3];
            float rho[3];
            float secondOrigin[3];
            float thirdOrigin[3];

            float imgPos[2];

            //random parameters
            float lensPos[2];
            float time;

            // input/output colors:
            float inputColors[3];
            float outputColors[3];
};

#endif //SAMPLE_DATA_H
