#include <thrust/sort.h>
#include <thrust/device_vector.h>
#include <thrust/copy.h>

#include "cuda_device.h"

struct compSD : public thrust::binary_function<SampleData,SampleData,bool>
{
__host__ __device__ bool operator()(const SampleData &sd1, const SampleData &sd2) {
	if (sd1.y == sd2.y)
		return sd1.x < sd2.x;
	else return sd1.y < sd2.y;
}
};

void sort_on_device(std::vector<SampleData>& h_vec)
{
    // transfer data to the device
    thrust::device_vector<SampleData> d_vec(h_vec);

    // sort data on the device
    thrust::sort(d_vec.begin(), d_vec.end());

    // transfer data back to host
    thrust::copy(d_vec.begin(), d_vec.end(), h_vec.begin());
}

struct plusSD : public thrust::binary_function<SampleData,SampleData,SampleData>
{
__host__ __device__ SampleData operator()(const SampleData &sd1, const SampleData &sd2) {
	//TODO
	SampleData s;
/*
	for(int f=0; f<25; f++) {
		(float*(&s))[f] = (float*(&sd1))[f] + (float*(&sd2))[f];
	}
	*/
	return s;
}
};

void normalize(std::vector<SampleData>& v) {
	// transfer data to the device
	thrust::device_vector<SampleData> d_v(v);

	SampleData mean, meanSquare;
	mean.reset(); meanSquare.reset();
	thrust::reduce(d_v.begin(), d_v.end(), mean, plusSD());

	thrust::copy(d_v.begin(), d_v.end(), v.begin());
}
