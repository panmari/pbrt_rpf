#include <thrust/sort.h>
#include <thrust/device_vector.h>
#include <thrust/copy.h>

#include "cuda_device.h"

void sort_on_device(std::vector<SampleData>& h_vec)
{
    // transfer data to the device
    thrust::device_vector<SampleData> d_vec(h_vec);

    // sort data on the device
    thrust::sort(d_vec.begin(), d_vec.end());

    // transfer data back to host
    thrust::copy(d_vec.begin(), d_vec.end(), h_vec.begin());
}

struct plusSquare : public thrust::binary_function<SampleData,SampleData,SampleData>
{
__device__ SampleData operator()(const SampleData &sd1, const SampleData &sd2) {
	SampleData s;
	for(int f=0; f<SampleData::getLastNormalizedOffset(); f++) {
		s[f] = sd1[f]*sd1[f] + sd2[f]*sd2[f];
	}
	return s;
}
};

struct minusMean : public thrust::unary_function
{
__device__ SampleData operator()(const SampleData &sd1, const SampleData &sd2) {
	SampleData s;
	for(int f=0; f<SampleData::getLastNormalizedOffset(); f++) {
		s[f] = sd1[f]*sd1[f] + sd2[f]*sd2[f];
	}
	return s;
}
};

void normalize(std::vector<SampleData>& v) {
	// transfer data to the device
	thrust::device_vector<SampleData> d_v(v);

	SampleData mean, meanSquare;
	mean.reset(); meanSquare.reset();
	mean = thrust::reduce(d_v.begin(), d_v.end(), mean);
	meanSquare = thrust::reduce(d_v.begin(), d_v.end(), meanSquare, plusSquare());
	thrust::device_vector<SampleData> temp(v.size());
	thrust::transform(d_v.begin(), d_v.end(), d_v, minusMean(mean));

	thrust::copy(d_v.begin(), d_v.end(), v.begin());
}
