#include <thrust/sort.h>
#include <thrust/device_vector.h>
#include <thrust/copy.h>
#include "cuda_device.h"
#include <cmath>

void sort_on_device(std::vector<SampleData>& h_vec)
{
    // transfer data to the device
    thrust::device_vector<SampleData> d_vec(h_vec);

    // sort data on the device
    thrust::sort(d_vec.begin(), d_vec.end());

    // transfer data back to host
    thrust::copy(d_vec.begin(), d_vec.end(), h_vec.begin());
}

struct square_functor : public thrust::unary_function<SampleData,SampleData>
{
	__device__
	SampleData operator()(const SampleData &sd1) {
	SampleData s;
	for(int f=0; f<SampleData::getLastNormalizedOffset(); f++) {
		s[f] = sd1[f]*sd1[f];
	}
	return s;
}
};

struct norm : public thrust::unary_function<SampleData, SampleData>
{
	const SampleData mean, std;

	norm(SampleData mean, SampleData std) : mean(mean), std(std) {}

	__device__
	SampleData operator()(const SampleData &sd1) {
		SampleData s;
		//s = sd1 - mean;
		//s /= std;
	return s;
}
};

void normalize(std::vector<SampleData>& v) {
	// transfer data to the device
	thrust::device_vector<SampleData> d_v(v);

	SampleData mean, meanSquare, std;
	mean.reset(); meanSquare.reset();
	mean = thrust::reduce(d_v.begin(), d_v.end(), mean);
	meanSquare = thrust::transform_reduce(d_v.begin(), d_v.end(), square_functor(), meanSquare, thrust::plus<SampleData>());
	mean.divide(v.size());
	meanSquare.divide(v.size());
	for (int f=0; f < SampleData::getLastNormalizedOffset(); f++) {
		std[f] = sqrt(meanSquare[f] - mean[f]*mean[f]);
	}
	//TODO?
	thrust::transform(d_v.begin(), d_v.end(), d_v.begin(), d_v.end(), norm(mean, std));

	thrust::copy(d_v.begin(), d_v.end(), v.begin());
}
