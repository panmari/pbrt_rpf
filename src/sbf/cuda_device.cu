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
    thrust::sort(d_vec.begin(), d_vec.end(), compSD());

    // transfer data back to host
    thrust::copy(d_vec.begin(), d_vec.end(), h_vec.begin());
}

