#include <iostream>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <complex.h>
#include <chrono>

#include "helper.h"
#include "myMath.h"

using namespace std;
using namespace std::chrono;

__global__
void calc(float* out, const float* in, const int N, const int n) {
    size_t ix = blockIdx.x * blockDim.x + threadIdx.x;
    float vm;
    float* jv  = new float[n+1];
    float* yv  = new float[n+1];
    float* djv = new float[n+1];
    float* dyv = new float[n+1];

    if (ix < N) {
        redefined::sphBessjyv<float>(n, in[ix], vm, &jv[0], &yv[0], &djv[0], &dyv[0]);
        out[ix] = jv[n];
    }
}


int main(int argc, char* argv[]) {
    // parse command line arguments
    if (argc != 2)
        return 1;
    string filename(argv[1]);

    // define simulation parameters.
    Parameters parameters("./" + filename);

    const int N = 10;
    thrust::host_vector<float> in(N);
    thrust::host_vector<float> out(N);

    thrust::device_vector<float> inGpu(N);
    float* inGpuPtr = thrust::raw_pointer_cast(inGpu.data());

    thrust::device_vector<float> outGpu(N);
    float* outGpuPtr = thrust::raw_pointer_cast(outGpu.data());

    for (int n = 0; n < 3; n++) {
        for (int i = 0; i < in.size(); i++) {
            in[i] = static_cast<float>(i + 0.01);
        }
        inGpu = in;
        calc<<<N, 1>>>(outGpuPtr, inGpuPtr, N, n);
        out = outGpu;

        for (int i = 0; i < in.size(); i++) {
            cout << out[i] - std::sph_bessel<float>(n, in[i]) << endl;
        }
        cout << endl;
    }


    // cout << "Legendre Polynomial" << endl;
    // for (unsigned int n = 0; n <= 1; n++) {
    //     for (float x = 0.0; x <= 10.0; x++) {
    //         cout << redefined::legendre<float>(n, x) - std::legendre<float>(n, x) << ", ";
    //     }
    //     cout << endl;
    // }
    // cout << endl;
    // cout << endl;

    // cout << "hankel function" << endl;

    return 0;
}
