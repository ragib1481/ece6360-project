#include <iostream>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <string>
#include <complex.h>
#include <chrono>

#ifndef Nl
#define Nl 30
#endif

#include "Scatter.cuh"
#include "helper.cuh"

using namespace std;
using namespace std::chrono;

template <typename T>
void simulateCpu(Parameters<T> parameters) {
    scatter::cpu::MieScatter<T> scatterer(parameters);

    // run the simulation
    auto start = high_resolution_clock::now();
    scatterer.scatter();
    auto end = high_resolution_clock::now();
    cout << "cpu elapsed time: " << (chrono::duration_cast<milliseconds >(end-start)).count() << "ms" << endl;

    // save simulation results
    scatterer.saveResult("Cpu");
}

template <typename T>
void simulateGpu(Parameters<T> parameters) {
    scatter::gpu::MieScatter<T> scatterer(parameters);

    // run the simulation
    auto start = high_resolution_clock::now();
    scatterer.scatter();
    auto end = high_resolution_clock::now();
    cout << "gpu elapsed time: " << (chrono::duration_cast<milliseconds >(end-start)).count() << "ms" << endl;

    // save simulation results
    scatterer.saveResult("Gpu");
}

int main(int argc, char* argv[]) {
    // parse command line arguments
    if (argc != 2)
        return 1;
    string filename(argv[1]);

    cout << filename << endl;
    // load all simulation parameters from file
    Parameters<float> parameters("./" + filename);
    parameters.print();

    // run the simulation on cpu
    simulateCpu<float>(parameters);

    // run the simulation of gpu
    simulateGpu<float>(parameters);

    return 0;
}
