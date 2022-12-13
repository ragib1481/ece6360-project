#include <iostream>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <string>
#include <complex.h>
#include <chrono>

#include "Scatter.cuh"
#include "helper.cuh"

using namespace std;
using namespace std::chrono;


int main(int argc, char* argv[]) {
    // parse command line arguments
    if (argc != 2)
        return 1;
    string filename(argv[1]);

    cout << filename << endl;
    // load all simulation parameters from file
    Parameters<float> parameters("./" + filename);
    parameters.print();

    // pass the simulation parameters to the MieScatter object for simulation
    scatter::cpu::MieScatter<float> scatterer(parameters);

    // run the simulation
    auto start = high_resolution_clock::now();
    scatterer.scatter();
    auto end = high_resolution_clock::now();
    cout << "cpu elapsed time: " << (chrono::duration_cast<milliseconds >(end-start)).count() << endl;

    // save simulation results
    scatterer.saveResult();



    return 0;
}
