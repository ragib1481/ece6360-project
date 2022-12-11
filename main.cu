#include <iostream>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <string>
#include <complex.h>
#include <chrono>

#include "helper.cuh"
#include "Scatter.cuh"
#include "Vec2.cuh"

using namespace std;
using namespace std::chrono;


int main(int argc, char* argv[]) {
    // parse command line arguments
    if (argc != 2)
        return 1;
    string filename(argv[1]);

    cout << filename << endl;
    // define simulation parameters.
    Parameters<float> parameters("./" + filename);

    scatter::cpu::MieScatter<float> scatterer(parameters);
    scatterer.scatter();

    // load all simulation parameters from file
    // pass the simulation parameters to the MieScatter object for simulation
    // save simulation results

    return 0;
}
