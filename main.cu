#include <iostream>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <complex.h>
#include <chrono>

#include "helper.h"
#include "Scatter.cuh"

using namespace std;
using namespace std::chrono;

int main(int argc, char* argv[]) {
    // parse command line arguments
    if (argc != 2)
        return 1;
    string filename(argv[1]);

    // define simulation parameters.
    Parameters parameters("./" + filename);

    // load all simulation parameters from file
    // pass the simulation parameters to the MieScatter object for simulation
    // save simulation results

    return 0;
}
