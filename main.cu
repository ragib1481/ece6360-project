#include <iostream>
#include <thrust/host_vector.h>
#include <complex.h>
#include <chrono>

#include "Beam.cuh"
#include "helper.h"
#include "Propagator.cuh"

using namespace std;
using namespace std::chrono;


int main(int argc, char* argv[]) {
    // parse command line arguments
    if (argc != 2)
        return 1;
    string filename(argv[1]);

    // define simulation parameters.
    Parameters parameters("./" + filename);

    // get beam at the focus with the simulation parameters.
    auto start = high_resolution_clock::now();
    beams::cpu::GaussianBeam beam(parameters);
    auto end = high_resolution_clock::now();
    cout << "(CPU) Beam At Focus :                   " << (duration_cast<milliseconds >(end-start)).count() << "ms" << endl;

    // plot beam at the focus
    start = high_resolution_clock::now();
    beams::cpu::plotBeamFromFourier<double>(beam.getbeamAtFocusTilde(),
                                            parameters.numSamples.y,
                                            parameters.numSamples.x,
                                            "Beam At Focus");
    end = high_resolution_clock::now();
    cout << "(CPU) Plot Beam :                       " << (duration_cast<milliseconds >(end-start)).count() << "ms" << endl;

    // declare fresnel propagator with simulation parameters.
    propagator::cpu::Propagator propagator(parameters, beam.getbeamAtFocusTilde());

    // perform propagation.
    start = high_resolution_clock::now();
    auto beamTilde = propagator.propagate(128);
    end = high_resolution_clock::now();
    cout << "(CPU) Propagation to single z plane :   " << (duration_cast<milliseconds >(end-start)).count() << "ms" << endl;

    // plot beam at the propagated distance
    beams::cpu::plotBeamFromFourier<std::complex<double>>(beamTilde,
                                                            parameters.numSamples.y,
                                                            parameters.numSamples.x,
                                                            "Beam At 100um");


    // compute z propagation profile
    start = high_resolution_clock::now();
    auto zProfile = propagator.propagationProfile(parameters.zStart, parameters.zEnd, parameters.zSampling);
    end = high_resolution_clock::now();
    cout << "(CPU) Propagation Profile:              " << (duration_cast<milliseconds >(end-start)).count() << "ms" << endl;
    beams::cpu::plotBeamFromSpatial(zProfile,
                                    ((parameters.zEnd-parameters.zStart)/parameters.zSampling),
                                    parameters.numSamples.x, "z profile");

    return 0;
}
