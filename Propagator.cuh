//
// Created by ragib1481 on 11/28/22.
//

#ifndef PROJECT_PROPAGATOR_CUH
#define PROJECT_PROPAGATOR_CUH

#include <fstream>
#include <string>
#include <thrust/host_vector.h>
#include <complex.h>
#include <fftw3.h>
#include <opencv2/opencv.hpp>

#include "helper.h"

template <typename T>
struct doublet {
    T x = 0; T y = 0;
};

class Parameters {
public:
    doublet<unsigned short> dim;                        // dimension of sample along x and y
    doublet<double> sampling;                            // sampling interval along x and Y
    doublet<unsigned short> numSamples;                 // samples along x and y
    double na;                                           // numerical aperture of focusing optics
    double k;                                            // wavenubmer for simulation
    double zStart;
    double zEnd;
    double zSampling;

    // default constructor
    Parameters() {
        dim.x           = 0;
        dim.y           = 0;
        sampling.x      = 0;
        sampling.y      = 0;
        numSamples.x    = 0;
        numSamples.y    = 0;
        na              = 0;
        k               = 0;
        zStart          = 0;
        zEnd            = 0;
        zSampling       = 0;
    }

    // construct object from specified parameter file
    Parameters(std::string filename) {
        std::string discard;
        std::ifstream infile;
        infile.open(filename);
        if (infile.is_open()) {
            infile >> discard;
            infile >> dim.x;
            infile >> dim.y;

            infile >> discard;
            infile >> sampling.x;
            infile >> sampling.y;

            numSamples.x = (unsigned short) ((double)dim.x / sampling.x);
            numSamples.y = (unsigned short) ((double)dim.y / sampling.y);

            infile >> discard;
            infile >> na;

            infile >> discard;
            infile >> k;                                                            // load wavelength
            k = 2 * M_PI / k;                                                       // wavelength to wave-number

            infile >> discard;
            infile >> zStart;

            infile >> discard;
            infile >> zEnd;

            infile >> discard;
            infile >> zSampling;

            infile.close();
        }
    };

    // copy constructor
    Parameters(const Parameters& params) {
        dim              = params.dim;
        sampling         = params.sampling;
        numSamples       = params.numSamples;
        na               = params.na;
        k                = params.k;
        zStart           = params.zStart;
        zEnd             = params.zEnd;
        zSampling        = params.zSampling;
    }

    // overload assignment operator
    Parameters& operator=(const Parameters& params) {
        dim                 = params.dim;
        sampling            = params.sampling;
        numSamples          = params.numSamples;
        na                  = params.na;
        k                   = params.k;
        zStart              = params.zStart;
        zEnd                = params.zEnd;
        zSampling           = params.zSampling;

        return *this;
    }
};


namespace propagator {

    namespace cpu {

        class Propagator {
            Parameters parameters;
            thrust::host_vector<double> kxx;             // spatial frequencies along x
            thrust::host_vector<double> kyy;             // spatial frequencies along y
            thrust::host_vector<double> beamAtFocusTilde;     // spatial beam profile at focus
            int nRows;
            int nCols;

        public:
            Propagator(Parameters parameters, thrust::host_vector<double> beamAtFocusTilde) {
                this->parameters = parameters;

                // generate spatial frequencies
                kxx = fftfreq(parameters.numSamples.x, parameters.sampling.x, 2 * M_PI);
                kyy = fftfreq(parameters.numSamples.y, parameters.sampling.y, 2 * M_PI);

                nCols = kxx.size();
                nRows = kyy.size();

                // set beam
                this->beamAtFocusTilde = beamAtFocusTilde;
            }

            thrust::host_vector<std::complex<double>> propagate(double z) {
                thrust::host_vector<std::complex<double>> beamTilde(kxx.size() * kyy.size());

                double kx;
                double ky;
                double kz;

                // multiply the beam at focus with the accumulated phase at z
                for (int i = 0; i < kyy.size(); i++) {
                    for (int j = 0; j < kxx.size(); j++) {
                        kx = kxx[j]; ky = kyy[i];
                        kz = parameters.k - ((kx * kx + ky * ky) / 2 / parameters.k);
                        beamTilde[i * kxx.size() + j] = beamAtFocusTilde[i * kxx.size() + j] *
                                                    exp(std::complex<double>(0.0, 1.0) *
                                                            kz * z);
                    }
                }

                return beamTilde;
            }

            thrust::host_vector<double> propagationProfile(double zStart, double zEnd, double zSampling) {
                // declare z sampled location
                auto zArray = arange(zStart, zEnd, zSampling);

                thrust::host_vector<std::complex<double>> beamTilde(kxx.size() * kyy.size());
                thrust::host_vector<std::complex<double>> beam(kxx.size() * kyy.size());

                thrust::host_vector<double> beamZProfile(kxx.size() * zArray.size());

                double kx, ky, kz;

                // declare fftw_plan
                fftw_plan plan = fftw_plan_dft_2d(nRows, nCols, reinterpret_cast<fftw_complex*>(beamTilde.data()),
                                                  reinterpret_cast<fftw_complex*>(beam.data()), FFTW_BACKWARD, FFTW_MEASURE);

                /*  for each z
                 *      compute the beam profile for present z in fourier domain
                 *      perform ffthsift2d
                 *      compute the spatial beam profile from the fourier profile
                 *      perform ffthift2d
                 *      take the absolute value of the x-line at the center along the y-dimension and copy into beamZProfile
                 */
                for (int zIndex = 0; zIndex < zArray.size(); zIndex++) {
                    // compute the beam profile at z
                    for (int i = 0; i < kyy.size(); i++) {
                        for (int j = 0; j < kxx.size(); j++) {
                            kx = kxx[j]; ky = kyy[i];
                            kz = parameters.k - ((kx * kx + ky * ky) / 2 / parameters.k);
                            beamTilde[i * kxx.size() + j] = beamAtFocusTilde[i * kxx.size() + j] *
                                                            exp(std::complex<double>(0.0, 1.0) *
                                                                kz * zArray[zIndex]);
                        }
                    }

                    // perfrom 2d fftshift
                    fftshift2d<std::complex<double>>(beamTilde, nCols, nRows);

                    // perform fourier transform
                    fftw_execute(plan);

                    // perfrom 2d fftshift
                    fftshift2d<std::complex<double>>(beam, nCols, nRows);
                    for (int i = 0; i < kxx.size(); i++) {
                        beamZProfile[zIndex * kxx.size() + i] = abs(beam[kyy.size() / 2 * kxx.size() + i]);
                    }
                }

                fftw_destroy_plan(plan);

                return beamZProfile;
            }

        };

    }//cpu

}

#endif //PROJECT_PROPAGATOR_CUH
