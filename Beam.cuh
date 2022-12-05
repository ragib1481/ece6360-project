//
// Created by ragib1481 on 12/3/22.
//

#ifndef PROJECT_BEAM_CUH
#define PROJECT_BEAM_CUH

#include <fstream>
#include <string>
#include <thrust/host_vector.h>
#include <complex.h>
#include <opencv2/opencv.hpp>

#include "Propagator.cuh"
#include "helper.h"

namespace beams {

    namespace cpu {
        template <typename T>
        void plotBeamFromFourier(const thrust::host_vector<T>& beamTilde, const int nRows, const int nCols, std::string plotName) {
            thrust::host_vector<std::complex<double>> beamTildeTemp;
            thrust::host_vector<std::complex<double>> beam(nRows * nCols);
            thrust::host_vector<double> beamAbs(nRows * nCols);

            beamTildeTemp = beamTilde;

            fftshift2d<std::complex<double>>(beamTildeTemp, nCols, nRows);

            fftw_plan plan = fftw_plan_dft_2d(nRows, nCols, reinterpret_cast<fftw_complex*>(beamTildeTemp.data()),
                                              reinterpret_cast<fftw_complex*>(beam.data()), FFTW_BACKWARD, FFTW_ESTIMATE);

            fftw_execute(plan);
            fftw_destroy_plan(plan);

            // return magnitude of the beam.
            for (int i = 0; i < beam.size(); i++) {
                beamAbs[i] = abs(beam[i]);
            }
            fftshift2d<double>(beamAbs, nCols, nRows);

            // convert to opencv mat for easy plotting
            cv::Mat img(nRows, nCols, CV_64FC1, &beamAbs[0]);
            cv::Mat normImage;
            cv::normalize(img, normImage, 0, 255, cv::NORM_MINMAX);

            normImage.convertTo(normImage, CV_8UC1);

            cv::Mat imgColor;
            cv::applyColorMap(normImage, imgColor, cv::COLORMAP_JET);
            cv::imwrite(plotName + ".png", imgColor);
        }

        void plotBeamFromSpatial(const thrust::host_vector<double> beam, const int nRows, const int nCols, std::string plotName) {
            thrust::host_vector<double> b = beam;

            // convert to opencv mat for easy plotting
            cv::Mat img(nRows, nCols, CV_64FC1, &b[0]);
            cv::Mat normImage;
            cv::normalize(img, normImage, 0, 255, cv::NORM_MINMAX);

            normImage.convertTo(normImage, CV_8UC1);

            cv::Mat imgColor;
            cv::applyColorMap(normImage, imgColor, cv::COLORMAP_JET);
            cv::imwrite(plotName + ".png", imgColor);
        }

        class GaussianBeam {
            Parameters parameters;
            thrust::host_vector<double> beamAtFocusTilde;

        public:
            GaussianBeam(Parameters params) {
                parameters = params;
                double kx, ky;
                double alpha = M_PI / parameters.na;

                // generate spatial frequencies
                auto kxx = fftfreq(parameters.numSamples.x, parameters.sampling.x, 2 * M_PI);
                auto kyy = fftfreq(parameters.numSamples.y, parameters.sampling.y, 2 * M_PI);

                auto nCols = kxx.size();
                auto nRows = kyy.size();

                // generate beam at the focus
                beamAtFocusTilde.resize(kxx.size() * kyy.size());

                for(int j = 0; j < kyy.size(); j++) {
                    for (int i = 0; i < kxx.size(); i++) {
                        kx = kxx[i]; ky = kyy[j];
                        beamAtFocusTilde[j * kxx.size() + i] = exp(-(kx * kx + ky * ky) * alpha * alpha
                                / (2 * parameters.k * parameters.k));
                    }
                }
            }

            thrust::host_vector<double> getbeamAtFocusTilde() { return beamAtFocusTilde;}
        };

        class GaussBesselBeam {
            Parameters parameters;
            thrust::host_vector<double> beamAtFocusTilde;

        public:
            GaussBesselBeam(Parameters params) {
                parameters = params;
                double kx, ky;

                double alpha = M_PI / parameters.na;

                const double gbratio = 1.0;
                const double ringWidth = 10.0;
                double kr = gbratio * 2 * parameters.k/alpha;
                double mask = 0.0;
                double ksq = 0.0;

                // generate spatial frequencies
                auto kxx = fftfreq(parameters.numSamples.x, parameters.sampling.x, 2 * M_PI);
                auto kyy = fftfreq(parameters.numSamples.y, parameters.sampling.y, 2 * M_PI);

                auto nCols = kxx.size();
                auto nRows = kyy.size();

                // generate beam at the focus
                beamAtFocusTilde.resize(kxx.size() * kyy.size());

                for(int j = 0; j < kyy.size(); j++) {
                    for (int i = 0; i < kxx.size(); i++) {
                        kx = kxx[i]; ky = kyy[j];
                        ksq = kx * kx + ky * ky;

                        if ((ksq > kr * kr) && (ksq < ((kr + ringWidth) * (kr + ringWidth)))) {
                            mask = 1.0;
                        }
                        else
                            mask = 0.0;

                        beamAtFocusTilde[j * kxx.size() + i] = exp(-(kx * kx + ky * ky) * alpha * alpha
                                                                   / (2 * parameters.k * parameters.k)) * mask;
                    }
                }
            }

            thrust::host_vector<double> getbeamAtFocusTilde() { return beamAtFocusTilde;}
        };
    }
}


#endif //PROJECT_BEAM_CUH
