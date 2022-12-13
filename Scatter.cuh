//
// Created by ragib1481 on 12/8/22.
//

#ifndef PROJECT_SCATTER_CUH
#define PROJECT_SCATTER_CUH

#include <opencv2/opencv.hpp>
#include <thrust/complex.h>

#include <math.h>
#include <limits>


#include "Mesh.cuh"
#include "Vec3.cuh"
#include "Vec2.cuh"
#include "helper.cuh"
#include "myMath.cuh"
#include "Sphere.cuh"
#include "Bessel.cuh"


namespace scatter {

    template <typename T>
    __host__ __device__
    thrust::complex<T> Es(const thrust::complex<T>* bl, int Nl, const T k, const T r, const T cosTheta) {
        thrust::complex<T> es(0, 0);

        // calculate hl_kr upto order Nl
        thrust::complex<T>* hl_ka = new thrust::complex<T>[Nl+1];
        stimLab::chankelva_sph<T>(Nl, thrust::complex<T>(k*r, 0), hl_ka);

        // for each order l compute es_l and add to es;
        for (int l = 0; l <= Nl; l++) {
            es += bl[l] * hl_ka[l] * thrust::complex<T>(redefined::legendre<T>(l, cosTheta), 0);
        }

        // delete temporary variables
        delete[] hl_ka;
        return es;
    }

    template <typename T>
    __host__ __device__
    thrust::complex<T> Ef(const int Nl, const T k, const T r, const T cosTheta,
                          const T cosAlpha1, const T cosAlpha2) {
        thrust::complex<T> ef;
        T cl;

        // calcualte jl_kr upto order Nl
        T vm;
        thrust::complex<T>* jl_kr = new thrust::complex<T>[Nl+2];
        thrust::complex<T>* yl_kr = new thrust::complex<T>[Nl+2];
        thrust::complex<T>* jlp_kr = new thrust::complex<T>[Nl+2];
        thrust::complex<T>* ylp_kr = new thrust::complex<T>[Nl+2];
        stimLab::cbessjyva_sph<T>(Nl, thrust::complex<T>(k * r, 0), vm,
                               jl_kr, yl_kr, jlp_kr, ylp_kr);

        for (int l = 0; l <= Nl; l++) {
            cl = redefined::legendre<T>(l+1, cosAlpha1) - redefined::legendre<T>(l+1, cosAlpha2)
                    - redefined::legendre<T>(l-1, cosAlpha1) + redefined::legendre<T>(l-1, cosAlpha2);
            ef += pow(thrust::complex<T> (0, 1), l) * jl_kr[l] *
                    thrust::complex<T>(redefined::legendre<T>(l, cosTheta), 0) *
                            thrust::complex<T>(cl, 0);
        }

        // delete temporary variables
        delete[] jl_kr; delete[] yl_kr; delete[] jlp_kr; delete[] ylp_kr;
        return thrust::complex<T>(2 * M_PI, 0) * ef;
    }

    template <typename T>
    __host__ __device__
    thrust::complex<T> Ei(const thrust::complex<T>* al, const int Nl, const T k, const T r, const T cosTheta,
                          const thrust::complex<T> n) {
        thrust::complex<T> ei;

        // calculate hl_kr upto order Nl
        T vm;
        thrust::complex<T>* jl_knr = new thrust::complex<T>[Nl+2];
        thrust::complex<T>* yl_knr = new thrust::complex<T>[Nl+2];
        thrust::complex<T>* jlp_knr = new thrust::complex<T>[Nl+2];
        thrust::complex<T>* ylp_knr = new thrust::complex<T>[Nl+2];
        stimLab::cbessjyva_sph<T>(Nl, thrust::complex<T>(k * r, 0) * n, vm,
                               jl_knr, yl_knr, jlp_knr, ylp_knr);

        // for each order l compute es_l and add to es;
        for (int l = 0; l <= Nl; l++) {
            ei += al[l] * jl_knr[l] * thrust::complex<T>(redefined::legendre<T>(l, cosTheta), 0);
        }

        // delete temporary variables
        delete[] jl_knr; delete[] yl_knr; delete[] jlp_knr; delete[] ylp_knr;
        return ei;
    }

    namespace cpu {
        template <typename T>
        class MieScatter {
            /* simulation algorithm
             *      from the parameters create a cartesian mesh grid.
             *      for each point in the mesh accumulate the field for each scatterer.
             *      save the absolute value of the field in the mesh as the result
             * */
            mesh::cpu::CartesianMesh<T> mesh;               // cartesian meshgrid to evaluate the simulation on
            thrust::host_vector<Sphere<T>> spheres;         // list of spheres in the mesh
            T k;                                            // wavenumber of interest
            Vec2<T> kUnitVec;                               // direction of k vector
            Vec2<T> focus;                                  // illumination focus
            thrust::host_vector<thrust::complex<T>> field;
            T cosAlpha1;
            T cosAlpha2;

        public:
            // define constructor with the specified parameters
            MieScatter(Parameters<T> parameters) {

                // initialize carteisan mesh
                mesh = mesh::cpu::CartesianMesh<T>(parameters.xStart, parameters.zStart,
                                                   parameters.xEnd, parameters.zEnd,
                                                   parameters.sampling.getElement(0),
                                                   parameters.sampling.getElement(1));

                // resize the field vector
                field.resize(mesh.getHeight() * mesh.getWidth());
                for (int i = 0; i < field.size(); i++) {
                    field[i] = thrust::complex<T>(0, 0);
                }

                // initialize other simulation parameters;
                k           = parameters.k;
                spheres     = parameters.spheres;
                kUnitVec    = parameters.kUnitVector;
                focus       = parameters.focus;
                cosAlpha2   = cos(asin(parameters.na));
                cosAlpha1   = cos(0);
            }

            // perform scattering simulation
            void scatter() {
                /* for each point i, j in the mesh:
                 *      set: totalField  = complex(0, 0)
                 *      for each sphere:
                 *          if the point is inside the radius of the sphere:
                 *              calculate: scatteredField = Ei
                 *          else:
                 *              calculate: scatteredField = Ef + Es
                 *          totalField += scatteredField
                 *      set: field at point i , j = totalField
                 * */
                thrust::complex<T> totalField;
                thrust::complex<T> ei_l_temp;
                thrust::complex<T> es_l_temp;
                thrust::complex<T> ef_l_temp;
                Vec2<T> r;
                T distance;
                T cosTheta;
                Vec2<T> c;
                thrust::complex<T> shiftedPhase;
                for (int j = 0; j < mesh.getHeight(); j++) {
                    for (int i = 0; i < mesh.getWidth(); i++) {
                        // reset temporary variables
                        totalField = thrust::complex<T> (0, 0);

                        // iterate through all spheres to calculate the field at the present point
                        for (auto sphere: spheres) {
                            // reset temporary variables for each sphere
                            ei_l_temp = thrust::complex<T> (0, 0);
                            es_l_temp = thrust::complex<T> (0, 0);
                            ef_l_temp = thrust::complex<T> (0, 0);

                            c = sphere.center() - focus;
                            shiftedPhase = exp(thrust::complex<T>(0, k * kUnitVec.dot(c)));
                            r = mesh.getPoint(j, i);
                            distance = r.mag();

                            cosTheta = r.getElement(0) / distance;

                            if (distance < sphere.radius()) {
                                ei_l_temp = Ei<T>(sphere.getAl(), sphere.getMaxOrder(), k, distance, cosTheta,
                                                  sphere.refractiveIndex());

                                totalField +=  shiftedPhase * ei_l_temp;
                            }
                            else {
                                es_l_temp = Es<T>(sphere.getBl(), sphere.getMaxOrder(), k, distance, cosTheta);
                                ef_l_temp = Ef<T>(sphere.getMaxOrder(), k, distance, cosTheta, cosAlpha1, cosAlpha2);

                                totalField += ef_l_temp + shiftedPhase * es_l_temp;
                            }
                        }
                        field[j * mesh.getWidth() + i] = totalField;
                    }
                }
            }

            void saveResult() {
                thrust::host_vector<T> mag(field.size());
                thrust::host_vector<T> real(field.size());
                thrust::host_vector<T> imag(field.size());

                for (int i = 0; i < field.size(); i++) {
                    mag[i] = abs(field[i]);
                    real[i] = field[i].real();
                    imag[i] = field[i].imag();
                }

                // save magnitude of the image
                cv::Mat magImg(mesh.getHeight(), mesh.getWidth(), CV_32FC1, &mag[0]);
                cv::Mat normMagImage;
                cv::normalize(magImg, normMagImage, 0, 255, cv::NORM_MINMAX);
                normMagImage.convertTo(normMagImage, CV_8UC1);
                cv::Mat magImgColor;
                cv::applyColorMap(normMagImage, magImgColor, cv::COLORMAP_JET);
                cv::imwrite("./magnitude.png", magImgColor);

                // save real part of the result
                cv::Mat realImg(mesh.getHeight(), mesh.getWidth(), CV_32FC1, &real[0]);
                cv::Mat normRealImage;
                cv::normalize(realImg, normRealImage, 0, 255, cv::NORM_MINMAX);
                normRealImage.convertTo(normRealImage, CV_8UC1);
                cv::Mat realImgColor;
                cv::applyColorMap(normRealImage, realImgColor, cv::COLORMAP_JET);
                cv::imwrite("./real.png", realImgColor);

                // save imag part of the result
                cv::Mat imagImg(mesh.getHeight(), mesh.getWidth(), CV_32FC1, &imag[0]);
                cv::Mat normImagImage;
                cv::normalize(imagImg, normImagImage, 0, 255, cv::NORM_MINMAX);
                normImagImage.convertTo(normImagImage, CV_8UC1);
                cv::Mat imagImgColor;
                cv::applyColorMap(normImagImage, imagImgColor, cv::COLORMAP_JET);
                cv::imwrite("./imag.png", imagImgColor);
            }
        };
    }
}


#endif //PROJECT_SCATTER_CUH
