//
// Created by ragib1481 on 12/8/22.
//

#ifndef PROJECT_SCATTER_CUH
#define PROJECT_SCATTER_CUH

#include <thrust/complex.h>

#include <math.h>
#include <opencv2/opencv_modules.hpp>
#include <opencv2/opencv.hpp>

#include "Mesh.cuh"
#include "Vec3.cuh"
#include "Vec2.cuh"
#include "helper.cuh"
#include "myMath.cuh"
#include "Sphere.cuh"


namespace scatter {

    template <typename T>
    __host__ __device__
    thrust::complex<T> Es_l(const thrust::complex<T> bl, const T k, const T r, const T cosTheta, const int l) {
        thrust::complex<T> val;
        val = bl * redefined::spHankel1Complex<T>(l, thrust::complex<T>(k * r, 0));
        val *= thrust::complex<T>(redefined::legendre<T>(l, cosTheta), 0);
        return val;
    }

    template <typename T>
    __host__ __device__
    thrust::complex<T> Ei_l(const thrust::complex<T> al, const thrust::complex<T> n,
                            const T k, const T r, const T cosTheta, const int l) {
        thrust::complex<T> val;
        val = al;
        val *= redefined::spBesselJComplex<T>(l, thrust::complex<T>(k * r, 0) * n);
        val *= thrust::complex<T>(redefined::legendre<T>(l, cosTheta), 0);
        return val;
    }

    template <typename T>
    __host__ __device__
    thrust::complex<T> Ef_l(const T k, const T r, const T cosTheta,
                            const int l, const T cosAlpha1, const T cosAlpha2) {
        thrust::complex<T> val;
        val = pow(thrust::complex<T>(0, 1), l);
        val *= redefined::spBesselJComplex<T>(l, thrust::complex<T>(k * r, 0));
        val *= thrust::complex<T>(redefined::legendre<T>(l, cosTheta), 0);
        val *= (redefined::legendre<T>(l+1, cosAlpha1) - redefined::legendre<T>(l+1, cosAlpha2) -
                redefined::legendre<T>(l-1, cosAlpha1) + redefined::legendre<T>(l-1, cosAlpha2));
        return val;
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
                            r = mesh.getPoint(j, i) - sphere.center();
                            distance = r.mag();

                            // if (distance < 0.001) {
                            //     cosTheta = 0;
                            // }
                            // else {
                            //     cosTheta = cos(r.getElement(0) / distance);
                            // }
                            cosTheta = r.getElement(1) / distance;

                            if (distance < sphere.radius()) {
                                for (unsigned int l = 0; l <= sphere.getMaxOrder(); l++) {
                                    ei_l_temp += Ei_l<T>(sphere.getAl(l), sphere.refractiveIndex(),
                                                         k, distance, cosTheta, l);
                                }
                                totalField +=  shiftedPhase * ei_l_temp;
                            }
                            else {
                                for (unsigned int l = 0; l <= sphere.getMaxOrder(); l++) {
                                    es_l_temp += Es_l<T>(sphere.getBl(l), k, distance, cosTheta, l);
                                    ef_l_temp += Ef_l<T>(k, distance, cosTheta, l, cosAlpha1, cosAlpha2);
                                }
                                // totalField += thrust::complex<T>(2 * M_PI, 0) * ef_l_temp + shiftedPhase * es_l_temp;
                                // totalField += shiftedPhase * es_l_temp;
                                // totalField += thrust::complex<T>(2 * M_PI, 0) * ef_l_temp;
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

                // save magnitude of the result
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
