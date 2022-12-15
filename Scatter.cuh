//
// Created by ragib1481 on 12/8/22.
//

#ifndef PROJECT_SCATTER_CUH
#define PROJECT_SCATTER_CUH

#include <opencv2/opencv.hpp>
#include <thrust/complex.h>
#include <string>

#include <math.h>
#include <limits>


#include "Mesh.cuh"
#include "Vec3.cuh"
#include "Vec2.cuh"
#include "helper.cuh"
#include "myMath.cuh"
#include "Sphere.cuh"
#include "Bessel.cuh"

#ifndef Nl
#define Nl 30
#endif

namespace scatter {

    template <typename T>
    __host__ __device__
    thrust::complex<T> Es(const thrust::complex<T>* bl, const T k, const T r, const T cosTheta) {
        thrust::complex<T> es(0, 0);

        // calculate hl_kr upto order Nl
        thrust::complex<T> hl_ka[Nl+1];
        redefined::chankelva_sph<T>(thrust::complex<T>(k*r, 0), hl_ka);

        // for each order l compute es_l and add to es;
        for (int l = 0; l <= Nl; l++) {
            es += bl[l] * hl_ka[l] * thrust::complex<T>(redefined::legendre<T>(l, cosTheta), 0);
        }

        // delete temporary variables
        return es;
    }

    template <typename T>
    __host__ __device__
    thrust::complex<T> Ef(const T k, const T r, const T cosTheta,
                          const T cosAlpha1, const T cosAlpha2) {
        thrust::complex<T> ef(0, 0);
        T cl;

        // calcualte jl_kr upto order Nl
        T vm;
        thrust::complex<T> jl_kr[Nl+2];
        thrust::complex<T> yl_kr[Nl+2];
        thrust::complex<T> jlp_kr[Nl+2];
        thrust::complex<T> ylp_kr[Nl+2];
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
        return thrust::complex<T>(2 * M_PI, 0) * ef;
    }

    template <typename T>
    __host__ __device__
    thrust::complex<T> Ei(const thrust::complex<T>* al, const T k, const T r, const T cosTheta,
                          const thrust::complex<T> n) {
        thrust::complex<T> ei(0, 0);

        // calculate hl_kr upto order Nl
        T vm;
        thrust::complex<T> jl_knr[Nl+2];
        thrust::complex<T> yl_knr[Nl+2];
        thrust::complex<T> jlp_knr[Nl+2];
        thrust::complex<T> ylp_knr[Nl+2];
        stimLab::cbessjyva_sph<T>(Nl, thrust::complex<T>(k * r, 0) * n, vm,
                               jl_knr, yl_knr, jlp_knr, ylp_knr);

        // for each order l compute es_l and add to es;
        for (int l = 0; l <= Nl; l++) {
            ei += al[l] * jl_knr[l] * thrust::complex<T>(redefined::legendre<T>(l, cosTheta), 0);
        }

        // delete temporary variables
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
            mesh::CartesianMesh<T> mesh;               // cartesian meshgrid to evaluate the simulation on
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
                mesh = mesh::CartesianMesh<T>(parameters.xStart, parameters.zStart,
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

                            c = sphere.center() - focus;
                            shiftedPhase = exp(thrust::complex<T>(0, k * kUnitVec.dot(c)));
                            r = mesh.getPoint(j, i) - sphere.center();
                            distance = r.mag();

                            cosTheta = kUnitVec.dot(r) / (r.mag() * kUnitVec.mag());

                            if (distance < sphere.radius()) {
                                totalField += shiftedPhase * Ei<T>(sphere.getAl(), k, distance, cosTheta,
                                                                    sphere.refractiveIndex());
                            }
                            else {
                                totalField += shiftedPhase * Es<T>(sphere.getBl(), k, distance, cosTheta);
                                totalField += Ef<T>(k, distance, cosTheta, cosAlpha1, cosAlpha2);
                            }
                        }
                        field[j * mesh.getWidth() + i] = totalField;
                    }
                }
            }

            void saveResult(std::string suffix) {
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
                cv::imwrite("./magnitude"+suffix+".png", magImgColor);

                // save real part of the result
                cv::Mat realImg(mesh.getHeight(), mesh.getWidth(), CV_32FC1, &real[0]);
                cv::Mat normRealImage;
                cv::normalize(realImg, normRealImage, 0, 255, cv::NORM_MINMAX);
                normRealImage.convertTo(normRealImage, CV_8UC1);
                cv::Mat realImgColor;
                cv::applyColorMap(normRealImage, realImgColor, cv::COLORMAP_JET);
                cv::imwrite("./real"+suffix+".png", realImgColor);

                // save imag part of the result
                cv::Mat imagImg(mesh.getHeight(), mesh.getWidth(), CV_32FC1, &imag[0]);
                cv::Mat normImagImage;
                cv::normalize(imagImg, normImagImage, 0, 255, cv::NORM_MINMAX);
                normImagImage.convertTo(normImagImage, CV_8UC1);
                cv::Mat imagImgColor;
                cv::applyColorMap(normImagImage, imagImgColor, cv::COLORMAP_JET);
                cv::imwrite("./imag"+suffix+".png", imagImgColor);
            }
        };
    };

    namespace gpu {
        // scatter kernel. compute the scattering for a single pixel for all spheres
        template <typename T>
        __global__
        void scatterkernel(thrust::complex<T>* field, const Sphere<T>* spheres, const Vec2<T>* mesh,
                                  const unsigned int width, const unsigned int height, const T k, const Vec2<T> focus,
                                  const Vec2<T> kUnitVec, const T cosAlpha1, const T cosAlpha2, const int numSpheres) {
            /* for present point in the mesh
             * set: totalField  = complex(0, 0)
             * for each sphere:
             *     if the point is inside the radius of the sphere:
             *         calculate: scatteredField = Ei
             *     else:
             *         calculate: scatteredField = Ef + Es
             *     totalField += scatteredField
             * set: field at the present point = totalField
             * */
            size_t ix = blockDim.x * blockIdx.x + threadIdx.x;
            if ( ix < (width * height) ) {
                thrust::complex<T> totalField(0, 0);                // field accumulator
                thrust::complex<T> shiftedPhase(0, 0);              // phase shift for the sphere center
                T cosTheta;
                T distance;
                Vec2<T> c;
                Vec2<T> r;

                for (int s = 0; s < numSpheres; s++) {
                    c = spheres[s].center() - focus;
                    shiftedPhase = exp(thrust::complex<T>(0, k * kUnitVec.dot(c)));
                    r = mesh[ix] - spheres[s].center();
                    distance = r.mag();

                    cosTheta = kUnitVec.dot(r) / (r.mag() * kUnitVec.mag());

                    if (distance < spheres[s].radius()) {
                        totalField += shiftedPhase * Ei<T>(spheres[s].getAl(), k, distance, cosTheta,
                                                           spheres[s].refractiveIndex());
                    }
                    else {
                        totalField += shiftedPhase * Es<T>(spheres[s].getBl(), k, distance, cosTheta);
                        totalField += Ef<T>(k, distance, cosTheta, cosAlpha1, cosAlpha2);
                    }
                }
                field[ix] = totalField;
            }
        }

        template <typename T>
        class MieScatter {
            /* simulation algorithm
             *      from the parameters create a cartesian mesh grid.
             *      for each point in the mesh accumulate the field for each scatterer.
             *      save the absolute value of the field in the mesh as the result
             * */
            mesh::CartesianMesh<T> mesh;               // cartesian meshgrid to evaluate the simulation on
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
                mesh = mesh::CartesianMesh<T>(parameters.xStart, parameters.zStart,
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


            // perform scattering simulation in parallel on gpu
            void scatter(bool verbose=false) {
                /* copy the mesh grid to the gpu.
                 * copy the spheres to the gpu.
                 * declare the field data as the output from the kernel.
                 * launch a thread for each point in the mesh grid.
                 * launch the threads in 2d blocks and grids.
                 * copy the result back to the host.
                 * */

                cudaError_t error;

                // copy mesh to the gpu and get a pointer
                thrust::device_vector<Vec2<T>> meshVectorGpu = mesh.getGrid();
                Vec2<T> *meshGpuPtr = thrust::raw_pointer_cast(meshVectorGpu.data());

                if (verbose){
                    cudaDeviceSynchronize();
                    std::cout << "copying mesh to gpu" << std::endl;
                    error = cudaGetLastError();
                    if (error != cudaSuccess) {
                        std::cout << cudaGetErrorString(error) << std::endl;
                    } else {
                        std::cout << "Success" << std::endl;
                    }
                }

                // copy mesh to the gpu and get a pointer
                thrust::device_vector<Sphere<T>> spheresGpu = spheres;
                Sphere<T>* spheresGpuPtr = thrust::raw_pointer_cast(spheresGpu.data());

                if (verbose) {
                    cudaDeviceSynchronize();
                    std::cout << "copying spheres to gpu" << std::endl;
                    error = cudaGetLastError();
                    if (error != cudaSuccess) {
                        std::cout << cudaGetErrorString(error) << std::endl;
                    } else {
                        std::cout << "Success" << std::endl;
                    }
                }

                // declare device_vector to contain the calculated field
                thrust::device_vector<thrust::complex<T>> fieldGpu(field.size());
                thrust::complex<T>* fieldGpuPtr = thrust::raw_pointer_cast(fieldGpu.data());

                unsigned int width = mesh.getWidth();
                unsigned int height= mesh.getHeight();

                // define blocks and grid
                uint threads = 256;
                uint blocks = width * height / threads + 1;
                scatterkernel<T><<<blocks, threads>>>(fieldGpuPtr, spheresGpuPtr, meshGpuPtr, width, height,
                                                      k, focus, kUnitVec, cosAlpha1, cosAlpha2,
                                                      spheresGpu.size());
                // scatterkernel(thrust::complex<T>* field, Sphere<T>* spheres, const Vec2<T>* mesh,
                //               const unsigned int width, const unsigned int height, const T k, const Vec2<T>& focus,
                //               Vec2<T>& kUnitVec, const T cosAlpha1, const T cosAlpha2, const int numSpheres) {

                if (verbose) {
                    cudaDeviceSynchronize();
                    std::cout << "performing computation" << std::endl;
                    error = cudaGetLastError();
                    if (error != cudaSuccess) {
                        std::cout << cudaGetErrorString(error) << std::endl;
                    } else {
                        std::cout << "Success" << std::endl;
                    }
                }

                // copy data back to the gpu
                field = fieldGpu;

                if (verbose) {
                    std::cout << "copying data back to host" << std::endl;
                    error = cudaGetLastError();
                    if (error != cudaSuccess) {
                        std::cout << cudaGetErrorString(error) << std::endl;
                    } else {
                        std::cout << "Success" << std::endl;
                    }
                }
            }

            void saveResult(std::string suffix) {
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
                cv::imwrite("./magnitude"+suffix+".png", magImgColor);

                // save real part of the result
                cv::Mat realImg(mesh.getHeight(), mesh.getWidth(), CV_32FC1, &real[0]);
                cv::Mat normRealImage;
                cv::normalize(realImg, normRealImage, 0, 255, cv::NORM_MINMAX);
                normRealImage.convertTo(normRealImage, CV_8UC1);
                cv::Mat realImgColor;
                cv::applyColorMap(normRealImage, realImgColor, cv::COLORMAP_JET);
                cv::imwrite("./real"+suffix+".png", realImgColor);

                // save imag part of the result
                cv::Mat imagImg(mesh.getHeight(), mesh.getWidth(), CV_32FC1, &imag[0]);
                cv::Mat normImagImage;
                cv::normalize(imagImg, normImagImage, 0, 255, cv::NORM_MINMAX);
                normImagImage.convertTo(normImagImage, CV_8UC1);
                cv::Mat imagImgColor;
                cv::applyColorMap(normImagImage, imagImgColor, cv::COLORMAP_JET);
                cv::imwrite("./imag"+suffix+".png", imagImgColor);
            }
        };

    };
}


#endif //PROJECT_SCATTER_CUH
