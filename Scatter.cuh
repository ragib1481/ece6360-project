//
// Created by ragib1481 on 12/8/22.
//

#ifndef PROJECT_SCATTER_CUH
#define PROJECT_SCATTER_CUH

#include <opencv2/opencv_modules.hpp>
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
    thrust::complex<T> Es_l(const thrust::complex<T> bl, const T k, const T r, const T cosTheta, const int l) {
        thrust::complex<T> val;
        val = bl;
        // val *= redefined::spHankel1Complex<T>(l, thrust::complex<T>(k * r, 0));
        val *= thrust::complex<T>(redefined::legendre<T>(l, cosTheta), 0);
        return val;
    }

    template <typename T>
    __host__ __device__
    thrust::complex<T> Ei_l(const thrust::complex<T> al, const thrust::complex<T> n,
                            const T k, const T r, const T cosTheta, const int l) {
        thrust::complex<T> val;
        val = al;
        // val *= redefined::spBesselJComplex<T>(l, thrust::complex<T>(k * r, 0) * n);
        val *= thrust::complex<T>(redefined::legendre<T>(l, cosTheta), 0);
        return val;
    }

    template <typename T>
    __host__ __device__
    thrust::complex<T> Ef_l(const T k, const T r, const T cosTheta,
                            const int l, const T cosAlpha1, const T cosAlpha2) {
        thrust::complex<T> val;
        val = pow(thrust::complex<T>(0, 1), l);
        // val *= redefined::spBesselJComplex<T>(l, thrust::complex<T>(k * r, 0));
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
                            r = mesh.getPoint(j, i);
                            distance = r.mag();

                            // if (distance < 0.001) {
                            //     cosTheta = 0;
                            // }
                            // else {
                            //     cosTheta = cos(r.getElement(0) / distance);
                            // }
                            cosTheta = r.getElement(0) / distance;

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
                                totalField += thrust::complex<T>(2 * M_PI, 0) * ef_l_temp + shiftedPhase * es_l_temp;
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

                T maxMag  = std::numeric_limits<T>::min();
                T maxReal = std::numeric_limits<T>::min();
                T maxImag = std::numeric_limits<T>::min();

                T minMag  = std::numeric_limits<T>::max();
                T minReal = std::numeric_limits<T>::max();
                T minImag = std::numeric_limits<T>::max();

                for (int i = 0; i < field.size(); i++) {
                    mag[i] = abs(field[i]);
                    if (maxMag > mag[i]) maxMag = mag[i];
                    if (minMag < mag[i]) minMag = mag[i];

                    real[i] = field[i].real();
                    if (maxReal > real[i]) maxReal = real[i];
                    if (minReal < real[i]) minReal = real[i];

                    imag[i] = field[i].imag();
                    if (maxImag > imag[i]) maxImag = imag[i];
                    if (minImag < imag[i]) minImag = imag[i];
                }

                thrust::host_vector<char> bytesMag(3 * field.size());
                thrust::host_vector<char> bytesReal(3 * field.size());
                thrust::host_vector<char> bytesImag(3 * field.size());
                for (int i = 0; i < field.size(); i++) {
                    bytesMag[3 * i] = (char) (255.0 * (mag[i] - minMag) / (maxMag - minMag));
                    bytesMag[3 * i + 1] = (char) (255.0 * (mag[i] - minMag) / (maxMag - minMag));
                    bytesMag[3 * i + 2] = (char) (255.0 * (mag[i] - minMag) / (maxMag - minMag));

                    bytesReal[3 * i] = (char) (255.0 * (real[i] - minReal) / (maxReal - minReal));
                    bytesReal[3 * i + 1] = (char) (255.0 * (real[i] - minReal) / (maxReal - minReal));
                    bytesReal[3 * i + 2] = (char) (255.0 * (real[i] - minReal) / (maxReal - minReal));

                    bytesImag[3 * i] = (char) (255.0 * (imag[i] - minImag) / (maxImag- minImag));
                    bytesImag[3 * i + 1] = (char) (255.0 * (imag[i] - minImag) / (maxImag- minImag));
                    bytesImag[3 * i + 2] = (char) (255.0 * (imag[i] - minImag) / (maxImag- minImag));
                }

                saveImage(bytesMag, "./magnitude.tga", mesh.getWidth(), mesh.getHeight());
                saveImage(bytesReal, "./real.tga", mesh.getWidth(), mesh.getHeight());
                saveImage(bytesImag, "./imag.tga", mesh.getWidth(), mesh.getHeight());
            }

            void saveImage(const thrust::host_vector<char>& bytes, std::string fileName,
                           const short& width, const short& height) {
                /* save the vector of char as an uncompressed tga file*/
                std::ofstream outfile;

                outfile.open(fileName, std::ios::binary | std::ios::out);	// open a binary file
                outfile.put(0);	// id length (field 1)
                outfile.put(0);	// color map type (field 2)
                outfile.put(2);	// image_type (field 3)
                outfile.put(0); outfile.put(0);	// color map field entry index (field 4)
                outfile.put(0); outfile.put(0);	// color map length (field 4)
                outfile.put(0);	// color map entry size (field 4)
                outfile.put(0); outfile.put(0);	// x origin (field 5)
                outfile.put(0); outfile.put(0);	// y origin (field 5)
                outfile.write((char*)&width, 2);	// image width (field 5)
                outfile.write((char*)&height, 2);	// image height (field 5)
                outfile.put(24);	// pixel depth (field 5)
                outfile.put(0);	// image descriptor (field 5)
                outfile.write(&bytes[0], width * height * 3);	// write the image data
                outfile.close();	// close the file
            }

        };
    }
}


#endif //PROJECT_SCATTER_CUH
