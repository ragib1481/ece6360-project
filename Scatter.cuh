//
// Created by ragib1481 on 12/8/22.
//

#ifndef PROJECT_SCATTER_CUH
#define PROJECT_SCATTER_CUH

#include <thrust/complex.h>

#include <math.h>

#include "Mesh.cuh"
#include "Vec3.cuh"
#include "Vec2.cuh"
#include "helper.cuh"
#include "myMath.cuh"
#include "Sphere.cuh"


namespace scatter {

    template <typename T>
    //__host__ __device__
    thrust::complex<T> Es_l(const thrust::complex<T> bl, const T k, const T r, const T cosTheta, const int l) {
        thrust::complex<T> val;
        val = bl * redefined::spHankel1Complex(l, thrust::complex<T>(k * r, 0));
        val *= thrust::complex<T>(redefined::legendre<T>(l, cosTheta), 0);
        return val;
    }

    template <typename T>
    //__host__ __device__
    thrust::complex<T> Ei_l(const thrust::complex<T> al, const thrust::complex<T> n,
                            const T k, const T r, const T cosTheta, const int l) {
        thrust::complex<T> val;
        val = al;
        val *= redefined::spBesselJComplex<T>(l, thrust::complex<T>(k * r) * n);
        val *= thrust::complex<T>(redefined::legendre<T>(l, cosTheta), 0);
        return val;
    }

    template <typename T>
    //__host__ __device__
    thrust::complex<T> Ef_l(const T k, const T r, const T cosTheta,
                            const int l, const T cosAlpha1, const T cosAlpha2) {
        thrust::complex<T> val;
        val = pow(thrust::complex<T>(0, 1), l);
        val *= redefined::spBesselJComplex(l, thrust::complex<T>(k * r, 0));
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
                for (int j = 0; j < mesh.getHeight(); j++) {
                    for (int i = 0; i < mesh.getWidth(); i++) {
                        // reset temporary variables
                        totalField = thrust::complex<T> (0, 0);
                        ei_l_temp = thrust::complex<T> (0, 0);
                        es_l_temp = thrust::complex<T> (0, 0);
                        ef_l_temp = thrust::complex<T> (0, 0);

                        // iterate through all spheres to calculate the field at the present point
                        for (auto sphere: spheres) {
                            Vec2<T> r = mesh.getPoint(j, i) - sphere.center();
                            T distance = r.mag();
                            T cosTheta;
                            if (distance < 0.001) {
                                cosTheta = 1;
                            }
                            else {
                                cosTheta = cos(r.getElement(0) / distance); // placeholder
                            }
                            if (distance < sphere.radius()) {
                                for (unsigned int l = 0; l <= sphere.getMaxOrder(); l++) {
                                    ei_l_temp += Ei_l<T>(sphere.getAl(l), sphere.refractiveIndex(),
                                                         k, sphere.radius(), cosTheta, l);
                                }
                                totalField += ei_l_temp;
                            }
                            else {
                                for (unsigned int l = 0; l <= sphere.getMaxOrder(); l++) {
                                    es_l_temp += Es_l<T>(sphere.getBl(l), k, sphere.radius(), cosTheta, l);
                                    ef_l_temp += Ef_l<T>(k, sphere.radius(), cosTheta, l, cosAlpha1, cosAlpha2);
                                }
                                totalField += thrust::complex<T>(2 * M_PI, 0) * ef_l_temp + es_l_temp;
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
                    mag[i] = field[i];
                    real[i] = field[i].real();
                    imag[i] = field[i].imag();
                }
            }
        };
    }
}


#endif //PROJECT_SCATTER_CUH
