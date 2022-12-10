//
// Created by ragib1481 on 12/8/22.
//

#ifndef PROJECT_SCATTER_CUH
#define PROJECT_SCATTER_CUH

#include <thrust/complex.h>

#include "Mesh.cuh"
#include "Vec3.cuh"
#include "Vec2.cuh"
#include "helper.cuh"
#include "myMath.cuh"
#include "Sphere.cuh"


namespace scatter {

    template <typename T>
    __host__ __device__
    thrust::complex<T> Es(const T r, const T theta, const T lambda,
                          const T k, const thrust::complex<T> n, const T a, const int Nl) {

        thrust::complex<T> val(0,0);
        for (int l = 0; l <= Nl; l++) {
            val += Bl(l, k, n , a) *
                   redefined::spHankel1Complex<T>(l, thrust::complex<T>(k*r, 0)) *
                   thrust::complex<T>(redefined::legendre<T>(l, cos(theta)), 0);
        }

        return val;
    }

    template <typename T>
    __host__ __device__
    thrust::complex<T> Ei(const T r, const T theta, const T lambda,
                          const T k, const thrust::complex<T> n, const T a, const int Nl) {
        thrust::complex<T> val(0,0);
        for (int l = 0; l <= Nl; l++) {
            val += Al(l, k, n , a) *
                    redefined::spBesselJComplex<T>(l, thrust::complex<T>(k*r, 0) * n) *
                    thrust::complex<T>(redefined::legendre<T>(l, cos(theta)), 0);
        }
    }


    namespace cpu {
        class MieScatter {

            // define constructor with the specified parameters

            /* simulation algorithm
             *      from the parameters create a cartesian mesh grid.
             *      for each point in the mesh accumulate the field for each scatterer.
             *      save the absolute value of the field in the mesh as the result
             * */

        };
    }
}


#endif //PROJECT_SCATTER_CUH
