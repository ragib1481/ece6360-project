//
// Created by ragib1481 on 12/8/22.
//

#ifndef PROJECT_SCATTER_CUH
#define PROJECT_SCATTER_CUH

#include <thrust/complex.h>

#include "Mesh.cuh"
#include "Vec3.cuh"
#include "helper.h"
#include "myMath.h"


namespace scatter {

    template <typename T>
    __host__ __device__
    thrust::complex<T> Bl(T jl_ka, T jlp_ka, T jl_kna, T jlp_kna, T hl_ka, T hlp_ka,
                          T n, int l) {
        thrust::complex<T> val(0,0);
    }

    template <typename T>
    __host__ __device__
    thrust::complex<T> Al(int l, T lambda, T n, T a) {
        std::complex<T> val(0, 0);
        return val;
    }

    template <typename T>
    __host__ __device__
    thrust::complex<T> Es(const T r, const T theta, const T lambda,
                          const T k, const T n, const T a, const int Nl) {
        thrust::complex<T> val;
        thrust::complex<T> hl_ka, hlp_ka;
        T vm;

        // calculate jl_kna(x), jlp_kna(x), for kna
        T* jl_kna   = new T[Nl+1];
        T* jlp_kna  = new T[Nl+1];
        T* yl_kna   = new T[Nl+1];
        T* ylp_kna  = new T[Nl+1];
        // redefined::sphBessjyv(Nl, k * n * a, vm, jl_kna, yl_kna, jlp_kna, ylp_kna);

        // calculate jl(x), jl_prime(x), jyl(x), y_prime(x) for ka
        T* jl_ka   = new T[Nl+1];
        T* jlp_ka  = new T[Nl+1];
        T* yl_ka   = new T[Nl+1];
        T* ylp_ka  = new T[Nl+1];
        // redefined::sphBessjyv(Nl, k * n * a, vm, jl_ka, yl_ka, jlp_ka, ylp_ka);

        // calculate  hl(x) and hl_prime(x) for ka
        for (int l = 0; l <= Nl; l++) {
            hl_ka   = thrust::complex<T>(jl_ka[l], yl_ka[l]);
            hlp_ka  = thrust::complex<T>(jlp_ka[l], ylp_ka[l]);
        }

        return val;
    }

    template <typename T>
    __host__ __device__
    thrust::complex<T> Ei(const T r, const T theta, const T lambda,
                          const T k, const T n, const T a, const int Nl) {
        thrust::complex<T> val;
        for (int l = Nl; l >=  0; l--) {
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
