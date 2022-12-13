//
// Created by ragib1481 on 12/9/22.
//

#ifndef PROJECT_SPHERE_CUH
#define PROJECT_SPHERE_CUH

#include <thrust/complex.h>

#include "Vec2.cuh"
#include "Vec3.cuh"
#include "myMath.cuh"
#include "Bessel.cuh"

// ************************ functions to calculate scattering co-efficients ****************************

// TODO: revise this code
template <typename T>
__host__ __device__
void AlBl(thrust::complex<T>* al, thrust::complex<T>* bl,
          int Nl, T k, thrust::complex<T> n, T a) {
    T vm;
    thrust::complex<T>* jl_kna  = new thrust::complex<T>[Nl + 2];
    thrust::complex<T>* jlp_kna = new thrust::complex<T>[Nl + 2];
    thrust::complex<T>* yl_kna  = new thrust::complex<T>[Nl + 2];
    thrust::complex<T>* ylp_kna = new thrust::complex<T>[Nl + 2];

    stim::cbessjyva_sph(Nl, thrust::complex<T>(k * a, 0) * n, vm,
                        jl_kna, yl_kna, jlp_kna, ylp_kna);

    thrust::complex<T>* jl_ka  = new thrust::complex<T>[Nl + 2];
    thrust::complex<T>* jlp_ka = new thrust::complex<T>[Nl + 2];
    thrust::complex<T>* yl_ka  = new thrust::complex<T>[Nl + 2];
    thrust::complex<T>* ylp_ka = new thrust::complex<T>[Nl + 2];

    stim::cbessjyva_sph(Nl, thrust::complex<T>(k * a, 0), vm,
                        jl_ka, yl_ka, jlp_ka, ylp_ka);

    thrust::complex<T>* hl_ka  = new thrust::complex<T>[Nl + 1];
    stim::chankelva_sph(Nl, thrust::complex<T>(k * a, 0), hl_ka);

    thrust::complex<T>* hlp_ka = new thrust::complex<T>[Nl + 1];
    stim::chankelvap_sph(Nl, thrust::complex<T>(k * a, 0), hl_ka);

    for (int l = 0; l <= Nl; l++) {
        al[l] = (jl_ka[l] * hlp_ka[l] - jlp_ka[l] * hl_ka[l]) /
                (jl_kna[l] * hlp_ka[l] - hl_ka[l] * jlp_kna[l] * n);

        bl[l] = (jl_ka[l] * jlp_kna[l] * n - jlp_kna[l] * jl_ka[l]) /
                (jl_kna[l] * hlp_ka[l] - hl_ka[l] * jlp_kna[l] * n);

        al[l] = thrust::complex<T>(2 * l + 1, 0) * pow(thrust::complex<T>(0, 1), l)* al[l];
        bl[l] = thrust::complex<T>(2 * l + 1, 0) * pow(thrust::complex<T>(0, 1), l)* bl[l];
    }
    delete[] jl_kna; delete[] jlp_kna; delete[] yl_kna; delete[] ylp_kna;
    delete[] jl_ka;  delete[] jlp_ka;  delete[] yl_ka;  delete[] ylp_ka;
}

// ************************ class declaration for a scattering sphere **********************************
template <typename T>
__host__ __device__
class Sphere {
    Vec2<T> c;                                  // center of the sphere
    T r;                                        // radius of the sphere in um
    thrust::complex<T> n;                       // complex refractive index of the sphere
    thrust::complex<T>* bl;                     // scattering co-efficient Bl
    thrust::complex<T>* al;                     // scattering co-efficient Al
    unsigned int Nl = 30;                                     // max order to compute the scattering co-efficients for.

public:
    __host__ __device__
    Sphere(Vec2<T> c, T r, thrust::complex<T>n, T k):
                    c(c), r(r), n(n){
        // Nl = static_cast<int>(ceil(2.0 + k * r + 4.0 * cbrt(k * r)));
        // Nl = 500;
        al = new thrust::complex<T>[Nl+1];
        bl = new thrust::complex<T>[Nl+1];

        AlBl<T>(al, bl, Nl, k, n, r);
    };

    __host__ __device__
    Vec2<T> center() { return c; }

    __host__ __device__
    T radius() { return r; }

    __host__ __device__
    thrust::complex<T> getAl(unsigned int i) {return al[i];}

    __host__ __device__
    thrust::complex<T> getBl(unsigned int i) {return bl[i];}

    __host__ __device__
    unsigned int getMaxOrder() {return Nl;}

    __host__ __device__
    thrust::complex<T> refractiveIndex() {return n;}
};

#endif //PROJECT_SPHERE_CUH
