//
// Created by ragib1481 on 12/9/22.
//

#ifndef PROJECT_SPHERE_CUH
#define PROJECT_SPHERE_CUH

#include <thrust/complex.h>

#include "Vec2.cuh"
#include "Vec3.cuh"
#include "myMath.cuh"

// ************************ functions to calculate scattering co-efficients ****************************
// TODO: revise this code
template <typename T>
__host__ __device__
thrust::complex<T> Bl(int l, T k, thrust::complex<T> n, T a) {
    thrust::complex<T> num;
    thrust::complex<T> den;
    thrust::complex<T> jl_kna = redefined::spBesselJComplex<T>(l, thrust::complex<T>(k * a, 0) * n);
    thrust::complex<T> jl_ka  = redefined::spBesselJComplex<T>(l, thrust::complex<T>(k * a, 0));
    thrust::complex<T> jlp_ka  = redefined::spBesselJPComplex<T>(l, thrust::complex<T>(k * a, 0));
    thrust::complex<T> jlp_kna = redefined::spBesselJPComplex<T>(l, thrust::complex<T>(k * a, 0) * n);

    num = jl_ka * jlp_kna * n - jl_kna * jlp_ka;
    den = jl_kna * redefined::spHankel1PComplex<T>(l, thrust::complex<T>(k * a, 0)) -
          redefined::spHankel1Complex<T>(l, thrust::complex<T>(k * a, 0)) * jlp_kna * n;

    return thrust::complex<T>(2*l+1,0) * pow(thrust::complex<T>(0,1), l) * num / den;

}

// TODO: revise this code
template <typename T>
__host__ __device__
thrust::complex<T> Al(int l, T k, thrust::complex<T> n, T a) {
    thrust::complex<T> jl_kna = redefined::spBesselJComplex<T>(l, thrust::complex<T>(k * a, 0) * n);
    thrust::complex<T> jl_ka  = redefined::spBesselJComplex<T>(l, thrust::complex<T>(k * a, 0));
    thrust::complex<T> jlp_ka  = redefined::spBesselJPComplex<T>(l, thrust::complex<T>(k * a, 0));
    thrust::complex<T> jlp_kna = redefined::spBesselJPComplex<T>(l, thrust::complex<T>(k * a, 0) * n);
    thrust::complex<T> hl_ka  = redefined::spHankel1Complex<T>(l, thrust::complex<T>(k * a, 0));
    thrust::complex<T> hlp_ka  = redefined::spHankel1PComplex<T>(l, thrust::complex<T>(k * a, 0));

    thrust::complex<T> num = jl_ka * hlp_ka - jlp_ka * hl_ka;
    thrust::complex<T> den = jl_kna * hlp_ka - hl_ka * jlp_kna * n;
    return thrust::complex<T>(2 * l + 1, 0) * pow(thrust::complex<T>(0,1), l) * num / den;
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

        for (int l = 0; l <=Nl; l++) {
            al[l] = Al<T>(l, k, n, r);
            bl[l] = Bl<T>(l, k, n, r);
        }
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
