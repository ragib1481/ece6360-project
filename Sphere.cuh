//
// Created by ragib1481 on 12/9/22.
//

#ifndef PROJECT_SPHERE_CUH
#define PROJECT_SPHERE_CUH

#include <thrust/complex.h>

#include "Vec2.cuh"
#include "Vec3.cuh"
#include "myMath.cuh"

// ************************ class declaration for a point in the cartesian co-ordinate system **********
template <typename T>
__host__ __device__
class Point {
    /* This data structure holds one point in the cartesian system.
     * Each Point is a Vec2 data structure with some added functionality,
     * e.g. distance and angle measurement.
     *
     * */
    Vec2<T> p;
public:
    __host__ __device__
    Point(const Vec2<T>& p1){
        p = p1;
    }

    __host__ __device__
    Point& operator=(const Vec2<T>& p1){
        p = p1;
        return *this;
    }

    __host__ __device__
    T distance(const Vec2<T>& p1) {
        return (p-p1).mag();
    }

    __host__ __device__
    T angle(const Vec2<T>& p1) {
        Vec2<T> p2 = p - p1;
        return p2.getElement(0) / p2.mag();
    }
};


// ************************ functions to calculate scattering co-efficients ****************************
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
          redefined::spHankel1PComplex<T>(l, thrust::complex<T>(k * a, 0)) * jlp_kna * n;

    return thrust::complex<T>(2*l+1,0) * pow(thrust::complex<T>(0,1), l) * num / den;

}

template <typename T>
__host__ __device__
thrust::complex<T> Al(int l, T k, thrust::complex<T> n, T a) {
    thrust::complex<T> jl_kna = redefined::spBesselJComplex<T>(l, thrust::complex<T>(k * a, 0) * n);
    thrust::complex<T> jl_ka  = redefined::spBesselJComplex<T>(l, thrust::complex<T>(k * a, 0));
    thrust::complex<T> jlp_ka  = redefined::spBesselJPComplex<T>(l, thrust::complex<T>(k * a, 0));
    thrust::complex<T> jlp_kna = redefined::spBesselJPComplex<T>(l, thrust::complex<T>(k * a, 0) * n);
    thrust::complex<T> hl_ka  = redefined::spHankel1Complex<T>(l, thrust::complex<T>(k * a, 0));
    thrust::complex<T> hlp_ka  = redefined::spHankel1PComplex<T>(l, thrust::complex<T>(k * a, 0));
    return hlp_ka;
}

// ************************ class declaration for a scattering sphere **********************************
template <typename T>
__host__ __device__
class Sphere {
    Point<T> c;                                 // center of the sphere
    T r;                                        // radius of the sphere in um
    thrust::complex<T> n;                       // complex refractive index of the sphere
    T* bl;                                      // scattering co-efficient Bl
    T* al;                                      // scattering co-efficient Al
    int Nl;                                     // max order to compute the scattering co-efficients for.

public:
    __host__ __device__
    Sphere(Point<T> c, T r, thrust::complex<T>n, T k):
                    c(c), r(r), n(n){
        Nl = static_cast<int>(ceil(2.0 + k * r + 4.0 * cbrt(k * r)));
        al = new T[Nl+1];
        bl = new T[Nl+1];

        for (int l = 0; l <=Nl; l++) {
            al[l] = Al<T>(l, k, n, r);
            bl[l] = Bl<T>(l, k, n, r);
        }
    };

    __host__ __device__
    Point<T> center() { return c; }

    __host__ __device__
    T radius() { return r; }

     __host__ __device__
     T getAl(unsigned int i) {return al[i];}

    __host__ __device__
    T getBl(unsigned int i) {return bl[i];}

    __host__ __device__
    int getMaxOrder() {return Nl;}

    __host__ __device__
    thrust::complex<T> refractiveIndex() {return n;}
};

#endif //PROJECT_SPHERE_CUH
