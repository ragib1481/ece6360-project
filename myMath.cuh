//
// Created by ragib1481 on 12/6/22.
//

#ifndef PROJECT_MYMATH_CUH
#define PROJECT_MYMATH_CUH

#include <math.h>
#include "helper.cuh"
#include "Bessel.cuh"

#ifndef Nl
#define Nl 30
#endif

#define el 0.5772156649015329

__host__ __device__
int factorial(int n) {
    int val = 1;
    for (int i = 2; i <= n; i++) {
        val *= i;
    }
    return val;
}

namespace redefined {

    template<typename P>
    __host__ __device__
    int chankelva_sph(const thrust::complex<P> z, thrust::complex<P>* chv) {
        P vm;
        int returnVal;
        thrust::complex<P> cjv  [Nl+2];
        thrust::complex<P> cyv  [Nl+2];
        thrust::complex<P> cjvp [Nl+2];
        thrust::complex<P> cyvp [Nl+2];

        returnVal = stimLab::cbessjyva_sph<P>(Nl, z, vm, cjv, cyv, cjvp, cyvp);

        for (int i = 0; i <= Nl; i++) {
            chv[i] = cjv[i] + thrust::complex<P> (0, 1) * cyv[i];
        }

        return returnVal;
    }

    template <typename T>
    __host__ __device__ T legendre(const int n, const T& x) {
        T pn;
        T pn_2 = x;
        T pn_1 = (static_cast<T>(3.0) * x * x - static_cast<T>(1.0)) * static_cast<T>(0.5);

        if (n == 0) return static_cast<T>(1.0);
        if (n == 1) return pn_2;
        if (n == 2) return pn_1;
        if (n == -1) return static_cast<T>(1.0);

        for (int l = 3; l <= n; l++) {
            pn = ((static_cast<T>(2.0) * static_cast<T>(l) - static_cast<T>(1.0)) * x * pn_1 - (static_cast<T>(l) - static_cast<T>(1.0)) * pn_2) / static_cast<T>(l);
            pn_2 = pn_1;
            pn_1 = pn;
        }

        return pn;
    }

    template <typename T>
    __host__ __device__
    thrust::complex<T> spBesselJComplex(int n, thrust::complex<T> x) {
        float temp;
        thrust::complex<float> val(0, 0);
        for (int s = 0; s < 20; s++) {
            temp = pow<T>(-1.0, (T)s) / ((T)factorial(s) * tgamma((T)(s + n) + 1.5));
            val += pow<T>(x/thrust::complex<T>(2.0f, 0.0f), 2.0*s+n+0.5) * thrust::complex<T>(temp, 0);
        }

        return val * sqrt(thrust::complex<T>(M_PI / 2.0, 0)/ x);
    }

    template <typename T>
    __host__ __device__
    thrust::complex<T> spBesselYComplex(int n, thrust::complex<T> x) {
        thrust::complex<T> mult(1, 0);
        if ((n % 2) == 0)  mult.real(-1.0);
        return mult * spBesselJComplex<T>(-n-1, x);
    }

    template <typename T>
    __host__ __device__
    thrust::complex<T> spHankel1Complex(int n, thrust::complex<T> x) {
        return spBesselJComplex<T>(n, x) + thrust::complex<T>(0.0, 1.0) * spBesselYComplex<T>(n, x);
    }

    template <typename T>
    __host__ __device__
    thrust::complex<T> spBesselJPComplex(int n, thrust::complex<T> x) {
        thrust::complex<T> val = thrust::complex<T>(n, 0) * spBesselJComplex<T>(n-1, x);
        val -= thrust::complex<T>(n+1, 0) * spBesselJComplex<T>(n+1, x);
        val /= thrust::complex<T>(2*n+1, 0);
        return val;
    }

    template <typename T>
    __host__ __device__
    thrust::complex<T> spBesselYPComplex(int n, thrust::complex<T> x) {
        thrust::complex<T> val = thrust::complex<T>(n, 0) * spBesselYComplex<T>(n-1, x);
        val -= thrust::complex<T>(n+1, 0) * spBesselYComplex<T>(n+1, x);
        val /= thrust::complex<T>(2*n+1, 0);
        return val;
    }

    template <typename T>
    __host__ __device__
    thrust::complex<T> spHankel1PComplex(int n, thrust::complex<T> x) {
        thrust::complex<T> val = thrust::complex<T>(n, 0) * spHankel1Complex<T>(n-1, x);
        val -= thrust::complex<T>(n+1, 0) * spHankel1Complex<T>(n+1, x);
        val /= thrust::complex<T>(2*n+1, 0);
        return val;
    }

}


#endif //PROJECT_MYMATH_CUH
