//
// Created by ragib1481 on 12/6/22.
//

#ifndef PROJECT_MYMATH_H
#define PROJECT_MYMATH_H

//#include <iostream>
#include <math.h>

namespace redefined {
    template <class T>
    __host__ __device__ inline T spBessel0(const T& x) {
        if (x == static_cast<T>(0.0)) return static_cast<T>(1.0);
        return sin(x) / x;
    }

    template <class T>
    __host__ __device__ inline T spBessel1(const T& x) {
        if (x == static_cast<T>(0.0)) return static_cast<T>(0.0);
        return (sin(x) - x * cos(x))/ (x * x);
    }

    template <class T>
    __host__ __device__ T spBesselN(const unsigned int n, const T& x) {
        if (n == 0) return spBessel0(x);
        if (n == 1) return spBessel1(x);
        if (x == static_cast<T>(0.0)) return static_cast<T>(0.0);

        T zm_2 = spBessel0(x);
        T zm_1 = spBessel1(x);
        T zm;

        for (unsigned int m = 2; m <= n; m++) {
            zm = (static_cast<T>(2.0) * static_cast<T>(m) - static_cast<T>(1.0)) * zm_1 / x - zm_2;
            zm_2 = zm_1;
            zm_1 = zm;
        }

        return zm;
    }

    template <class T>
    __host__ __device__ T spBesselDer1(const unsigned int n, const T& x) {
    }

    __host__ __device__ float spHankel(const int n, const float x) {
        return x;
    }

    __host__ __device__ float spHankel1(const int n, const float x) {
        return x;
    }

    template <class T>
    __host__ __device__ T legendre(const unsigned int n, const T& x) {
        T pn;
        T pn_2 = x;
        T pn_1 = (static_cast<T>(3.0) * x * x - static_cast<T>(1.0)) * static_cast<T>(0.5);

        if (n == 0) return static_cast<T>(1.0);
        if (n == 1) return pn_2;
        if (n == 2) return pn_1;

        for (unsigned int l = 3; l <= n; l++) {
            pn = ((static_cast<T>(2.0) * static_cast<T>(l) - static_cast<T>(1.0)) * x * pn_1 - (static_cast<T>(l) - static_cast<T>(1.0)) * pn_2) / static_cast<T>(l);
            pn_2 = pn_1;
            pn_1 = pn;
        }

        return pn;
    }
}


#endif //PROJECT_MYMATH_H
