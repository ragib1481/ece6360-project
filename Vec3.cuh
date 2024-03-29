//
// Created by ragib1481 on 12/7/22.
//

#ifndef PROJECT_VEC3_CUH
#define PROJECT_VEC3_CUH

#include <iostream>
#include <ostream>
#include <math.h>

template <typename T>
__host__ __device__
class Vec3 {
    T vec[3];

public:
    __host__ __device__
    Vec3() {
        vec[0] = static_cast<T>(0);
        vec[1] = static_cast<T>(0);
        vec[2] = static_cast<T>(0);
    }

    __host__ __device__
    Vec3(T i, T j, T k) {
        vec[0] = i;
        vec[1] = j;
        vec[2] = k;
    }

    __host__ __device__
    Vec3& operator=(const Vec3& r) {
        vec[0] = r.vec[0];
        vec[1] = r.vec[1];
        vec[2] = r.vec[2];
        return *this;
    }

    __host__ __device__
    Vec3& operator=(const T& val) {
        vec[0] = val;
        vec[1] = val;
        vec[2] = val;
        return *this;
    }

    __host__ __device__
    Vec3 operator+(const Vec3& r) {
        Vec3 v(0, 0, 0);
        v.vec[0] = vec[0] + r.vec[0];
        v.vec[1] = vec[1] + r.vec[1];
        v.vec[2] = vec[2] + r.vec[2];
        return v;
    }

    __host__ __device__
    Vec3 operator+(const T& val) {
        Vec3 v(0, 0, 0);
        v.vec[0] = vec[0] + val;
        v.vec[1] = vec[1] + val;
        v.vec[2] = vec[2] + val;
        return v;
    }

    __host__ __device__
    Vec3 operator-(const Vec3& r) {
        Vec3 v(0, 0, 0);
        v.vec[0] = vec[0] - r.vec[0];
        v.vec[1] = vec[1] - r.vec[1];
        v.vec[2] = vec[2] - r.vec[2];
        return v;
    }

    __host__ __device__
    Vec3 operator-(const T& val) {
        Vec3 v(0, 0, 0);
        v.vec[0] = vec[0] - val;
        v.vec[1] = vec[1] - val;
        v.vec[2] = vec[2] - val;
        return v;
    }

    __host__ __device__
    Vec3 operator/(const Vec3& r) {
        Vec3 v(0, 0, 0);
        v.vec[0] = vec[0] / r.vec[0];
        v.vec[1] = vec[1] / r.vec[1];
        v.vec[2] = vec[2] / r.vec[2];
        return v;
    }

    __host__ __device__
    Vec3 operator/(const T& val) {
        Vec3 v(0, 0, 0);
        v.vec[0] = vec[0] / val;
        v.vec[1] = vec[1] / val;
        v.vec[2] = vec[2] / val;
        return v;
    }

    __host__ __device__
    Vec3 operator*(const Vec3& r) {
        Vec3 v(0, 0, 0);
        v.vec[0] = vec[0] * r.vec[0];
        v.vec[1] = vec[1] * r.vec[1];
        v.vec[2] = vec[2] * r.vec[2];
        return v;
    }

    __host__ __device__
    Vec3 operator*(const T& val) {
        Vec3 v(0, 0, 0);
        v.vec[0] = vec[0] * val;
        v.vec[1] = vec[1] * val;
        v.vec[2] = vec[2] * val;
        return v;
    }

    __host__ __device__
    T dot(const Vec3& r) {
        T x = vec[0] * r.vec[0] + vec[1] * r.vec[1] + vec[2] * r.vec[2];
        return x;
    }

    __host__ __device__
    T mag() {
        return sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);
    }

    __host__ __device__
    T getElement(unsigned short i) {
        return vec[i];
    }

    __host__ __device__
    void setElement(unsigned short i, T val) {
        vec[i] = val;
    }

    friend std::ostream& operator<<(std::ostream& output, const Vec3& v) {
        output << "(" << v.vec[0] << ", " << v.vec[1] << ", " << v.vec[2] << ")";
        return output;
    }
};

#endif //PROJECT_VEC3_CUH
