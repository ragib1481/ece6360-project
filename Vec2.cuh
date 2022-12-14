//
// Created by ragib1481 on 12/9/22.
//

#ifndef PROJECT_VEC2_CUH
#define PROJECT_VEC2_CUH

#include <iostream>
#include <ostream>
#include <math.h>

template <typename T>
__host__ __device__
class Vec2 {
    T vec[2];

public:
    __host__ __device__
    Vec2() {
        vec[0] = static_cast<T>(0);
        vec[1] = static_cast<T>(0);
    }

    __host__ __device__
    Vec2(T i, T j) {
        vec[0] = i;
        vec[1] = j;
    }

    __host__ __device__
    Vec2& operator=(const Vec2& r) {
        vec[0] = r.vec[0];
        vec[1] = r.vec[1];
        return *this;
    }

    __host__ __device__
    Vec2& operator=(const T& val) const{
        vec[0] = val;
        vec[1] = val;
        return *this;
    }

    __host__ __device__
    Vec2 operator+(const Vec2& r) const{
        Vec2 v(0, 0);
        v.vec[0] = vec[0] + r.vec[0];
        v.vec[1] = vec[1] + r.vec[1];
        return v;
    }

    __host__ __device__
    Vec2 operator+(const T& val) const{
        Vec2 v(0, 0);
        v.vec[0] = vec[0] + val;
        v.vec[1] = vec[1] + val;
        return v;
    }

    __host__ __device__
    Vec2 operator-(const Vec2& r) const{
        Vec2 v(0, 0);
        v.vec[0] = vec[0] - r.vec[0];
        v.vec[1] = vec[1] - r.vec[1];
        return v;
    }

    __host__ __device__
    Vec2 operator-(const T& val) const{
        Vec2 v(0, 0);
        v.vec[0] = vec[0] - val;
        v.vec[1] = vec[1] - val;
        return v;
    }

    __host__ __device__
    Vec2 operator/(const Vec2& r) const{
        Vec2 v(0, 0);
        v.vec[0] = vec[0] / r.vec[0];
        v.vec[1] = vec[1] / r.vec[1];
        return v;
    }

    __host__ __device__
    Vec2 operator/(const T& val) const{
        Vec2 v(0, 0);
        v.vec[0] = vec[0] / val;
        v.vec[1] = vec[1] / val;
        return v;
    }

    __host__ __device__
    Vec2 operator*(const Vec2& r) const{
        Vec2 v(0, 0);
        v.vec[0] = vec[0] * r.vec[0];
        v.vec[1] = vec[1] * r.vec[1];
        return v;
    }

    __host__ __device__
    Vec2 operator*(const T& val) const{
        Vec2 v(0, 0);
        v.vec[0] = vec[0] * val;
        v.vec[1] = vec[1] * val;
        return v;
    }

    __host__ __device__
    T dot(const Vec2& r) const{
        T x = vec[0] * r.vec[0] + vec[1] * r.vec[1];
        return x;
    }

    __host__ __device__
    T mag() {
        return sqrt(vec[0] * vec[0] + vec[1] * vec[1]);
    }

    __host__ __device__
    T getElement(unsigned short i) const{
        return vec[i];
    }

    __host__ __device__
    void setElement(unsigned short i, T val) {
        vec[i] = val;
    }

    friend std::ostream& operator<<(std::ostream& output, const Vec2& v) {
        output << "(" << v.vec[0] << ", " << v.vec[1] << ")";
        return output;
    }
};

#endif //PROJECT_VEC2_CUH
