//
// Copyright 2013 Yangbin Lin. All Rights Reserved.
//
// Author: yblin@jmu.edu.cn (Yangbin Lin)
//
// This file is part of the Code Library.
//

#ifndef GEOMETRY_KERNEL_VECTOR_3D_H_
#define GEOMETRY_KERNEL_VECTOR_3D_H_

#include <cassert>
#include <cfloat>
#include <cmath>

#include "codelibrary/base/equal.h"

namespace cl {

/**
 * 3D vector.
 */
template<typename T>
class Vector3D {
public:
    using value_type = T;

    Vector3D() = default;

    Vector3D(const T& x1, const T& y1, const T& z1)
        : x(x1), y(y1), z(z1) {}

    bool operator ==(const Vector3D& rhs) const {
        return x == rhs.x && y == rhs.y && z == rhs.z;
    }

    bool operator !=(const Vector3D& rhs) const {
        return !(*this == rhs);
    }

    bool operator < (const Vector3D& rhs) const {
        return x == rhs.x ? (y == rhs.y ? z < rhs.z : y < rhs.y) : x < rhs.x;
    }

    bool operator <=(const Vector3D& rhs) const {
        return !(rhs < *this);
    }

    bool operator > (const Vector3D& rhs) const {
        return rhs < *this;
    }

    bool operator >=(const Vector3D& rhs) const {
        return !(*this < rhs);
    }

    const Vector3D& operator +=(const Vector3D& rhs) {
        x += rhs.x;
        y += rhs.y;
        z += rhs.z;
        return *this;
    }

    const Vector3D& operator -=(const Vector3D& rhs) {
        x -= rhs.x;
        y -= rhs.y;
        z -= rhs.z;
        return *this;
    }

    const Vector3D& operator *=(const T& rhs) {
        x *= rhs;
        y *= rhs;
        z *= rhs;
        return *this;
    }

    /**
     * Return the squared euclidean norm of the vector.
     */
    double squared_norm() const {
        return static_cast<double>(x) * x + static_cast<double>(y) * y +
               static_cast<double>(z) * z;
    }

    /**
     * Return the euclidean norm of the vector.
     */
    double norm() const {
        return std::sqrt(squared_norm());
    }

    /**
     * Return the i-th component value of vector.
     */
    const T& operator[] (int i) const {
        return (i == 0) ? x : (i == 1 ? y : z);
    }

    /**
     * Return the reference value of the i-th component of vector.
     */
    T& operator[] (int i) {
        return (i == 0) ? x : (i == 1 ? y : z);
    }

    int size() const { return 3; }

    
    /**
     * Normalize the vector.
     * The length of the normalized vector is always one.
     */
    void Normalize() {
        static_assert(std::is_floating_point<T>::value,
                      "T must be a floating point.");

        double len = this->norm();
        if (Equal(len, 0.0)) {
            x = 0;
            y = 0;
            z = 1;
        } else {
            T t = static_cast<T>(1.0 / len);
            x *= t;
            y *= t;
            z *= t;
        }
    }

    friend Vector3D operator +(const Vector3D& lhs, const Vector3D& rhs) {
        return Vector3D(lhs.x + rhs.x, lhs.y + rhs.y, lhs.z + rhs.z);
    }

    friend Vector3D operator -(const Vector3D& lhs, const Vector3D& rhs) {
        return Vector3D(lhs.x - rhs.x, lhs.y - rhs.y, lhs.z - rhs.z);
    }

    friend Vector3D operator -(const Vector3D& rhs) {
        return Vector3D(-rhs.x, -rhs.y, -rhs.z);
    }

    friend Vector3D operator *(const T& lhs, const Vector3D& rhs) {
        return Vector3D(lhs * rhs.x, lhs * rhs.y, lhs * rhs.z);
    }

    friend Vector3D operator *(const Vector3D& lhs, const T& rhs) {
        return Vector3D(lhs.x * rhs, lhs.y * rhs, lhs.z * rhs);
    }

    friend double operator *(const Vector3D& lhs, const Vector3D& rhs) {
        return static_cast<double>(lhs.x) * rhs.x +
               static_cast<double>(lhs.y) * rhs.y +
               static_cast<double>(lhs.z) * rhs.z;
    }

    /**
     * Return the cross product of two vectors.
     */
    template <typename T1>
    friend Vector3D<T1> CrossProduct(const Vector3D<T1>& v1, 
                                     const Vector3D<T1>& v2);

    friend std::ostream& operator <<(std::ostream& os, const Vector3D& v) {
        os << fmt::format("({}, {}, {})", v.x, v.y, v.z);
        return os;
    }

    T x = 0; // X component.
    T y = 0; // Y component.
    T z = 0; // Z component.
};

using IVector3D = Vector3D<int>;
using FVector3D = Vector3D<float>;
using RVector3D = Vector3D<double>;

/**
 * Return the cross product of two vectors.
 */
template <typename T>
Vector3D<T> CrossProduct(const Vector3D<T>& v1, const Vector3D<T>& v2) {
    Vector3D<T> v;
    v.x = v1.y * v2.z - v1.z * v2.y;
    v.y = v1.z * v2.x - v1.x * v2.z;
    v.z = v1.x * v2.y - v1.y * v2.x;
    return v;
}

} // namespace cl

#endif // MATH_VECTOR_VECTOR_3D_H_
