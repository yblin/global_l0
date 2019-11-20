//
// Copyright 2011 Yangbin Lin. All Rights Reserved.
//
// Author: yblin@jmu.edu.cn (Yangbin Lin)
//
// This file is part of the Code Library.
//

#ifndef GEOMETRY_KERNEL_VECTOR_2D_H_
#define GEOMETRY_KERNEL_VECTOR_2D_H_

#include <cassert>
#include <cmath>

namespace cl {

/**
 * 2D vector.
 */
template<typename T>
class Vector2D {
public:
    using value_type = T;

    Vector2D() = default;

    Vector2D(const T& x1, const T& y1)
        : x(x1), y(y1) {}

    bool operator ==(const Vector2D& rhs) const {
        return x == rhs.x && y == rhs.y;
    }

    bool operator !=(const Vector2D& rhs) const {
        return !(*this == rhs);
    }

    bool operator < (const Vector2D& rhs) const {
        return x < rhs.x || (x == rhs.x && y < rhs.y);
    }

    bool operator <=(const Vector2D& rhs) const {
        return !(rhs < *this);
    }

    bool operator > (const Vector2D& rhs) const {
        return rhs < *this;
    }

    bool operator >=(const Vector2D& rhs) const {
        return !(*this < rhs);
    }

    const Vector2D& operator +=(const Vector2D& rhs) {
        x += rhs.x;
        y += rhs.y;
        return *this;
    }

    const Vector2D& operator -=(const Vector2D& rhs) {
        x -= rhs.x;
        y -= rhs.y;
        return *this;
    }

    const Vector2D& operator *=(const T& rhs) {
        x *= rhs;
        y *= rhs;
        return *this;
    }

    /**
     * Return the squared euclidean norm of the vector.
     */
    double squared_norm() const {
        return static_cast<double>(x) * x + static_cast<double>(y) * y;
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
        return i == 0 ? x : y;
    }

    /**
     * Return the reference value of the i-th component of vector.
     */
    T& operator[] (int i) {
        return i == 0 ? x : y;
    }

    int size() const { return 2; }

    friend Vector2D operator +(const Vector2D& lhs, const Vector2D& rhs) {
        return Vector2D(lhs.x + rhs.x, lhs.y + rhs.y);
    }

    friend Vector2D operator -(const Vector2D& lhs, const Vector2D& rhs) {
        return Vector2D(lhs.x - rhs.x, lhs.y - rhs.y);
    }

    friend Vector2D operator -(const Vector2D& rhs) {
        return Vector2D(-rhs.x, -rhs.y);
    }

    friend Vector2D operator *(const T& lhs, const Vector2D& rhs) {
        return Vector2D(lhs * rhs.x, lhs * rhs.y);
    }

    friend Vector2D operator *(const Vector2D& lhs, const T& rhs) {
        return Vector2D(lhs.x * rhs, lhs.y * rhs);
    }

    friend double operator *(const Vector2D& lhs, const Vector2D& rhs) {
        return static_cast<double>(lhs.x) * rhs.x +
               static_cast<double>(lhs.y) * rhs.y;
    }
    
    friend std::ostream& operator <<(std::ostream& os, const Vector2D& v) {
        os << fmt::format("({}, {})", v.x, v.y);
        return os;
    }

    T x = 0; // X component.
    T y = 0; // Y component.
};

using IVector2D = Vector2D<int>;
using FVector2D = Vector2D<float>;
using RVector2D = Vector2D<double>;

} // namespace cl

#endif // GEOMETRY_KERNEL_VECTOR_2D_H_
