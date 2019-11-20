//
// Copyright 2013 Yangbin Lin. All Rights Reserved.
//
// Author: yblin@jmu.edu.cn (Yangbin Lin)
//
// This file is part of the Code Library.
//

#ifndef GEOMETRY_KERNEL_POINT_3D_H_
#define GEOMETRY_KERNEL_POINT_3D_H_

#include <cassert>
#include <ostream>

#include "codelibrary/geometry/kernel/box_3d.h"
#include "codelibrary/geometry/kernel/vector_3d.h"

namespace cl {

/**
 * 3D point.
 */
template<typename T>
class Point3D {
public:
    using value_type = T;

    Point3D() = default;

    Point3D(const T& x1, const T& y1, const T& z1)
        : x(x1), y(y1), z(z1) {}

    bool operator ==(const Point3D& rhs) const {
        return x == rhs.x && y == rhs.y && z == rhs.z;
    }

    bool operator !=(const Point3D& rhs) const {
        return !(*this == rhs);
    }

    bool operator < (const Point3D& rhs) const {
        return x == rhs.x ? (y == rhs.y ? z < rhs.z : y < rhs.y) : x < rhs.x;
    }

    bool operator <=(const Point3D& rhs) const {
        return !(rhs < *this);
    }

    bool operator > (const Point3D& rhs) const {
        return rhs < *this;
    }

    bool operator >=(const Point3D& rhs) const {
        return !(*this < rhs);
    }

    const Point3D& operator +=(const Vector3D<T>& rhs) {
        x += rhs.x;
        y += rhs.y;
        z += rhs.z;
        return *this;
    }

    const Point3D& operator -=(const Vector3D<T>& rhs) {
        x -= rhs.x;
        y -= rhs.y;
        z -= rhs.z;
        return *this;
    }

    /**
     * Return the i-th component value of point.
     */
    const T& operator[] (int i) const {
        return (i == 0) ? x : (i == 1 ? y : z);
    }

    /**
     * Return the reference value of the i-th component of point.
     */
    T& operator[] (int i) {
        return (i == 0) ? x : (i == 1 ? y : z);
    }

    /**
     * Return the bounding box of this point.
     */
    const Box3D<T> bounding_box() const {
        return Box3D<T>(x, x, y, y, z, z);
    }

    /**
     * Return the dimension.
     */
    int size() const {
        return 3;
    }

    Vector3D<T> ToVector() const {
        return Vector3D<T>(x, y, z);
    }

    friend Point3D operator +(const Point3D& lhs, const Vector3D<T>& rhs) {
        return Point3D(lhs.x + rhs.x, lhs.y + rhs.y, lhs.z + rhs.z);
    }

    friend Point3D operator -(const Point3D<T>& lhs, const Vector3D<T>& rhs) {
        return Point3D(lhs.x - rhs.x, lhs.y - rhs.y, lhs.z - rhs.z);
    }

    friend Vector3D<T> operator -(const Point3D& lhs, const Point3D& rhs) {
        return Vector3D<T>(lhs.x - rhs.x, lhs.y - rhs.y, lhs.z - rhs.z);
    }

    T x = 0; // X coordinate.
    T y = 0; // Y coordinate.
    T z = 0; // Z coordinate.
};

/**
 * Compare two 3D points along the given direction.
 */
class PointDirectionCompare3D {
public:
    template <typename T>
    explicit PointDirectionCompare3D(const Vector3D<T>& direction)
        : direction_(direction.x, direction.y, direction.z) {}

    template <typename T>
    bool operator() (const Point3D<T>& p1, const Point3D<T>& p2) const {
        RVector3D v1(p1.x, p1.y, p1.z);
        RVector3D v2(p2.x, p2.y, p2.z);
        return v1 * direction_ < v2 * direction_;
    }

private:
    RVector3D direction_;
};

using IPoint3D = Point3D<int>;
using FPoint3D = Point3D<float>;
using RPoint3D = Point3D<double>;

} // namespace cl

#endif // GEOMETRY_KERNEL_POINT_3D_H_
