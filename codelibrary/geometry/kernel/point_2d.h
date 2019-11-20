//
// Copyright 2011 Yangbin Lin. All Rights Reserved.
//
// Author: yblin@jmu.edu.cn (Yangbin Lin)
//
// This file is part of the Code Library.
//

#ifndef GEOMETRY_KERNEL_POINT_2D_H_
#define GEOMETRY_KERNEL_POINT_2D_H_

#include <cassert>
#include <ostream>

#include "codelibrary/geometry/kernel/box_2d.h"
#include "codelibrary/geometry/kernel/vector_2d.h"

namespace cl {

/**
 * 2D point.
 */
template<typename T>
class Point2D {
public:
    using value_type = T;

    Point2D() = default;

    Point2D(const T& x1, const T& y1)
        : x(x1), y(y1) {}

    bool operator ==(const Point2D& rhs) const {
        return x == rhs.x && y == rhs.y;
    }

    bool operator !=(const Point2D& rhs) const {
        return !(*this == rhs);
    }

    bool operator <(const Point2D& rhs) const {
        return x < rhs.x || (x == rhs.x && y < rhs.y);
    }

    bool operator <=(const Point2D& rhs) const {
        return !(rhs < *this);
    }

    bool operator >(const Point2D& rhs) const {
        return rhs < *this;
    }

    bool operator >=(const Point2D& rhs) const {
        return !(*this < rhs);
    }

    const Point2D& operator +=(const Vector2D<T>& rhs) {
        x += rhs.x;
        y += rhs.y;
        return *this;
    }

    const Point2D& operator -=(const Vector2D<T>& rhs) {
        x -= rhs.x;
        y -= rhs.y;
        return *this;
    }

    /**
     * Return the i-th component value of point.
     */
    const T& operator[] (int i) const {
        return i == 0 ? x : y;
    }

    /**
     * Return the reference value of the i-th component of point.
     */
    T& operator[] (int i) {
        return i == 0 ? x : y;
    }

    /**
     * Return the bounding box of this point.
     */
    Box2D<T> bounding_box() const {
        return Box2D<T>(x, x, y, y);
    }

    /**
     * Return the dimension.
     */
    int size() const {
        return 2;
    }

    Vector2D<T> ToVector() const {
        return Vector2D<T>(x, y);
    }

    friend Point2D operator +(const Point2D& lhs, const Vector2D<T>& rhs) {
        return Point2D<T>(lhs.x + rhs.x, lhs.y + rhs.y);
    }

    friend Point2D operator -(const Point2D& lhs, const Vector2D<T>& rhs) {
        return Point2D<T>(lhs.x - rhs.x, lhs.y - rhs.y);
    }

    friend Vector2D<T> operator -(const Point2D& lhs, const Point2D& rhs) {
        return Vector2D<T>(lhs.x - rhs.x, lhs.y - rhs.y);
    }

    T x = 0; // X coordinate.
    T y = 0; // Y coordinate.
};

/**
 * Compare two 2D points along the given direction.
 */
template <typename T>
class PointDirectionCompare2D {
public:
    explicit PointDirectionCompare2D(const Vector2D<T>& direction)
        : direction_(direction) {}

    bool operator() (const Point2D<T>& p1, const Point2D<T>& p2) const {
        return p1.ToVector() * direction_ < p2.ToVector() * direction_;
    }

private:
    Vector2D<T> direction_;
};

using IPoint2D = Point2D<int>;
using FPoint2D = Point2D<float>;
using RPoint2D = Point2D<double>;

} // namespace cl

#endif // GEOMETRY_KERNEL_POINT_2D_H_
