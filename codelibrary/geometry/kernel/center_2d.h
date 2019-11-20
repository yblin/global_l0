//
// Copyright 2014 Yangbin Lin. All Rights Reserved.
//
// Author: yblin@jmu.edu.cn (Yangbin Lin)
//
// This file is part of the Code Library.
//
// Computer the center of the geometric objects.
//

#ifndef GEOMETRY_KERNEL_CENTER_2D_H_
#define GEOMETRY_KERNEL_CENTER_2D_H_

#include "codelibrary/base/array.h"
#include "codelibrary/geometry/kernel/point_2d.h"
#include "codelibrary/geometry/kernel/segment_2d.h"

namespace cl {

/**
 * Compute the centroid of two 2D points.
 */
template <typename T>
Point2D<T> Centroid(const Point2D<T>& p1, const Point2D<T>& p2) {
    static_assert(std::is_floating_point<T>::value,
                  "T must be a floating point.");

    T x = (p1.x + p2.x) / T(2);
    T y = (p1.y + p2.y) / T(2);
    return Point2D<T>(x, y);
}

/**
 * Compute the centroid of given 2D line segment.
 */
template <typename T>
Point2D<T> Centroid(const Segment2D<T>& line) {
    return Centroid(line.lower_point(), line.upper_point());
}

/**
 * Compute the centroid of given 2D box.
 */
template <typename T>
Point2D<T> Centroid(const Box2D<T>& box) {
    return Centroid(Point2D<T>(box.x_min(), box.y_min()),
                    Point2D<T>(box.x_max(), box.y_max()));
}

/**
 * Compute the centroid of the given 2D points.
 */
template <typename T>
Point2D<T> Centroid(const Array<Point2D<T>>& points) {
    static_assert(std::is_floating_point<T>::value,
                  "T must be a floating point.");

    assert(!points.empty());

    T x = 0, y = 0;
    for (const Point2D<T>& p : points) {
        x += p.x;
        y += p.y;
    }
    return Point2D<T>(x / points.size(), y / points.size());
}

/**
 * Compute weighted centroid of given 2D points.
 */
template <typename T1, typename T2>
Point2D<T2> Centroid(const Array<Point2D<T1>>& points,
                     const Array<T2>& weights) {
    static_assert(std::is_floating_point<T2>::value,
                  "T must be a floating point.");

    assert(!points.empty());
    assert(points.size() == weights.size());

    T2 sum = 0, x = 0, y = 0;
    for (int i = 0; i < points.size(); ++i) {
        const Point2D<T2>& p = points[i];
        const T2& w = weights[i];
        assert(w >= 0);

        x += w * p.x;
        y += w * p.y;
        sum += w;
    }
    assert(sum != 0);

    return Point2D<T2>(x / sum, y / sum);
}

} // namespace cl

#endif // GEOMETRY_KERNEL_CENTER_2D_H_
