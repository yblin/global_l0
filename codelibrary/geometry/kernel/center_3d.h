//
// Copyright 2014 Yangbin Lin. All Rights Reserved.
//
// Author: yblin@jmu.edu.cn (Yangbin Lin)
//
// This file is part of the Code Library.
//
// Computer the center of geometric objects.
//

#ifndef GEOMETRY_KERNEL_CENTER_3D_H_
#define GEOMETRY_KERNEL_CENTER_3D_H_

#include "codelibrary/base/array.h"
#include "codelibrary/geometry/kernel/point_3d.h"
#include "codelibrary/geometry/kernel/segment_3d.h"

namespace cl {

/**
 * Compute the centroid of two 3D points.
 */
template <typename T>
Point3D<T> Centroid(const Point3D<T>& p1, const Point3D<T>& p2) {
    static_assert(std::is_floating_point<T>::value,
                  "T must be a floating point.");

    return Point3D<T>((p1.x + p2.x) / T(2), (p1.y + p2.y) / T(2),
                      (p1.z + p2.z) / T(2));
}

/**
 * Compute the centroid of given 3D bounding box.
 */
template <typename T>
Point3D<T> Centroid(const Box3D<T>& box) {
    return Centroid(Point3D<T>(box.x_min(), box.y_min(), box.z_min()),
                    Point3D<T>(box.x_max(), box.y_max(), box.z_max()));
}

/**
 * Compute the centroid of given 3D line segment.
 */
template <typename T>
Point3D<T> Centroid(const Segment3D<T>& line) {
    return Centroid(line.lower_point(), line.upper_point());
}

/**
 * Compute the centroid of the given 3D points.
 */
template <typename T>
Point3D<T> Centroid(const Array<Point3D<T>>& points) {
    static_assert(std::is_floating_point<T>::value,
                  "T must be a floating point.");

    assert(!points.empty());

    T x = 0, y = 0, z = 0;
    for (const Point3D<T>& p : points) {
        x += p.x;
        y += p.y;
        z += p.z;
    }
    return Point3D<T>(x / points.size(), y / points.size(), z / points.size());
}

/**
 * Compute weighted centroid of given 3D points.
 */
template <typename T1, typename T2>
Point3D<T2> Centroid(const Array<Point3D<T1>>& points,
                     const Array<T2>& weights) {
    static_assert(std::is_floating_point<T2>::value,
                  "T2 must be a floating point.");

    assert(!points.empty());
    assert(points.size() == weights.size());

    T2 sum = 0, x = 0, y = 0, z = 0;
    for (int i = 0; i < points.size(); ++i) {
        const Point3D<T2>& p = points[i];
        const T2& w = weights[i];
        assert(w >= 0);

        x += w * p.x;
        y += w * p.y;
        z += w * p.z;
        sum += w;
    }
    assert(sum != 0);

    return Point3D<T2>(x / sum, y / sum, z / sum);
}

} // namespace cl

#endif // GEOMETRY_KERNEL_CENTER_3D_H_
