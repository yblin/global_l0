//
// Copyright 2016 Yangbin Lin. All Rights Reserved.
//
// Author: yblin@jmu.edu.cn (Yangbin Lin)
//
// This file is part of the Code Library.
//

#ifndef GEOMETRY_KERNEL_RAY_3D_H_
#define GEOMETRY_KERNEL_RAY_3D_H_

#include "codelibrary/geometry/kernel/segment_3d.h"
#include "codelibrary/math/vector.h"

namespace cl {

/**
 * A 3D ray can be uniquely defined by its origin point and direction vector.
 */
template <typename T>
class Ray3D {
public:
    using value_type = T;

    Ray3D() = default;

    /**
     * Construct from a origin point and a direction vector.
     */
    Ray3D(const Point3D<T>& origin, const Vector3D<T>& direction)
        : origin_(origin), direction_(direction) {}

    /**
     * Construct from two distinct points.
     */
    Ray3D(const Point3D<T>& p1, const Point3D<T>& p2)
        : Ray3D(p1, p2 - p1) {}

    const Vector3D<T>& direction() const {
        return direction_;
    }

    const Point3D<T>& origin() const {
        return origin_;
    }

private:
    Point3D<T> origin_;     // The origin point of the ray.
    Vector3D<T> direction_; // The direction vector of line.
};

using IRay3D = Ray3D<int>;
using FRay3D = Ray3D<float>;
using RRay3D = Ray3D<double>;

} // namespace cl

#endif // GEOMETRY_KERNEL_RAY_3D_H_
