//
// Copyright 2013 Yangbin Lin. All Rights Reserved.
//
// Author: yblin@jmu.edu.cn (Yangbin Lin)
//
// This file is part of the Code Library.
//

#ifndef GEOMETRY_KERNEL_PLANE_3D_H_
#define GEOMETRY_KERNEL_PLANE_3D_H_

#include "codelibrary/geometry/kernel/line_3d.h"

namespace cl {

/**
 * 3D plane class.
 *
 * A 3D plane can be uniquely defined by a point on plane and a normal vector.
 */
template <typename T>
class Plane3D {
public:
    using value_type = T;

    Plane3D() = default;

    Plane3D(const Point3D<T>& point, const Vector3D<T>& normal)
        : point_(point), normal_(normal) {}

    /**
     * Construct plane by three non collinear points.
     */
    Plane3D(const Point3D<T>& a, const Point3D<T>& b, const Point3D<T>& c) {
        Construct(a, b, c);
    }

    const Vector3D<T>& normal() const { return normal_; }
    const Point3D<T>& point()   const { return point_;  }

private:
    /**
     * Construct 3D plane from three points.
     */
    void Construct(const Point3D<T>& a, const Point3D<T>& b,
                   const Point3D<T>& c) {
        Vector3D<T> v1(b - a), v2(c - a);
        normal_ = CrossProduct(v1, v2);
        point_ = a;
    }

    Point3D<T> point_;   // A point on the plane.
    Vector3D<T> normal_; // The normal vector of plane.
};

using IPlane3D = Plane3D<int>;
using FPlane3D = Plane3D<float>;
using RPlane3D = Plane3D<double>;

} // namespace cl

#endif // GEOMETRY_KERNEL_PLANE_3D_H_
