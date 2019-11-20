//
// Copyright 2019 Yangbin Lin. All Rights Reserved.
//
// Author: yblin@jmu.edu.cn (Yangbin Lin)
//
// This file is part of the Code Library.
//

#ifndef GEOMETRY_KERNEL_TRANSFORM_3D_H_
#define GEOMETRY_KERNEL_TRANSFORM_3D_H_

#include "codelibrary/geometry/kernel/line_3d.h"
#include "codelibrary/geometry/kernel/plane_3d.h"
#include "codelibrary/geometry/kernel/point_3d.h"
#include "codelibrary/geometry/kernel/quaternion.h"
#include "codelibrary/geometry/kernel/sphere_3d.h"

namespace cl {

/**
 * Return the orthogonal projection of a point 'p' on a line.
 */
template <typename T1, typename T2>
RPoint3D Project(const Point3D<T1>& p, const Line3D<T2>& line) {
    const Point3D<T2>& q = line.point1();
    const Vector3D<T2>& v = line.direction();

    RVector3D v0(v.x, v.y, v.z);
    double norm = v0 * v0;
    if (norm == 0.0) return line.point1();

    RVector3D v1(static_cast<double>(p.x) - q.x,
                 static_cast<double>(p.y) - q.y,
                 static_cast<double>(p.z) - q.z);

    double b = (v0 * v1) / norm;
    return RPoint3D(b * v0.x + q.x, b * v0.y + q.y, b * v0.z + q.z);
}

/**
 * Return the orthogonal projection of a point 'p' on a plane.
 */
template <typename T1, typename T2>
RPoint3D Project(const Point3D<T1>& p, const Plane3D<T2>& plane) {
    const Vector3D<T2>& normal = plane.normal();
    RVector3D direction(normal.x, normal.y, normal.z);
    RLine3D line(RPoint3D(p.x, p.y, p.z), direction);
    return Project(plane.point(), line);
}

/**
 * Return the projection of a point 'p' on a sphere.
 */
template <typename T1, typename T2>
RPoint3D Project(const Point3D<T1>& p, const Sphere3D<T2>& sphere) {
    double radius = sphere.radius();
    RPoint3D c(sphere.center().x, sphere.center().y, sphere.center().z);
    RPoint3D q(p.x, p.y, p.z);

    if (q == c) {
        // If p is the center of sphere, return the point that has with maximal
        // Z value.
        return c + RVector3D(0.0, 0.0, radius);
    }

    RVector3D v = q - c;
    v *= radius / v.norm();

    return c + v;
}

/**
 * Rotate the given point.
 */
template <typename T>
RPoint3D Rotate(const Point3D<T>& point, const Quaternion& rotation) {
    Quaternion q(point.x, point.y, point.z, 0.0);
    q = rotation * q * rotation.Reverse();

    return RPoint3D(q.x(), q.y(), q.z());
}

/**
 * Rotate the given plane.
 */
template <typename T>
RPlane3D Rotate(const Plane3D<T>& plane, const Quaternion& rotation) {
    return RPlane3D(Rotate(plane.point(), rotation),
                    rotation.Rotate(plane.normal()));
}

/**
 * Rotate the given segment.
 */
template <typename T>
RSegment3D Rotate(const Segment3D<T>& seg, const Quaternion& rotation) {
    return RSegment3D(Rotate(seg.lower_point(), rotation),
                      Rotate(seg.upper_point(), rotation));
}

/**
 * Rotate the given line.
 */
template <typename T>
RLine3D Rotate(const Line3D<T>& line, const Quaternion& rotation) {
    return RLine3D(Rotate(line.point(), rotation),
                   rotation.Rotate(line.direction()));
}

} // namespace cl

#endif // GEOMETRY_KERNEL_TRANSFORM_3D_H_
