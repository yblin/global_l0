//
// Copyright 2014 Yangbin Lin. All Rights Reserved.
//
// Author: yblin@jmu.edu.cn (Yangbin Lin)
//
// This file is part of the Code Library.
//
// Robust intersection test between 3D kernel geometric objects.
//

#ifndef GEOMETRY_KERNEL_INTERSECT_3D_H_
#define GEOMETRY_KERNEL_INTERSECT_3D_H_

#include <algorithm>
#include <cfloat>
#include <utility>

#include "codelibrary/geometry/kernel/box_3d.h"
#include "codelibrary/geometry/kernel/line_3d.h"
#include "codelibrary/geometry/kernel/plane_3d.h"
#include "codelibrary/geometry/kernel/point_3d.h"
#include "codelibrary/geometry/kernel/ray_3d.h"
#include "codelibrary/geometry/kernel/segment_3d.h"
#include "codelibrary/geometry/kernel/sphere_3d.h"

namespace cl {
namespace geometry {

template <typename T>
bool Intersect(const Box3D<T>& box1, const Box3D<T>& box2) {
    return box1.x_max() >= box2.x_min() && box1.x_min() <= box2.x_max() &&
           box1.y_max() >= box2.y_min() && box1.y_min() <= box2.y_max() &&
           box1.z_max() >= box2.z_min() && box1.z_min() <= box2.z_max();
}

template <typename T>
bool Intersect(const Point3D<T>& point, const Box3D<T>& box) {
    return (point.x >= box.x_min() && point.x <= box.x_max() &&
            point.y >= box.y_min() && point.y <= box.y_max() &&
            point.z >= box.z_min() && point.z <= box.z_max());
}

template <typename T>
bool Intersect(const Box3D<T>& box, const Point3D<T>& point) {
    return Intersect(point, box);
}

template <typename T>
bool Cross(const Box3D<T>& box, const Ray3D<T>& ray,
           std::pair<double, double>* intersection = nullptr) {
    if (!box.is_valid()) return false;

    double t_min = 0, t_max = DBL_MAX;
    for (int i = 0; i < 3; ++i) {
        // It handles the inf cases correct.
        double inverse_direction = 1.0 / ray.direction()[i];
        double t1 = (box.min(i) - ray.origin()[i]) * inverse_direction;
        double t2 = (box.max(i) - ray.origin()[i]) * inverse_direction;
        if (inverse_direction < 0.0) {
            std::swap(t1, t2);
        }
        t_min = std::max(t1, t_min);
        t_max = std::max(t2, t_max);
        if (t_max < t_min) return false;
    }
    if (intersection) {
        *intersection = std::make_pair(t_min, t_max);
    }
    return true;
}

template <typename T>
bool Cross(const Ray3D<T>& ray, const Box3D<T>& box,
           std::pair<double, double>* intersection = nullptr) {
    return Cross(box, ray, intersection);
}

template <typename T>
bool Cross(const Box3D<T>& box, const Line3D<T>& line,
           std::pair<double, double>* intersection = nullptr) {
    if (!box.is_valid()) return false;

    double t_min = -DBL_MAX, t_max = DBL_MAX;
    RVector2D v(line.direction().x, line.direction().y, line.direction().z);
    for (int i = 0; i < 3; ++i) {
        // It handles the inf case correctly.
        double inverse_direction = 1.0 / v[i];
        double t1 = (box.min(i) - line.point1()[i]) * inverse_direction;
        double t2 = (box.max(i) - line.point1()[i]) * inverse_direction;
        if (inverse_direction < 0.0) {
            std::swap(t1, t2);
        }
        t_min = std::max(t1, t_min);
        t_max = std::min(t2, t_max);
        if (t_max < t_min) return false;
    }
    if (intersection) {
        *intersection = std::make_pair(t_min, t_max);
    }
    return true;
}

template <typename T>
bool Cross(const Line3D<T>& line, const Box3D<T>& box,
           std::pair<double, double>* intersection = nullptr) {
    return Cross(box, line, intersection);
}

template <typename T>
bool Cross(const Box3D<T>& box, const Segment3D<T>& segment,
           std::pair<double, double>* intersection = nullptr) {
    if (!box.is_valid()) return false;

    RVector3D direction(segment.direction().x, segment.direction().y,
                        segment.direction().z);
    double t_min = 0.0, t_max = 1.0;
    for (int i = 0; i < 3; ++i) {
        // It handles the inf cases correct.
        double inverse_direction = 1.0 / direction[i];
        double t1 = (box.min(i) - segment.lower_point()[i]) * inverse_direction;
        double t2 = (box.max(i) - segment.lower_point()[i]) * inverse_direction;
        if (inverse_direction < 0.0) {
            std::swap(t1, t2);
        }
        t_min = std::max(t_min, t1);
        t_max = std::min(t_max, t2);
        if (t_max < t_min) return false;
    }
    if (intersection) {
        *intersection = std::make_pair(t_min, t_max);
    }
    return true;
}

template <typename T>
bool Cross(const Segment3D<T>& segment, const Box3D<T>& box,
           std::pair<double, double>* intersection = nullptr) {
    return Cross(box, segment, intersection);
}

template <typename T>
bool Cross(const Line3D<T>& line, const Sphere3D<T>& sphere,
           std::pair<double, double>* intersection = nullptr) {
    double a = line.direction() * line.direction();
    if (a == 0) return false;

    Vector3D<T> diff = sphere.center() - line.point1();
    double b = line.direction() * diff;
    double c = diff * diff - sphere.radius() * sphere.radius();

    double discriminant = b * b - a * c;
    if (discriminant >= 0.0) {
        if (intersection) {
            double root = std::sqrt(discriminant);
            *intersection = std::make_pair((b - root) / a, (b + root) / a);
        }
        return true;
    }

    return false;
}

template <typename T>
bool Cross(const Sphere3D<T>& sphere, const Line3D<T>& line,
           std::pair<double, double>* intersection = nullptr) {
    return Cross(line, sphere, intersection);
}

template <typename T>
bool Cross(const Line3D<T>& line, const Plane3D<T>& plane,
           RPoint3D* intersection_point = nullptr) {
    const Vector3D<T>& u = line.direction();
    const Vector3D<T> w = line.point1() - plane.point();

    double d =  plane.normal() * u;
    double n = -(plane.normal() * w);

    // They are not parallel.
    double t = n / d;
    if (!std::isfinite(t)) return false;

    if (intersection_point) {
        intersection_point->x = line.point1().x + t * u.x;
        intersection_point->y = line.point1().y + t * u.y;
        intersection_point->z = line.point1().z + t * u.z;
    }

    return true;
}
template <typename T>
bool Cross(const Plane3D<T>& plane, const Line3D<T>& line,
           RPoint3D* intersection_point = nullptr) {
    return Cross(line, plane, intersection_point);
}

template <typename T>
bool Cross(const Ray3D<T>& ray, const Sphere3D<T>& sphere,
           std::pair<double, double>* intersection = nullptr) {
    double a = ray.direction() * ray.direction();
    if (a == 0.0) return false;

    Vector3D<T> diff = sphere.center() - ray.origin();
    double b = ray.direction() * diff;
    double c = diff * diff - sphere.radius() * sphere.radius();

    double discriminant = b * b - a * c;
    if (discriminant >= 0.0) {
        double root = std::sqrt(discriminant);
        if (b + root < 0.0) return false;

        if (intersection) {
            if (b - root <= 0.0) {
                *intersection = std::make_pair(0.0, (b + root) / a);
            } else {
                *intersection = std::make_pair((b - root) / a, (b + root) / a);
            }
        }
        return true;
    }

    return false;
}

template <typename T>
bool Cross(const Sphere3D<T>& sphere, const Ray3D<T>& ray,
           std::pair<double, double>* intersection = nullptr) {
    return Cross(ray, sphere, intersection);
}

template <typename T>
bool Cross(const Plane3D<T>& plane1, const Plane3D<T>& plane2,
           RLine3D* intersection = nullptr) {
    RVector3D v1(plane1.normal().x, plane1.normal().y, plane1.normal().z);
    RVector3D v2(plane2.normal().x, plane2.normal().y, plane2.normal().z);
    const RVector3D u = CrossProduct(v1, v2);
    double as[3] = { std::fabs(u.x), std::fabs(u.y), std::fabs(u.z) };

    // First, determine the maximum absolute value of cross production.
    int maxc = static_cast<int>(std::max_element(as + 0, as + 3) - as);

    // Next, compute a point on the intersection line.
    RPoint3D p;
    double d1 = -(v1 * plane1.point().ToVector());
    double d2 = -(v2 * plane2.point().ToVector());

    switch (maxc) {
    case 0: // Intersected with plane x = 0.
        p.x = 0.0;
        p.y = (d2 * v1.z - d1 * v2.z) / u.x;
        p.z = (d1 * v2.y - d2 * v1.y) / u.x;
        break;
    case 1: // Intersected with plane y = 0.
        p.x = (d1 * v2.z - d2 * v1.z) / u.y;
        p.y = 0.0;
        p.z = (d2 * v1.x - d1 * v2.x) / u.y;
        break;
    case 2: // Intersected with plane z = 0.
        p.x = (d2 * v1.y - d1 * v2.y) / u.z;
        p.y = (d1 * v2.x - d2 * v1.x) / u.z;
        p.z = 0.0;
    }

    if (intersection) *intersection = RLine3D(p, u);

    return std::isfinite(p.x) && std::isfinite(p.y) && std::isfinite(p.z);
}

} // namespace geometry
} // namespace cl

#endif // GEOMETRY_KERNEL_INTERSECT_3D_H_
