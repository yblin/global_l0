//
// Copyright 2012 Yangbin Lin. All Rights Reserved.
//
// Author: yblin@jmu.edu.cn (Yangbin Lin)
//
// This file is part of the Code Library.
//
// Distance between 2D geometric kernel objects.
//

#ifndef GEOMETRY_KERNEL_DISTANCE_2D_H_
#define GEOMETRY_KERNEL_DISTANCE_2D_H_

#include <algorithm>
#include <cfloat>
#include <cmath>

#include "codelibrary/geometry/kernel/line_2d.h"
#include "codelibrary/geometry/kernel/multi_polygon_2d.h"
#include "codelibrary/geometry/kernel/point_2d.h"
#include "codelibrary/geometry/kernel/segment_2d.h"

namespace cl {

/**
 * Return the squared Euclidean distance between two points.
 */
template <typename T>
double SquaredDistance(const Point2D<T>& p1, const Point2D<T>& p2) {
    double t1 = static_cast<double>(p1.x) - p2.x;
    double t2 = static_cast<double>(p1.y) - p2.y;
    return t1 * t1 + t2 * t2;
}

/**
 * Return the Euclidean distance between two points.
 */
template <typename T>
double Distance(const Point2D<T>& p1, const Point2D<T>& p2) {
    return std::sqrt(SquaredDistance(p1, p2));
}

/**
 * Return the signed Euclidean distance from point to line.
 */
template <typename T>
double SignedDistance(const Point2D<T>& p, const Line2D<T>& line) {
    double a = line.direction().y;
    double b = -line.direction().x;
    double c = -a * line.point1().x - b * line.point1().y;

    double norm = std::sqrt(a * a + b * b);
    if (norm == 0.0) {
        return Distance(p, line.point1());
    }

    return (a * p.x + b * p.y + c) / norm;
}
template <typename T>
double SignedDistance(const Line2D<T>& line, const Point2D<T>& p) {
    return SignedDistance(p, line);
}

/**
 * Return the Euclidean distance from point to line.
 */
template <typename T>
double Distance(const Point2D<T>& p, const Line2D<T>& line) {
    return std::fabs(SignedDistance(p, line));
}
template <typename T>
double Distance(const Line2D<T>& line, const Point2D<T>& p) {
    return Distance(p, line);
}

/**
 * Return the Euclidean distance from point to segment.
 */
template <typename T>
double Distance(const Point2D<T>& p, const Segment2D<T>& segment) {
    if (segment.lower_point() == segment.upper_point()) {
        return Distance(p, segment.lower_point());
    }

    Vector2D<T> v0 = segment.upper_point() - segment.lower_point();
    Vector2D<T> v1 = p - segment.lower_point();

    double t1 = v0 * v1;
    if (t1 < 0.0) {
        return Distance(p, segment.lower_point());
    }

    double t2 = v0 * v0;
    if (t2 <= t1) {
        return Distance(p, segment.upper_point());
    }

    double b = t1 / t2;
    RPoint2D pb(b * v0.x + segment.lower_point().x,
                b * v0.y + segment.lower_point().y);
    RPoint2D pt(p.x, p.y);

    return Distance(pt, pb);
}
template <typename T>
double Distance(const Segment2D<T>& segment, const Point2D<T>& p) {
    return Distance(p, segment);
}

/**
 * Return the Euclidean distance from point to polygon.
 */
template <typename T>
double Distance(const Point2D<T>& p, const Polygon2D<T>& polygon) {
    double dis = DBL_MAX;

    for (int i = 0; i < polygon.size(); ++i) {
        dis = std::min(dis, Distance(p, polygon.edge(i)));
    }
    return dis;
}
template <typename T>
double Distance(const Polygon2D<T>& polygon, const Point2D<T>& p) {
    return Distance(p, polygon);
}

/**
 * Return the Euclidean distance from point to multi-polygon.
 */
template <typename T>
double Distance(const Point2D<T>& p, const MultiPolygon2D<T>& polygon) {
    double dis = DBL_MAX;

    for (const auto& boundary : polygon.boundaries()) {
        dis = std::min(dis, Distance(p, boundary.polygon));
    }
    return dis;
}
template <typename T>
double Distance(const MultiPolygon2D<T>& polygon, const Point2D<T>& p) {
    return Distance(p, polygon);
}

} // namespace cl

#endif // GEOMETRY_KERNEL_DISTANCE_2D_H_
