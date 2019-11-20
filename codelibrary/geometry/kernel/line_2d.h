//
// Copyright 2011 Yangbin Lin. All Rights Reserved.
//
// Author: yblin@jmu.edu.cn (Yangbin Lin)
//
// This file is part of the Code Library.
//

#ifndef GEOMETRY_KERNEL_LINE_2D_H_
#define GEOMETRY_KERNEL_LINE_2D_H_

#include "codelibrary/geometry/kernel/segment_2d.h"
#include "codelibrary/geometry/kernel/vector_2d.h"

namespace cl {

/**
 * 2D line class.
 *
 * A 2D line can be uniquely defined by a point and a direction vector.
 */
template <typename T>
class Line2D {
public:
    using value_type = T;

    Line2D() = default;

    /**
     * Construct from a point on the line and direction of line.
     */
    Line2D(const Point2D<T>& p, const Vector2D<T>& direction)
        : point1_(p), point2_(p + direction) {}

    /**
     * Construct from two distinct points.
     */
    Line2D(const Point2D<T>& p1, const Point2D<T>& p2)
        : point1_(p1), point2_(p2) {}

    /**
     * Construct from segment.
     */
    explicit Line2D(const Segment2D<T>& segment)
        : Line2D(segment.lower_point(), segment.upper_point()) {}

    Vector2D<T> direction() const {
        return point2_ - point1_;
    }

    const Point2D<T>& point1() const {
        return point1_;
    }

    const Point2D<T>& point2() const {
        return point2_;
    }

protected:
    Point2D<T> point1_, point2_; // Two points on the line.
};

using ILine2D = Line2D<int>;
using FLine2D = Line2D<float>;
using RLine2D = Line2D<double>;

} // namespace cl

#endif // GEOMETRY_KERNEL_LINE_2D_H_
