//
// Copyright 2013 Yangbin Lin. All Rights Reserved.
//
// Author: yblin@jmu.edu.cn (Yangbin Lin)
//
// This file is part of the Code Library.
//

#ifndef GEOMETRY_KERNEL_SEGMENT_3D_H_
#define GEOMETRY_KERNEL_SEGMENT_3D_H_

#include <algorithm>

#include "codelibrary/geometry/kernel/box_3d.h"
#include "codelibrary/geometry/kernel/point_3d.h"

namespace cl {

/**
 * 3D line segment.
 */
template <typename T>
class Segment3D {
public:
    using Point = Point3D<T>;
    using value_type = T;

    Segment3D() = default;

    Segment3D(const Point& p1, const Point& p2)
        : lower_point_(std::min(p1, p2)),
          upper_point_(std::max(p1, p2)),
          bounding_box_(std::min(p1.x, p2.x), std::max(p1.x, p2.x),
                        std::min(p1.y, p2.y), std::max(p1.y, p2.y),
                        std::min(p1.z, p2.z), std::max(p1.z, p2.z)) {}

    /**
     * The left bottom endpoint of line segment.
     */
    const Point& lower_point() const {
        return lower_point_;
    }

    /**
     * The right top endpoint of line segment.
     */
    const Point& upper_point() const {
        return upper_point_;
    }

    /**
     * The bounding box of line segment.
     */
    const Box3D<T>& bounding_box() const {
        return bounding_box_;
    }

    /**
     * Return the direction of segment.
     */
    Vector3D<T> direction() const {
        return upper_point_ - lower_point_;
    }

    /**
     * Return the length of segment.
     */
    double length() const {
        return direction().norm();
    }

    bool operator ==(const Segment3D& rhs) const {
        return lower_point_ == rhs.lower_point_ &&
               upper_point_ == rhs.upper_point_;
    }

    bool operator !=(const Segment3D& rhs) const {
        return !(*this == rhs);
    }

    bool operator <(const Segment3D& rhs) const {
        return lower_point_ < rhs.lower_point_ ||
               (lower_point_ == rhs.lower_point_ &&
                upper_point_ < rhs.upper_point_);
    }

private:
    Point lower_point_;     // The left bottom endpoint of line segment.
    Point upper_point_;     // The right top endpoint of line segment.
    Box3D<T> bounding_box_; // The bounding box of line segment.
};

using ISegment3D = Segment3D<int>;
using FSegment3D = Segment3D<float>;
using RSegment3D = Segment3D<double>;

} // namespace cl

#endif // GEOMETRY_KERNEL_SEGMENT_3D_H_
