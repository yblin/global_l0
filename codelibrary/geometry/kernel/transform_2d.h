//
// Copyright 2019 Yangbin Lin. All Rights Reserved.
//
// Author: yblin@jmu.edu.cn (Yangbin Lin)
//
// This file is part of the Code Library.
//

#ifndef GEOMETRY_KERNEL_TRANSFORM_2D_H_
#define GEOMETRY_KERNEL_TRANSFORM_2D_H_

#include "codelibrary/geometry/kernel/line_2d.h"
#include "codelibrary/geometry/kernel/point_2d.h"

namespace cl {

/**
 * Project point 'p' on a line and return the projection point.
 */
template <typename T1, typename T2>
RPoint2D Project(const Point2D<T1>& p, const Line2D<T2>& line) {
    const RPoint2D q(line.point1().x, line.point1().y);
    const RVector2D v0(line.direction().x, line.direction().y);
    RVector2D v1 = p - q;

    double norm = v0 * v0;
    if (norm == 0.0) return q;

    double b = (v0 * v1) / norm;
    return RPoint2D(b * v0.x + q.x, b * v0.y + q.y);
}

} // namespace cl

#endif // GEOMETRY_KERNEL_TRANSFORM_2D_H_
