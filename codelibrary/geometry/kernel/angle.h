//
// Copyright 2012 Yangbin Lin. All Rights Reserved.
//
// Author: yblin@jmu.edu.cn (Yangbin Lin)
//
// This file is part of the Code Library.
//

#ifndef GEOMETRY_KERNEL_ANGLE_H_
#define GEOMETRY_KERNEL_ANGLE_H_

#include <cmath>

#include "codelibrary/base/algorithm.h"
#include "codelibrary/base/units.h"
#include "codelibrary/geometry/kernel/vector_2d.h"

namespace cl {

/**
 * Get the radian angle of 2D vector.
 */
template <typename T>
inline double Radian(const Vector2D<T>& v) {
    double radian = std::atan2(v.y, v.x);
    return radian < 0.0 ? radian + M_PI  + M_PI : radian;
}

/**
 * Get the degree angle of 2D vector.
 */
template <typename T>
inline double Degree(const Vector2D<T>& v) {
    return unit::angle::RadianToDegree(Radian(v));
}

/**
 * Get the radian angle between two vectors.
 *
 * The return value is in the range: [0.0, PI].
 */
template <class Vector>
double Radian(const Vector& v1, const Vector& v2) {
    double l1 = v1.norm();
    double l2 = v2.norm();

    if (l1 == 0.0 || l2 == 0.0) return 0.0;
    return std::acos(Clamp(v1 * v2 / l1 / l2, -1.0, 1.0));
}

/**
 * Get the degree between two vectors.
 */
template <typename Vector>
double Degree(const Vector& v1, const Vector& v2) {
    return unit::angle::RadianToDegree(Radian(v1, v2));
}

} // namespace cl

#endif // GEOMETRY_KERNEL_ANGLE_H_
