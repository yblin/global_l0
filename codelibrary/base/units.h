//
// Copyright 2019 Yangbin Lin. All Rights Reserved.
//
// Author: yblin@jmu.edu.cn (Yangbin Lin)
//
// This file is part of the Code Library.
//

#ifndef BASE_UNITS_H_
#define BASE_UNITS_H_

#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#include <math.h>
#else
#include <cmath>
#endif // _USE_MATH_DEFINES

#include "codelibrary/base/format.h"

namespace cl {
namespace unit {

namespace time {

/**
 * Convert a time in seconds to human readable string.
 */
inline std::string ToString(double time) {
    if (time < 1e-6) {
        return std::to_string(static_cast<int>(time * 1e9)) + " ns";
    }

    if (time < 1e-3) {
        return std::to_string(static_cast<int>(time * 1e6)) + " us";
    }

    if (time < 1.0) {
        return std::to_string(static_cast<int>(time * 1e3)) + " ms";
    }

    if (time < 60.0) {
        return fmt::format("{:0.3g} s", time);
    }

    time /= 60.0;
    if (time < 60.0) {
        return fmt::format("{:0.3g} min", time);
    }

    time /= 60.0;
    if (time < 24.0) {
        return fmt::format("{:0.3g} h", time);
    }

    time /= 24.0;
    if (time < 30.0) {
        return fmt::format("{:0.3g} days", time);
    }

    if (time < 365.2425) {
        return fmt::format("{:0.3g} months", time / 30.43687);
    }

    return fmt::format("{:0.3g} years", time / 365.2425);
}

} // namespace time

namespace angle {

/**
 * Convert degree angle to radian angle.
 */
inline constexpr double DegreeToRadian(double degree) {
    return degree * M_PI / 180.0;
}

/**
 * Convert radian angle to degree angle.
 */
inline constexpr double RadianToDegree(double radian) {
    return radian * 180.0 / M_PI;
}

} // namespace angle

} // namespace unit
} // namespace cl

#endif // BASE_UNITS_H_
