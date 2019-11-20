//
// Copyright 2014 Yangbin Lin. All Rights Reserved.
//
// Author: yblin@jmu.edu.cn (Yangbin Lin)
//
// This file is part of the Code Library.
//

#ifndef UTIL_METRIC_ANGULAR_H_
#define UTIL_METRIC_ANGULAR_H_

#include <cmath>

#include "codelibrary/util/metric/cosine.h"

namespace cl {
namespace metric {

/**
 * The metric to measure the angle between two vectors.
 *
 * The range of angular metric is [0, 1].
 */
class Angular {
public:
    Angular() = default;

    Angular(const Angular&) = delete;

    Angular& operator=(const Angular&) = delete;

    template <typename T>
    double operator() (const T& a, const T& b) const {
        return std::acos(cosine_(a, b)) / M_PI;
    }

private:
    Cosine cosine_;
};

} // namespace metric
} // namespace cl

#endif // UTIL_METRIC_ANGULAR_H_
