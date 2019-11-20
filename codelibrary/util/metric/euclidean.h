//
// Copyright 2014 Yangbin Lin. All Rights Reserved.
//
// Author: yblin@jmu.edu.cn (Yangbin Lin)
//
// This file is part of the Code Library.
//

#ifndef UTIL_METRIC_EUCLIDEAN_H_
#define UTIL_METRIC_EUCLIDEAN_H_

#include <cassert>
#include <cmath>

namespace cl {
namespace metric {

/**
 * Euclidean distance for two n-dimensional points.
 */
class Euclidean {
public:
    Euclidean() = default;

    Euclidean(const Euclidean&) = delete;

    Euclidean& operator=(const Euclidean&) = delete;

    template <typename T>
    double operator() (const T& a, const T& b) const {
        int size1 = a.size();
        int size2 = b.size();
        assert(size1 == size2);

        double t = 0.0;
        for (int i = 0; i < size1; ++i) {
            t += static_cast<double>(a[i] - b[i]) * (a[i] - b[i]);
        }

        return std::sqrt(t);
    }
};

} // namespace metric
} // namespace cl

#endif // UTIL_METRIC_EUCLIDEAN_H_
