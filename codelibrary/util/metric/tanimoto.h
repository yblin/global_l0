//
// Copyright 2018 Yangbin Lin. All Rights Reserved.
//
// Author: yblin@jmu.edu.cn (Yangbin Lin)
//
// This file is part of the Code Library.
//

#ifndef UTIL_METRIC_TANIMOTO_H_
#define UTIL_METRIC_TANIMOTO_H_

#include <algorithm>
#include <cassert>
#include <cmath>

namespace cl {
namespace metric {

/**
 * Tanimoto distance can be written as:
 *
 *                         AB
 * f(A, B) = 1 -  --------------------
 *                 |A|^2 + |B|^2 - AB
 */
class Tanimoto {
public:
    Tanimoto() = default;

    Tanimoto(const Tanimoto&) = delete;

    Tanimoto& operator=(const Tanimoto&) = delete;

    template <typename T>
    double operator() (const T& a, const T& b) const {
        assert(a.size() == b.size());

        double t1 = a * b;
        double t2 = a * a + b * b - t1;
        return 1.0 - t1 / t2;
    }
};

} // namespace metric
} // namespace cl

#endif // UTIL_METRIC_TANIMOTO_H_
