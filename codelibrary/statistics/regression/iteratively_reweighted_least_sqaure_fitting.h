//
// Copyright 2015 Yangbin Lin. All Rights Reserved.
//
// Author: yblin@jmu.edu.cn (Yangbin Lin)
//
// This file is part of the Code Library.
//

#ifndef STATISTICS_REGRESSION_ITERATIVELY_REWEIGHTED_LEAST_SQAURE_FITTING_H_
#define STATISTICS_REGRESSION_ITERATIVELY_REWEIGHTED_LEAST_SQAURE_FITTING_H_

#include <cassert>

#include "codelibrary/base/algorithm.h"
#include "codelibrary/geometry/util/distance_3d.h"
#include "codelibrary/statistics/regression/linear_least_squares_fitting.h"

namespace cl {
namespace statistics {

/**
 * Iteratively linear least square fitting for geometric object from input
 * points.
 *
 * Parameters:
 *  [first, last) - the input points.
 *  n_iterations  - the number of iterations (an optional value is 3).
 *  scale         - distance where point gets zero weight.
 *  object        - fitted geometric object (plane or line).
 *
 * Return:
 *  the fitting quality: 0 is best, 1 is worst.
 */
template <typename Iterator, typename FittingObject>
double IterativelyReweightedLeastSqauresFitting(Iterator first, Iterator last,
                                                int n_iterations,
                                                double scale,
                                                FittingObject* object) {
    int n = CountElements(first, last);
    assert(n > 0);
    assert(n_iterations >= 1);
    assert(scale > 0.0);
    assert(object);

    double metric = 1.0;

    int iteration = 0;
    do {
        Array<double> weights(n, 1.0); // Least squares weights.

        if (iteration > 0) {
            // Compute the weight for each point according to its distance to
            // the current fitted plane.
            int i = 0;
            for (Iterator p = first; p != last; ++p, ++i) {
                weights[i] = geometry::Distance(*p, *object);
            }

            // True if all the distances from point to object is larger than
            // scale.
            bool no_fitted_points = true;

            // The farther away, the smaller the weight:
            //   w = (1 - (distance / scale)^2)^2.
            for (double& w : weights) {
                if (w < scale) {
                    w /= scale;
                    w = 1.0 - w * w;
                    w *= w;
                    no_fitted_points = false;
                } else {
                    w = 0.0;
                }
            }

            if (no_fitted_points) return 1.0;
        }

        metric = LinearLeastSquaresFitting(first, last, weights, object);
    } while (++iteration < n_iterations);

    return metric;
}

} // namespace statistics
} // namespace cl

#endif // STATISTICS_REGRESSION_ITERATIVELY_REWEIGHTED_LEAST_SQAURE_FITTING_H_

