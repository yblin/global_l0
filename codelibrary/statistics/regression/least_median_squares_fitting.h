//
// Copyright 2013 Yangbin Lin. All Rights Reserved.
//
// Author: yblin@jmu.edu.cn (Yangbin Lin)
//
// This file is part of the Code Library.
//

#ifndef STATISTICS_REGRESSION_LEAST_MEDIAN_SQUARES_FITTING_H_
#define STATISTICS_REGRESSION_LEAST_MEDIAN_SQUARES_FITTING_H_

#include <cassert>
#include <cfloat>
#include <random>

#include "codelibrary/base/algorithm.h"
#include "codelibrary/geometry/kernel/distance_3d.h"
#include "codelibrary/statistics/kernel/median.h"

namespace cl {
namespace statistics {

/**
 * Least median square fitting for 3D plane of 3D points [first, last).
 *
 * Return the median distance from points to the plane.
 */
template <typename Iterator>
double LeastMedianSquaresFitting(Iterator first, Iterator last,
                                 RPlane3D* plane) {
    assert(first != last);
    assert(plane);

    int size = CountElements(first, last);

    // A default horizontal plane that goes through the first input point.
    *plane = RPlane3D(*first, RVector3D(0.0, 0.0, 1.0));

    // The iterations number m can be determined by following equation:
    //   P = 1 - (1 - e^k)^m,
    // where e is the fraction of outliers, k is the least number of points to
    // define the model and P is the probability that at least one of the m
    // sub-samples is good.
    // Impose e = 0.5, k = 3 and m = 50, we have P = 99.87%.
    const int n_iterations = 50;

    double min_distance = DBL_MAX;

    std::mt19937 random_engine;
    std::uniform_int_distribution<int> uniform(0, size - 1);

    RPoint3D pa, pb, pc;
    for (int itr = 0; itr < n_iterations; ++itr) {
        int a = uniform(random_engine);
        int b = uniform(random_engine);
        int c = uniform(random_engine);
        if (a == b || b == c || a == c) continue;

        for (int j = 0; j < 3; ++j) {
            pa[j] = first[a][j];
            pb[j] = first[b][j];
            pc[j] = first[c][j];
        }

        if (pa == pb || pa == pc || pb == pc) continue;

        RPlane3D plane_of_3_points(pa, pb, pc);

        Array<double> distances(size);
        for (int i = 0; i < size; ++i) {
            distances[i] = Distance(first[i], plane_of_3_points);
        }

        double distance = Median(distances.begin(), distances.end());
        if (distance < min_distance) {
            min_distance = distance;
            *plane = plane_of_3_points;
        }
    }

    return min_distance;
}

/**
 * Least median square fitting for 3D line.
 *
 * Return the median distance from points to the line.
 */
template <typename Iterator>
static double LeastMedianSquaresFitting(Iterator first, Iterator last,
                                        RLine3D* line) {
    assert(line);
    assert(first != last);

    // A default line along X axis.
    *line = RLine3D(*first, RVector3D(1.0, 0.0, 0.0));

    // The iterations number m can be determined by following equation:
    //   P = 1 - (1 - e^k)^m,
    // where e is the fraction of outliers, k is the least number of points to
    // define the model and P is the probability that at least one of the m
    // sub-samples is good.
    // Impose e = 0.5, k = 2 and m = 30, we have P = 99.98%.
    const int n_iterations = 30;

    double min_distance = DBL_MAX;

    int size = CountElements(first, last);
    std::mt19937 random_engine;
    std::uniform_int_distribution<int> uniform(0, size - 1);

    for (int itr = 0; itr < n_iterations; ++itr) {
        const RPoint3D& pa = first[uniform(random_engine)];
        const RPoint3D& pb = first[uniform(random_engine)];

        if (pa == pb) continue;

        RLine3D l(pa, pb);

        Array<double> distances(size);
        for (int i = 0; i < size; ++i) {
            distances[i] = Distance(first[i], l);
        }

        double distance = Median(distances.begin(), distances.end());
        if (distance < min_distance) {
            min_distance = distance;
            *line = l;
        }
    }

    return min_distance;
}

} // namespace statistics
} // namespace cl

#endif // STATISTICS_REGRESSION_LEAST_MEDIAN_SQUARES_FITTING_H_
