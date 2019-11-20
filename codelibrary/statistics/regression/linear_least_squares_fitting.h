//
// Copyright 2017 Yangbin Lin. All Rights Reserved.
//
// Author: yblin@jmu.edu.cn (Yangbin Lin)
//
// This file is part of the Code Library.
//

#ifndef STATISTICS_REGRESSION_LINEAR_LEAST_SQUARES_FITTING_H_
#define STATISTICS_REGRESSION_LINEAR_LEAST_SQUARES_FITTING_H_

#include <algorithm>
#include <cassert>

#include "codelibrary/base/algorithm.h"
#include "codelibrary/geometry/kernel/line_2d.h"
#include "codelibrary/geometry/kernel/line_3d.h"
#include "codelibrary/geometry/kernel/plane_3d.h"
#include "codelibrary/geometry/kernel/transform_2d.h"
#include "codelibrary/geometry/kernel/transform_3d.h"
#include "codelibrary/math/matrix/matrix.h"
#include "codelibrary/statistics/regression/principal_component_analysis_2d.h"
#include "codelibrary/statistics/regression/principal_component_analysis_3d.h"

namespace cl {
namespace statistics {

/**
 * Linear Least Square Fitting for 2D line.
 *
 * Return the fitting quality: 0 is best, 1 is worst.
 */
template <typename T>
double LinearLeastSquaresFitting(const Array<Point2D<T>>& points, 
                                 RLine2D* line) {    
    assert(line);
    assert(!points.empty());

    if (points.size() == 1) {
        *line = RLine2D(RPoint2D(points.back().x, points.back().y), 
                        RVector2D(0.0, 1.0));
        return 1.0;
    } else if (points.size() == 2) {
        RPoint2D p1(points[0].x, points[0].y);
        RPoint2D p2(points[1].x, points[1].y);
        *line = RLine2D(p1, p2);
        return 0.0;
    } else {
        PrincipalComponentAnalysis2D pca(points);
        return pca.ToLine(line);
    }
}

/**
 * Weighted Linear Least Square Fitting for 2D line.
 *
 * Return the fitting quality: 0 is best, 1 is worst.
 */
template <typename T>
double LinearLeastSquaresFitting(const Array<Point2D<T>>& points,
                                 const Array<double>& weights,
                                 RLine2D* line) {
    assert(!points.empty());

    if (points.size() <= 2) {
        return LinearLeastSquaresFitting(points, line);
    }

    PrincipalComponentAnalysis2D pca(points, weights);
    return pca.ToLine(line);
}

/**
 * Linear Least Square Fitting for 2D line segment.
 *
 * Return the fitting quality: 0 is best, 1 is worst.
 */
template <typename T>
double LinearLeastSquaresFitting(const Array<Point2D<T>>& points,
                                 RSegment2D* seg) {
    assert(seg);
    assert(!points.empty());

    RLine2D l;
    double metric = LinearLeastSquaresFitting(points, &l);

    PointDirectionCompare2D<T> compare(l.direction());
    RPoint2D p1 = Project(*std::min_element(points.begin(), points.end(),
                                            compare), l);
    RPoint2D p2 = Project(*std::max_element(points.begin(), points.end(),
                                            compare), l);
    *seg = RSegment2D(p1, p2);

    return metric;
}

/**
 * Weighted Linear Least Square Fitting for 2D line segment.
 *
 * Return the fitting quality: 0 is best, 1 is worst.
 */
template <typename T>
double LinearLeastSquaresFitting(const Array<Point2D<T>>& points,
                                 const Array<double>& weights,
                                 RSegment2D* seg) {
    assert(seg);
    assert(!points.empty());

    RLine2D l;
    double metric = LinearLeastSquaresFitting(points, weights, &l);

    PointDirectionCompare2D<T> compare(l.direction());
    RPoint2D p1 = Project(*std::min_element(points, compare), l);
    RPoint2D p2 = Project(*std::max_element(points, compare), l);
    *seg = RSegment2D(p1, p2);

    return metric;
}

/**
 * Linear Least Square Fitting for 3D line.
 *
 * Return the fitting quality: 0 is best, 1 is worst.
 */
template <typename T>
double LinearLeastSquaresFitting(const Array<Point3D<T>>& points, 
                                 RLine3D* line) {
    assert(line);
    assert(!points.empty());

    if (points.size() == 1) {
        *line = RLine3D(RPoint3D(points[0].x, points[0].y, points[0].z),
                        RVector3D(0.0, 0.0, 1.0));
        return 1.0;
    } else if (points.size() == 2) {
        RPoint3D p1(points[0].x, points[0].y, points[0].z);
        RPoint3D p2(points[1].x, points[1].y, points[1].z);
        *line = RLine3D(p1, p2);
        return 0.0;
    } else {
        PrincipalComponentAnalysis3D pca(points);
        return pca.ToLine(line);
    }
}

/**
 * Linear Least Square Fitting for 3D segment.
 *
 * Return the fitting quality: 0 is best, 1 is worst.
 */
template <typename T>
double LinearLeastSquaresFitting(const Array<Point3D<T>>& points,
                                 RSegment3D* seg) {
    assert(seg);
    assert(!points.empty());

    RLine3D l;
    double metric = LinearLeastSquaresFitting(points, &l);

    PointDirectionCompare3D compare(l.direction());
    RPoint3D p1 = Project(*std::min_element(points.begin(), points.end(),
                                            compare), l);
    RPoint3D p2 = Project(*std::max_element(points.begin(), points.end(),
                                            compare), l);
    *seg = RSegment3D(p1, p2);

    return metric;
}

/**
 * Linear Least Square Fitting for 3D plane.
 *
 * Return the fitting quality: 0 is best, 1 is worst.
 */
template <typename T>
double LinearLeastSquaresFitting(const Array<Point3D<T>>& points,
                                 RPlane3D* plane) {
    assert(plane);
    assert(!points.empty());

    if (points.size() == 1) {
        RPoint3D p(points[0].x, points[0].y, points[0].z);
        *plane = RPlane3D(p, RVector3D(0.0, 0.0, 1.0));
        return 0.0;
    } else if (points.size() == 2) {
        RPoint3D p(points[0].x, points[0].y, points[0].z);
        *plane = RPlane3D(p, RVector3D(0.0, 0.0, 1.0));
        return 0.0;
    } else if (points.size() == 3) {
        RPoint3D p1(points[0].x, points[0].y, points[0].z);
        RPoint3D p2(points[1].x, points[1].y, points[1].z);
        RPoint3D p3(points[2].x, points[2].y, points[2].z);

        *plane = RPlane3D(p1, p2, p3);
        return 1.0;
    }

    PrincipalComponentAnalysis3D pca(points);
    return pca.ToPlane(plane);
}

/**
 * Weighted Linear Least Square Fitting for 3D plane.
 *
 * Return the fitting quality: 0 is best, 1 is worst.
 */
template <typename T>
double LinearLeastSquaresFitting(const Array<Point3D<T>>& points,
                                 const Array<double>& weights,
                                 RPlane3D* plane) {
    if (points.size() <= 3) {
        return LinearLeastSquaresFitting(points, plane);
    }

    PrincipalComponentAnalysis3D pca(points, weights);
    return pca.ToPlane(plane);
}

} // namespace statistics
} // namespace cl

#endif // STATISTICS_REGRESSION_LINEAR_LEAST_SQUARES_FITTING_H_
