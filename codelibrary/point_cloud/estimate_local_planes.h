//
// Copyright 2018 Yangbin Lin. All Rights Reserved.
//
// Author: yblin@jmu.edu.cn (Yangbin Lin)
//
// This file is part of the Code Library.
//

#ifndef POINT_CLOUD_ESTIMATE_LOCAL_PLANE_H_
#define POINT_CLOUD_ESTIMATE_LOCAL_PLANE_H_

#include <cassert>

#include "codelibrary/base/array.h"
#include "codelibrary/geometry/kernel/angle.h"
#include "codelibrary/geometry/kernel/distance_3d.h"
#include "codelibrary/geometry/kernel/plane_3d.h"
#include "codelibrary/geometry/kernel/point_3d.h"
#include "codelibrary/statistics/regression/linear_least_squares_fitting.h"
#include "codelibrary/util/tree/kd_tree.h"

namespace cl {
namespace point_cloud {

/**
 * Estimate a local plane for each point. It optional outputs a fitting score
 * for each plane.
 *
 * Reference:
 *   Chauve A L, Labatut P, Pons J P. Robust piecewise-planar 3D reconstruction
 *   and completion from large-scale unstructured point data[C]. IEEE Conference
 *   on Computer Vision and Pattern Recognition (CVPR), 2010: 1261-1268.
 *
 * Parameters:
 *  kd_tree    - the input points are stored in KD tree.
 *  resolution - resolution to determine a local plane.
 *  planes     - the estimated local plane for each point.
 *  smoothness - the fitting quality of each plane.
 */
template <typename Point>
void EstimateLocalPlanes(const KDTree<Point>& kd_tree,
                         double resolution,
                         Array<RPlane3D>* local_planes,
                         Array<double>* smoothness = nullptr) {
    assert(resolution > 0.0);
    assert(local_planes);

    int n_points = kd_tree.size();

    local_planes->resize(n_points);
    if (smoothness) smoothness->assign(n_points, 1.0);

    const int max_iters = 5;

    Array<Point> neighbor_points;
    double radius = 4.0 * resolution * resolution;

    const Array<Point>& points = kd_tree.points();
    for (int i = 0; i < n_points; ++i) {
        kd_tree.FindRadiusNeighbors(points[i], radius, &neighbor_points);

        int iter = 0;
        do {
            double t = statistics::LinearLeastSquaresFitting(
                        neighbor_points, &(*local_planes)[i]);
            if (smoothness) (*smoothness)[i] = t;

            // Delete the points whose distance to the tangent plane is
            // smaller than resolution / 2.
            int size = 0;
            for (const Point& p : neighbor_points) {
                if (Distance(p, (*local_planes)[i]) <= resolution * 0.5) {
                    neighbor_points[size++] = p;
                }
            }
            if (size == neighbor_points.size()) break;
            if (size < 3) {
                if (smoothness) (*smoothness)[i] = 1.0;
                break;
            }
            neighbor_points.resize(size);
        } while (iter++ < max_iters);
    }
}

} // namespace point_cloud
} // namespace cl

#endif // POINT_CLOUD_ESTIMATE_LOCAL_PLANE_H_
