//
// Copyright 2017 Yangbin Lin. All Rights Reserved.
//
// Author: yblin@jmu.edu.cn (Yangbin Lin)
//
// This file is part of the Code Library.
//

#ifndef POINT_CLOUD_BILATERAL_FILTER_H_
#define POINT_CLOUD_BILATERAL_FILTER_H_

#include <cassert>
#include <cmath>

#include "codelibrary/base/algorithm.h"
#include "codelibrary/base/array.h"
#include "codelibrary/geometry/kernel/box_3d.h"
#include "codelibrary/geometry/kernel/distance_3d.h"
#include "codelibrary/geometry/kernel/point_3d.h"
#include "codelibrary/geometry/pca_estimate_normals_3d.h"
#include "codelibrary/point_cloud/point_octree.h"
#include "codelibrary/util/tree/octree.h"

namespace cl {
namespace point_cloud {

/**
 * Perform bilateral filter on the given points. After filtering, the input 
 * points will be replaced by the filtered points.
 *
 * The normals of filtered points are also outputted.
 *
 * Parameters:
 *  resolution - the resolution of the filter.
 *  sigma_n    - Gaussian weight for the distance to the tangent plane.
 *  points     - the input points. It will be replaced by the filtered points.
 *  normals    - the normal vectors of the filtered points.
 *
 * Reference:
 *   Julie Digne, and Carlo de Franchis, The Bilateral Filter for Point Clouds,
 *   Image Processing On Line, 7 (2017), pp. 278–287.
 */
void BilateralFilter(double resolution, double sigma_n, Array<RPoint3D>* points,
                     Array<RVector3D>* normals) {
    assert(resolution > 0.0);
    assert(sigma_n > 0.0);
    assert(points);
    assert(normals);

    double sigma_d = resolution / 3.0;
    double t1 = -0.5 / sigma_d / sigma_d;
    double t2 = -0.5 / sigma_n / sigma_n;

    int n = points->size();
    normals->resize(n);

    PointOctree octree(*points, resolution);

    Array<RPoint3D> neighbors;
    for (int i = 0; i < n; ++i) {
        RPoint3D& p = (*points)[i];
        octree.FindNeighbors(p, resolution * resolution, &neighbors);
        if (neighbors.size() < 3) continue;

        RVector3D np;
        geometry::PCAEstimateNormal(neighbors, &np);
        if ((*normals)[i] * np < 0.0) {
            np = -np;
        }
        (*normals)[i] = np;

        double sum_w = 0.0;
        double delta_p = 0.0;
        for (const RPoint3D& q : neighbors) {
            double dis = SquaredDistance(p, q);
            double dn = (q - p) * np;
            double w = std::exp(t1 * dis) * std::exp(t2 * dn * dn);
            delta_p += w * dn;
            sum_w += w;
        }

        if (sum_w > 0.0) p += delta_p / sum_w * np;
    }
}

} // namespace point_cloud
} // namespace cl

#endif // POINT_CLOUD_BILATERAL_FILTER_H_
