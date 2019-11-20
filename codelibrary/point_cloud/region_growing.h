//
// Copyright 2018 Yangbin Lin. All Rights Reserved.
//
// Author: yblin@jmu.edu.cn (Yangbin Lin)
//
// This file is part of the Code Library.
//

#ifndef POINT_CLOUD_REGION_GROWING_H_
#define POINT_CLOUD_REGION_GROWING_H_

#include <cassert>

#include "codelibrary/base/array.h"
#include "codelibrary/geometry/kernel/angle.h"
#include "codelibrary/geometry/kernel/distance_3d.h"
#include "codelibrary/geometry/kernel/plane_3d.h"
#include "codelibrary/geometry/kernel/point_3d.h"
#include "codelibrary/point_cloud/estimate_local_planes.h"
#include "codelibrary/point_cloud/k_nearest_neighbors.h"
#include "codelibrary/statistics/regression/linear_least_squares_fitting.h"
#include "codelibrary/util/tree/kd_tree.h"

namespace cl {
namespace point_cloud {

/**
 * RegionGrowing finds all 3D planes present in the input cloud, and outputs a
 * vector of planes, as well as a vector of labels to determine each point
 * belong to which plane.
 *
 * Only planes with more than 'min_support_points' points are detected.
 *
 * Reference:
 *   Chauve A L, Labatut P, Pons J P. Robust piecewise-planar 3D reconstruction
 *   and completion from large-scale unstructured point data[C]. IEEE Conference
 *   on Computer Vision and Pattern Recognition (CVPR), 2010: 1261-1268.
 *
 * Parameters:
 *   kd_tree - the input points are stored in KD tree.
 *   resolution         - distance resolution.
 *   k_neighbors        - number of neighbors in growing.
 *   angle_threshold    - the maximal angle (in degree) between tangent plane
 *                        and fitting plane.
 *   min_support_points - the minimal number of points that can support a plane.
 *   labels             - determine which plane each point belongs to.
 *   planes             - the extracted planes.
 */
template <typename Point>
void RegionGrowing(const KDTree<Point>& kd_tree,
                   double resolution,
                   int k_neighbors,
                   double angle_threshold,
                   int min_support_points,
                   Array<int>* labels,
                   Array<RPlane3D>* planes) {
    assert(resolution > 0.0);
    assert(k_neighbors > 0);
    assert(angle_threshold > 0.0 && angle_threshold < 180.0);
    assert(min_support_points > 0);
    assert(labels);
    assert(planes);

    int n_points = kd_tree.size();
    labels->assign(n_points, -1);
    planes->clear();

    if (n_points < k_neighbors) return;

    // Then growing each region according to the sequence of smoothness.
    Array<RPlane3D> local_planes;
    Array<double> smoothness;
    EstimateLocalPlanes(kd_tree, resolution, &local_planes, &smoothness);
    Array<int> seq;
    IndexSort(smoothness.begin(), smoothness.end(), &seq);

    Array<Array<int>> neighbors;
    MutualKNearestNeighbors(kd_tree, k_neighbors, &neighbors);

    const Array<Point>& points = kd_tree.points();
    int n_planes = 0;
    for (int index : seq) {
        if ((*labels)[index] != -1) continue;

        RPlane3D plane = local_planes[index];

        double best = DBL_MAX;
        Array<int> queue;

        while (true) {
            int front = 0;
            queue.push_back(index);
            (*labels)[index] = n_planes;

            while (front < queue.size()) {
                int cur = queue[front++];
                for (int j : neighbors[cur]) {
                    if ((*labels)[j] != -1) continue;

                    double angle = Degree(local_planes[j].normal(),
                                          plane.normal());
                    if (angle > angle_threshold &&
                        angle < 180.0 - angle_threshold) continue;

                    if (Distance(points[j], plane) <= resolution) {
                        (*labels)[j] = (*labels)[cur];
                        queue.push_back(j);
                    }
                }
            }

            // If the cluster is too small, drop it.
            if (queue.size() < 3) {
                for (int v : queue) {
                    (*labels)[v] = -1;
                }
                queue.clear();
                break;
            }

            Array<Point> ps(queue.size());
            for (int i = 0; i < queue.size(); ++i) {
                ps[i] = points[queue[i]];
            }

            double metrics = statistics::LinearLeastSquaresFitting(ps, &plane);
            if (metrics + DBL_EPSILON < best) {
                best = metrics;
            } else {
                if (queue.size() < min_support_points) {
                    for (int v : queue) {
                        (*labels)[v] = -1;
                    }
                    queue.clear();
                }
                break;
            }

            for (int v : queue) {
                (*labels)[v] = -1;
            }
            queue.clear();
        }

        if (!queue.empty()) {
            planes->push_back(plane);
            ++n_planes;
        }
    }

    for (int& label : *labels) {
        if (label == -2) label = -1;
    }
}

} // namespace point_cloud
} // namespace cl

#endif // POINT_CLOUD_REGION_GROWING_H_
