//
// Copyright 2014 Yangbin Lin. All Rights Reserved.
//
// Author: yblin@jmu.edu.cn (Yangbin Lin)
//
// This file is part of the Code Library.
//

#ifndef POINT_CLOUD_MST_ORIENT_NORMALS_H_
#define POINT_CLOUD_MST_ORIENT_NORMALS_H_

#include <algorithm>
#include <cassert>
#include <functional>
#include <queue>
#include <utility>

#include "codelibrary/base/array.h"
#include "codelibrary/geometry/kernel/distance_3d.h"
#include "codelibrary/util/metric/cosine.h"
#include "codelibrary/util/metric/squared_euclidean.h"
#include "codelibrary/util/tree/kd_tree.h"

namespace cl {
namespace point_cloud {

namespace internal {

/**
 * This class is used for MST-based normal orientation, it computes the distance
 * between two points with normals. This distance is the weight of MST edge.
 *
 * Reference:
 *   Huang H, Li D, Zhang H, et al. Consolidation of unorganized point clouds
 *   for surface reconstruction[J]. ACM Transactions on Graphics (TOG). ACM,
 *   2009, 28(5): 176.
 */
class NormalDistance {
public:
    template <typename Point>
    double operator() (const Point& p1, const RVector3D& v1,
                       const Point& p2, const RVector3D& v2) {
        if (p1 == p2) return 0.0;

        RPoint3D p11(p1.x + v1.x, p1.y + v1.y, p1.z + v1.z);
        RPoint3D p12(p1.x - v1.x, p1.y - v1.y, p1.z - v1.z);
        RPoint3D p21(p2.x + v2.x, p2.y + v2.y, p2.z + v2.z);
        RPoint3D p22(p2.x - v2.x, p2.y - v2.y, p2.z - v2.z);

        RPoint3D m11 = p11 + 0.5 * (p21 - p11);
        RPoint3D m12 = p11 + 0.5 * (p22 - p11);
        RPoint3D m21 = p12 + 0.5 * (p21 - p12);
        RPoint3D m22 = p12 + 0.5 * (p22 - p12);

        RLine3D line(p1, p2);
        double d11 = Distance(m11, line);
        double d12 = Distance(m12, line);
        double d21 = Distance(m21, line);
        double d22 = Distance(m22, line);

        double max_d = std::max(std::max(d11, d12), std::max(d21, d22));

        return 1.0 - std::fabs(v1 * v2) * max_d / (1.0 + Distance(p1, p2));
    }
};

} // namespace internal

/**
 * Orient the normals of the input points using the propagation of a seed
 * orientation through a minimum spanning tree of the K-Nearest Neighbors graph.
 *
 * The input point cloud is stored in the KD tree.
 *
 * Note that, it requires the normals direction are already computed.
 *
 * Reference:
 * [1] Hoppe, H., DeRose, T., Duchamp, T., et al. Surface reconstruction from
 *     unorganized points[J]. ACM Transactions on Graphics (TOG). ACM,
 *     1992, 26(2): 71-78.
 *
 * [2] Huang H, Li D, Zhang H, et al. Consolidation of unorganized point clouds
 *     for surface reconstruction[J]. ACM Transactions on Graphics (TOG). ACM,
 *     2009, 28(5): 176.
 */
template <typename NormalDistance = internal::NormalDistance>
void MSTOrientNormals(const KDTree<RPoint3D>& kd_tree, int k,
                      Array<RVector3D>* normals) {
    assert(!kd_tree.empty());
    assert(k > 0);
    assert(normals);
    assert(kd_tree.size() == normals->size());

    int n = kd_tree.size();
    k = std::min(k, n);

    const Array<RPoint3D>& points = kd_tree.points();

    // Store neighbors.
    Array<int> neighbors;

    // Heap for minimum spanning tree.
    using Edge = std::pair<double, std::pair<int, int> >;
    std::priority_queue<Edge, std::vector<Edge>, std::greater<Edge> > heap;

    // Visited label.
    Array<bool> is_visit(n, false);

    // Get the sequence of the input points according to their Z values. 
    Array<int> seq;
    auto compare = [&](const RPoint3D& a, const RPoint3D& b) {
        return a.z > b.z;
    };
    IndexSort(points.begin(), points.end(), compare, &seq);

    // For cosine distance metric.
    metric::Cosine cosine;

    for (int s : seq) {
        if (is_visit[s]) continue;

        // Orient the normal of the point with maximum Z towards +Z axis.
        int start = s;
        (*normals)[start] = RVector3D(0.0, 0.0, 1.0);

        heap.push(Edge(0, std::pair<int, int>(start, start)));

        while (!heap.empty()) {
            int source = heap.top().second.first;
            int target = heap.top().second.second;
            heap.pop();

            if (is_visit[target]) continue;
            is_visit[target] = true;

            if ((*normals)[source] * (*normals)[target] < 0.0) {
                (*normals)[target] = -(*normals)[target];
            }

            RVector3D v(points[target].x - points[source].x,
                        points[target].y - points[source].y,
                        points[target].z - points[source].z);
            if (std::fabs(cosine(v, (*normals)[source])) > 0.8 &&
                std::fabs(cosine(v, (*normals)[target])) > 0.8) {
                if ((*normals)[source] * (*normals)[target] > 0.0) {
                    (*normals)[target] = -(*normals)[target];
                }
            }

            // Add the adjacent KNN edges into heap.
            kd_tree.FindKNearestNeighbors(points[target], k, &neighbors);
            for (int neighbor : neighbors) {
                if (is_visit[neighbor]) continue;

                // Compute the distance measure.
                double dis = NormalDistance()(
                                points[target], (*normals)[target],
                                points[neighbor], (*normals)[neighbor]);
                heap.push(Edge(dis, std::pair<int, int>(target, neighbor)));
            }
        }
    }
}

/**
 * Similar to the previous one.
 */
template <typename Iterator, typename NormalDistance = internal::NormalDistance>
void MSTOrientNormals(Iterator first, Iterator last, int k,
                      Array<RVector3D>* normals) {
    using Point = typename std::iterator_traits<Iterator>::value_type;

    KDTree<Point> kd_tree(first, last);
    MSTOrientNormals<Point, NormalDistance>(kd_tree, k, normals);
}

} // namespace point_cloud
} // namespace cl

#endif // POINT_CLOUD_MST_ORIENT_NORMALS_H_
