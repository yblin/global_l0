//
// Copyright 2014 Yangbin Lin. All Rights Reserved.
//
// Author: yblin@jmu.edu.cn (Yangbin Lin)
//
// This file is part of the Code Library.
//

#ifndef GEOMETRY_PCA_ESTIMATE_NORMALS_3D_H_
#define GEOMETRY_PCA_ESTIMATE_NORMALS_3D_H_

#include <cassert>

#include "codelibrary/base/units.h"
#include "codelibrary/geometry/kernel/center_3d.h"
#include "codelibrary/geometry/kernel/point_3d.h"
#include "codelibrary/statistics/regression/principal_component_analysis_3d.h"
#include "codelibrary/util/tree/kd_tree.h"

namespace cl {
namespace geometry {

/**
 * Estimates normal direction by weighted PCA method.
 *
 * The output normal is normalize to the unit length. Its orientation is
 * randomly assigned.
 *
 * Note that, we only need to compute the least eigenvector of the covariance
 * matrix.
 */
template <typename T>
void PCAEstimateNormal(const Array<Point3D<T>>& points,
                       const Array<T>& weights,
                       RVector3D* normal) {
    static_assert(std::is_floating_point<T>::value, "");

    statistics::PrincipalComponentAnalysis3D pca(points, weights);
    *normal = pca.eigenvector(0);
}

/**
 * Estimates normal vector by PCA method.
 *
 * The output normal is normalize to the unit length, and its orientation is
 * randomly assigned.
 */
template <typename T>
void PCAEstimateNormal(const Array<Point3D<T>>& points, RVector3D* normal) {
    static_assert(std::is_floating_point<T>::value, "");

    Array<T> weights(points.size(), 1);
    PCAEstimateNormal(points, weights, normal);
}

/**
 * Estimates normal directions of the [first, last) range of 2D/3D points PCA
 * method over the k nearest neighbors.
 *
 * The output normals are normalize to the unit length, and their orientation
 * are randomly assigned.
 *
 *   kd_tree  - the input points are stored in the KD tree.
 *   k        - used to define the k-nearest neighbors.
 *   normals  - the output normals.
 */
template <typename T>
void PCAEstimateNormals(const KDTree<Point3D<T>>& kd_tree, int k,
                        Array<RVector3D>* normals) {
    assert(!kd_tree.empty());
    assert(k > 0);
    assert(normals);

    int n = kd_tree.size();
    k = std::min(k, n);

    normals->resize(n);

    const Array<Point3D<T>>& points = kd_tree.points();
    Array<Point3D<T>> neighbors;

    for (int i = 0; i < n; ++i) {
        kd_tree.FindKNearestNeighbors(points[i], k, &neighbors);
        PCAEstimateNormal(neighbors, &(*normals)[i]);
    }
}

/**
 * Estimate a set of normal vectors.
 */
template <typename Iterator, typename T>
void PCAEstimateNormals(Iterator first, Iterator last, int k,
                        Array<RVector3D>* normals) {
    static_assert(std::is_floating_point<T>::value, "");

    using Point = typename std::iterator_traits<Iterator>::value_type;
    KDTree<Point> kd_tree(first, last);
    PCAEstimateNormals(kd_tree, k, normals);
}

/**
 * Orientation aware PCA normal estimation.
 *
 * This function will re-orient the normal vector for each point according to
 * the neighbors whose normal vector has the same orientation.
 *
 * Parameters:
 *   kd_tree - the input points are stored in the KD tree.
 *   k       - used to define the k-nearest neighbors.
 *   normals - the output normals.
 */
template <typename T>
void OrientationAwarePCAEstimateNormals(const KDTree<Point3D<T>>& kd_tree,
                                        int k,
                                        Array<RVector3D>* normals) {
    static_assert(std::is_floating_point<T>::value, "");

    assert(!kd_tree.empty());
    assert(k > 0);
    assert(normals);
    assert(normals->size() == kd_tree.size());

    int n = kd_tree.size();
    k = std::min(k, n);

    const Array<Point3D<T>>& points = kd_tree.points();
    Array<int> neighbors;

    for (int i = 0; i < n; ++i) {
        kd_tree.FindKNearestNeighbors(points[i], k, &neighbors);

        Array<Point3D<T>> neighbor_points;
        neighbor_points.reserve(k);
        for (int j = 0; j < k; ++j) {
            if ((*normals)[i] * (*normals)[neighbors[j]] >= 0.0) {
                neighbor_points.push_back(points[neighbors[j]]);
            }
        }

        RVector3D normal;
        PCAEstimateNormal(neighbor_points, &normal);
        if (normal * (*normals)[i] < 0) {
            (*normals)[i] = -normal;
        } else {
            (*normals)[i] =  normal;
        }
    }
}

/**
 * Similar to the previous one.
 */
template <typename Iterator, typename T>
void OrientationAwarePCAEstimateNormals(Iterator first, Iterator last, int k,
                                        Array<RVector3D>* normals) {
    static_assert(std::is_floating_point<T>::value, "");

    using Point = typename std::iterator_traits<Iterator>::value_type;

    KDTree<Point> kd_tree(first, last);
    OrientationAwarePCAEstimateNormals(kd_tree, k, normals);
}

} // namespace geometry
} // namespace cl

#endif // GEOMETRY_PCA_ESTIMATE_NORMALS_3D_H_
