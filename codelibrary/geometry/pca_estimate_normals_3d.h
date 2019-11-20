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
                       Vector3D<T>* normal) {
    static_assert(std::is_floating_point<T>::value,
                  "T must be a floating point.");

    assert(!points.empty());
    assert(points.size() == weights.size());
    assert(normal);

    Point3D<T> centroid = Centroid(points, weights);

    T a00 = 0, a01 = 0, a02 = 0, a11 = 0, a12 = 0, a22 = 0;
    int i = 0;
    T sum = 0;
    for (const Point3D<T>& p : points) {
        T x = p.x - centroid.x;
        T y = p.y - centroid.y;
        T z = p.z - centroid.z;
        T w = weights[i++];

        a00 += w * x * x;
        a01 += w * x * y;
        a02 += w * x * z;
        a11 += w * y * y;
        a12 += w * y * z;
        a22 += w * z * z;
        sum += w;
    }
    if (sum == 0) {
        *normal = Vector3D<T>(0, 0, 1);
        return;
    }

    T t = 1 / sum;
    a00 = a00 * t;
    a01 = a01 * t;
    a02 = a02 * t;
    a11 = a11 * t;
    a12 = a12 * t;
    a22 = a22 * t;

    // Computing the least eigenvalue of the covariance matrix.
    T q = (a00 + a11 + a22) / 3;
    T pq = (a00 - q) * (a00 - q) + (a11 - q) * (a11 - q) +
            (a22 - q) * (a22 - q) + 2 * (a01 * a01 + a02 * a02 + a12 * a12);
    pq = std::sqrt(pq / 6);
    T mpq = std::pow(1 / pq, 3);
    T det_b = mpq * ((a00 - q) * ((a11 - q) * (a22 - q) - a12 * a12) -
               a01 * (a01 * (a22 - q) - a12 * a02) +
               a02 * (a01 * a12 - (a11 - q) * a02));
    T r = det_b / 2;
    T phi = 0;
    if (r <= -1)
        phi = M_PI / 3;
    else if (r >= 1)
        phi = 0;
    else
        phi = std::acos(r) / 3;
    T eig = q + 2 * pq * std::cos(phi + M_PI * (T(2) / T(3)));

    // Computing the corresponding eigenvector.
    normal->x = a01 * a12 - a02 * (a11 - eig);
    normal->y = a01 * a02 - a12 * (a00 - eig);
    normal->z = (a00 - eig) * (a11 - eig) - a01 * a01;

    // Normalize.
    double norm = normal->norm();
    if (norm == 0.0) {
        *normal = Vector3D<T>(0, 0, 1);
    } else {
        *normal *= static_cast<T>(1.0 / norm);
    }
}

/**
 * Estimates normal vector by PCA method.
 *
 * The output normal is normalize to the unit length, and its orientation is
 * randomly assigned.
 */
template <typename T>
void PCAEstimateNormal(const Array<Point3D<T>>& points, Vector3D<T>* normal) {
    static_assert(std::is_floating_point<T>::value,
                  "T must be a floating point.");

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
                        Array<Vector3D<T>>* normals) {
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
                        Array<Vector3D<T>>* normals) {
    static_assert(std::is_floating_point<T>::value,
                  "T must be a floating point.");

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
                                        Array<Vector3D<T>>* normals) {
    static_assert(std::is_floating_point<T>::value,
                  "T must be a floating point.");

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

        Vector3D<T> normal;
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
                                        Array<Vector3D<T>>* normals) {
    static_assert(std::is_floating_point<T>::value,
                  "T must be a floating point.");

    using Point = typename std::iterator_traits<Iterator>::value_type;

    KDTree<Point> kd_tree(first, last);
    OrientationAwarePCAEstimateNormals(kd_tree, k, normals);
}

} // namespace geometry
} // namespace cl

#endif // GEOMETRY_PCA_ESTIMATE_NORMALS_3D_H_
