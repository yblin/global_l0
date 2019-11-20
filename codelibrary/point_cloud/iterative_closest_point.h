//
// Copyright 2015 Yangbin Lin. All Rights Reserved.
//
// Author: yblin@jmu.edu.cn (Yangbin Lin)
//
// This file is part of the Code Library.
//

#ifndef GEOMETRY_POINT_CLOUD_ITERATIVE_CLOSEST_POINT_H_
#define GEOMETRY_POINT_CLOUD_ITERATIVE_CLOSEST_POINT_H_

#include <cassert>
#include <cfloat>
#include <cmath>
#include <memory>
#include <utility>

#include "codelibrary/base/algorithm.h"
#include "codelibrary/base/array.h"
#include "codelibrary/base/equal.h"
#include "codelibrary/geometry/kernel/distance_3d.h"
#include "codelibrary/geometry/kernel/quaternion.h"
#include "codelibrary/geometry/kernel/transform_3d.h"
#include "codelibrary/util/tree/kd_tree.h"

namespace cl {
namespace geometry {
namespace point_cloud {

/**
 * Iterative Closest Point (ICP) method is an algorithm to estimate the rigid
 * transformation between two point clouds. The algorithm iteratively revises
 * the transformation needed to minimize the distance from the source to the
 * reference point cloud.
 *
 * Sample usage:
 *
 *   typedef IterativeClosestPoint ICP;
 *
 *   // Define callback.
 *   class ICPCallback : public ICP::Callback {
 *   public:
 *       ICPCallback() {
 *           printf("ICP begin\n");
 *       }
 *
 *       void operator() (const ICP::Result& result) {
 *           printf("iteration: %d, correspondences: %d, mse: %lf\n",
 *                  result.iteration,
 *                  result.correspondences.size(),
 *                  result.mean_square_error);
 *       }
 *   };
 *
 *   ICP icp(target.begin(), target.end());
 *   ICPCallback callback;
 *   icp.set_callback(std::make_shared<ICPCallback>(&callback));
 *   ICP::Result result;
 *   icp.Align(source, 10, 0.1, &result);
 *
 * Reference:
 *   Mckay N D, Besl P J. A Method for Registration of 3-D Shapes [J].
 *   IEEE Transactions on Pattern Analysis & Machine Intelligence, 1992,
 *   14(11):239-256.
 */
class IterativeClosestPoint {
public:
    /**
     * Result of ICP.
     */
    struct Result {
        // Current iteration number.
        int iteration;

        // Rotation from source point cloud to target point cloud.
        Quaternion rotation;

        // Translation from source point cloud to target point cloud.
        RVector3D translation;

        // Current aligned points.
        Array<RPoint3D> aligned_points;

        // Correspondences between two point clouds.
        Array<std::pair<int, int> > correspondences;

        // The mean square error between correspondences points.
        double mean_square_error;
    };

    /**
     * Callback for Iteration results.
     */
    class Callback {
    public:
        virtual ~Callback() = default;

        virtual void operator() (const Result& result) = 0;
    };

    /**
     * Construct ICP by the target point cloud.
     */
    template <typename Iterator>
    IterativeClosestPoint(Iterator first, Iterator last)
        : kd_tree_(first, last), callback_(nullptr) {}

    IterativeClosestPoint(const IterativeClosestPoint&) = delete;

    IterativeClosestPoint& operator =(const IterativeClosestPoint&) = delete;

    /**
     * Set callback.
     */
    void set_callback(const std::shared_ptr<Callback>& callback) {
        callback_ = callback;
    }

    /**
     * Align the given points to target points.
     *
     *  [first, last)                 - points to be aligned.
     *  max_iterations                - maximal number of iterations.
     *  max_correspondences_distances - maximum distance threshold between two
     *                                  correspondent points.
     *  result                        - result of ICP.
     */
    template <typename Iterator>
    void Align(Iterator first, Iterator last,
               int max_iterations,
               double max_correspondences_distance,
               Result* result) const {
        assert(max_iterations > 0);
        assert(max_correspondences_distance > 0.0);
        assert(result);

        int n = CountElements(first, last);
        assert(n == kd_tree_.size());

        // Initialize the result.
        result->aligned_points.reserve(n);
        for (Iterator p = first; p != last; ++p) {
            result->aligned_points.emplace_back(p->x, p->y, p->z);
        }
        result->iteration = 0;
        result->mean_square_error = DBL_MAX;
        result->rotation = Quaternion(0.0, 0.0, 0.0, 1.0);
        result->translation = RVector3D(0.0, 0.0, 0.0);

        Array<int> neighbors;
        double r = max_correspondences_distance * max_correspondences_distance;

        for (; result->iteration < max_iterations; ++result->iteration) {
            result->correspondences.clear();
            Array<RPoint3D> x, y;

            double error = 0.0;

            for (int i = 0; i < n; ++i) {
                const RPoint3D& p = result->aligned_points[i];
                kd_tree_.FindKNearestInRadiusNeighbors(p, 1, r, &neighbors);

                if (!neighbors.empty()) {
                    int j = neighbors.front();
                    const RPoint3D& q = kd_tree_.points()[j];
                    result->correspondences.emplace_back(i, j);
                    x.emplace_back(p);
                    y.emplace_back(q);
                    error += SquaredDistance(p, q);
                }
            }

            if (result->correspondences.empty()) {
                if (callback_) (*callback_)(*result);
                return;
            }

            Rigid(x, y, &(result->rotation), &(result->translation));

            for (int i = 0; i < n; ++i) {
                RPoint3D& p = result->aligned_points[i];
                p = Rotate(p, result->rotation);
                p = Translate(p, result->translation);
            }

            error /= result->correspondences.size();
            if (Equal(error, result->mean_square_error)) {
                result->mean_square_error = error;
                if (callback_) (*callback_)(*result);
                return;
            }

            result->mean_square_error = error;
            if (callback_) (*callback_)(*result);
        }
    }

private:
    // KDtree stores the target point cloud.
    KDTree<RPoint3D> kd_tree_;

    // Iteartion callback.
    std::shared_ptr<Callback> callback_;
};

} // namespace point_cloud
} // namespace geometry
} // namespace cl

#endif // GEOMETRY_POINT_CLOUD_ITERATIVE_CLOSEST_POINT_H_
