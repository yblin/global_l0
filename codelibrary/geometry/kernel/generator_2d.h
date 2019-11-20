//
// Copyright 2019 Yangbin Lin. All Rights Reserved.
//
// Author: yblin@jmu.edu.cn (Yangbin Lin)
//
// This file is part of the Code Library.
//
// 2D geometric objects generators.
//

#ifndef GEOMETRY_KERNEL_GENERATOR_2D_H_
#define GEOMETRY_KERNEL_GENERATOR_2D_H_

#include <algorithm>

#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#include <math.h>
#else
#include <cmath>
#endif // _USE_MATH_DEFINES

#include <queue>
#include <random>

#include "codelibrary/base/array.h"
#include "codelibrary/base/array_nd.h"
#include "codelibrary/base/units.h"
#include "codelibrary/geometry/kernel/circle_2d.h"
#include "codelibrary/geometry/kernel/distance_2d.h"
#include "codelibrary/geometry/kernel/intersect_2d.h"
#include "codelibrary/geometry/kernel/segment_2d.h"
#include "codelibrary/geometry/kernel/triangle_2d.h"
#include "codelibrary/geometry/kernel/vector_2d.h"

namespace cl {
namespace geometry {

/**
 * Randomly generate 2D points in the box.
 */
template <typename T>
class RandomPointInBox2D {
    static_assert(std::is_floating_point<T>::value);

public:
    explicit RandomPointInBox2D(const Box2D<T>& box)
        : uniform_x_(box.x_min(), box.x_max()),
          uniform_y_(box.y_min(), box.y_max()) {}

    /**
     * Generate a 2D point.
     */
    template <typename RandomEngine>
    Point2D<T> operator()(RandomEngine* engine) {
        return Point2D<T>(uniform_x_(*engine), uniform_y_(*engine));
    }

private:
    // For generating uniform 2D points.
    std::uniform_real_distribution<T> uniform_x_, uniform_y_;
};

/**
 * Randomly generate integer value points in the 2D box.
 */
template <>
class RandomPointInBox2D<int> {
public:
    explicit RandomPointInBox2D(const IBox2D& box)
        : uniform_x_(box.x_min(), box.x_max()),
          uniform_y_(box.y_min(), box.y_max()) {}

    /**
     * Generate a 2D point.
     */
    template <typename RandomEngine>
    IPoint2D operator()(RandomEngine* engine) {
        return IPoint2D(uniform_x_(*engine), uniform_y_(*engine));
    }

private:
    // For generating uniform 2D points.
    std::uniform_int_distribution<int> uniform_x_, uniform_y_;
};

/**
 * Randomly generate real point in 2D circle.
 */
template <typename T>
class RandomPointInCircle2D {
    static_assert(std::is_floating_point<T>::value);

public:
    explicit RandomPointInCircle2D(const Circle2D<T>& circle)
        : uniform_angle_(0, T(2) * M_PI),
          uniform_r_(0, T(1)),
          radius_(circle.radius()),
          center_x_(circle.center().x),
          center_y_(circle.center().y) {}

    /**
     * Generate a 2D point in the circle.
     */
    template <typename RandomEngine>
    Point2D<T> operator()(RandomEngine* engine) {
        T angle = uniform_angle_(*engine);
        T sqrt_r = std::sqrt(uniform_r_(*engine));
        T x = sqrt_r * std::cos(angle) * radius_;
        T y = sqrt_r * std::sin(angle) * radius_;

        return Point2D<T>(x + center_x_, y + center_y_);
    }

private:
    // For generating uniform angles in range [0, 2PI).
    std::uniform_real_distribution<T> uniform_angle_;

    // For generating uniform value in range [0, 1].
    std::uniform_real_distribution<T> uniform_r_;

    // The circle that is used to generate the random points in it.
    T radius_, center_x_, center_y_;
};

/**
 * Randomly generate 2D points on 2D circle.
 * Note that, the generated points just almost on the circle, not exact.
 */
template <typename T>
class RandomPointOnCircle2D {
    static_assert(std::is_floating_point<T>::value);

public:
    explicit RandomPointOnCircle2D(const Circle2D<T>& circle)
        : uniform_angle_(0, T(2) * M_PI),
          radius_(circle.radius()),
          center_x_(circle.center().x),
          center_y_(circle.center().y) {}

    /**
     * Generate a 2D point.
     */
    template <typename RandomEngine>
    Point2D<T> operator()(RandomEngine* engine) {
        // Generate the random point in the unite circle.
        T angle = uniform_angle_(*engine);
        T x = std::cos(angle) * radius_ + center_x_;
        T y = std::sin(angle) * radius_ + center_y_;
        return RPoint2D(x, y);
    }

private:
    // For generating uniform angles in range [0, 2PI).
    std::uniform_real_distribution<T> uniform_angle_;

    // The circle that is used to generate the random points in it.
    T radius_, center_x_, center_y_;
};

/**
 * Randomly generate 2D point in the 2D triangle.
 *
 * According to [1], a random point, P, uniformly from within triangle ABC,
 * can be computed by the following convex combination of the vertices:
 *
 *   P = (1 - \sqrt(r1)) A + (\sqrt(r1)(1-r2)) B + (r2 \sqrt(r1)) C
 *
 * where r1, r2 ~ U[0, 1]
 *
 * Note that the generated points just almost on the triangle, not exact.
 *
 * Reference:
 * [1] Osada R, Funkhouser T, Chazelle B, et al. Shape distributions[J]. ACM
 *     Transactions on Graphics, 2002, 21(4):807-832.
 */
template <typename T>
class RandomPointInTriangle2D {
    static_assert(std::is_floating_point<T>::value);

public:
    explicit RandomPointInTriangle2D(const Triangle2D<T>& triangle)
        : x1_(triangle.p1().x),
          x2_(triangle.p2().x),
          x3_(triangle.p3().x),
          y1_(triangle.p1().y),
          y2_(triangle.p2().y),
          y3_(triangle.p3().y),
          uniform_(0, 1) {}

    /**
     * Random generate point on the given triangle.
     */
    template <typename RandomEngine>
    Point2D<T> operator()(RandomEngine* engine) {
        T r1 = uniform_(*engine);
        T r2 = uniform_(*engine);
        T r1_sqrt = std::sqrt(r1);
        T t1 = 1 - r1_sqrt;
        T t2 = r1_sqrt * (1 - r2);
        T t3 = r2 * r1_sqrt;

        return Point2D<T>(t1 * x1_ + t2 * x2_ + t3 * x3_,
                          t1 * y1_ + t2 * y2_ + t3 * y3_);
    }

private:
    T x1_, x2_, x3_, y1_, y2_, y3_;
    std::uniform_real_distribution<T> uniform_;
};

/**
 * Generate random point on a 2D line segment.
 */
template <typename T>
class RandomPointOnSegment2D {
    static_assert(std::is_floating_point<T>::value);

public:
    explicit RandomPointOnSegment2D(const Segment2D<T>& line)
        : uniform_(0, 1),
          direction_(line.direction().x, line.direction().y),
          x_(line.lower_point().x),
          y_(line.lower_point().y) {}

    /**
     * Generate a 2D point.
     */
    template <typename RandomEngine>
    Point2D<T> operator()(RandomEngine* engine) {
        T t = uniform_(*engine);
        return Point2D<T>(t * direction_.x + x_, t * direction_.y + y_);
    }

private:
    // Uniform generator.
    std::uniform_real_distribution<T> uniform_;

    // The input segment.
    Vector2D<T> direction_;
    T x_, y_;
};

/**
 * Linear time Poisson-disk sampling in a box.
 *
 * Reference:
 *   Bridson, Robert. Fast Poisson disk sampling in arbitrary dimensions[J].
 *   2007, 49(2):22.
 *
 * Parameters:
 *  box           - the sample domain.
 *  min_distance  - the minimum distance between samples.
 *  random_engine - random engine for generate random values.
 *  points        - the output sample points.
 */
template <typename T, typename RandomEngine>
void PoissonDiskSampling(const Box2D<T>& box, double min_distance,
                         RandomEngine* random_engine,
                         Array<Point2D<T>>* points) {
    static_assert(std::is_floating_point<T>::value);

    assert(min_distance > 0.0);
    assert(random_engine);
    assert(points);

    points->clear();

    // Initialize grid for storing samples and accelerating spatial searches.
    // The cell size is bounded by min_distances / sqrt(dimension).
    double cell_size = min_distance / std::sqrt(2.0);
    assert(box.x_length() / cell_size < INT_MAX);
    assert(box.y_length() / cell_size < INT_MAX);
    int size1 = static_cast<int>(box.x_length() / cell_size) + 1;
    int size2 = static_cast<int>(box.y_length() / cell_size) + 1;
    assert(size1 <= INT_MAX / size2);
    
    ArrayND<bool> grid(size1, size2);
    ArrayND<Point2D<T>> grid_point(size1, size2);

    RandomPointInBox2D<T> random(box);
    std::queue<Point2D<T>> processing_list;
    Point2D<T> p = random(random_engine);
    processing_list.push(p);
    int x = static_cast<int>((p.x - box.x_min()) / cell_size);
    int y = static_cast<int>((p.y - box.y_min()) / cell_size);
    x = Clamp(x, 0, size1 - 1);
    y = Clamp(y, 0, size2 - 1);
    grid(x, y) = true;
    grid_point(x, y) = p;

    std::uniform_real_distribution<T> uniform_angle(0, T(2) * M_PI);
    std::uniform_real_distribution<T> uniform_r(1, 2);

    const double resolution = min_distance * min_distance;

    while (!processing_list.empty()) {
        const Point2D<T> p = processing_list.front();
        processing_list.pop();
        points->push_back(p);

        // k is the limit of samples to choose before rejection.
        // Typically k = 30 (see the reference).
        const int k = 30;
        for (int i = 0; i < k; ++i) {
            // Randomly generate new point.
            T sqrt_r = std::sqrt(uniform_r(*random_engine));
            T radius = sqrt_r * min_distance;
            T angle  = uniform_angle(*random_engine);
            T x1 = p.x + sqrt_r * std::cos(angle) * radius;
            T y1 = p.y + sqrt_r * std::sin(angle) * radius;

            Point2D<T> p1(x1, y1);
            if (!Intersect(box, p1)) continue;

            int x = static_cast<int>((x1 - box.x_min()) / cell_size);
            int y = static_cast<int>((y1 - box.y_min()) / cell_size);
            x = Clamp(x, 0, size1 - 1);
            y = Clamp(y, 0, size2 - 1);
            if (grid(x, y)) continue;

            // Check if there is a point whose distance to p1 is smaller then
            // min_distance.
            bool check_neighbors = true;
            for (int k1 = -2; k1 <= 2; ++k1) {
                int x2 = x + k1;
                if (x2 < 0 || x2 >= size1) continue;
                for (int k2 = -2; k2 <= 2; ++k2) {
                    int y2 = y + k2;
                    if (y2 < 0 || y2 >= size2) continue;

                    // It already tested.
                    if (k1 == 0 && k2 == 0) continue;

                    // Do not need to test the corners.
                    if (std::abs(k1) + std::abs(k2) == 4) continue;

                    if (grid(x2, y2) &&
                        SquaredDistance(grid_point(x2, y2), p1) < resolution) {
                        check_neighbors = false;
                        break;
                    }
                }
                if (!check_neighbors) break;
            }

            if (check_neighbors) {
                grid(x, y) = true;
                grid_point(x, y) = p1;
                processing_list.push(p1);
            }
        }
    }
}

} // namespace geometry
} // namespace cl

#endif // GEOMETRY_KERNEL_GENERATOR_2D_H_
