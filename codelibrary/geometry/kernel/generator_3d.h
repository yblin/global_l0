//
// Copyright 2019 Yangbin Lin. All Rights Reserved.
//
// Author: yblin@jmu.edu.cn (Yangbin Lin)
//
// This file is part of the Code Library.
//
// 3D geometric objects generators.
//

#ifndef GEOMETRY_KERNEL_GENERATOR_3D_H_
#define GEOMETRY_KERNEL_GENERATOR_3D_H_

#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#include <math.h>
#else
#include <cmath>
#endif // _USE_MATH_DEFINES

#include <memory>
#include <random>

#include "codelibrary/geometry/kernel/box_3d.h"
#include "codelibrary/geometry/kernel/point_3d.h"
#include "codelibrary/geometry/kernel/sphere_3d.h"
#include "codelibrary/geometry/kernel/triangle_3d.h"
#include "codelibrary/geometry/kernel/vector_3d.h"

namespace cl {
namespace geometry {

/**
 * Randomly generate real value points in the 3D box.
 */
template <typename T>
class RandomPointInBox3D {
    static_assert(std::is_floating_point<T>::value);

public:
    explicit RandomPointInBox3D(const Box3D<T>& box)
        : uniform_x_(box.x_min(), box.x_max()),
          uniform_y_(box.y_min(), box.y_max()),
          uniform_z_(box.z_min(), box.z_max()) {}

    /**
     * Generate a 2D point.
     */
    template <typename RandomEngine>
    Point3D<T> operator()(RandomEngine* engine) {
        return Point3D<T>(uniform_x_(*engine), uniform_y_(*engine),
                          uniform_z_(*engine));
    }

private:
    // For generating uniform 3D points.
    std::uniform_real_distribution<T> uniform_x_, uniform_y_, uniform_z_;
};

/**
 * Randomly generate integer value points in the 3D box.
 */
template <>
class RandomPointInBox3D<int> {
public:
    explicit RandomPointInBox3D(const IBox3D& box)
        : uniform_x_(box.x_min(), box.x_max()),
          uniform_y_(box.y_min(), box.y_max()),
          uniform_z_(box.z_min(), box.z_max()) {}

    /**
     * Generate a 3D point.
     */
    template <typename RandomEngine>
    IPoint3D operator()(RandomEngine* engine) {
        return IPoint3D(uniform_x_(*engine), uniform_y_(*engine), 
                        uniform_z_(*engine));
    }

private:
    // For generating uniform 3D points.
    std::uniform_int_distribution<int> uniform_x_, uniform_y_, uniform_z_;
};

/**
 * Randomly generate 3D points on the sphere.
 * Note that the generated points just almost on the sphere, not exact.
 */
template <typename T>
class RandomPointOnSphere3D {
    static_assert(std::is_floating_point<T>::value);

public:
    explicit RandomPointOnSphere3D(const Sphere3D<T>& sphere)
        : radius_(sphere.radius()),
          center_x_(sphere.center().x),
          center_y_(sphere.center().y),
          center_z_(sphere.center().z),
          uniform_(0, 1) {}

    /**
     * Random generate point on the given sphere.
     */
    template <typename RandomEngine>
    Point3D<T> operator()(RandomEngine* engine) {
        T u = uniform_(*engine);
        T v = uniform_(*engine);
        T theta = u * T(2) * M_PI;
        T phi = std::acos(v * 2 - 1);

        T x = std::cos(theta) * std::sin(phi) * radius_ + center_x_;
        T y = std::sin(theta) * std::sin(phi) * radius_ + center_y_;
        T z = std::cos(phi) * radius_ + center_z_;
        return Point3D<T>(x, y, z);
    }

private:
    T radius_;                         // Radius of sphere.
    T center_x_, center_y_, center_z_; // Center of sphere.
    std::uniform_real_distribution<T> uniform_;
};

/**
 * Randomly generate 3D point in the 3D triangle.
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
class RandomPointInTriangle3D {
    static_assert(std::is_floating_point<T>::value);

public:
    explicit RandomPointInTriangle3D(const Triangle3D<T>& triangle)
        : x1_(triangle.p1().x),
          x2_(triangle.p2().x),
          x3_(triangle.p3().x),
          y1_(triangle.p1().y),
          y2_(triangle.p2().y),
          y3_(triangle.p3().y),
          z1_(triangle.p1().z),
          z2_(triangle.p2().z),
          z3_(triangle.p3().z),
          uniform_(0, 1) {}

    /**
     * Random generate point on the given triangle.
     */
    template <typename RandomEngine>
    Point3D<T> operator()(RandomEngine* engine) {
        T r1 = uniform_(*engine);
        T r2 = uniform_(*engine);
        T r1_sqrt = std::sqrt(r1);
        T t1 = 1 - r1_sqrt;
        T t2 = r1_sqrt * (1 - r2);
        T t3 = r2 * r1_sqrt;

        return Point3D<T>(t1 * x1_ + t2 * x2_ + t3 * x3_,
                          t1 * y1_ + t2 * y2_ + t3 * y3_,
                          t1 * z1_ + t2 * z2_ + t3 * z3_);
    }

private:
    T x1_, x2_, x3_, y1_, y2_, y3_, z1_, z2_, z3_;
    std::uniform_real_distribution<T> uniform_;
};

} // namespace geometry
} // namespace cl

#endif // GEOMETRY_KERNEL_GENERATOR_3D_H_
