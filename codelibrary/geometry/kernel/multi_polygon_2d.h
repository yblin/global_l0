//
// Copyright 2012 Yangbin Lin. All Rights Reserved.
//
// Author: yblin@jmu.edu.cn (Yangbin Lin)
//
// This file is part of the Code Library.
//

#ifndef GEOMETRY_KERNEL_MULTI_POLYGON_2D_H_
#define GEOMETRY_KERNEL_MULTI_POLYGON_2D_H_

#include <cassert>

#include "codelibrary/base/array.h"
#include "codelibrary/geometry/kernel/polygon_2d.h"

namespace cl {

/**
 * Multi-polygon is a set of polygon with boundary flag (inner or outer).
 *
 * It must uphold the following conditions:
 *   1. Every boundaries must be a simple polygon (allow coincide edges).
 *   2. All boundaries must be disjoint from others (except on vertices).
 */
template <typename T>
class MultiPolygon2D {
    using Polygon = Polygon2D<T>;

public:
    using value_type = T;

    // The boundary for multi-polygon
    struct Boundary {
        Boundary() = default;

        Boundary(const Polygon& poly, bool outer)
            : polygon(poly), is_outer(outer) {}

        Polygon polygon;      // The polygon of the boundary.
        bool is_outer = true; // True if this boundary is outer boundary.
    };

    MultiPolygon2D() = default;

    /**
     * Construct multi-polygon with given polygon as its outer boundary.
     */
    explicit MultiPolygon2D(const Polygon& polygon) {
        if (!polygon.empty()) Insert(polygon, true);
    }

    /**
     * Return the area of the polygon.
     */
    double Area() const {
        double area = 0.0;
        for (const Boundary& b : boundaries_) {
            if (b.is_outer) {
                area += b.polygon.Area();
            } else {
                area -= b.polygon.Area();
            }
        }

        return area;
    }

    /**
     * Insert a new boundary.
     */
    void Insert(const Polygon& polygon, bool is_outer) {
        assert(!polygon.empty());

        boundaries_.emplace_back(polygon, is_outer);
        bounding_box_.Join(polygon.bounding_box());
    }

    /**
     * Clear the polygon.
     */
    void clear() {
        bounding_box_ = Box2D<T>();
        boundaries_.clear();
    }

    bool empty() const {
        return boundaries_.empty();
    }

    const Array<Boundary>& boundaries() const {
        return boundaries_;
    }

    /**
     * Return the bounding box of multi-polygon.
     */
    const Box2D<T>& bounding_box() const {
        return bounding_box_;
    }

    /**
     * Return the number of boundaries.
     */
    int n_boundaries() const {
        return boundaries_.size();
    }

protected:
    Box2D<T> bounding_box_;      // The bounding box of multi-polygon.
    Array<Boundary> boundaries_; // The boundaries of multi-polygon.
};

using IMultiPolygon2D = MultiPolygon2D<int>;
using FMultiPolygon2D = MultiPolygon2D<float>;
using RMultiPolygon2D = MultiPolygon2D<double>;

} // namespace cl

#endif // GEOMETRY_KERNEL_MULTI_POLYGON_2D_H_
