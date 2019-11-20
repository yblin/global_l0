//
// Copyright 2011 Yangbin Lin. All Rights Reserved.
//
// Author: yblin@jmu.edu.cn (Yangbin Lin)
//
// This file is part of the Code Library.
//

#ifndef GEOMETRY_KERNEL_POLYGON_2D_H_
#define GEOMETRY_KERNEL_POLYGON_2D_H_

#include <algorithm>
#include <cassert>

#include "codelibrary/base/array.h"
#include "codelibrary/geometry/kernel/box_2d.h"
#include "codelibrary/geometry/kernel/point_2d.h"
#include "codelibrary/geometry/kernel/predicate_2d.h"
#include "codelibrary/geometry/kernel/segment_2d.h"

namespace cl {

/**
 * 2D Simple Polygon2D.
 *
 * A polygon is a closed chain of edges. It can be used as a container of
 * ordered points.
 *
 * The templated parameter, 'T', is the value type of the polygon vertices.
 *
 * Note that, if the number of vertices smaller than 3, the polygon will set to
 * be empty.
 */
template <typename T>
class Polygon2D {
public:
    using Point = Point2D<T>;
    using value_type = T;
    using Iterator = typename Array<Point>::const_iterator;

    Polygon2D() = default;

    explicit Polygon2D(const Array<Point>& vertices)
        : Polygon2D(vertices.begin(), vertices.end()) {}

    template <typename InputIterator>
    Polygon2D(InputIterator first, InputIterator last)
        : vertices_(first, last), bounding_box_(first, last) {
        Initialize();
    }

    /**
     * Construct a axis aligned rectangle polygon from a box.
     */
    explicit Polygon2D(const Box2D<T>& rectangle)
        : bounding_box_(rectangle) {
        assert(!rectangle.empty());

        vertices_.emplace_back(rectangle.x_min(), rectangle.y_min());
        vertices_.emplace_back(rectangle.x_max(), rectangle.y_min());
        vertices_.emplace_back(rectangle.x_max(), rectangle.y_max());
        vertices_.emplace_back(rectangle.x_min(), rectangle.y_max());
        Initialize();
    }

    /**
     * Return true if the polygon does not have any vertices.
     */
    bool empty() const {
        return size_ == 0;
    }

    /**
     * Clear the polygon.
     */
    void clear() {
        bounding_box_ = Box2D<T>();
        vertices_.clear();
        size_ = 0;
    }

    /**
     * Return the size of vertices.
     */
    int size() const {
        return size_;
    }

    /**
     * Get the i-th edge of polygon.
     */
    Segment2D<T> edge(int i) const {
        return (i + 1 == size_) ? Segment2D<T>(vertices_[i], vertices_[0])
                                : Segment2D<T>(vertices_[i], vertices_[i + 1]);
    }

    /**
     * Get the i-th vertex of polygon.
     */
    const Point& vertex(int i) const {
        return vertices_[i];
    }

    /**
     * Return the vertices of polygon.
     */
    const Array<Point>& vertices() const {
        return vertices_;
    }

    /**
     * Return the bounding box of this polygon.
     */
    const Box2D<T>& bounding_box() const {
        return bounding_box_;
    }

    /**
     * Return the iterator to the first vertex.
     */
    Iterator begin() const {
        return vertices_.begin();
    }

    /**
     * Return the Iterator to the last vertex + 1.
     */
    Iterator end() const {
        return vertices_.end();
    }

    /**
     * Return the area of polygon.
     */
    double Area() const {
        double s = 0.0;
        for (int i = 1; i < size_; ++i) {
            s += static_cast<double>(vertices_[i - 1].x) * vertices_[i].y -
                 static_cast<double>(vertices_[i].x) * vertices_[i - 1].y;
        }
        if (size_ > 0) {
            s += static_cast<double>(vertices_.back().x) * vertices_[0].y -
                 static_cast<double>(vertices_[0].x) * vertices_.back().y;
        }
        return 0.5 * std::fabs(s);
    }

    /**
     * Check if two polygons are equal.
     */
    bool operator ==(const Polygon2D& rhs) const {
        if (size_ != rhs.size_) return false;
        if (size_ == 0) return true;

        Iterator s1 = std::min_element(begin(), end());
        Iterator s2 = std::min_element(rhs.begin(), rhs.end());
        Iterator t1 = s1, t2 = s2;

        bool is_equal = true;
        for (int i = 0; i < size_; ++i) {
            if (*t1 != *t2) {
                is_equal = false;
                break;
            }

            if (++t1 == end()) t1 = begin();
            if (++t2 == rhs.end()) t2 = rhs.begin();
        }

        if (!is_equal) {
            is_equal = true;
            t1 = s1;
            t2 = s2;
            for (int i = 0; i < size_; ++i) {
                if (*t1 != *t2) {
                    is_equal = false;
                    break;
                }

                if (++t1 == end()) t1 = begin();
                if (t2 == rhs.begin()) t2 = rhs.end();
                --t2;
            }
        }

        return is_equal;
    }

    bool operator !=(const Polygon2D& rhs) const {
        return !(*this == rhs);
    }

    /**
     * Check if this polygon is in clockwise order.
     */
    bool IsClockwise() const {
        if (size_ == 0) return true;

        Iterator left_most = std::min_element(vertices_.begin(),
                                              vertices_.end());
        Iterator next = left_most;
        if (++next == vertices_.end()) {
            next = vertices_.begin();
        }

        Iterator prev = left_most == vertices_.begin() ? vertices_.end()
                                                       : left_most;
        --prev;

        return geometry::Orientation(*prev, *left_most, *next) < 0;
    }

    /**
     * Check if this polygon is in anticlockwise order.
     */
    bool IsAnticlockwise() const {
        if (size_ == 0) return true;

        Iterator left_most = std::min_element(vertices_.begin(),
                                              vertices_.end());
        Iterator next = left_most;
        if (++next == vertices_.end()) {
            next = vertices_.begin();
        }

        Iterator prev = left_most == vertices_.begin() ? vertices_.end()
                                                       : left_most;
        --prev;

        return geometry::Orientation(*prev, *left_most, *next) > 0;
    }

protected:
    /**
     * Initialize the polygon.
     */
    void Initialize() {
        if (vertices_.empty()) return;

        // Erase the duplicate points.
        int n = static_cast<int>(std::unique(vertices_.begin(), 
                                             vertices_.end()) - 
                                 vertices_.begin());
        vertices_.resize(n);
        if (vertices_.back() == vertices_.front()) {
            vertices_.pop_back();
        }
        size_ = vertices_.size();

        if (size_ < 3) {
            vertices_.clear();
            size_ = 0;
        }
    }

    // Number of vertices.
    int size_ = 0;

    // Vertices of polygon.
    Array<Point> vertices_;

    // Bounding box of polygon.
    Box2D<T> bounding_box_;
};

using IPolygon2D = Polygon2D<int>;
using FPolygon2D = Polygon2D<float>;
using RPolygon2D = Polygon2D<double>;

} // namespace cl

#endif // GEOMETRY_KERNEL_POLYGON_2D_H_
