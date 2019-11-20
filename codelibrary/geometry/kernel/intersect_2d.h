//
// Copyright 2014 Yangbin Lin. All Rights Reserved.
//
// Author: yblin@jmu.edu.cn (Yangbin Lin)
//
// This file is part of the Code Library.
//
// Robust intersection test between 2D kernel geometric objects.
//

#ifndef GEOMETRY_KERNEL_INTERSECT_2D_H_
#define GEOMETRY_KERNEL_INTERSECT_2D_H_

#include <cfloat>

#include "codelibrary/geometry/kernel/circle_2d.h"
#include "codelibrary/geometry/kernel/line_2d.h"
#include "codelibrary/geometry/kernel/multi_polygon_2d.h"
#include "codelibrary/geometry/kernel/point_2d.h"
#include "codelibrary/geometry/kernel/polygon_2d.h"
#include "codelibrary/geometry/kernel/predicate_2d.h"
#include "codelibrary/geometry/kernel/segment_2d.h"

namespace cl {
namespace geometry {

// All Intersect() functions are used to test weather two objects intersect
// (including touch and inside). The computation is strictly exact.
template <typename T>
bool Intersect(const Box2D<T>& box1, const Box2D<T>& box2) {
    return box1.x_max() >= box2.x_min() && box1.x_min() <= box2.x_max() &&
           box1.y_max() >= box2.y_min() && box1.y_min() <= box2.y_max();
}

template <typename T>
bool Intersect(const Point2D<T>& point, const Box2D<T>& box) {
    return (point.x >= box.x_min() && point.x <= box.x_max() &&
            point.y >= box.y_min() && point.y <= box.y_max());
}

template <typename T>
bool Intersect(const Box2D<T>& box, const Point2D<T>& point) {
    return Intersect(point, box);
}

template <typename T>
bool Intersect(const Point2D<T>& point, const Segment2D<T>& seg) {
    return point >= seg.lower_point() && point <= seg.upper_point() &&
           Orientation(seg.lower_point(), seg.upper_point(), point) == 0;
}

template <typename T>
bool Intersect(const Segment2D<T>& seg, const Point2D<T>& point) {
    return Intersect(point, seg);
}

template <typename T>
bool Intersect(const Point2D<T>& point, const Line2D<T>& line) {
    return Orientation(line.point1(), line.point2(), point) == 0;
}

template <typename T>
bool Intersect(const Line2D<T>& line, const Point2D<T>& point) {
    return Intersect(point, line);
}

template <typename T>
bool Intersect(const Point2D<T>& p, const Polygon2D<T>& polygon) {
    if (!Intersect(p, polygon.bounding_box())) {
        return false;
    }

    int count = 0;
    for (int i = 0; i < polygon.size(); ++i) {
        const Segment2D<T> seg = polygon.edge(i);
        if (seg.lower_point().x > p.x || seg.upper_point().x <= p.x) continue;

        int o = Orientation(seg.lower_point(), seg.upper_point(), p);
        if (o == 0) return true;
        if (o < 0) ++count;
    }
    return count % 2 != 0;
}

template <typename T>
bool Intersect(const Polygon2D<T>& polygon, const Point2D<T>& point) {
    return Intersect(point, polygon);
}

template <typename T>
bool Intersect(const Point2D<T>& p, const MultiPolygon2D<T>& multi_polygon) {
    if (!Intersect(p, multi_polygon.bounding_box())) {
        return false;
    }

    int count = 0;
    for (const auto& b : multi_polygon.boundaries()) {
        const Polygon2D<T>& polygon = b.polygon;
        if (!Intersect(p, polygon.bounding_box())) continue;

        for (int j = 0; j < polygon.size(); ++j) {
            Segment2D<T> seg = polygon.edge(j);
            if (seg.lower_point().x > p.x || seg.upper_point().x <= p.x) {
                continue;
            }

            int o = Orientation(seg.lower_point(), seg.upper_point(), p);
            if (o == 0) return true;
            if (o < 0)  ++count;
        }
    }

    return (count % 2 != 0);
}

template <typename T>
bool Intersect(const MultiPolygon2D<T>& multi_polygon, const Point2D<T>& p) {
    return Intersect(p, multi_polygon);
}

template <typename T>
bool Intersect(const Segment2D<T>& s1, const Segment2D<T>& s2) {
    if (!Intersect(s1.bounding_box(), s2.bounding_box())) return false;

    // Cross test.
    int o1 = Orientation(s1.lower_point(), s1.upper_point(), s2.lower_point());
    int o2 = Orientation(s1.lower_point(), s1.upper_point(), s2.upper_point());
    if (o1 == o2 && o1 != 0) {
        // Two endpoints of s2 lay on the same side (left or right) of s1, thus
        // two segments never be intersected.
        return false;
    }

    int o3 = Orientation(s2.lower_point(), s2.upper_point(), s1.lower_point());
    int o4 = Orientation(s2.lower_point(), s2.upper_point(), s1.upper_point());

    // Two endpoints of s1 lay on the same side (left or right) of s2, thus
    // two segments never be intersected.
    return !(o3 == o4 && o3 != 0);
}

template <typename T>
bool Intersect(const Segment2D<T>& seg, const Box2D<T>& box) {
    const Box2D<T>& box1 = seg.bounding_box();
    if (!Intersect(box1, box)) return false;

    if (Intersect(seg.lower_point(), box) ||
        Intersect(seg.upper_point(), box)) return true;

    Point2D<T> p1(box.x_min(), box.y_min());
    Point2D<T> p2(box.x_max(), box.y_min());
    Point2D<T> p3(box.x_max(), box.y_max());
    Point2D<T> p4(box.x_min(), box.y_max());

    int o1 = Orientation(seg.lower_point(), seg.upper_point(), p1);
    int o2 = Orientation(seg.lower_point(), seg.upper_point(), p2);
    int o3 = Orientation(seg.lower_point(), seg.upper_point(), p3);
    int o4 = Orientation(seg.lower_point(), seg.upper_point(), p4);

    if (box.x_min() > box1.x_min() && box.x_min() < box1.x_max()) {
        if (o1 == 0) return true;
        if (o4 == 0) return true;
        if (o1 != o4) return true;
    }

    if (box.x_max() > box1.x_min() && box.x_max() < box1.x_max()) {
        if (o2 == 0) return true;
        if (o3 == 0) return true;
        if (o2 != o3) return true;
    }

    if (box.y_min() > box1.y_min() && box.y_min() < box1.y_max()) {
        if (o1 == 0) return true;
        if (o2 == 0) return true;
        if (o1 != o2) return true;
    }

    if (box.y_max() > box1.y_min() && box.y_max() < box1.y_max()) {
        if (o3 == 0) return true;
        if (o4 == 0) return true;
        if (o3 != o4) return true;
    }

    return false;
}

template <typename T>
bool Intersect(const Box2D<T>& box, const Segment2D<T>& seg) {
    return Intersect(seg, box);
}

template <typename T>
bool Intersect(const Segment2D<T>& seg, const Line2D<T>& line) {
    int o1 = Orientation(line.point1(), line.point2(), seg.lower_point());
    if (o1 == 0) return true;

    int o2 = Orientation(line.point1(), line.point2(), seg.upper_point());
    if (o2 == 0) return true;

    return o1 != o2;
}

template <typename T>
bool Intersect(const Line2D<T>& line, const Segment2D<T>& seg) {
    return Intersect(seg, line);
}

template <typename T>
bool Intersect(const Segment2D<T>& seg, const Polygon2D<T>& polygon) {
    if (!Intersect(seg, polygon.bounding_box())) return false;

    if (Intersect(seg.lower_point(), polygon) ||
        Intersect(seg.upper_point(), polygon)) return true;

    for (int i = 0; i < polygon.size(); ++i) {
        if (Intersect(seg, polygon.edge(i))) return true;
    }

    return false;
}

template <typename T>
bool Intersect(const Polygon2D<T>& polygon, const Segment2D<T>& seg) {
    return Intersect(seg, polygon);
}

template <typename T>
bool Intersect(const Segment2D<T>& seg, const MultiPolygon2D<T>& polygon) {
    if (!Intersect(seg, polygon.bounding_box())) return false;

    if (Intersect(seg.lower_point(), polygon) ||
        Intersect(seg.upper_point(), polygon)) return true;

    for (int i = 0; i < polygon.size_boundaries(); ++i) {
        const Polygon2D<T>& poly = polygon.boundaries()[i].polygon;
        for (int j = 0; j < poly.size(); ++j) {
            if (Intersect(seg, poly.edge(i))) return true;
        }
    }

    return false;
}

template <typename T>
bool Intersect(const MultiPolygon2D<T>& polygon, const Segment2D<T>& seg) {
    return Intersect(seg, polygon);
}

template <typename T>
bool Intersect(const Polygon2D<T>& polygon1, const Polygon2D<T>& polygon2) {
    if (!Intersect(polygon1.bounding_box(), polygon2.bounding_box()))
        return false;

    for (int i = 0; i < polygon1.size(); ++i) {
        if (Intersect(polygon1.points()[i], polygon2)) return true;
        if (Intersect(polygon1.edge(i), polygon2)) return true;
    }

    return false;
}

/**
 * Compute the intersection point between line(p1, p2) and line(p3, p4).
 * 
 * This function do not check if two lines are intersect. Therefore, the user
 * should make sure the two lines are intersect. Otherwise, an INF value may be
 * returned.
 */
inline void IntersectionPoint(double x1, double y1, double x2, double y2,
                              double x3, double y3, double x4, double y4,
                              RPoint2D* intersection) {
    double t1 = x1 * y2 - y1 * x2;
    double t2 = x3 * y4 - y3 * x4;
    double t = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4);
    intersection->x = (t1 * (x3 - x4) - (x1 - x2) * t2) / t;
    intersection->y = (t1 * (y3 - y4) - (y1 - y2) * t2) / t;
}

/**
 * Check if two lines are intersect (excluding touch). It also optional computes
 * the intersection point if two lines are intersect. We choose double as
 * coordinate type of the intersection to avoid the overflow. Also note that we
 * do not use the exact predication here.
 */
template <typename T>
bool Cross(const Line2D<T>& line1, const Line2D<T>& line2,
           RPoint2D* intersection = nullptr) {
    RPoint2D p;
    IntersectionPoint(line1.point1().x, line1.point1().y,
                      line1.point2().x, line1.point2().y,
                      line2.point1().x, line2.point1().y,
                      line2.point2().x, line2.point2().y,
                      &p);
    if (intersection) *intersection = p;
    return std::isfinite(p.x) && std::isfinite(p.y);
}

/**
 * Check if two line segments are intersect (excluding touch). It also optional
 * computes the intersection point if two line segments are intersect.
 *
 * The predication here is exact.
 */
template <typename T>
bool Cross(const Segment2D<T>& s1, const Segment2D<T>& s2,
           RPoint2D* intersection = nullptr) {
    if (!Intersect(s1.bounding_box(), s2.bounding_box())) return false;

    // Cross test.
    int o1 = Orientation(s1.lower_point(), s1.upper_point(), s2.lower_point());
    int o2 = Orientation(s1.lower_point(), s1.upper_point(), s2.upper_point());
    if (o1 == o2 || o1 == 0 || o2 == 0) {
        // Two endpoints of s2 lay on the same side (left or right) of s1, thus
        // two segments never be intersected.
        return false;
    }

    int o3 = Orientation(s2.lower_point(), s2.upper_point(), s1.lower_point());
    int o4 = Orientation(s2.lower_point(), s2.upper_point(), s1.upper_point());
    if (o3 == o4 || o3 == 0 || o4 == 0) {
        // Two endpoints of s1 lay on the same side (left or right) of s1, thus
        // two segments never be intersected.
        return false;
    }

    if (intersection) {
        // There must be an intersection between two segments.
        RPoint2D p;
        IntersectionPoint(s1.lower_point().x, s1.lower_point().y, 
                          s1.upper_point().x, s1.upper_point().y,
                          s2.lower_point().x, s2.lower_point().y,
                          s2.upper_point().x, s2.upper_point().y,
                          &p);

        // Ensure the intersection is inside the bounding box.
        double x_min = std::min(s1.bounding_box().x_min(),
                                s2.bounding_box().x_min());
        double x_max = std::max(s1.bounding_box().x_max(),
                                s2.bounding_box().x_max());
        double y_min = std::min(s1.bounding_box().y_min(),
                                s2.bounding_box().y_min());
        double y_max = std::max(s1.bounding_box().y_max(),
                                s2.bounding_box().y_max());

        intersection->x = Clamp(p.x, x_min, x_max);
        intersection->y = Clamp(p.y, y_min, y_max);
    }

    return true;
}

/**
 * Check if a line (or line segment) is intersect with a box (excluding touch).
 * It optimal computes the intersection range (t1, t2), which means the
 * intersection line is:
 *   (line.p + t1 * line.direction, line.p + t2 * line.direction).
 */
template <typename T>
bool Cross(const Box2D<T>& box, const Segment2D<T>& segment,
           std::pair<double, double>* intersection = nullptr) {
    if (!box.is_valid()) return false;

    const Vector2D v = segment.direction();
    RVector2D direction(v.x, v.y);
    double t_min = 0.0, t_max = 1.0;
    for (int i = 0; i < 2; ++i) {
        // It handles the inf cases correct.
        double inverse_direction = 1.0 / direction[i];
        double t1 = (box.min(i) - segment.lower_point()[i]) * inverse_direction;
        double t2 = (box.max(i) - segment.lower_point()[i]) * inverse_direction;
        if (inverse_direction < 0.0) {
            std::swap(t1, t2);
        }

        // min(a, NAN) -> a, min(NAN, a) -> NAN
        t_min = std::max(t_min, t1);
        t_max = std::min(t_max, t2);

        if (t_max < t_min) return false;
    }
    if (intersection) {
        *intersection = std::make_pair(t_min, t_max);
    }
    return true;
}
template <typename T>
bool Cross(const Segment2D<T>& segment, const Box2D<T>& box,
           std::pair<double, double>* intersection = nullptr) {
    return Cross(box, segment, intersection);
}

template <typename T>
bool Cross(const Box2D<T>& box, const Line2D<T>& line,
           std::pair<double, double>* intersection = nullptr) {
    if (!box.is_valid()) return false;

    double t_min = -DBL_MAX, t_max = DBL_MAX;
    const Vector2D v = line.direction();
    RVector2D direction(v.x, v.y);
    for (int i = 0; i < 2; ++i) {
        // It handles the inf cases correct.
        double inverse_direction = 1.0 / direction[i];
        double t1 = (box.min(i) - line.point1()[i]) * inverse_direction;
        double t2 = (box.max(i) - line.point1()[i]) * inverse_direction;
        if (inverse_direction < 0.0) {
            std::swap(t1, t2);
        }
        t_min = std::max(t1, t_min);
        t_max = std::min(t2, t_max);
        if (t_max < t_min) return false;
    }

    if (intersection) {
        *intersection = std::make_pair(t_min, t_max);
    }
    return true;
}
template <typename T>
bool Cross(const Line2D<T>& line, const Box2D<T>& box,
           std::pair<double, double>* intersection = nullptr) {
    return Cross(box, line, intersection);
}

template <typename T>
bool Cross(const Box2D<T>& box, const Line2D<T>& line,
           RSegment2D* intersection = nullptr) {
    std::pair<double, double> p;
    if (!Cross(box, line, &p)) {
        return false;
    }

    if (intersection) {
        RPoint2D p1, p2;
        RVector2D direction(line.direction().x, line.direction().y);
        p1.x = line.point1().x + direction.x * p.first;
        p1.y = line.point1().y + direction.y * p.first;
        p2.x = line.point1().x + direction.x * p.second;
        p2.y = line.point1().y + direction.y * p.second;
        *intersection = RSegment2D(p1, p2);
    }
    return true;
}
template <typename T>
bool Cross(const Line2D<T>& line, const Box2D<T>& box,
           RSegment2D* intersection = nullptr) {
    return Cross(box, line, intersection);
}

/**
 * Check if the line is cross with a polygon. It optional computes all
 * intersection points. And these intersection points are ordered along the
 * line.
 * 
 * The predication here is exact.
 */
template <typename T>
bool Cross(const Line2D<T>& line, const Polygon2D<T>& polygon,
           Array<RPoint2D>* intersection = nullptr) {
    using Point = Point2D<T>;

    // Make sure the polygon is in anti-clockwise.
    int n = polygon.size();

    // Generate 1 to n.
    Array<int> seq(n);
    std::iota(seq.begin(), seq.end(), 0);
    if (polygon.IsClockwise()) std::reverse(seq.begin(), seq.end());

    Point p1 = std::min(line.point1(), line.point2());
    Point p2 = std::max(line.point1(), line.point2());

    // Store the intersection points.
    int cur_o = 0, next_o = 0, prev_o = 0;
    int n_intersection = 0;
    for (int i = 0; i < n; ++i) {
        const Point& cur_p  = polygon.vertices()[seq[i]];
        const Point& next_p = i + 1 < n ? polygon.vertex(seq[i + 1])
                                        : polygon.vertex(seq[0]);
        const Point& prev_p = i == 0 ? polygon.vertex(seq.back())
                                     : polygon.vertex(seq[i - 1]);

        // Cross test.
        if (i == 0) {
            prev_o = Orientation(p1, p2, prev_p);
            cur_o  = Orientation(p1, p2, cur_p);
        } else {
            prev_o = cur_o;
            cur_o = next_o;
        }
        next_o = Orientation(p1, p2, next_p);

        if (cur_o == next_o && cur_o != 0) continue;
        if (cur_o == 0) {
            if (prev_o == next_o) continue;

            if (prev_o == 0) {
                if (cur_p < prev_p && next_o < 0) {
                    ++n_intersection;
                    if (intersection)
                        intersection->emplace_back(cur_p.x, cur_p.y);
                }
                if (prev_p < cur_p && next_o > 0) {
                    ++n_intersection;
                    if (intersection)
                        intersection->emplace_back(cur_p.x, cur_p.y);
                }
            } else if (next_o == 0) {
                if (next_p < cur_p && prev_o < 0) {
                    ++n_intersection;
                    if (intersection)
                        intersection->emplace_back(cur_p.x, cur_p.y);
                }
                if (cur_p < next_p && prev_o > 0) {
                    ++n_intersection;
                    if (intersection)
                        intersection->emplace_back(cur_p.x, cur_p.y);
                }
            } else {
                ++n_intersection;
                if (intersection)
                    intersection->emplace_back(cur_p.x, cur_p.y);
            }
        } else if (cur_o != 0 && next_o != 0) {
            ++n_intersection;

            if (intersection) {
                RPoint2D p;
                IntersectionPoint(p1.x, p1.y, p2.x, p2.y,
                                  cur_p.x, cur_p.y, next_p.x, next_p.y,
                                  &p);

                // Ensure the intersection is inside the bounding box.
                double x_min = std::min(cur_p.x, next_p.x);
                double x_max = std::max(cur_p.x, next_p.x);
                double y_min = std::min(cur_p.y, next_p.y);
                double y_max = std::max(cur_p.y, next_p.y);

                p.x = Clamp(p.x, x_min, x_max);
                p.y = Clamp(p.y, y_min, y_max);
                intersection->push_back(p);
            }
        }
    }
    assert(n_intersection % 2 == 0);

    return n_intersection != 0;
}
template <typename T>
bool Cross(const Polygon2D<T>& polygon, const Line2D<T>& line,
           Array<RPoint2D>* intersection) {
    return Cross(line, polygon, intersection);
}

/**
 * Check if circle1 is cross with circle2. It optional computes two intersection
 * points of two circles.
 * 
 * Note that, the predication here is not exact.
 */
template <typename T>
bool Cross(const Circle2D<T>& circle1, const Circle2D<T>& circle2,
           std::pair<RPoint2D, RPoint2D>* p = nullptr) {
    assert(p);

    double dx = circle2.center().x - circle1.center().x;
    double dy = circle2.center().y - circle1.center().y;

    double dis = std::sqrt(dx * dx + dy * dy);
    if (dis == 0.0) return false;

    double r1 = circle1.radius();
    double r2 = circle2.radius();

    // No intersection.
    if (dis > r1 + r2) return false;

    // One circle is contained in the other.
    if (dis < std::fabs(r1 - r2)) {
        return false;
    }

    if (p) {
        double a = (r1 * r1 - r2 * r2 + dis * dis) / (2.0 * dis);
        double x2 = dx * a / dis + circle1.center().x;
        double y2 = dy * a / dis + circle1.center().y;

        double dis2 = std::sqrt((r1 * r1) - (a * a));

        double rx = -dy * (dis2 / dis);
        double ry =  dx * (dis2 / dis);

        p->first.x  = x2 + rx;
        p->second.x = x2 - rx;
        p->first.y  = y2 + ry;
        p->second.y = y2 - ry;
    }

    return true;
}

template <typename T>
bool Cross(const Line2D<T>& line, const Circle2D<T>& circle,
           std::pair<RPoint2D, RPoint2D>* p = nullptr) {
    double dx = line.direction().x;
    double dy = line.direction().y;
    double dr = std::sqrt(dx * dx + dy * dy);
    if (dr == 0.0) return false;

    double cx = circle.center().x;
    double cy = circle.center().y;
    double d = dy * line.point1().x - dy * cx - dx * line.point1().y + dx * cy;

    double r = circle.radius();
    double delt = r * r * dr * dr - d * d;

    // No intersection.
    if (delt < 0.0) return false;

    if (p) {
        double tmp1 = std::sqrt(delt);
        double tmp2 = dr * dr;
        p->first.x  = ( d * dy + (dy > 0 ? dx : -dx) * tmp1) / tmp2 + cx;
        p->first.y  = (-d * dx + std::fabs(dy) * tmp1) / tmp2 + cy;
        p->second.x = ( d * dy - (dy > 0 ? dx : -dx) * tmp1) / tmp2 + cx;
        p->second.y = (-d * dx - std::fabs(dy) * tmp1) / tmp2 + cy;
    }

    return true;
}

template <typename T>
bool Cross(const Circle2D<T>& circle, const Line2D<T>& line,
           std::pair<RPoint2D, RPoint2D>* p = nullptr) {
    return Cross(circle, line, p);
}

} // namespace geometry
} // namespace cl

#endif // GEOMETRY_KERNEL_INTERSECT_2D_H_
