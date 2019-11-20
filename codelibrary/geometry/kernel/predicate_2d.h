//
// Copyright 2017 Yangbin Lin. All Rights Reserved.
//
// Author: yblin@jmu.edu.cn (Yangbin Lin)
//
// This file is part of the Code Library.
//

#ifndef GEOMETRY_KERNEL_PREDICATE_2D_H_
#define GEOMETRY_KERNEL_PREDICATE_2D_H_

#include <algorithm>
#include <cfenv>

#include "codelibrary/geometry/kernel/point_2d.h"
#include "codelibrary/math/number/determinant.h"
#include "codelibrary/math/number/exact_float.h"
#include "codelibrary/math/number/interval_float.h"

namespace cl {
namespace geometry {

/**
 * Compute the determinant for orientation test, i.e.,
 *
 * |px, py, 1|   |qx - px, qy - py|
 * |qx, qy, 1| = |                |
 * |rx, ry, 1|   |rx - px, ry - py|
 *
 * The return value (stored in 'det') is twice of the signed area of (p, q, r).
 */
template <typename T>
void OrientationDeterminant(const T& px, const T& py, const T& qx, const T& qy,
                            const T& rx, const T& ry, T* det) {
    // Optimized for ExactFloat and IntervalFloat.
    T t = qy - py;
    t *= rx - px;
    *det = qx - px;
    *det *= ry - py;
    *det -= t;
}

/**
 * Check the orientation from pq to pr.
 *
 * Return +1 if three points occur in counter-clockwise;
 *        -1 if they occur in clockwise;
 *         0 if they are collinear.
 */
template <typename T>
int Orientation(const Point2D<T>& p, const Point2D<T>& q, const Point2D<T>& r) {
    // Static filter, adapted from CGAL.
    double pqx = static_cast<double>(q.x) - p.x;
    double pqy = static_cast<double>(q.y) - p.y;
    double prx = static_cast<double>(r.x) - p.x;
    double pry = static_cast<double>(r.y) - p.y;

    // Then semi-static filter.
    double maxx = std::fabs(pqx);
    double maxy = std::fabs(pqy);
    double aprx = std::fabs(prx);
    double apry = std::fabs(pry);

    if (maxx < aprx) maxx = aprx;
    if (maxy < apry) maxy = apry;
    if (maxx > maxy) std::swap(maxx, maxy);

    // Protect against underflow in the computation of eps.
    if (maxx < 1e-146) {
        if (maxx == 0.0) return 0;
    // Protect against overflow in the computation of det.
    } else if (maxy < 1e153) {
        double eps = 8.8872057372592798e-16 * maxx * maxy;
        double det = pqx * pry - pqy * prx;
        if (det > eps)  return +1;
        if (det < -eps) return -1;
    }

    // Interval filter.
    if (fegetround() != FE_UPWARD) fesetround(FE_UPWARD);
    IntervalFloat det1;
    OrientationDeterminant(IntervalFloat(p.x), IntervalFloat(p.y),
                           IntervalFloat(q.x), IntervalFloat(q.y),
                           IntervalFloat(r.x), IntervalFloat(r.y),
                           &det1);
    if (det1.lower() > 0.0) return +1;
    if (det1.upper() < 0.0) return -1;
    if (det1.lower() == 0.0 && det1.upper() == 0.0) return 0;

    // Exact computation.
    ExactFloat det2;
    OrientationDeterminant(ExactFloat(p.x), ExactFloat(p.y),
                           ExactFloat(q.x), ExactFloat(q.y),
                           ExactFloat(r.x), ExactFloat(r.y),
                           &det2);
    return det2.sign();
}

/**
 * Compute the determinant for in circle test, i.e.,
 *
 * |px, py, px^2 + py^2, 1|
 * |qx, qy, qx^2 + qy^2, 1|
 * |rx, ry, rx^2 + ry^2, 1|
 * |tx, ty, tx^2 + ty^2, 1|
 */
template <typename T>
void IncircleDeterminant(const T& px, const T& py, const T& qx, const T& qy,
                         const T& rx, const T& ry, const T& tx, const T& ty,
                         T* det) {
    T pqx = qx - px, pqy = qy - py, prx = rx - px, pry = ry - py, ptx = tx - px;
    T pty = ty - py, qtx = tx - qx, qty = ty - qy, qrx = rx - qx, qry = ry - qy;

    T t1, t2;
    Determinant(pqx, pqy, ptx, pty, &t1);
    ptx *= qtx;
    pty *= qty;
    ptx += pty;

    Determinant(pqx, pqy, prx, pry, &t2);
    prx *= qrx;
    pry *= qry;
    prx += pry;

    Determinant(t1, ptx, t2, prx, det);
}

/**
 * In circle test.
 *
 * Check if point t is inside a circle passing through p, q and r.
 * Require p, q, r appear in counter-clockwise order.
 *
 * Return +1 if the point pd lies inside the circle;
 *        -1 if it lies outside;
 *         0 if the four points are co-circular.
 */
template <typename T>
int InCircle(const Point2D<T>& p, const Point2D<T>& q,
             const Point2D<T>& r, const Point2D<T>& t) {
    if (p == t || q == t || r == t) return 0;

    // Static filter, adapted from CGAL.
    double qpx = static_cast<double>(q.x) - p.x;
    double qpy = static_cast<double>(q.y) - p.y;
    double rpx = static_cast<double>(r.x) - p.x;
    double rpy = static_cast<double>(r.y) - p.y;
    double tpx = static_cast<double>(t.x) - p.x;
    double tpy = static_cast<double>(t.y) - p.y;
    double tqx = static_cast<double>(t.x) - q.x;
    double tqy = static_cast<double>(t.y) - q.y;
    double rqx = static_cast<double>(r.x) - q.x;
    double rqy = static_cast<double>(r.y) - q.y;

    // We compute the semi-static bound.
    double maxx = std::fabs(qpx);
    double maxy = std::fabs(qpy);
    double arpx = std::fabs(rpx);
    double arpy = std::fabs(rpy);
    double atqx = std::fabs(tqx);
    double atqy = std::fabs(tqy);
    double atpx = std::fabs(tpx);
    double atpy = std::fabs(tpy);
    double arqx = std::fabs(rqx);
    double arqy = std::fabs(rqy);

    if (maxx < arpx) maxx = arpx;
    if (maxx < atpx) maxx = atpx;
    if (maxx < atqx) maxx = atqx;
    if (maxx < arqx) maxx = arqx;
    if (maxy < arpy) maxy = arpy;
    if (maxy < atpy) maxy = atpy;
    if (maxy < atqy) maxy = atqy;
    if (maxy < arqy) maxy = arqy;

    if (maxx > maxy) std::swap(maxx, maxy);

    // Protect against underflow in the computation of eps.
    if (maxx < 1e-73) {
        if (maxx == 0.0) return 0;
    } else if (maxy < 1e76) {
        double eps = 8.8878565762001373e-15 * maxx * maxy * (maxy * maxy);
        double det = (qpx * tpy - qpy * tpx) * (rpx * rqx + rpy * rqy) -
                     (tpx * tqx + tpy * tqy) * (qpx * rpy - qpy * rpx);
        if (det > eps)  return +1;
        if (det < -eps) return -1;
    }

    // Interval filter.
    if (fegetround() != FE_UPWARD) fesetround(FE_UPWARD);
    IntervalFloat det1;
    IncircleDeterminant(IntervalFloat(p.x), IntervalFloat(p.y),
                        IntervalFloat(q.x), IntervalFloat(q.y),
                        IntervalFloat(r.x), IntervalFloat(r.y),
                        IntervalFloat(t.x), IntervalFloat(t.y),
                        &det1);
    if (det1.lower() > 0.0) return +1;
    if (det1.upper() < 0.0) return -1;
    if (det1.lower() == 0.0 && det1.upper() == 0.0) return 0;

    // Exact computation.
    ExactFloat det;
    IncircleDeterminant(ExactFloat(p.x), ExactFloat(p.y),
                        ExactFloat(q.x), ExactFloat(q.y),
                        ExactFloat(r.x), ExactFloat(r.y),
                        ExactFloat(t.x), ExactFloat(t.y),
                        &det);
    return det.sign();
}

} // namespace geometry
} // namespace cl

#endif // GEOMETRY_KERNEL_PREDICATE_2D_H_
