//
// Copyright 2017 Yangbin Lin. All Rights Reserved.
//
// Author: yblin@jmu.edu.cn (Yangbin Lin)
//
// This file is part of the Code Library.
//

#ifndef GEOMETRY_KERNEL_PREDICATE_3D_H_
#define GEOMETRY_KERNEL_PREDICATE_3D_H_

#include <algorithm>
#include <cfenv>

#include "codelibrary/geometry/kernel/predicate_2d.h"
#include "codelibrary/geometry/kernel/point_3d.h"
#include "codelibrary/math/number/determinant.h"
#include "codelibrary/math/number/exact_float.h"
#include "codelibrary/math/number/interval_float.h"

namespace cl {
namespace geometry {

/**
 * Compute the determinant for 3D orientation test, i.e.,
 *
 * |px, py, pz, 1|
 * |qx, qy, qz, 1|
 * |rx, ry, rz, 1|
 * |sx, sy, sz, 1|
 */
template <typename T>
void OrientationDeterminant(const T& px, const T& py, const T& pz,
                            const T& qx, const T& qy, const T& qz,
                            const T& rx, const T& ry, const T& rz,
                            const T& sx, const T& sy, const T& sz,
                            T* det) {
    T pqx = qx - px, pqy = qy - py, pqz = qz - pz, prx = rx - px, pry = ry - py,
      prz = rz - pz, psx = sx - px, psy = sy - py, psz = sz - pz;
    Determinant(pqx, pqy, pqz, prx, pry, prz, psx, psy, psz, det);
}

/**
 * Compute the orientation from s lies above the plane to pr.
 *
 * Return +1 if s lies above the plane defined by p, q, and r.
 *           "Above" is defined so that p, q, and r appear in counter-clockwise
 *           order when viewed from above the plane.
 *        -1 if s lies below the plane.
 *         0 if they are coplanar.
 */
template <typename T>
int Orientation(const Point3D<T>& p, const Point3D<T>& q,
                const Point3D<T>& r, const Point3D<T>& s) {
    // Static filter, adapted from CGAL.
    double px = static_cast<double>(p.x);
    double py = static_cast<double>(p.y);
    double pz = static_cast<double>(p.z);
    double qx = static_cast<double>(q.x);
    double qy = static_cast<double>(q.y);
    double qz = static_cast<double>(q.z);
    double rx = static_cast<double>(r.x);
    double ry = static_cast<double>(r.y);
    double rz = static_cast<double>(r.z);
    double sx = static_cast<double>(s.x);
    double sy = static_cast<double>(s.y);
    double sz = static_cast<double>(s.z);

    double pqx = qx - px;
    double pqy = qy - py;
    double pqz = qz - pz;
    double prx = rx - px;
    double pry = ry - py;
    double prz = rz - pz;
    double psx = sx - px;
    double psy = sy - py;
    double psz = sz - pz;

    double maxx = std::fabs(pqx);
    double maxy = std::fabs(pqy);
    double maxz = std::fabs(pqz);

    double aprx = std::fabs(prx);
    double apsx = std::fabs(psx);
    double apry = std::fabs(pry);
    double apsy = std::fabs(psy);
    double aprz = std::fabs(prz);
    double apsz = std::fabs(psz);

    if (maxx < aprx) maxx = aprx;
    if (maxx < apsx) maxx = apsx;
    if (maxy < apry) maxy = apry;
    if (maxy < apsy) maxy = apsy;
    if (maxz < aprz) maxz = aprz;
    if (maxz < apsz) maxz = apsz;

    // Sort maxx < maxy < maxz.
    if (maxx > maxz) std::swap(maxx, maxz);
    if (maxy > maxz)
        std::swap(maxy, maxz);
    else if (maxy < maxx)
        std::swap(maxx, maxy);

    // Protect against underflow in the computation of eps.
    if (maxx < 1e-97) {
        if (maxx == 0.0) return 0;
    // Protect against overflow in the computation of det.
    } else if (maxz < 1e102) {
        double eps = 5.1107127829973299e-15 * maxx * maxy * maxz;

        double det = 0.0;
        Determinant(pqx, pqy, pqz, prx, pry, prz, psx, psy, psz, &det);

        if (det > eps)  return +1;
        if (det < -eps) return -1;
    }

    // Interval filter.
    if (fegetround() != FE_UPWARD) fesetround(FE_UPWARD);
    IntervalFloat det1;
    OrientationDeterminant(IntervalFloat(p.x), IntervalFloat(p.y),
                           IntervalFloat(p.z),
                           IntervalFloat(q.x), IntervalFloat(q.y),
                           IntervalFloat(q.z),
                           IntervalFloat(r.x), IntervalFloat(r.y),
                           IntervalFloat(r.z),
                           IntervalFloat(s.x), IntervalFloat(s.y),
                           IntervalFloat(s.z),
                           &det1);
    if (det1.lower() > 0.0) return 1;
    if (det1.upper() < 0.0) return -1;
    if (det1.upper() == 0.0 && det1.lower() == 0.0) return 0;

    // Exact computation.
    ExactFloat det2;
    OrientationDeterminant(ExactFloat(p.x), ExactFloat(p.y), ExactFloat(p.z),
                           ExactFloat(q.x), ExactFloat(q.y), ExactFloat(q.z),
                           ExactFloat(r.x), ExactFloat(r.y), ExactFloat(r.z),
                           ExactFloat(s.x), ExactFloat(s.y), ExactFloat(s.z),
                           &det2);
    return det2.sign();
}

/**
 * Compute the determinant for in sphere test, i.e.,
 *
 * |px, py, pz, px^2 + py^2 + pz^2, 1|
 * |qx, qy, qz, qx^2 + qy^2 + qz^2, 1|
 * |rx, ry, rz, rx^2 + ry^2 + rz^2, 1|
 * |sx, sy, sz, sx^2 + sy^2 + sz^2, 1|
 * |tx, ty, tz, tx^2 + ty^2 + tz^2, 1|
 */
template <typename T>
void InSphereDeterminant(const T& px, const T& py, const T& pz,
                         const T& qx, const T& qy, const T& qz,
                         const T& rx, const T& ry, const T& rz,
                         const T& sx, const T& sy, const T& sz,
                         const T& tx, const T& ty, const T& tz,
                         T* det) {
    T ptx = px - tx, pty = py - ty, ptz = pz - tz, qtx = qx - tx, qty = qy - ty,
      qtz = qz - tz, rtx = rx - tx, rty = ry - ty, rtz = rz - tz, stx = sx - tx,
      sty = sy - ty, stz = sz - tz;

    // pt2 = ptx * ptx + pty * pty + ptz * ptz.
    // qt2 = qtx * qtx + qty * qty + qtz * qtz.
    // rt2 = rtx * rtx + rty * rty + rtz * rtz.
    // st2 = stx * stx + sty * sty + stz * stz.
    T pt2 = ptx * ptx, qt2 = qtx * qtx, rt2 = rtx * rtx, st2 = stx * stx;
    pt2 += pty * pty, qt2 += qty * qty, rt2 += rty * rty, st2 += sty * sty;
    pt2 += ptz * ptz, qt2 += qtz * qtz, rt2 += rtz * rtz, st2 += stz * stz;

    Determinant(ptx, pty, ptz, pt2, rtx, rty, rtz, rt2,
                qtx, qty, qtz, qt2, stx, sty, stz, st2, &det);
}

/**
 * In 3D sphere test.
 *
 * The points p, q, r, and s must be ordered so that they have a positive
 * orientation (as defined by Orientation()), or the sign of the result
 * will be reversed.
 *
 * Return +1 if the point t lies inside the sphere passing through p, q, r, s.
 *        -1 if s lies outside;
 *         0 if five points are co-spherical.
 */
template <typename T>
int InShpere(const Point3D<T>& p, const Point3D<T>& q,
             const Point3D<T>& r, const Point3D<T>& s,
             const Point3D<T>& t) {
    // Static filter, adapted from CGAL.
    double px = static_cast<double>(p.x);
    double py = static_cast<double>(p.y);
    double pz = static_cast<double>(p.z);
    double qx = static_cast<double>(q.x);
    double qy = static_cast<double>(q.y);
    double qz = static_cast<double>(q.z);
    double rx = static_cast<double>(r.x);
    double ry = static_cast<double>(r.y);
    double rz = static_cast<double>(r.z);
    double sx = static_cast<double>(s.x);
    double sy = static_cast<double>(s.y);
    double sz = static_cast<double>(s.z);
    double tx = static_cast<double>(t.x);
    double ty = static_cast<double>(t.y);
    double tz = static_cast<double>(t.z);

    double ptx = px - tx;
    double pty = py - ty;
    double ptz = pz - tz;
    double pt2 = ptx * ptx + pty * pty + ptz * ptz;
    double qtx = qx - tx;
    double qty = qy - ty;
    double qtz = qz - tz;
    double qt2 = qtx * qtx + qty * qty + qtz * qtz;
    double rtx = rx - tx;
    double rty = ry - ty;
    double rtz = rz - tz;
    double rt2 = rtx * rtx + rty * rty + rtz * rtz;
    double stx = sx - tx;
    double sty = sy - ty;
    double stz = sz - tz;
    double st2 = stx * stx + sty * sty + stz * stz;

    // Compute the semi-static bound.
    double maxx = std::fabs(ptx);
    double maxy = std::fabs(pty);
    double maxz = std::fabs(ptz);
    double aqtx = std::fabs(qtx);
    double artx = std::fabs(rtx);
    double astx = std::fabs(stx);
    double aqty = std::fabs(qty);
    double arty = std::fabs(rty);
    double asty = std::fabs(sty);
    double aqtz = std::fabs(qtz);
    double artz = std::fabs(rtz);
    double astz = std::fabs(stz);

    if (maxx < aqtx) maxx = aqtx;
    if (maxx < artx) maxx = artx;
    if (maxx < astx) maxx = astx;

    if (maxy < aqty) maxy = aqty;
    if (maxy < arty) maxy = arty;
    if (maxy < asty) maxy = asty;

    if (maxz < aqtz) maxz = aqtz;
    if (maxz < artz) maxz = artz;
    if (maxz < astz) maxz = astz;

    // Sort maxx < maxy < maxz.
    if (maxx > maxz) std::swap(maxx, maxz);
    if (maxy > maxz)
        std::swap(maxy, maxz);
    else if (maxy < maxx)
        std::swap(maxx, maxy);

    // Protect against underflow in the computation of eps.
    if (maxx < 1e-58) {
        if (maxx == 0.0) return 0;
    // Protect against overflow in the computation of det.
    } else if (maxz < 1e61) {
        double eps = 1.2466136531027298e-13 * maxx * maxy * maxz * maxz * maxz;
        double det = 0.0;
        Determinant(ptx, pty, ptz, pt2, rtx, rty, rtz, rt2,
                    qtx, qty, qtz, qt2, stx, sty, stz, st2, &det);
        if (det > eps)  return +1;
        if (det < -eps) return -1;
    }

    // Interval filter.
    if (fegetround() != FE_UPWARD) fesetround(FE_UPWARD);
    IntervalFloat det1;
    InSphereDeterminant(IntervalFloat(p.x), IntervalFloat(p.y),
                        IntervalFloat(p.z),
                        IntervalFloat(q.x), IntervalFloat(q.y),
                        IntervalFloat(q.z),
                        IntervalFloat(r.x), IntervalFloat(r.y),
                        IntervalFloat(r.z),
                        IntervalFloat(s.x), IntervalFloat(s.y),
                        IntervalFloat(s.z),
                        IntervalFloat(t.x), IntervalFloat(t.y),
                        IntervalFloat(t.z),
                        &det1);
    if (det1.lower() > 0.0) return 1;
    if (det1.upper() < 0.0) return -1;
    if (det1.upper() == 0.0 && det1.lower() == 0.0) return 0;

    // Exact computation.
    ExactFloat det;
    InSphereDeterminant(ExactFloat(p.x), ExactFloat(p.y), ExactFloat(p.z),
                        ExactFloat(q.x), ExactFloat(q.y), ExactFloat(q.z),
                        ExactFloat(r.x), ExactFloat(r.y), ExactFloat(r.z),
                        ExactFloat(s.x), ExactFloat(s.y), ExactFloat(s.z),
                        ExactFloat(t.x), ExactFloat(t.y), ExactFloat(t.z),
                        &det);
    return det.sign();
}

} // namespace geometry
} // namespace cl

#endif // GEOMETRY_KERNEL_PREDICATE_3D_H_
