//
// Copyright 2013 Yangbin Lin. All Rights Reserved.
//
// Author: yblin.x_mu@qq.com (Yangbin Lin)
//
// This file is part of the Code Library.
//

#ifndef GEOMETRY_KERNEL_QUATERNION_H_
#define GEOMETRY_KERNEL_QUATERNION_H_

#include <cassert>
#include <cfloat>

#include "codelibrary/base/equal.h"
#include "codelibrary/geometry/kernel/angle.h"
#include "codelibrary/geometry/kernel/vector_3d.h"
#include "codelibrary/math/matrix/matrix.h"

namespace cl {

/**
 * Quaternions provide a convenient mathematical notation for representing
 * orientations and rotations of objects in three dimensions.
 *
 * Note that, Quaternion q and -q represent the same rotation.
 */
class Quaternion {
public:
    Quaternion() = default;

    Quaternion(double x, double y, double z, double w)
        : x_(x), y_(y), z_(z), w_(w) {}

    /**
     * Construct quaternion from a rotation axis a rotation angle (in radian).
     */
    Quaternion(const RVector3D& rotation_axis, double rotation_angle) {
        double l = rotation_axis.norm();
        if (Equal(l, 0.0)) return;

        double c = std::cos(0.5 * rotation_angle);
        double s = std::sin(0.5 * rotation_angle);
        double norm = s / l;

        x_ = rotation_axis.x * norm;
        y_ = rotation_axis.y * norm;
        z_ = rotation_axis.z * norm;
        w_ = c;
    }

    /**
     * Construct quaternion from two vectors.
     * Represent rotation from v1 to v2.
     */
    Quaternion(const RVector3D& v1, const RVector3D& v2) {
        RVector3D rotation_axis = CrossProduct(v1, v2);
        double phi = Radian(v1, v2);
        *this = Quaternion(rotation_axis, phi);
    }

    /**
     * Construct quaternion from 3 X 3 rotation matrix.
     */
    explicit Quaternion(const RMatrix& rotation) {
        assert(rotation.n_rows() == 3 && rotation.n_columns() == 3);

        double trace = rotation(0, 0) + rotation(1, 1) + rotation(2, 2);

        if (trace > 0.0) {
            double t = 1.0 + trace;
            double s = 0.5 / std::sqrt(t);
            w_ = s * t;
            x_ = (rotation(2, 1) - rotation(1, 2)) * s;
            y_ = (rotation(0, 2) - rotation(2, 0)) * s;
            z_ = (rotation(1, 0) - rotation(0, 1)) * s;
        } else {
            int next[3] = { 1, 2, 0 };

            int i = 0;
            if (rotation(1, 1) > rotation(0, 0)) {
                i = 1;
            }
            if (rotation(2, 2) > rotation(i, i)) {
                i = 2;
            }
            int j = next[i];
            int k = next[j];
            double t = (rotation(i, i) - rotation(j, j) - rotation(k, k)) + 1.0;
            double s = 0.5 / std::sqrt(t);
            double q[4];
            q[i] = s * t;
            q[3] = (rotation(k, j) - rotation(j, k)) * s;
            q[j] = (rotation(j, i) + rotation(i, j)) * s;
            q[k] = (rotation(k, i) + rotation(i, k)) * s;
            x_ = q[0];
            y_ = q[1];
            z_ = q[2];
            w_ = q[3];
        }
    }

    /**
     * Construct Quaternion from Euler angles.
     */
    Quaternion(double phi, double theta, double psi) {
        double c1 = std::cos(0.5 * phi);
        double c2 = std::cos(0.5 * theta);
        double c3 = std::cos(0.5 * psi);
        double s1 = std::sin(0.5 * phi);
        double s2 = std::sin(0.5 * theta);
        double s3 = std::sin(0.5 * psi);

        w_ = c1 * c2 * c3 - s1 * s2 * s3;
        x_ = s1 * s2 * c3 + c1 * c2 * s3;
        y_ = s1 * c2 * c3 + c1 * s2 * s3;
        z_ = c1 * s2 * c3 - s1 * c2 * s3;
    }

    /**
     * Return the reverse of this quaternion.
     */
    Quaternion Reverse() const {
        return {-x_, -y_, -z_, w_};
    }

    /**
     * Normalize the quaternion into the unit quaternion.
     */
    void Normalize() {
        double norm1 = norm();
        if (norm1 == 0.0) {
            x_ = y_ = z_ = 0.0;
            w_ = 1.0;
            return;
        }

        norm1 = 1.0 / norm1;
        x_ *= norm1;
        y_ *= norm1;
        z_ *= norm1;
        w_ *= norm1;
    }

    /**
     * Norm of quaternion.
     */
    double norm() const {
        return std::sqrt(x_ * x_ + y_ * y_ + z_ * z_ + w_ * w_);
    }

    /**
     * Convert quaternion to a rotation angle(in radian) and rotation axis.
     */
    void ToRotation(double* rotation_angle, RVector3D* rotation_axis) const {
        assert(rotation_angle && rotation_axis);

        *rotation_angle = 2.0 * std::acos(w_);

        if (*rotation_angle == 0.0) {
            rotation_axis->x = rotation_axis->y = 0.0;
            rotation_axis->z = 1.0;
            return;
        }

        rotation_axis->x = x_;
        rotation_axis->y = y_;
        rotation_axis->z = z_;
        double norm = 1.0 / rotation_axis->norm();
        *rotation_axis *= norm;
    }

    /**
     * Convert quaternion to rotation matrix.
     */
    void ToRotationMatrix(RMatrix* rotation) const {
        Quaternion q = *this;
        q.Normalize();

        double sqw = q.w_ * q.w_;
        double sqx = q.x_ * q.x_;
        double sqy = q.y_ * q.y_;
        double sqz = q.z_ * q.z_;

        rotation->resize(3, 3);
        (*rotation)(0, 0) =  sqx - sqy - sqz + sqw;
        (*rotation)(1, 1) = -sqx + sqy - sqz + sqw;
        (*rotation)(2, 2) = -sqx - sqy + sqz + sqw;

        double tmp1 = q.x_ * q.y_;
        double tmp2 = q.z_ * q.w_;
        (*rotation)(1, 0) = 2.0 * (tmp1 + tmp2);
        (*rotation)(0, 1) = 2.0 * (tmp1 - tmp2);

        tmp1 = q.x_ * q.z_;
        tmp2 = q.y_ * q.w_;
        (*rotation)(2, 0) = 2.0 * (tmp1 - tmp2);
        (*rotation)(0, 2) = 2.0 * (tmp1 + tmp2);

        tmp1 = q.y_ * q.z_;
        tmp2 = q.x_ * q.w_;
        (*rotation)(2, 1) = 2.0 * (tmp1 + tmp2);
        (*rotation)(1, 2) = 2.0 * (tmp1 - tmp2);
    }

    /**
     * Convert quaternion to Euler angles (in radian form).
     */
    void ToEulerAngles(double* phi, double* theta, double* psi) {
        assert(phi && theta && psi);

        double sw = w_ * w_;
        double sx = x_ * x_;
        double sy = y_ * y_;
        double sz = z_ * z_;

        double unit = sw + sx + sy + sz;
        double test = x_ * y_ + z_ * w_;
        if (test > 0.499 * unit) {
            // Singularity at north pole.
            *phi = 2.0 * std::atan2(x_, w_);
            *theta = 0.5 * M_PI;
            *psi = 0.0;
            return;
        }

        if (test < -0.499 * unit) {
            // Singularity at south pole.
            *phi = -2.0 * std::atan2(x_, w_);
            *theta = -0.5 * M_PI;
            *psi = 0.0;
            return;
        }

        *phi   = std::atan2(2.0 * y_ * w_ - 2.0 * x_ * z_, sx - sy - sz + sw);
        *theta = std::asin(2.0 * test / unit);
        *psi   = std::atan2(2.0 * x_ * w_ - 2.0 * y_ * z_, -sx + sy - sz + sw);
    }

    /**
     * this *= lhs.
     */
    Quaternion& operator *=(const Quaternion& rhs) {
        return Multiply(*this, rhs);
    }

    /**
     * Set this quaternion to two quaternions' product, and return the product.
     * Two rotation quaternions can be combined into one equivalent quaternion
     * by multiplication.
     */
    Quaternion& Multiply(const Quaternion& a, const Quaternion& b) {
        double x = a.w_ * b.x_ + a.x_ * b.w_ + a.y_ * b.z_ - a.z_ * b.y_;
        double y = a.w_ * b.y_ + a.y_ * b.w_ + a.z_ * b.x_ - a.x_ * b.z_;
        double z = a.w_ * b.z_ + a.z_ * b.w_ + a.x_ * b.y_ - a.y_ * b.x_;
        double w = a.w_ * b.w_ - a.x_ * b.x_ - a.y_ * b.y_ - a.z_ * b.z_;
        x_ = x;
        y_ = y;
        z_ = z;
        w_ = w;
        return *this;
    }

    /**
     * Rotate the given vector by this quaternion.
     */
    template <typename T>
    RVector3D Rotate(const Vector3D<T>& v) const {
        Quaternion q(v.x, v.y, v.z, 0.0);
        q = *this * q * this->Reverse();

        return {q.x_, q.y_, q.z_};
    }

    /**
     * Return c = lhs * rhs.
     */
    friend Quaternion operator *(const Quaternion& lhs, const Quaternion& rhs) {
        return Quaternion().Multiply(lhs, rhs);
    }

    double x() const { return x_; }
    double y() const { return y_; }
    double z() const { return z_; }
    double w() const { return w_; }

private:
    double x_ = 0.0, y_ = 0.0, z_ = 0.0, w_ = 1.0;
};

} // namespace cl

#endif // GEOMETRY_KERNEL_QUATERNION_H_
