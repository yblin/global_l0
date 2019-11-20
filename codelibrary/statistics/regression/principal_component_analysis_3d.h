//
// Copyright 2015 Yangbin Lin. All Rights Reserved.
//
// Author: yblin@jmu.edu.cn (Yangbin Lin)
//
// This file is part of the Code Library.
//

#ifndef STATISTICS_REGRESSION_PRINCIPAL_COMPONENT_ANALYSIS_3D_H_
#define STATISTICS_REGRESSION_PRINCIPAL_COMPONENT_ANALYSIS_3D_H_

#include <cassert>

#include "codelibrary/geometry/kernel/center_3d.h"
#include "codelibrary/geometry/kernel/plane_3d.h"
#include "codelibrary/geometry/kernel/point_3d.h"
#include "codelibrary/math/matrix/eigen/symmetric_eigen.h"

namespace cl {
namespace statistics {

/**
 * Principal Component Analysis for 3D data.
 *
 * Principal Component Analysis (PCA) is an orthogonal linear transformation
 * that transforms the data to a new coordinate system such that the greatest
 * variance by some projection of the data comes to lie on the first coordinate
 * (called the first principal component), the second greatest variance on the
 * second coordinate, and so on.
 */
class PrincipalComponentAnalysis3D {
public:
    /**
     * Compute PCA by covariance method.
     */
    template <typename T>
    PrincipalComponentAnalysis3D(const Array<Point3D<T>>& points) {
        assert(!points.empty());

        // Get the covariance matrix.
        centroid_ = Centroid(points);
        double a00 = 0.0, a01 = 0.0, a02 = 0.0, a11 = 0.0, a12 = 0.0, a22 = 0.0;
        for (const Point3D<T>& p : points) {
            double x = p.x - centroid_.x;
            double y = p.y - centroid_.y;
            double z = p.z - centroid_.z;

            a00 += x * x;
            a01 += x * y;
            a02 += x * z;
            a11 += y * y;
            a12 += y * z;
            a22 += z * z;
        }

        RMatrix mat(3, 3);
        mat(0, 0) = a00;
        mat(0, 1) = mat(1, 0) = a01;
        mat(0, 2) = mat(2, 0) = a02;
        mat(1, 1) = a11;
        mat(1, 2) = mat(2, 1) = a12;
        mat(2, 2) = a22;

        // Get the eigenvectors of covariance matrix.
        matrix::SymmetricEigen3(mat, &eigenvalues_, &eigenvectors_);
    }

    /**
     * Compute weighted PCA by covariance method.
     */
    template <typename T>
    PrincipalComponentAnalysis3D(const Array<Point3D<T>>& points,
                                 const Array<double>& weights) {
        assert(!points.empty());
        assert(points.size() == weights.size());

        // Get the covariance matrix.
        centroid_ = Centroid(points, weights);
        double a00 = 0.0, a01 = 0.0, a02 = 0.0, a11 = 0.0, a12 = 0.0, a22 = 0.0;
        int i = 0;
        for (const Point3D<T>& p : points) {
            double x = p.x - centroid_.x;
            double y = p.y - centroid_.y;
            double z = p.z - centroid_.z;
            double w = weights[i++];

            a00 += w * x * x;
            a01 += w * x * y;
            a02 += w * x * z;
            a11 += w * y * y;
            a12 += w * y * z;
            a22 += w * z * z;
        }

        RMatrix mat(3, 3);
        mat(0, 0) = a00;
        mat(0, 1) = mat(1, 0) = a01;
        mat(0, 2) = mat(2, 0) = a02;
        mat(1, 1) = a11;
        mat(1, 2) = mat(2, 1) = a12;
        mat(2, 2) = a22;

        // Get the eigenvectors of covariance matrix.
        matrix::SymmetricEigen3(mat, &eigenvalues_, &eigenvectors_);
    }

    PrincipalComponentAnalysis3D(const PrincipalComponentAnalysis3D&) = delete;

    PrincipalComponentAnalysis3D&
    operator=(const PrincipalComponentAnalysis3D&) = delete;

    const RMatrix& covariance_matrix() const {
        return covariance_matrix_;
    }

    const RVector& eigenvalues() const {
        return eigenvalues_;
    }

    const RMatrix& eigenvectors() const {
        return eigenvectors_;
    }

    const RPoint3D& centroid() const {
        return centroid_;
    }

    /**
     * Convert 3D PCA to 3D line.
     *
     * Return the fitting quality: 0 is best, 1 is worst.
     */
    double ToLine(RLine3D* line) const {
        assert(line);

        if (eigenvalues_[2] <= eigenvalues_[1]) {
            // Return a degenerated line.
            *line = RLine3D(centroid_, RVector3D(0.0, 0.0, 0.0));
            return 1.0;
        }

        RVector3D direction(eigenvectors_(0, 2), eigenvectors_(1, 2),
                            eigenvectors_(2, 2));
        *line = RLine3D(centroid_, direction);
        return eigenvalues_[1] / eigenvalues_[2];
    }

    /**
     * Convert 3D PCA to 3D plane.
     *
     * Return the fitting quality: 0 is best, 1 is worst.
     */
    double ToPlane(RPlane3D* plane) const {
        assert(plane);

        if (eigenvalues_[2] <= eigenvalues_[1] &&
            eigenvalues_[1] <= eigenvalues_[0]) {
            // Return a degenerated plane.
            *plane = RPlane3D(centroid_, RVector3D(0.0, 0.0, 0.0));
            return 1.0;
        }

        if (eigenvectors_.get_column(0).squared_norm() == 0.0) {
            RVector3D v1(eigenvectors_(0, 1), eigenvectors_(1, 1),
                         eigenvectors_(2, 1));
            RVector3D v2(eigenvectors_(0, 2), eigenvectors_(1, 2),
                         eigenvectors_(2, 2));
            RVector3D v = CrossProduct(v1, v2);
            *plane = RPlane3D(centroid_, v);
        } else {
            RVector3D v(eigenvectors_(0, 0), eigenvectors_(1, 0), 
                        eigenvectors_(2, 0));
            *plane = RPlane3D(centroid_, v);
        }
        return eigenvalues_[0] / eigenvalues_[1];
    }

private:
    // The centroid point of the input data.
    RPoint3D centroid_;

    // Covariance matrix of the input data.
    RMatrix covariance_matrix_;

    // Eigenvalues of the covariance matrix, sorted in ascending order. 
    RVector eigenvalues_;

    // Corresponding eigenvectors.
    RMatrix eigenvectors_;
};

} // namespace statistics
} // namespace cl

#endif // STATISTICS_REGRESSION_PRINCIPAL_COMPONENT_ANALYSIS_3D_H_
