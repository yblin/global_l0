//
// Copyright 2015 Yangbin Lin. All Rights Reserved.
//
// Author: yblin@jmu.edu.cn (Yangbin Lin)
//
// This file is part of the Code Library.
//

#ifndef STATISTICS_REGRESSION_PRINCIPAL_COMPONENT_ANALYSIS_2D_H_
#define STATISTICS_REGRESSION_PRINCIPAL_COMPONENT_ANALYSIS_2D_H_

#include <cassert>

#include "codelibrary/geometry/kernel/center_2d.h"
#include "codelibrary/geometry/kernel/line_2d.h"
#include "codelibrary/geometry/kernel/point_2d.h"
#include "codelibrary/math/matrix/eigen/symmetric_eigen.h"

namespace cl {
namespace statistics {

/**
 * Principal component analysis for 2D points.
 *
 * Principal Component Analysis (PCA) is an orthogonal linear transformation
 * that transforms the data to a new coordinate system such that the greatest
 * variance by some projection of the data comes to lie on the first coordinate
 * (called the first principal component), the second greatest variance on the
 * second coordinate, and so on.
 */
class PrincipalComponentAnalysis2D {
public:
    /**
     * Compute PCA by covariance method.
     */
    template <typename T>
    PrincipalComponentAnalysis2D(const Array<Point2D<T>>& points) {
        assert(!points.empty());

        // Get the covariance matrix.
        centroid_ = Centroid(points);
        double a00 = 0.0, a01 = 0.0, a11 = 0.0;
        for (const Point2D<T>& p : points) {
            double x = p.x - centroid_.x;
            double y = p.y - centroid_.y;

            a00 += x * x;
            a01 += x * y;
            a11 += y * y;
        }

        RMatrix mat(2, 2);
        mat(0, 0) = a00;
        mat(0, 1) = mat(1, 0) = a01;
        mat(1, 1) = a11;

        // Get the eigenvectors of covariance matrix.
        matrix::SymmetricEigen2(mat, &eigenvalues_, &eigenvectors_);
    }

    /**
     * Compute weighted PCA by covariance method.
     */
    template <typename T>
    PrincipalComponentAnalysis2D(const Array<Point2D<T>>& points,
                                 const Array<double>& weights) {
        assert(!points.empty());
        assert(points.size() == weights.size());

        // Get the covariance matrix.
        centroid_ = Centroid(points, weights);
        double a00 = 0.0, a01 = 0.0, a11 = 0.0;
        int i = 0;
        for (const Point2D<T>& p : points) {
            double x = p.x - centroid_.x;
            double y = p.y - centroid_.y;
            double w = weights[i++];

            a00 += w * x * x;
            a01 += w * x * y;
            a11 += w * y * y;
        }

        RMatrix mat(2, 2);
        mat(0, 0) = a00;
        mat(0, 1) = mat(1, 0) = a01;
        mat(1, 1) = a11;

        // Get the eigenvectors of covariance matrix.
        matrix::SymmetricEigen2(mat, &eigenvalues_, &eigenvectors_);
    }

    PrincipalComponentAnalysis2D(const PrincipalComponentAnalysis2D&) = delete;

    PrincipalComponentAnalysis2D&
    operator=(const PrincipalComponentAnalysis2D&) = delete;

    const RMatrix& covariance_matrix() const {
        return covariance_matrix_;
    }

    const RVector& eigenvalues() const {
        return eigenvalues_;
    }

    const RMatrix& eigenvectors() const {
        return eigenvectors_;
    }

    const RPoint2D& centroid() const {
        return centroid_;
    }

    /**
     * Convert 2D PCA to 2D line.
     *
     * Return fitting quality: 0 is best, 1 is worst.
     */
    double ToLine(RLine2D* line) const {
        assert(line);

        if (eigenvalues_[1] <= eigenvalues_[0]) {
            // Return a degenerated line.
            *line = RLine2D(centroid_, RVector2D(0.0, 0.0));
            return 1.0;
        }

        RVector2D direction(eigenvectors_(0, 1), eigenvectors_(1, 1));
        *line = RLine2D(centroid_, direction);
        return eigenvalues_[0] / eigenvalues_[1];
    }

private:
    // The centroid point of the input data.
    RPoint2D centroid_;

    // Covariance matrix of the input data.
    RMatrix covariance_matrix_;

    // Eigenvalues of the covariance matrix, sorted in ascending order. 
    RVector eigenvalues_;

    // Corresponding vectors.
    RMatrix eigenvectors_;
};

} // namespace statistics
} // namespace cl

#endif // STATISTICS_REGRESSION_PRINCIPAL_COMPONENT_ANALYSIS_2D_H_
