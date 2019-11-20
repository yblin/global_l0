//
// Copyright 2017 Yangbin Lin. All Rights Reserved.
//
// Author: yblin@jmu.edu.cn (Yangbin Lin)
//
// This file is part of the Code Library.
//

#ifndef STATISTICS_REGRESSION_RIDGE_REGRESSION_H_
#define STATISTICS_REGRESSION_RIDGE_REGRESSION_H_

#include "codelibrary/math/matrix/matrix.h"
#include "codelibrary/math/matrix/solver/cgnr_solver.h"

namespace cl {
namespace statistics {

/**
 * To solve the problem Ax = b. The standard approach is ordinary least squares
 * linear regression. However, if no x satisfies the equation or more than one x
 * does, that is, the solution is not unique. The problem is said to be ill 
 * posed.
 * 
 * In such cases, ordinary least squares estimation leads to an overdetermined
 * (over-fitted), or more often an underdetermined (under-fitted) system of
 * equations.
 *
 * Ridge regression addresses some of the problems of Linear Least Squares by
 * imposing a penalty on the size of coefficients. The ridge coefficients
 * minimize a penalized residual sum of squares:
 *
 *   min ||Xw - y||^2 + lambda ||w||^2
 *    w
 *
 * Here, lambda >= 0 is a complexity parameter that controls the amount of 
 * shrinkage: the larger the value of lambda, the greater the amount of 
 * shrinkage and thus the coefficients become more robust to collinearity.
 *
 * According to [1], there always exists a lambda such that the MSE of ridge 
 * linear model is less than the MSE of ordinary least squares linear model.
 *
 * Reference:
 *  [1] Arthur E. Hoerl, Robert W. Kennard. Ridge Regression: Biased Estimation
 *      for Nonorthogonal Problems. Technometric, 1970, 12(1):55-67.
 */
class RidgeRegression {
public:
    using TrainingResult = matrix::LinearSolver::Result;

    RidgeRegression() = default;

    RidgeRegression(const RidgeRegression&) = delete;

    RidgeRegression& operator= (const RidgeRegression&) = delete;

    /**
     * Set options for linear solver.
     */
    void set_options(const matrix::LinearSolver::Options& option) {
        options_ = option;
    }

    /**
     * Fit a linear model to the input data.
     */
    template <typename Matrix>
    TrainingResult Train(const Matrix& data, const RVector& y, double lambda) {
        int n = data.n_rows();
        assert(y.size() == n);
        assert(lambda >= 0.0);

        RVector v(n, 1.0);

        Matrix a = Concatenate(v, data);
        optimization::CGNRSolver cgnr_solver(options_);
        RVector d(data.n_columns() + 1, lambda);
        return cgnr_solver.Solve(a, d, y, &linear_model_);
    }

    /**
     * Predicate the y value of the given variable.
     */
    double Predicate(const RVector& x) const {
        assert(x.size() + 1 == linear_model().size());

        RVector b(x.size() + 1);
        b[0] = 1.0;
        for (int i = 0; i < x.size(); ++i) {
            b[i + 1] = x[i];
        }

        return b * linear_model_;
    }

    /**
     * Return the fitted linear model.
     */
    const RVector& linear_model() const {
        return linear_model_;
    }

private:
    // The fitted linear model: y = w0 + w1 * x1 + w2 * x2 + ... + wn * xn.
    RVector linear_model_;

    // Options for linear solver.
    matrix::LinearSolver::Options options_;
};

} // namespace statistics
} // namespace cl

#endif // STATISTICS_REGRESSION_RIDGE_REGRESSION_H_
