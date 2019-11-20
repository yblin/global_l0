//
// Copyright 2017 Yangbin Lin. All Rights Reserved.
//
// Author: yblin@jmu.edu.cn (Yangbin Lin)
//
// This file is part of the Code Library.
//

#ifndef STATISTICS_REGRESSION_LINEAR_REGRESSION_H_
#define STATISTICS_REGRESSION_LINEAR_REGRESSION_H_

#include "codelibrary/math/matrix/matrix.h"
#include "codelibrary/math/matrix/solver/cgnr_solver.h"

namespace cl {
namespace statistics {

/**
 * Given a data set X = {x_i1, x_i2, ..., x_ip} of n predictors and a response
 * vector Y = {y_i} as inputs, a linear regression model assumes that the
 * relationship between y_i and x_i is linear.
 *
 * The linear regression model takes the form:
 *
 *   y_i = w_0 + w_1 * x_1 + ... + w_n * x_n + e_i
 *       = [1, x_i] * w + e_i
 *
 * where e_i is error term.
 *
 * We minimize the above model by linear least squares:
 *
 *   min ||Xw - y||^2 = ||e||^2
 *    w
 */
class LinearRegression {
public:
    using TrainingResult = matrix::LinearSolver::Result;

    LinearRegression() = default;

    LinearRegression(const LinearRegression&) = delete;

    LinearRegression& operator= (const LinearRegression&) = delete;

    /**
     * Set options for linear solver.
     */
    void set_options(const matrix::LinearSolver::Options& option) {
        options_ = option;
    }

    /**
     * Fit a linear model to the input data.
     *
     * Return the training result (equal to the result of LinearSolver).
     */
    template <typename Matrix>
    TrainingResult Train(const Matrix& data, const RVector& y) {
        int n = data.n_rows();
        assert(y.size() == n);

        RVector v(n, 1.0);

        Matrix a = Concatenate(v, data);
        optimization::CGNRSolver cgnr_solver(options_);
        TrainingResult res = cgnr_solver.Solve(a, y, &linear_model_);

        error_ = y - a * linear_model_;
        return res;
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
     * Return fitted linear model.
     */
    const RVector& linear_model() const {
        return linear_model_;
    }

    /**
     * Return the error term.
     */
    const RVector& error() const {
        return error_;
    }

private:
    // The fitted linear model: y = w0 + w1 * x1 + w2 * x2 + ... + wn * xn.
    RVector linear_model_;

    // Error term: e_i = y_i - x_i * w.
    RVector error_;

    // Options for linear solver.
    matrix::LinearSolver::Options options_;
};

} // namespace statistics
} // namespace cl

#endif // STATISTICS_REGRESSION_LINEAR_REGRESSION_H_
