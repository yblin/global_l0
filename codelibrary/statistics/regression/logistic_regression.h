//
// Copyright 2019 Yangbin Lin. All Rights Reserved.
//
// Author: yblin@jmu.edu.cn (Yangbin Lin)
//
// This file is part of the Code Library.
//

#ifndef STATISTICS_REGRESSION_LOGISTIC_REGRESSION_H_
#define STATISTICS_REGRESSION_LOGISTIC_REGRESSION_H_

#include <cmath>
#include <memory>

#include "codelibrary/base/array.h"
#include "codelibrary/optimization/unconstrained/minimizer.h"
#include "codelibrary/optimization/unconstrained/separable_problem.h"

namespace cl {
namespace statistics {

/**
 * In statistics, the logistic model (or logit model) is a widely used
 * statistical model that, uses a logistic function to model a binary dependent
 * variable.
 *
 * The logistic function is a sigmoid function, which takes any real input t,
 * (t \in R), and outputs a value between zero and one. The logistic function
 * g(t) is defined as follows:
 *
 *   g(t) = e^t / (e^t + 1) = 1 / (1 + e^-t)
 *
 * Logistic regression is estimating the parameters of a logistic model; it is a
 * form of binomial regression. Let us assume that t is a linear function of x.
 * We can then express t as follows:
 *
 *   t = b_0 + bx
 *
 * And the logistic function can now be written as:
 *
 *   g(x) = 1 / (1 + e^-(b_0 + bx))
 *
 * Note that g(x) is interpreted as the probability of the dependent variable
 * equaling a "success" or "case" rather than a failure or non-case.
 *
 * Cost function:
 *
 * Instead of Mean Squared Error, we use a cost function called Cross-Entropy,
 * also known as Log Loss:
 *
 *   C(x) = -1/m (y * log(g(x)) + (1 - y) * log(1 - g(x)))
 *
 * where, m is the number of training data.
 *
 * The derivative of cost:
 *
 *   C'(x) = -1/m (y * 1/g(x) g'(x) - (1 - y) * (1 / (1 - g(x)) g'(x)))
 *
 * and the derivative of g(t) equal to:
 *
 *   g'(x) = g(x)(1 - g(x)) * x
 *
 * Thus,
 *
 *   C'(x) = -1/m (y * (1 - g(x)) * x - (1 - y) * g(x) * x)
 *         = -1/m (y * x - g(x) * x)
 *         =  1/m (g(x) - y) * x
 *
 * which is beautiful and convenient.
 */
class LogisticRegression {
    /**
     * Return 1 / (1 + e^-x)
     */
    static double Logistic(double x) {
        return 1.0 / (1.0 + std::exp(-x));
    }

    class CostFunction : public optimization::SeparableProblem {
        friend class LogisticRegression;

    public:
        CostFunction(int n_parameters, int n_functions)
            : SeparableProblem(n_parameters, n_functions) {}

        /**
         * Evaluate the i-th objective function value and the i-th gradient.
         *
         * f = -y * log(g(x)) - (1 - y) * log(1 - g(x))
         * g = (g(x) - y) * x
         */
        virtual void Evaluate(const RVector& t, int i, double* f_i,
                              RVector* g_i) override {
            double tt = x_[i] * t;
            double gt = Logistic(tt);

            // Avoid overflow.
            double t1 = (gt < 1e-6) ? tt : std::log(gt);
            double t2 = (gt > 1 - 1e-6) ? -tt : std::log(1.0 - gt);

            *f_i = -y_[i] * t1 - (1.0 - y_[i]) * t2;
            *g_i = (gt - y_[i]) * x_[i];
        }

    private:
        Array<RVector> x_;
        Array<double> y_;
    };

public:
    LogisticRegression() = default;

    ~LogisticRegression() = default;

    LogisticRegression(const LogisticRegression&) = delete;

    LogisticRegression& operator=(const LogisticRegression&) = delete;

    /**
     * Training logistic regression.
     *
     * Parameters:
     *  x      - the input features.
     *  y      - the label for each data.
     *  solver - the solver to minimize the problem.
     * 
     * Return:
     *  the convergence result of the minimizer.
     */
    template <class Solver>
    typename Solver::Result
    Train(const Array<RVector>& x, const Array<int>& y, Solver* solver) {
        assert(x.size() == y.size());
        assert(!x.empty());
        for (int t : y)
            assert(t == 0 || t == 1);

        // Note that the parameters is n + 1, with one intercept term.
        int n_parameters = x.front().size() + 1;
        int n_data = x.size();

        cost_function_.reset(new CostFunction(n_parameters, n_data));
        for (RVector t : x) {
            // 'tx + t0' can be simplified to 'tx' by adding an extra '1' to x.
            t.resize(t.size() + 1, 1.0);
            cost_function_->x_.push_back(t);
        }
        cost_function_->y_.insert(y.begin(), y.end());

        model_.resize(n_parameters, 1.0);
        return solver->Minimize(cost_function_.get(), &model_);
    }

    /**
     * Predicate the y value of the given data.
     *
     * Return the probability of 'x' equaling a "success" label.
     */
    double Predicate(const RVector& x) const {
        assert(cost_function_.get());
        assert(x.size() == cost_function_->n_parameters() - 1);

        RVector t(x.size() + 1);
        for (int i = 0; i < x.size(); ++i) {
            t[i] = x[i];
        }
        t[x.size()] = 1.0;

        return Logistic(model_ * t);
    }

    /**
     * Set an guessed or pre-trained model.
     */
    void InitializeModel(const RVector& model) {
        model_ = model;
    }

    /**
     * Return the trained model.
     */
    const RVector& model() const {
        return model_;
    }

private:
    // Cost function of the logistic regression problem.
    std::unique_ptr<CostFunction> cost_function_;

    // The parameters of trained model.
    RVector model_;
};

} // namespace statistics
} // namespace cl

#endif // STATISTICS_REGRESSION_LOGISTIC_REGRESSION_H_
