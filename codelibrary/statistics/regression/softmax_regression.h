//
// Copyright 2019 Yangbin Lin. All Rights Reserved.
//
// Author: yblin@jmu.edu.cn (Yangbin Lin)
//
// This file is part of the Code Library.
//

#ifndef STATISTICS_REGRESSION_SOFTMAX_REGRESSION_H_
#define STATISTICS_REGRESSION_SOFTMAX_REGRESSION_H_

#include <cmath>
#include <memory>
#include <unordered_map>

#include "codelibrary/base/algorithm.h"
#include "codelibrary/optimization/unconstrained/lbfgs.h"
#include "codelibrary/optimization/unconstrained/separable_problem.h"

namespace cl {
namespace statistics {

/**
 * Softmax regression model generalizes logistic regression to classification
 * problems where the class label y can take on more than two possible values.
 * This will be useful for such problems as MNIST digit classification, where
 * the goal is to distinguish between 10 different numerical digits.
 *
 * The cost function of Softmax regression is:
 *                 m    k
 *   C(x) = -1/m  Sum [Sum 1{y^(i) = j} log p(y^(i) = j | x;t)]
 *                i=1  j=1
 *
 * Where, 1{.} is the indicator function, so that 1{a true statement} = 1, and
 * 1{a false statement} = 0. And, we have:
 *
 *                           e^(t_j x^(i))
 *   p(y^(i) = j | x;t) = --------------------
 *                         Sum_l(e^(t_l x^(i))
 *
 * Its gradient is:
 *
 *                  m
 *   g_j(x) = -1/m Sum [x^(i)(1{y^(i)=j} - p(y^(i)=j | x; t)]
 *                 i=1
 *
 * Note that, the above softmax model is overparameterized. Thus, We always
 * replace t_1 with the vector of all 0's.
 */
class SoftmaxRegression {
    /**
     * Softmax function to compute the probability for each label.
     */
    static void Softmax(RVector* p) {
        assert(p);

        double max_x = *std::max_element(p->begin(), p->end());
        double s = 0.0;
        for (int j = 0; j < p->size(); ++j) {
            (*p)[j] = std::exp((*p)[j] - max_x);
            s += (*p)[j];
        }
        *p *= 1.0 / s;
    }

public:
    class CostFunction : public optimization::SeparableProblem {
        friend class SoftmaxRegression;

    public:
        CostFunction(int n_parameters, int n_functions, int n_labels)
            : SeparableProblem(n_parameters, n_functions),
              n_labels_(n_labels) {}

        /**
         * Evaluate the i-th objective function value.
         */
        virtual void Evaluate(const RVector& parameter, int i, double* f_i,
                              RVector* g_i) override {
            int n = x_[i].size();

            // Note that, p[0] is always 0.
            RVector tx(n_labels_, 0.0);
            for (int j = 1; j < n_labels_; ++j) {
                for (int k = 0; k < n; ++k) {
                    tx[j] += parameter[(j - 1) * n + k] * x_[i][k];
                }
            }

            double max_tx = *std::max_element(tx.begin(), tx.end());
            double sum = 0.0;
            RVector p(n_labels_);
            for (int i = 0; i < n_labels_; ++i) {
                p[i] = std::exp(tx[i] - max_tx);
                sum += p[i];
            }
            p *= 1.0 / sum;

            // More accurate than -log(p/sum).
            *f_i = std::log(sum) - (tx[y_[i]] - max_tx);

            g_i->assign(n_parameters_, 0.0);
            for (int j = 1; j < n_labels_; ++j) {
                for (int k = 0; k < n; ++k) {
                    if (y_[i] == j) {
                        (*g_i)[(j - 1) * n + k] = (p[j] - 1.0) * x_[i][k];
                    } else {
                        (*g_i)[(j - 1) * n + k] = p[j] * x_[i][k];
                    }
                }
            }
        }

    private:
        Array<RVector> x_;
        Array<int> y_;
        int n_labels_;
    };

public:
    /**
     * Construct logistic regression by the given feature numbers.
     */
    SoftmaxRegression() = default;

    ~SoftmaxRegression() = default;

    SoftmaxRegression(const SoftmaxRegression&) = delete;

    SoftmaxRegression& operator=(const SoftmaxRegression&) = delete;

    /**
     * Training Softmax regression.
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

        // Normalize labels.
        Array<int> labels;
        NormalizeLabels(y, &labels);

        // Note that the parameters is (k - 1) * (n + 1).
        int n_parameters = (x.front().size() + 1) * (n_labels_ - 1);
        int n_data = x.size();

        cost_function_.reset(new CostFunction(n_parameters, n_data, n_labels_));
        for (RVector t : x) {
            // 'tx + t0' can be simplified to 'tx' by adding an extra '1' to x.
            t.resize(x.size() + 1, 1.0);
            cost_function_->x_.push_back(t);
        }
        cost_function_->y_.swap(labels);
        model_.resize(n_parameters, 0.0);
        return solver->Minimize(cost_function_.get(), &model_);
    }

    /**
     * Predicate the label of the given data.
     */
    int Predicate(const RVector& x) const {
        assert(cost_function_.get());
        assert(x.size() + 1 ==
               cost_function_->n_parameters() / (n_labels_ - 1));

        int n = x.size() + 1;

        cl::RVector t = x;
        t.resize(t.size() + 1, 1.0);

        double max = 0.0;
        int best = 0;
        for (int j = 1; j < n_labels_; ++j) {
            double tx = 0.0;
            for (int k = 0; k < n; ++k) {
                tx += model_[(j - 1) * n + k] * t[k];
            }
            if (tx > max) {
                max = tx;
                best = j;
            }
        }

        return hash_.at(best);
    }

    /**
     * Predicate the most likely k labels with probabilities of the given data.
     */
    void Predicate(const RVector& x, int k,
                   Array<std::pair<int, double> >* labels) const {
        assert(cost_function_.get());
        assert(x.size() + 1 ==
               cost_function_->n_parameters() / (n_labels_ - 1));
        assert(k > 0);
        assert(k <= n_labels_);
        assert(labels);

        int n = x.size() + 1;

        RVector t(x.size() + 1);
        for (int i = 0; i < x.size(); ++i) {
            t[i] = x[i];
        }
        t[x.size()] = 1.0;

        RVector p(n_labels_, 0.0);
        for (int j = 1; j < n_labels_; ++j) {
            for (int k = 0; k < n; ++k) {
                p[j] += model_[(j - 1) * n + k] * t[k];
            }
        }
        Softmax(&p);

        Array<int> seq;
        IndexSort(p.begin(), p.end(), &seq);
        std::reverse(seq.begin(), seq.end());
        labels->clear();
        for (int i = 0; i < k; ++i) {
            int j = seq[i];
            labels->emplace_back(hash_.at(j), p[j]);
        }
    }

    /**
     * Set an guessed or pretrained model.
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
    /**
     * The input Y label may be not continuous. We need to normalize them.
     */
    void NormalizeLabels(const Array<int>& y, Array<int>* labels) {
        std::unordered_map<int, int> map;
        hash_.clear();
        n_labels_ = 0;
        for (int i : y) {
            if (map.find(i) == map.end()) {
                map[i] = n_labels_;
                hash_[n_labels_] = i;
                labels->push_back(n_labels_);
                ++n_labels_;
            } else {
                labels->push_back(map[i]);
            }
        }
    }

    // Number of labels.
    int n_labels_;

    // Cost function of Logistic regression.
    std::unique_ptr<CostFunction> cost_function_;

    // The parameters of trained model.
    RVector model_;

    // Hash for input labels.
    std::unordered_map<int, int> hash_;
};

} // namespace statistics
} // namespace cl

#endif // STATISTICS_REGRESSION_SOFTMAX_REGRESSION_H_
