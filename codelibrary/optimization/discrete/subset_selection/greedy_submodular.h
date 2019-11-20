//
// Copyright 2018 Yangbin Lin. All Rights Reserved.
//
// Author: yblin@jmu.edu.cn (Yangbin Lin)
//
// This file is part of the Code Library.
//

#ifndef OPTIMIZATION_DISCRETE_SUBSET_SELECTION_GREEDY_SUBMODULAR_H_
#define OPTIMIZATION_DISCRETE_SUBSET_SELECTION_GREEDY_SUBMODULAR_H_

#include <algorithm>
#include <cfloat>
#include <random>

#include "codelibrary/base/algorithm.h"
#include "codelibrary/base/array.h"
#include "codelibrary/math/matrix/matrix.h"

namespace cl {
namespace optimization {

/**
 * Find the representives by deterministic submodular method. It ensures that
 * the energy of the obtained solution is atmost 3 times than optimal solution.
 *
 * Reference:
 *   Elhamifar E, Kaluza M C D P. Online Summarization via Submodular and Convex
 *   Optimization[C]. IEEE Conference on Computer Vision and Pattern
 *   Recognition. IEEE, 2017:1818-1826.
 *
 * Parameters:
 *  distances       - the pair-wise dissimilarity distances.
 *  lambda          - parameter for regularization.
 *  representatives - the output indices of representatives.
 */
void GreedySubmodularSubsetSelection(const RMatrix& distances, double lambda,
                                     Array<int>* representatives) {
    assert(!distances.empty());
    assert(distances.n_rows() == distances.n_columns());
    assert(lambda >= 0.0);
    assert(representatives);

    representatives->clear();

    int n = distances.n_rows();

    Array<int> x;
    int index = 0;
    x.push_back(index);

    Array<double> x_dis(n);
    for (int i = 0; i < n; ++i) {
        x_dis[i] = distances(i, index);
    }
    Array<double> x_dis1 = x_dis;

    Array<Array<int>> y_dis(n);
    for (int i = 0; i < n; ++i) {
        Array<int> seq;
        IndexSort(distances.begin() + i * n,
                  distances.begin() + (i + 1) * n,
                  &seq);
        for (int j = n - 1; j >= 0; --j) {
            y_dis[i].push_back(seq[j]);
        }
    }
    Array<Array<int>> y_dis1 = y_dis;

    Array<bool> is_marked(n, false);
    for (int index = 1; index < n; ++index) {
        is_marked[index] = true;

        // Compute a = F(X U {i}) - F(X)
        double s1 = std::accumulate(x_dis.begin(), x_dis.end(), 0.0);
        for (int i = 0; i < n; ++i) {
            x_dis[i] = std::min(distances(i, index), x_dis[i]);
        }
        double s2 = std::accumulate(x_dis.begin(), x_dis.end(), 0.0);
        double a = -lambda + (s1 - s2);

        // Compute b = F(Y - {i}) - F(Y)
        double d1 = 0.0;
        for (int i = 0; i < n; ++i) {
            d1 += distances(i, y_dis[i].back());
        }
        for (int i = 0; i < n; ++i) {
            while (is_marked[y_dis[i].back()]) {
                y_dis[i].pop_back();
            }
        }
        double d2 = 0.0;
        for (int i = 0; i < n; ++i) {
            d2 += distances(i, y_dis[i].back());
        }
        double b = (d1 - d2) + lambda;

        if (a >= b) {
            // X = X U {x_i}, Y = Y
            x.push_back(index);
            x_dis1 = x_dis;
            y_dis = y_dis1;
            is_marked[index] = false;
        } else {
            // Y = Y / {y_i}, X = X
            x_dis = x_dis1;
            y_dis1 = y_dis;
        }
    }

    for (int p : x) {
        representatives->push_back(p);
    }
    if (representatives->empty()) representatives->push_back(index);
}

} // namespace optimization
} // namespace cl

#endif // OPTIMIZATION_DISCRETE_SUBSET_SELECTION_GREEDY_SUBMODULAR_H_
