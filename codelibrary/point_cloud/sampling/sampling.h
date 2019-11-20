//
// Copyright 2019 Yangbin Lin. All Rights Reserved.
//
// Author: yblin@jmu.edu.cn (Yangbin Lin)
//
// This file is part of the Code Library.
//

#ifndef POINT_CLOUD_SAMPLING_H_
#define POINT_CLOUD_SAMPLING_H_

#include <algorithm>
#include <cassert>
#include <numeric>
#include <random>

#include "codelibrary/base/algorithm.h"
#include "codelibrary/base/array.h"

namespace cl {
namespace point_cloud {

/**
 * The basic sampling class for point clouds.
 */
class Sampling {
public:
    Sampling(int random_seed = 0)
        : random_(random_seed) {}

    void set_random_seed(int random_seed) {
        random_.seed(random_seed);
    }

    /**
     * Sampling the point clouds [first, last).
     * We store the indices of the sampled points.
     *
     * The default output samples is equal to the input with shuffled sequence.
     */
    template <typename Iterator>
    void Sample(Iterator first, Iterator last, Array<int>* samples) {
        Generate(first, last, samples);
    }

    void Sample(const Array<RPoint3D>& points, Array<int>* samples) {
        Sample(points.begin(), points.end(), samples);
    }

protected:
    /**
     * Generate a random sequence with (last - first) elements.
     */
    template <typename Iterator>
    void Generate(Iterator first, Iterator last, Array<int>* sequence) {
        assert(sequence);

        int n = CountElements(first, last);
        sequence->resize(n);
        std::iota(sequence->begin(), sequence->end(), 0);
        std::shuffle(sequence->begin(), sequence->end(), random_);
    }

    std::mt19937 random_;
};

} // namespace point_cloud
} // namespace cl

#endif // POINT_CLOUD_SAMPLING_H_
