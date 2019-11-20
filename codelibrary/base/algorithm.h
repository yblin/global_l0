//
// Copyright 2015 Yangbin Lin. All Rights Reserved.
//
// Author: yblin@jmu.edu.cn (Yangbin Lin)
//
// This file is part of the Code Library.
//

#ifndef BASE_ALGORITHM_H_
#define BASE_ALGORITHM_H_

#include <algorithm>
#include <numeric>

#include "codelibrary/base/array.h"

namespace cl {

template <typename T>
inline T Clamp(const T& v, const T& low, const T& high) {
    return v < low ? low : (high < v ? high : v);
}

/**
 * Count the elements of [first, last).
 */
template <typename Iterator>
int CountElements(Iterator first, Iterator last) {
    auto n = std::distance(first, last);
    assert(n >= 0);
    assert(n <= INT_MAX && "We only support INT_MAX elements at most.");

    return static_cast<int>(n);
}

/**
 * Return the index of the minimal element in [first, last).
 */
template <typename Iterator>
inline int MinIndex(Iterator first, Iterator last) {
    assert(first != last);

    return static_cast<int>(std::distance(first, 
                                          std::min_element(first, last)));
}

/**
 * Return the index of the maximal element in [first, last).
 */
template <typename Iterator>
inline int MaxIndex(Iterator first, Iterator last) {
    assert(first != last);

    return static_cast<int>(std::distance(first,
                                          std::max_element(first, last)));
}

/**
 * Get the sorted indices of [first, last) according to 'compare'.
 */
template <typename Iterator, typename Compare>
void IndexSort(Iterator first, Iterator last, Compare compare, 
               Array<int>* indices) {
    assert(indices);

    int n = CountElements(first, last);
    indices->resize(n);
    std::iota(indices->begin(), indices->end(), 0);

    std::sort(indices->begin(), indices->end(), [&](int a, int b) {
        return compare(first[a], first[b]);
    });
}

/**
 * Get the sorted indices of [first, last).
 */
template <typename Iterator>
void IndexSort(Iterator first, Iterator last, Array<int>* indices) {
    using T = typename std::iterator_traits<Iterator>::value_type;
    IndexSort(first, last, std::less<T>(), indices);
}

} // namespace cl

#endif // BASE_ALGORITHM_H_
