//
// Copyright 2014 Yangbin Lin. All Rights Reserved.
//
// Author: yblin@jmu.edu.cn (Yangbin Lin)
//
// This file is part of the Code Library.
//

#ifndef BASE_EQUAL_H_
#define BASE_EQUAL_H_

#include <cassert>
#include <cmath>
#include <limits>
#include <type_traits>

namespace cl {

/**
 * Return true if 'a' is at most 'ulps' away from 'b'.
 *
 * In particular, this function:
 *  - returns false if either number is (or both are) NAN;
 *  - returns false if either number is (or both are) infinity;
 *  - thinks + 0.0 and -0.0 are 0 DLP's apart.
 */
template <typename Float>
bool ULPEqual(const Float& a, const Float& b, int ulps = 4) {
    static_assert(std::is_floating_point<Float>::value,
        "template argument is not a numeric type");
    assert(ulps >= 0);

    if (std::isnan(a) || std::isnan(b)) return false;

    Float prev_a = std::nextafter(a, std::numeric_limits<Float>::lowest());
    Float next_a = std::nextafter(a, std::numeric_limits<Float>::max());

    Float min_a = a - (a - prev_a) * ulps;
    Float max_a = a + (next_a - a) * ulps;

    return min_a <= b && b <= max_a;
}

/**
 * Check if two objects are equal.
 */
template <typename T1, typename T2>
inline bool Equal(const T1& lhs, const T2& rhs) {
    return lhs == rhs;
}
inline bool Equal(float a, float b) {
    return ULPEqual(a, b);
}
inline bool Equal(double a, double b) {
    return ULPEqual(a, b);
}

/**
 * Compare two numbers with tolerance.
 *
 * It adopts a combination of absolute and relative tolerances method.
 * This function will fail when a and b near to 0.
 *
 * Here, we assume absolute tolerance is equal to the relative tolerance.
 */
template <typename T>
bool Equal(const T& lhs, const T& rhs, const T& tolerance) {
    static_assert(std::is_arithmetic<T>::value,
        "template argument is not a numeric type");

    return std::abs(lhs - rhs) <= tolerance * std::max({ lhs, rhs, T(1) });
}

/**
 * Check if two ranges are equal.
 */
template <typename T1, typename T2>
bool Equal(T1 first1, T1 last1, T2 first2, T2 last2) {
    for (; first1 != last1 || first2 != last2; ++first1, ++first2) {
        if ((first1 == last1 || first2 == last2)) {
            return false;
        }

        if (!Equal(*first1, *first2)) return false;
    }

    return true;
}

} // namespace cl

#endif // BASE_EQUAL_H_
