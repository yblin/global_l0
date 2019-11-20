//
// Copyright 2015 Yangbin Lin. All Rights Reserved.
//
// Author: yblin@jmu.edu.cn (Yangbin Lin)
//
// This file is part of the Code Library.
//
// This file defines some common bit utilities.
//

#ifndef BASE_BITS_H_
#define BASE_BITS_H_

#include <cassert>
#include <cstdint>
#include <limits>
#include <type_traits>

namespace cl {
namespace bits {

/**
 * Return the integer i such as 2^i <= n < 2^(i+1).
 *
 * Equal to n_bits - __builtin_clz(n)
 */
template <typename IntType>
int Log2Floor(IntType n) {
    static_assert(std::is_integral<IntType>::value, "n should be an integer.");
    assert(n > 0);

    int start = 0;
    switch (std::numeric_limits<IntType>::digits) {
    case 7:
    case 8:
        start = 2;
        break;
    case 15:
    case 16:
        start = 3;
        break;
    case 31:
    case 32:
        start = 4;
        break;
    case 63:
    case 64:
        start = 5;
        break;
    default:
        assert(false && "Over-sized type.");
    }

    int log = 0;
    for (int i = start; i >= 0; --i) {
        int shift = (1 << i);
        IntType x = n >> shift;
        if (x != 0) {
            n = x;
            log += shift;
        }
    }

    return log;
}

/**
 * Return the integer i such as 2^(i-1) < n <= 2^i.
 */
template <typename IntType>
inline IntType Log2Ceiling(IntType n) {
    assert(n >= 0);
    if (n == 0 || n == 1) return 0;

    return IntType(1) + Log2Floor(n - 1);
}

/**
 * Return the integer i such as 2^i <= n < 2^(i+1).
 */
template <typename IntType>
inline IntType Power2Floor(const IntType& n) {
     return IntType(1) << Log2Floor(n);
}

/**
 * Return the integer 2^i such as 2^(i-1) < n <= 2^i.
 */
template <typename IntType>
inline IntType Power2Ceiling(const IntType& n) {
    int n_bits = Log2Ceiling(n);
    assert(n_bits < std::numeric_limits<IntType>::digits);

    return IntType(1) << n_bits;
}

} // namespace bits
} // namespace cl

#endif // BASE_BITS_H_
