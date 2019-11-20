//
// Copyright 2017 Yangbin Lin. All Rights Reserved.
//
// Author: yblin@jmu.edu.cn (Yangbin Lin)
//
// This file is part of the Code Library.
//

#ifndef MATH_NUMBER_EXACT_FLOAT_H_
#define MATH_NUMBER_EXACT_FLOAT_H_

#include <cassert>
#include <cmath>
#include <cstdint>
#include <ostream>

#include "codelibrary/math/number/bigint.h"

namespace cl {

/**
 * ExactFloat is a multiple-precision floating point type. It supports exact
 * arithmetic including addition, substraction, and multiplication. Note that,
 * it does not support division, since the quotient of two floating numbers may
 * not be represented exactly.
 *
 * ExactFloat is useful for geometric algorithms, especially for disambiguating
 * cases where ordinary double-precision arithmetic yields an uncertain result.
 */
class ExactFloat {
    // Significand precision of double is 53.
    static const int DOUBLE_MANTISSA_BITS = 53;

    // Significand precision of float is 53.
    static const int FLOAT_MANTISSA_BITS = 24;

public:
    ExactFloat() = default;

    explicit ExactFloat(double n) {
        assert(!std::isnan(n));

        if (n != 0.0) {
            // x = significand * 2 ^ exponent
            double x = std::frexp(n, &exponent_);
            int64_t m = static_cast<int64_t>(std::ldexp(x,
                                                        DOUBLE_MANTISSA_BITS));
            exponent_ -= DOUBLE_MANTISSA_BITS;
            significand_.Assign(m);
        }
    }

    explicit ExactFloat(float n) {
        assert(!std::isnan(n));

        if (n != 0.0f) {
            // x = significand * 2 ^ exponent
            float x = std::frexp(n, &exponent_);
            int32_t m = static_cast<int32_t>(std::ldexp(x,
                                                        FLOAT_MANTISSA_BITS));
            exponent_ -= FLOAT_MANTISSA_BITS;
            significand_.Assign(m);
        }
    }

    explicit ExactFloat(int32_t n)  { significand_.Assign(n); }
    explicit ExactFloat(uint32_t n) { significand_.Assign(n); }
    explicit ExactFloat(int64_t n)  { significand_.Assign(n); }
    explicit ExactFloat(uint64_t n) { significand_.Assign(n); }

    /**
     * Return the sign of ExactFloat: -1 for negative, 0 for zero, and 1 for
     * positive.
     */
    int sign() const {
        return significand_.sign();
    }

    int exponent() const {
        return exponent_;
    }

    const BigInt& significand() const {
        return significand_;
    }

    /**
     * This = a + b.
     */
    ExactFloat& Add(const ExactFloat& a, const ExactFloat& b) {
        if (a.exponent_ > b.exponent_) {
            BigInt r;
            r.ShiftLeft(a.significand_, a.exponent_ - b.exponent_);
            significand_.Add(r, b.significand_);
            exponent_ = b.exponent_;
        } else {
            BigInt r;
            r.ShiftLeft(b.significand_, b.exponent_ - a.exponent_);
            significand_.Add(r, a.significand_);
            exponent_ = a.exponent_;
        }

        if (significand_.IsZero()) exponent_ = 0;

        return *this;
    }

    /**
     * This = a - b.
     */
    ExactFloat& Subtract(const ExactFloat& a, const ExactFloat& b) {
        if (a.exponent_ > b.exponent_) {
            BigInt r;
            r.ShiftLeft(a.significand_, a.exponent_ - b.exponent_);
            significand_.Subtract(r, b.significand_);
            exponent_ = b.exponent_;
        } else {
            BigInt r;
            r.ShiftLeft(b.significand_, b.exponent_ - a.exponent_);
            significand_.Subtract(a.significand_, r);
            exponent_ = a.exponent_;
        }

        if (significand_.IsZero()) exponent_ = 0;

        return *this;
    }

    /**
     * This = a * b.
     */
    ExactFloat& Multiply(const ExactFloat& a, const ExactFloat& b) {
        if (a.exponent_ < 0 && b.exponent_ < 0) {
            assert(a.exponent_ + b.exponent_ < a.exponent_ &&
                   "Exponenet overflow.");
        } else if (a.exponent_ > 0 && b.exponent_ > 0) {
            assert(a.exponent_ + b.exponent_ > a.exponent_ &&
                   "Exponenet overflow.");
        }

        significand_.Multiply(a.significand_, b.significand_);
        exponent_ = significand_.IsZero() ? 0 : a.exponent_ + b.exponent_;

        return *this;
    }

    ExactFloat& operator +=(const ExactFloat& rhs) {
        return Add(*this, rhs);
    }

    ExactFloat& operator -=(const ExactFloat& rhs) {
        return Subtract(*this, rhs);
    }

    ExactFloat& operator *=(const ExactFloat& rhs) {
        return Multiply(*this, rhs);
    }

    friend ExactFloat operator -(const ExactFloat& a) {
        ExactFloat r;
        r.significand_ = -a.significand_;
        r.exponent_ = a.exponent_;
        return r;
    }

    friend ExactFloat operator +(const ExactFloat& a, const ExactFloat& b) {
        ExactFloat r;
        r.Add(a, b);
        return r;
    }

    friend ExactFloat operator -(const ExactFloat& a, const ExactFloat& b) {
        ExactFloat r;
        r.Subtract(a, b);
        return r;
    }

    friend ExactFloat operator *(const ExactFloat& a, const ExactFloat& b) {
        ExactFloat r;
        r.Multiply(a, b);
        return r;
    }

private:
    // Base 2 exponent of ExactFloat.
    int exponent_ = 0;

    // Significant data.
    // ExactFloat = significand * 2^exponent
    BigInt significand_;
};

} // namespace cl

#endif // MATH_NUMBER_EXACT_FLOAT_H_
