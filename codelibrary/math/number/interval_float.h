//
// Copyright 2017 Yangbin Lin. All Rights Reserved.
//
// Author: yblin@jmu.edu.cn (Yangbin Lin)
//
// This file is part of the Code Library.
//

#ifndef MATH_NUMBER_INTERVAL_FLOAT_H_
#define MATH_NUMBER_INTERVAL_FLOAT_H_

#include <algorithm>
#include <cassert>
#include <cstdint>

namespace cl {

/**
 * A floating number type defined in interval form. It allows interval
 * arithmetics which is useful in the geometric predicates.
 *
 * It requires FP rounding set to FE_UPWARD.
 *
 * NOTE that, for VS complier, you need to set FP-mode to strict (/fp:strict).
 */
class IntervalFloat {
public:
    IntervalFloat() = default;
    explicit IntervalFloat(int32_t n) : lower_(-n), upper_(n) {}
    explicit IntervalFloat(float n)   : lower_(-n), upper_(n) {}
    explicit IntervalFloat(double n)  : lower_(-n), upper_(n) {}

    IntervalFloat(double lower, double upper)
        : lower_(-lower), upper_(upper) {
        assert(lower <= upper);
    }

    IntervalFloat(const IntervalFloat&) = default;

    IntervalFloat(IntervalFloat&&) = default;

    IntervalFloat& operator=(const IntervalFloat&) = default;

    IntervalFloat& operator=(IntervalFloat&&) = default;

    double lower() const {
        return -lower_;
    }

    double upper() const {
        return upper_;
    }

    /**
     * This = a + b.
     */
    IntervalFloat& Add(const IntervalFloat& a, const IntervalFloat& b) {
        lower_ = a.lower_ + b.lower_;
        upper_ = a.upper_ + b.upper_;

        return *this;
    }

    /**
     * This = a - b.
     */
    IntervalFloat& Subtract(const IntervalFloat& a, const IntervalFloat& b) {
        double l = a.lower_ + b.upper_; // Deal the case: this == &b.
        upper_ = a.upper_ + b.lower_;
        lower_ = l;

        return *this;
    }

    /**
     * This = a * b.
     */
    IntervalFloat& Multiply(const IntervalFloat& a, const IntervalFloat& b) {
        double l = 0.0, u = 0.0; // Deal the case: this == &a, or this == &b.
        if (a.lower_ <= 0.0) {
            l = ((b.lower_ >= 0.0) ?  a.upper_ : -a.lower_) * b.lower_;
            u = ((b.upper_ <= 0.0) ? -a.lower_ :  a.upper_) * b.upper_;
        } else if (a.upper_ <= 0.0) {
            l = ((b.upper_ <= 0.0) ? -a.upper_ :  a.lower_) * b.upper_;
            if (b.lower_ <= 0.0) {
                u = a.upper_ * -b.lower_;
            } else {
                u = a.lower_ * b.lower_;
            }
         } else {
            if (b.lower_ <= 0.0) {
                l = a.lower_ * b.upper_;
                u = a.upper_ * b.upper_;
            } else if (b.upper_ <= 0.0) {
                l = a.upper_ * b.lower_;
                u = a.lower_ * b.lower_;
            } else {
                l = std::max(a.lower_ * b.upper_, a.upper_ * b.lower_);
                u = std::max(a.upper_ * b.upper_, a.lower_ * b.lower_);
            }
        }

        lower_ = l;
        upper_ = u;

        return *this;
    }

    IntervalFloat& operator +=(const IntervalFloat& a) {
        return Add(*this, a);
    }

    IntervalFloat& operator -=(const IntervalFloat& a) {
        return Subtract(*this, a);
    }

    IntervalFloat& operator *=(const IntervalFloat& a) {
        return Multiply(*this, a);
    }

    friend IntervalFloat operator -(const IntervalFloat& a) {
        IntervalFloat r;
        r.lower_ = -a.upper_;
        r.upper_ = -a.lower_;
        return r;
    }

    friend IntervalFloat operator +(const IntervalFloat& a,
                                    const IntervalFloat& b) {
        IntervalFloat r;
        r.Add(a, b);
        return r;
    }

    friend IntervalFloat operator -(const IntervalFloat& a,
                                    const IntervalFloat& b) {
        IntervalFloat r;
        r.Subtract(a, b);
        return r;
    }

    friend IntervalFloat operator *(const IntervalFloat& a,
                                    const IntervalFloat& b) {
        IntervalFloat r;
        r.Multiply(a, b);
        return r;
    }

private:
    double lower_ = 0.0, upper_ = 0.0; // Two bounds of interval.
};

} // namespace cl

#endif // MATH_NUMBER_INTERVAL_FLOAT_H_
