//
// Copyright 2011 Yangbin Lin. All Rights Reserved.
//
// Author: yblin@jmu.edu.cn (Yangbin Lin)
//
// This file is part of the Code Library.
//

#ifndef MATH_NUMBER_BIGINT_H_
#define MATH_NUMBER_BIGINT_H_

#include <algorithm>
#include <cassert>
#include <climits>
#include <cstdint>
#include <istream>
#include <limits>
#include <ostream>
#include <random>
#include <string>
#include <vector>

namespace cl {

/**
 * Multiple precision integer.
 *
 * The radix of BigInt is 2^32 and its size can be up to INT_MAX / 32.
 * It is to ensure that the length of binary string can be stored in int.
 *
 * BigInt can be construct by int or string with given base.
 *
 * For example, the decimal -10 can be describe below:
 *   - BigInt a(-10);
 *   - BigInt a("-1010", 2);
 *   - BigInt a("-A", 16);
 */
class BigInt {
    // Threshold for Karatsuba multiplication.
    static const int KARATSUBA_MULTIPLY_THRESHOLD = 64;

    static const uint32_t RADIX_MASK = 0xffffffffU;

public:
    // Bits of radix.
    static const int RADIX_BITS = 32;

    // The maximum number of bytes of BigInt.
    // So that the maximum bits of BigInt is smaller than INT_MAX.
    static const int MAX_SIZE = INT_MAX / RADIX_BITS;

    BigInt()           : data_(2) {}
    BigInt(int32_t n)  : data_(2) { Assign(n); }
    BigInt(uint32_t n) : data_(2) { Assign(n); }
    BigInt(int64_t n)  : data_(2) { Assign(n); }
    BigInt(uint64_t n) : data_(2) { Assign(n); }

    explicit BigInt(const std::string& str, int radix = 10) {
        Assign(str, radix);
    }

    explicit BigInt(const char* str, int radix = 10)   {
        Assign(std::string(str), radix);
    }

    BigInt(const BigInt& n) = default;

    BigInt(BigInt&& n) noexcept = default;

    // Compare operators.
    bool operator ==(const BigInt& rhs) const { return Compare(rhs) == 0; }
    bool operator !=(const BigInt& rhs) const { return Compare(rhs) != 0; }
    bool operator < (const BigInt& rhs) const { return Compare(rhs) < 0;  }
    bool operator > (const BigInt& rhs) const { return Compare(rhs) > 0;  }
    bool operator <=(const BigInt& rhs) const { return Compare(rhs) <= 0; }
    bool operator >=(const BigInt& rhs) const { return Compare(rhs) >= 0; }
    bool operator ==(int32_t rhs)       const { return Compare(rhs) == 0; }
    bool operator !=(int32_t rhs)       const { return Compare(rhs) != 0; }
    bool operator < (int32_t rhs)       const { return Compare(rhs) < 0;  }
    bool operator > (int32_t rhs)       const { return Compare(rhs) > 0;  }
    bool operator <=(int32_t rhs)       const { return Compare(rhs) <= 0; }
    bool operator >=(int32_t rhs)       const { return Compare(rhs) >= 0; }
    bool operator ==(uint32_t rhs)      const { return Compare(rhs) == 0; }
    bool operator !=(uint32_t rhs)      const { return Compare(rhs) != 0; }
    bool operator < (uint32_t rhs)      const { return Compare(rhs) < 0;  }
    bool operator > (uint32_t rhs)      const { return Compare(rhs) > 0;  }
    bool operator <=(uint32_t rhs)      const { return Compare(rhs) <= 0; }
    bool operator >=(uint32_t rhs)      const { return Compare(rhs) >= 0; }

    // Assign operators.
    BigInt& operator =(int32_t rhs)            { return Assign(rhs);     }
    BigInt& operator =(uint32_t rhs)           { return Assign(rhs);     }
    BigInt& operator =(int64_t rhs)            { return Assign(rhs);     }
    BigInt& operator =(uint64_t rhs)           { return Assign(rhs);     }
    BigInt& operator =(const std::string& rhs) { return Assign(rhs, 10); }

    BigInt& operator =(const BigInt& n)     = default;
    BigInt& operator =(BigInt&& n) noexcept = default;

    // Arithmetic operators.
    BigInt& operator +=(const BigInt& rhs) { return Add(*this, rhs);        }
    BigInt& operator -=(const BigInt& rhs) { return Subtract(*this, rhs);   }
    BigInt& operator *=(const BigInt& rhs) { return Multiply(*this, rhs);   }
    BigInt& operator *=(int32_t rhs)       { return Multiply(*this, rhs);   }
    BigInt& operator *=(uint32_t rhs)      { return Multiply(*this, rhs);   }
    BigInt& operator /=(const BigInt& rhs) { return Divide(*this, rhs);     }
    BigInt& operator /=(int32_t rhs)       { return Divide(*this, rhs);     }
    BigInt& operator /=(uint32_t rhs)      { return Divide(*this, rhs);     }
    BigInt& operator %=(const BigInt& rhs) { return Mod(*this, rhs);        }
    BigInt& operator %=(int32_t rhs)       { return Mod(*this, rhs);        }
    BigInt& operator %=(uint32_t rhs)      { return Mod(*this, rhs);        }
    BigInt& operator <<=(int rhs)          { return ShiftLeft(*this, rhs);  }
    BigInt& operator >>=(int rhs)          { return ShiftRight(*this, rhs); }

    /**
     * Convert the BigInt to string.
     *
     * The radix must between [2, 36], the default is 10.
     * Return string based on radix, where 'A' to 'Z' means number 10 to 35.
     */
    std::string ToString(int radix = 10) const {
        assert(radix >= 2 && radix <= 36);

        if (IsZero()) return "0";

        std::string str;

        // If radix is 16 or 2, special treat.
        if (radix == 16) {
            for (int i = size_ - 1; i >= 0; --i) {
                str += UIntToHexString(data_[i]);
            }
            DeleteLeadingZero(&str);
        } else if (radix == 2) {
            for (int i = size_ - 1; i >= 0; --i) {
                str += UIntToBinaryString(data_[i]);
            }
            DeleteLeadingZero(&str);
        } else if (radix == 10) {
            // If radix is 10, expand the radix to 1000000000 for speed.
            const int decimal_base  = 1000000000;
            const int decimal_digit = 9;
            int length = 1 + size_ * RADIX_BITS / (3 * decimal_digit);
            std::vector<uint32_t> decimal(length);

            int size = 0;
            for (int i = size_ - 1; i >= 0; --i) {
                uint32_t hi = data_[i];
                uint64_t z = 0;

                for (int j = 0; j < size; ++j) {
                    z = (static_cast<uint64_t>(decimal[j]) << RADIX_BITS) | hi;
                    hi = static_cast<uint32_t>(z / decimal_base);
                    decimal[j] = static_cast<uint32_t>(z - hi * decimal_base);
                }
                while (hi != 0) {
                    decimal[size++] = hi % decimal_base;
                    hi /= decimal_base;
                }
            }

            for (int i = size - 1; i >= 0; --i) {
                std::string t(decimal_digit, '0');
                int j = 0;
                while (decimal[i] != 0) {
                    t[j++] = decimal[i] % 10 + '0';
                    decimal[i] /= 10;
                }
                std::reverse(t.begin(), t.end());
                str += t;
            }

            DeleteLeadingZero(&str);
        } else {
            BigInt a = *this;
            uint32_t remainder;
            while (a.size_ != 0) {
                // a = a / radix, remainder = a % radix.
                a.DivideData(a, static_cast<uint32_t>(radix), &remainder);
                str += UIntToChar(remainder);
            }
            std::reverse(str.begin(), str.end());
        }

        return sign_ > 0 ? str : std::string("-") + str;
    }

    /**
     * This = a << number.
     *
     * The number can not be negative.
     */
    BigInt& ShiftLeft(const BigInt& a, int number) {
        ShiftLeftData(a, number);
        sign_ = a.sign_;
        return *this;
    }

    /**
     * This = a >> number.
     *
     * The number can not be negative.
     */
    BigInt& ShiftRight(const BigInt& a, int number) {
        ShiftRightData(a, number);
        sign_ = size_ == 0 ? 0 : a.sign_;
        return *this;
    }

    /**
     * This = -BigInt.
     */
    BigInt& Negate() {
        sign_ = -sign_;
        return *this;
    }

    /**
     * This = -a.
     */
    BigInt& Negate(const BigInt& a) {
        size_ = a.size_;
        data_ = a.data_;
        sign_ = -a.sign_;
        return *this;
    }

    /**
     * This = a + b.
     */
    BigInt& Add(const BigInt& a, const BigInt& b) {
        if (a.IsZero()) return *this = b;
        if (b.IsZero()) return *this = a;

        if (a.sign_ == b.sign_) {
            sign_ = a.sign_;
            AddData(a, b);
        } else {
            int t = CompareData(a, b);
            if (t == 0) {
                this->SetZero();
            } else if (t > 0) {
                sign_ = a.sign_;
                SubtractData(a, b);
            } else {
                sign_ = -a.sign_;
                SubtractData(b, a);
            }
        }

        return *this;
    }

    /**
     * This = a - b.
     */
    BigInt& Subtract(const BigInt& a, const BigInt& b) {
        if (a.IsZero()) return this->Negate(b);
        if (b.IsZero()) return *this = a;

        if (a.sign_ != b.sign_) {
            sign_ = a.sign_;
            AddData(a, b);
        } else {
            int t = CompareData(a, b);
            if (t == 0) {
                SetZero();
            } else if (t > 0) {
                sign_ = a.sign_;
                SubtractData(a, b);
            } else {
                sign_ = -a.sign_;
                SubtractData(b, a);
            }
        }

        return *this;
    }

    /**
     * This = a * b.
     *
     * Auto choose multiplication algorithm by two BigInt's sizes.
     */
    BigInt& Multiply(const BigInt& a, const BigInt& b) {
        MultiplyData(a, b);
        sign_ = a.sign_ * b.sign_;
        return *this;
    }

    /**
     * This = a * b.
     */
    BigInt& Multiply(const BigInt &a, int32_t b) {
        if (b < 0) {
            MultiplyData(a, uint32_t(-b));
        } else {
            MultiplyData(a, uint32_t(b));
        }

        sign_ = a.sign_ * GetSign(b);
        return *this;
    }

    /**
     * This = a * b.
     */
    BigInt& Multiply(const BigInt &a, uint32_t b) {
        MultiplyData(a, b);
        sign_ = a.sign_;
        return *this;
    }

    /**
     * This = a / b, remainder = a % b.
     */
    BigInt& Divide(const BigInt& a, uint32_t b, uint32_t* remainder = nullptr) {
        assert(b != 0);

        uint32_t remainder_t;
        DivideData(a, b, &remainder_t);
        sign_ = a.sign_;
        if (remainder != nullptr) {
            *remainder = remainder_t;
        }
        return *this;
    }

    /**
     * This = a / b, remainder = a % b.
     */
    BigInt& Divide(const BigInt& a, int32_t b, int32_t* remainder = nullptr) {
        assert(b != 0);

        uint32_t t_b = (b < 0) ? uint32_t(-b) : uint32_t(b);
        uint32_t r;
        this->DivideData(a, t_b, &r);
        if (remainder != nullptr) {
            *remainder = static_cast<int32_t>(r) * a.sign_;
        }

        sign_ = a.sign_ * GetSign(b);
        if (size_ == 0) sign_ = 0;

        return *this;
    }

    /**
     * This = a / b, remainder = a % b.
     */
    BigInt& Divide(const BigInt& a, const BigInt& b,
                   BigInt* remainder = nullptr) {
        assert(!b.IsZero());

        BigInt remainder_t;
        DivideData(a, b, &remainder_t);

        if (remainder != nullptr) {
            remainder->size_ = remainder_t.size_;
            remainder->data_ = remainder_t.data_;
            remainder->sign_ = a.sign_ * (remainder_t.size_ != 0 ? 1 : 0);
        }
        sign_ = a.sign_ * b.sign_;
        if (size_ == 0) sign_ = 0;

        return *this;
    }

    /**
     * This = a % b.
     */
    BigInt& Mod(const BigInt& a, int32_t b) {
        assert(b != 0);

        int32_t remainder;
        this->Divide(a, b, &remainder);
        return this->Assign(remainder);
    }

    /**
     * This = a % b.
     */
    BigInt& Mod(const BigInt& a, uint32_t b) {
        assert(b != 0);

        uint32_t remainder;
        this->Divide(a, b, &remainder);
        return this->Assign(remainder);
    }

    /**
     * This = a % b.
     */
    BigInt& Mod(const BigInt& a, const BigInt& b) {
        assert(!b.IsZero());

        BigInt remainder;
        this->Divide(a, b, &remainder);
        return *this = remainder;
    }

    /**
     * Return the integer part of square root.
     */
    BigInt Sqrt() const {
        assert(sign_ >= 0);

        if (sign_ == 0) return BigInt(0);

        int bits = Bits();
        BigInt x(1);
        x.ShiftLeft(x, bits / 2 + bits % 2);

        BigInt y;
        while (true) {
            y = (x + *this / x) >> 1;
            if (y >= x) return x;
            x = y;
        }

        return BigInt(-1);
    }

    /**
     * Return the number of bits of the BigInt.
     */
    int Bits() const {
        if (this->IsZero()) return 0;
        return (size_ - 1) * RADIX_BITS + Log2(data_[size_ - 1]) + 1;
    }

    /**
     * Return true if BigInt is zero.
     */
    bool IsZero() const {
        return sign_ == 0;
    }

    /**
     * Set BigInt to zero.
     */
    BigInt& SetZero() {
        size_ = 0;
        sign_ = 0;
        return *this;
    }

    /**
     * Return the sign of BigInt.
     */
    int sign() const {
        return sign_;
    }

    /**
     * Return the current size.
     *
     * Note that, if the BigInt is zero, the size will be zero too.
     */
    int size() const {
        return size_;
    }

    /**
     * The internel data vector.
     */
    const std::vector<uint32_t>& data() const {
        return data_;
    }


    /**
     * Assign BigInt from string.
     *
     *  str   - the string that compose of letter or digit. Where, 'A' to 'Z' or
     *          ('a' to 'z' means number 10 to 35.
     *  radix - between [2, 36], the default is 10.
     */
    BigInt& Assign(const std::string& str, int radix = 10) {
        assert(radix >= 2 && radix <= 36);
        assert(!str.empty() && "The input string is empty.");

        this->SetZero();

        // Compute sign;
        int sign = 0;
        size_t begin = 0;
        if (str[0] == '+') {
            sign = 1;
            begin = 1;
            assert(str.size() > 1 && "The input string is invalid.");
        } else if (str[0] == '-') {
            sign = -1;
            begin = 1;
            assert(str.size() > 1 && "The input string is invalid.");
        } else {
            sign = 1;
        }

        // If radix is 16 or 2, special treat.
        if (radix == 2) {
            assert(str.size() / RADIX_BITS + 1 <= INT_MAX);
            data_.resize(str.size() / RADIX_BITS + 1);
            size_t i;
            for (i = str.size(); ; i -= RADIX_BITS) {
                if (i < RADIX_BITS) break;
                data_[size_++] = StringToUInt(str.substr(i - RADIX_BITS,
                                                         RADIX_BITS), 2);
            }
            if (i != begin) {
                data_[size_++] = StringToUInt(str.substr(begin, i), 2);
            }
            this->DeleteLeadingZero();
        } else if (radix == 16) {
            size_t offset = RADIX_BITS / 4;
            data_.resize(str.size() / offset + 1);

            size_t i;
            for (i = str.size(); ; i -= offset) {
                if (i < offset) break;
                data_[size_++] = StringToUInt(str.substr(i - offset, offset),
                                              16);
            }
            if (i != begin) {
                data_[size_++] = StringToUInt(str.substr(begin, i), 16);
            }
            this->DeleteLeadingZero();
        } else {
            size_t n = 0;
            uint32_t pow[32];
            pow[n] = 1U;
            while (pow[n] <= UINT32_MAX / radix) {
                pow[n + 1] = pow[n] * radix;
                ++n;
            }

            uint32_t number;
            for (size_t i = begin; i < str.size(); i += n) {
                number = StringToUInt(str.substr(i, n), radix);
                size_t index = std::min(n, str.size() - i);
                MultiplyData(*this, pow[index]);
                AddData(*this, BigInt(number));
            }
        }

        assert(size_ <= MAX_SIZE);

        sign_ = (size_ == 0) ? 0 : sign;

        return *this;
    }

    /**
     * Assign BigInt from int.
     */
    BigInt& Assign(int32_t n) {
        if (n == 0) return this->SetZero();

        size_ = 1;
        data_[0] = static_cast<uint32_t>(n < 0 ? -n : n);
        sign_ = GetSign(n);

        return *this;
    }
    BigInt& Assign(uint32_t n) {
        if (n == 0) return this->SetZero();

        size_ = 1;
        data_[0] = n;
        sign_ = 1;

        return *this;
    }
    BigInt& Assign(int64_t n) {
        data_[0] = static_cast<uint32_t>(n < 0 ? ((-n) & RADIX_MASK)
                                               : (n & RADIX_MASK));
        data_[1] = static_cast<uint32_t>(n < 0 ? ((-n) >> RADIX_BITS)
                                               : (n >> RADIX_BITS));

        size_ = data_[1] == 0 ? (data_[0] == 0 ? 0 : 1) : 2;
        sign_ = GetSign(n);

        return *this;
    }
    BigInt& Assign(uint64_t n) {
        data_[0] = n & RADIX_MASK;
        data_[1] = n >> RADIX_BITS;

        size_ = data_[1] == 0 ? (data_[0] == 0 ? 0 : 1) : 2;
        sign_ = GetSign(n);

        return *this;
    }

    /**
     * Convert this bigint to uint32.
     */
    uint32_t ToUInt32() const {
        assert(size_ <= 1);

        if (size_ == 0) return 0;
        return data_[0];
    }

    /**
     * Convert this bigint to uint64.
     */
    uint64_t ToUInt64() const {
        assert(size_ <= 2);

        if (size_ == 0) return 0;
        if (size_ == 1) return data_[0];
        return (static_cast<uint64_t>(data_[1]) << 32) + data_[0];
    }

    friend BigInt operator -(const BigInt& rhs) {
        BigInt t;
        t.Negate(rhs);
        return t;
    }

    friend BigInt operator +(const BigInt& lhs, const BigInt& rhs) {
        BigInt t;
        t.Add(lhs, rhs);
        return t;
    }

    friend BigInt operator -(const BigInt& lhs, const BigInt& rhs) {
        BigInt t;
        t.Subtract(lhs, rhs);
        return t;
    }

    friend BigInt operator *(const BigInt& lhs, int32_t rhs) {
        BigInt t;
        t.Multiply(lhs, rhs);
        return t;
    }

    friend BigInt operator *(const BigInt& lhs, uint32_t rhs) {
        BigInt t;
        t.Multiply(lhs, rhs);
        return t;
    }

    friend BigInt operator *(int32_t lhs, const BigInt& rhs) {
        BigInt t;
        t.Multiply(rhs, lhs);
        return t;
    }

    friend BigInt operator *(uint32_t lhs, const BigInt& rhs) {
        BigInt t;
        t.Multiply(rhs, lhs);
        return t;
    }

    friend BigInt operator *(const BigInt& lhs, const BigInt& rhs) {
        BigInt t;
        t.Multiply(lhs, rhs);
        return t;
    }

    friend BigInt operator /(const BigInt& lhs, int32_t rhs) {
        BigInt t;
        t.Divide(lhs, rhs);
        return t;
    }

    friend BigInt operator /(const BigInt& lhs, uint32_t rhs) {
        BigInt t;
        t.Divide(lhs, rhs);
        return t;
    }

    friend BigInt operator /(const BigInt& lhs, const BigInt& rhs) {
        BigInt t;
        t.Divide(lhs, rhs);
        return t;
    }

    friend int32_t operator %(const BigInt& lhs, int32_t rhs) {
        int32_t remainder;
        BigInt t;
        t.Divide(lhs, rhs, &remainder);
        return remainder;
    }

    friend uint32_t operator %(const BigInt& lhs, uint32_t rhs) {
        uint32_t remainder;
        BigInt t;
        t.Divide(lhs, rhs, &remainder);
        return remainder;
    }

    friend BigInt operator %(const BigInt& lhs, const BigInt& rhs) {
        BigInt t;
        t.Mod(lhs, rhs);
        return t;
    }

    friend BigInt operator <<(const BigInt& lhs, int rhs) {
        BigInt t;
        t.ShiftLeft(lhs, rhs);
        return t;
    }

    friend BigInt operator >>(const BigInt& lhs, int rhs) {
        BigInt t;
        t.ShiftRight(lhs, rhs);
        return t;
    }

    /**
     * Overloaded input stream, the radix is 10.
     */
    friend std::istream& operator >>(std::istream& is, BigInt& rhs) {
        std::string string;
        if (!(is >> string)) return is;
        rhs.Assign(string, 10);
        return is;
    }

    /**
     * Overloaded input stream, the default radix is 10.
     */
    friend std::ostream& operator <<(std::ostream& os, const BigInt& rhs) {
        if ((os.flags() & std::ostream::hex) != 0) {
            os << rhs.ToString(16);
        } else {
            os << rhs.ToString(10);
        }
        return os;
    }

private:
    /**
     * Get the sign of number.
     */
    template <typename T>
    static int GetSign(const T& number) {
        if (number == 0) return 0;
        return (number > 0) ? 1 : -1;
    }

    /**
     * Compare this BigInt to BigInt a.
     * Return  +1, if this > a;
     *         -1, if this < a;
     *          0, otherwise.
     */
    int Compare(const BigInt& a) const {
        if (sign_ == a.sign_) {
            return sign_ * CompareData(*this, a);
        }

        return sign_ == 0 ? a.sign_ * -1 : sign_;
    }

    /**
     * Compare this BigInt to uint32.
     */
    int Compare(uint32_t a) const {
        if (sign_ == 0) {
            return a > 0 ? -1 : 0;
        }

        if (sign_ == 1) {
            if (size_ > 1) return 1;
            return data_[0] == a ? 0 : (data_[0] > a ? 1 : -1);
        }
        return -1;
    }

    /**
     * Compare this BigInt to int32.
     */
    int Compare(int32_t a) const {
        if (sign_ == 0) {
            return a > 0 ? -1 : (a == 0 ? 0 : 1);
        }

        if (sign_ == 1) {
            if (size_ > 1 || a <= 0) return 1;
            return data_[0] > static_cast<uint32_t>(a) ?  1 :
                   data_[0] < static_cast<uint32_t>(a) ? -1 : 0;
        }

        if (size_ > 1 || a >= 0) return -1;
        if (data_[0] == static_cast<uint32_t>(std::abs(a))) return 0;
        return data_[0] < static_cast<uint32_t>(std::abs(a)) ? 1 : -1;
    }

    /**
     * Delete the unnecessary leading zero of the ouput string.
     */
    static void DeleteLeadingZero(std::string* str) {
        // Delete leading zero.
        size_t begin = 0;
        for (begin = 0; begin < str->size(); ++begin) {
            if ((*str)[begin] != '0') break;
        }
        *str = str->substr(begin);
    }

    /**
     * Convert a character to int.
     *
     * The character only allows letter or digit.
     */
    static uint32_t CharToUInt(char character, int radix) {
        assert((isdigit(character) || isalpha(character)) &&
               "Unknown input character.");

        int res = 0;

        if (isdigit(character) != 0) {
            res = character - '0';
        } else if (islower(character) != 0) {
            res = character - 'a' + 10;
        } else {
            res = character - 'A' + 10;
        }

        assert(res < radix && "Invalid input string.");

        return static_cast<uint32_t>(res);
    }

    /**
     * Convert int to char.
     *
     * The number must between [0, 36).
     */
    static char UIntToChar(uint32_t number) {
        assert(number < 36);

        return number >= 10 ? char(number + 'A' - 10) : char(number + '0');
    }

    /**
     * Return hexadecimal string of int, the string may contain leading zero.
     */
    static std::string UIntToHexString(uint32_t number) {
        std::string hex_string(RADIX_BITS / 4, '0');
        int index = RADIX_BITS / 4;
        while (number != 0) {
            hex_string[--index] = UIntToChar(number & 0x0000000F);
            number >>= 4;
        }
        return hex_string;
    }

    /**
     * Return binary string of int, the string may contain leading zero.
     */
    static std::string UIntToBinaryString(uint32_t number) {
        std::string binary_string(RADIX_BITS, '0');
        int index = RADIX_BITS;
        while (number != 0) {
            binary_string[--index] = (number & 0x00000001) + '0';
            number >>= 1;
        }
        return binary_string;
    }

    /**
     * Convert the string to unsigned int number.
     *
     * The radix must between [2, 36].
     */
    static uint32_t StringToUInt(const std::string& str, int radix) {
        uint32_t number = 0;

        for (char i : str) {
            uint32_t tmp = number * uint32_t(radix) + CharToUInt(i, radix);
            assert(tmp >= number);
            number = tmp;
        }

        return number;
    }

    /**
     * Find the log base 2 of an 32-bit integer in O(lg(N)) operations.
     *
     * n is the 32-bit unsigned integer to find the log2 of.
     */
    static int Log2(uint32_t n) {
        assert(n > 0);

        static const int MULTIPLY_DE_BRUIJN_BIT_POSITION[32] = {
            0, 9, 1, 10, 13, 21, 2, 29, 11, 14, 16, 18, 22, 25, 3, 30,
            8, 12, 20, 28, 15, 17, 24, 7, 19, 27, 23, 6, 26, 5, 4, 31
        };

        n |= n >> 1; // first round down to one less than a power of 2.
        n |= n >> 2;
        n |= n >> 4;
        n |= n >> 8;
        n |= n >> 16;

        return MULTIPLY_DE_BRUIJN_BIT_POSITION[(n * 0x07C4ACDDU) >> 27];
    }

    /**
     * Decreases size to eliminate any leading zero blocks.
     */
    BigInt& DeleteLeadingZero() {
        while (size_ > 0 && data_[size_ - 1] == 0) {
            --size_;
        }
        if (size_ == 0) sign_ = 0;
        return *this;
    }

    /**
     * Compare the data between two BigInts, ignore the sign.
     *
     * Return +1, if a.data > b.data;
     *        -1, if a.data < b.data;
     *         0, otherwise.
     */
    int CompareData(const BigInt& a, const BigInt& b) const {
        if (a.size_ != b.size_) {
            return a.size_ < b.size_ ? -1 : 1;
        }
        for (int i = a.size_ - 1; i >= 0; --i) {
            if (a.data_[i] == b.data_[i]) continue;
            return (a.data_[i] < b.data_[i]) ? -1 : 1;
        }
        return 0;
    }

    /**
     * this->data = a.data + b.data.
     */
    void AddData(const BigInt& a, const BigInt& b) {
        // a1 points to the longer input, b1 points to the shorter.
        const std::vector<uint32_t>*a1 = &a.data_, *b1 = &b.data_;
        int a_size = a.size_;
        int b_size = b.size_;
        if (a_size < b_size) {
            std::swap(a1, b1);
            std::swap(a_size, b_size);
        }

        data_.resize(a_size + 1);

        bool carry_in = false, carry_out;
        uint32_t tmp;
        int i;

        for (i = 0; i < b_size; ++i) {
            tmp = (*a1)[i] + (*b1)[i];
            if (!carry_in) {
                carry_in = tmp < (*a1)[i];
                data_[i] = tmp;
            } else {
                carry_out = tmp < (*a1)[i] || tmp == UINT32_MAX;
                data_[i] = tmp + static_cast<uint32_t>(carry_in);
                carry_in = carry_out;
            }
        }

        for (i = b_size; i < a_size; ++i) {
            tmp = (*a1)[i] + static_cast<uint32_t>(carry_in);
            carry_in = (tmp < (*a1)[i]);
            data_[i] = tmp;
        }

        if (carry_in) {
            data_[i] = 1;
            size_ = i + 1;
        } else {
            size_ = i;
        }

        assert(size_ <= MAX_SIZE && "Overflow detected.");
    }

    /**
     * this->data = a.data - b.data.
     *
     * Note that, a MUST be greater than b.
     */
    void SubtractData(const BigInt& a, const BigInt& b) {
        int a_size = a.size_;
        int b_size = b.size_;
        data_.resize(a_size);

        bool carry_in = false, carry_out;
        int i;
        uint32_t tmp;

        for (i = 0; i < b_size; ++i) {
            if (a.data_[i] == b.data_[i]) {
                data_[i] = 0 - static_cast<uint32_t>(carry_in);
            } else {
                tmp = a.data_[i] - b.data_[i];
                carry_out = tmp > a.data_[i];
                data_[i] = tmp - static_cast<uint32_t>(carry_in);
                carry_in = carry_out;
            }
        }

        for (i = b_size; i < a_size; ++i) {
            if (carry_in) {
                carry_in = (a.data_[i] == 0);
                data_[i] = a.data_[i] - 1;
            } else {
                data_[i] = a.data_[i];
            }
        }

        size_ = a_size;
        DeleteLeadingZero();
    }

    /**
     * this.data = a.data << number.
     *
     * The number can not be negative.
     */
    void ShiftLeftData(const BigInt& a, int number) {
        assert(number >= 0);

        if (a.size_ == 0) {
            size_ = 0;
            return;
        }
        if (number == 0) {
            size_ = a.size_;
            data_ = a.data_;
            return;
        }

        int t1 = number / RADIX_BITS;
        int t2 = number % RADIX_BITS;

        assert(MAX_SIZE - a.size_ >= t1);

        data_.resize(a.size_ + t1 + 1);

        std::copy_backward(a.data_.begin(), a.data_.begin() + a.size_,
                           data_.begin() + a.size_ + t1);
        std::fill_n(data_.begin(), t1, 0);

        size_ = a.size_ + t1 + 1;
        data_[size_ - 1] = 0;
        if (t2 != 0) {
            for (int i = size_ - 1; i > 0; --i) {
                data_[i] = (data_[i] << t2) |
                           (data_[i - 1] >> (RADIX_BITS - t2));
            }
            data_[0] <<= t2;
        }
        this->DeleteLeadingZero();

        assert(size_ <= MAX_SIZE && "Overflow detected.");
    }

    /**
     * this.data = a.data >> number.
     *
     * The number can not be negative.
     */
    void ShiftRightData(const BigInt& a, int number) {
        assert(number >= 0);

        if (a.size_ == 0) {
            size_ = 0;
            return;
        }
        if (number == 0) {
            size_ = a.size_;
            data_ = a.data_;
            return;
        }

        int t1 = number / RADIX_BITS;
        int t2 = number % RADIX_BITS;

        if (t1 >= a.size_) {
            size_ = 0;
            return;
        }

        data_.resize(a.size_);
        int size = a.size_ - t1;
        std::copy(a.data_.begin() + t1, a.data_.begin() + t1 + size,
                  data_.begin());

        if (t2 != 0) {
            for (int i = 0; i < size - 1; ++i) {
                data_[i] = (data_[i] >> t2) |
                           (data_[i + 1] << (RADIX_BITS - t2));
            }
            data_[size - 1] >>= t2;
            size_ = size;
            this->DeleteLeadingZero();
        } else {
            size_ = size;
        }

        sign_ = size_ == 0 ? 0 : a.sign_;
    }

    /**
     * this->data = a.data * b.
     */
    void MultiplyData(const BigInt& a, uint32_t b) {
        // Quick out.
        if (a.size_ == 0 || b == 0) {
            size_ = 0;
            return;
        }
        if (a.size_ == 1 && a.data_[0] == 1) {
            size_ = 1;
            data_[0] = b;
            return;
        }
        if (b == 1) {
            size_ = a.size_;
            data_ = a.data_;
            return;
        }
        if ((b & (b - 1)) == 0) {
            ShiftLeftData(a, Log2(b));
            return;
        }

        data_.resize(a.size_ + 1);

        uint32_t carry = 0;
        uint64_t t1, t2, t3, t4;
        int i = 0;
        for (; i + 4 < a.size_; i += 4) {
            t1 = static_cast<uint64_t>(a.data_[i])     * b;
            t2 = static_cast<uint64_t>(a.data_[i + 1]) * b;
            t3 = static_cast<uint64_t>(a.data_[i + 2]) * b;
            t4 = static_cast<uint64_t>(a.data_[i + 3]) * b;

            data_[i]     = static_cast<uint32_t>(t1 + carry);
            carry = static_cast<uint32_t>((t1 + carry) >> RADIX_BITS);
            data_[i + 1] = static_cast<uint32_t>(t2 + carry);
            carry = static_cast<uint32_t>((t2 + carry) >> RADIX_BITS);
            data_[i + 2] = static_cast<uint32_t>(t3 + carry);
            carry = static_cast<uint32_t>((t3 + carry) >> RADIX_BITS);
            data_[i + 3] = static_cast<uint32_t>(t4 + carry);
            carry = static_cast<uint32_t>((t4 + carry) >> RADIX_BITS);
        }
        for (; i < a.size_; ++i) {
            t1 = static_cast<uint64_t>(a.data_[i]) * b + carry;
            data_[i] = static_cast<uint32_t>(t1);
            carry = static_cast<uint32_t>(t1 >> RADIX_BITS);
        }
        if (carry != 0) {
            data_[a.size_] = carry;
            size_ = a.size_ + 1;
        } else {
            size_ = a.size_;
        }

        assert(size_ <= MAX_SIZE && "Overflow detected.");
    }

    /**
     * this.data = a.data * b.data.
     *
     * Auto choose multiplication algorithm by two BigInt's sizes.
     */
    void MultiplyData(const BigInt& a, const BigInt& b) {
        if (a.size_ > KARATSUBA_MULTIPLY_THRESHOLD &&
            b.size_ > KARATSUBA_MULTIPLY_THRESHOLD) {
            MultiplyKaratsubaData(a, b);
        } else {
            MultiplySimpleData(a, b);
        }
    }

    /**
     * this.data = a.data * b.data.
     *
     * Using baseline/comba algorithm.
     */
    void MultiplySimpleData(const BigInt& a, const BigInt& b) {
        // Quick out.
        if (a.size_ == 0 || b.size_ == 0) {
            size_ = 0;
            return;
        }
        if (a.size_ == 1 && a.data_[0] == 1) {
            size_ = b.size_;
            data_ = b.data_;
            return;
        }
        if (b.size_ == 1 && b.data_[0] == 1) {
            size_ = a.size_;
            data_ = a.data_;
            return;
        }

        int a_size = a.size_;
        int b_size = b.size_;
        int size = a_size + b_size;
        std::vector<uint32_t> data(size, 0);

        // If a == b, special treat.
        if (CompareData(a, b) == 0) {
            uint64_t tmp = 0;
            for (int i = 0; i < size; ++i) {
                int ty = std::min(a_size - 1, i);
                int tx = i - ty;
                int j = std::min(a.size_ - tx, ty + 1);
                j = std::min(j, (ty - tx + 1) >> 1);

                // Use 128-bit to store the result.
                uint64_t low = 0, high = 0;
                uint64_t t = 0;
                for (int k = 0; k < j; ++k) {
                    t = static_cast<uint64_t>(a.data_[tx++]) * a.data_[ty--];
                    if (low + t < low) ++high;
                    low += t;
                }

                high += high;
                if (low + low < low) {
                    ++high;
                }
                low += low;
                if (low + tmp < low) {
                    ++high;
                }
                low += tmp;

                if ((i & 1) == 0) {
                    t = static_cast<uint64_t>(a.data_[i >> 1]) *
                        a.data_[i >> 1];
                    if (low + t < low) {
                        ++high;
                    }
                    low += t;
                }

                data[i] = static_cast<uint32_t>(low);
                tmp = (low >> RADIX_BITS) + (high << RADIX_BITS);
            }
        } else {
            uint64_t tmp = 0;
            for (int i = 0; i < a_size; ++i) {
                uint32_t carry = 0;
                if (a.data_[i] == 0) continue;

                for (int j = 0; j < b_size; ++j) {
                    tmp = static_cast<uint64_t>(a.data_[i]) * b.data_[j] +
                          carry + data[i + j];
                    data[i + j] = static_cast<uint32_t>(tmp);
                    carry = static_cast<uint32_t>(tmp >> RADIX_BITS);
                }
                data[i + b_size] = carry;
            }
        }

        data_.swap(data);
        size_ = size;
        DeleteLeadingZero();

        assert(size_ <= MAX_SIZE && "Overflow detected.");
    }

    /**
     * this.data = a.data * b.data.
     *
     * Using Karatsuba algorithm.
     *
     * This is known as divide-and-conquer and leads to the famous O(N ^ log(3))
     * or O(N ^ 1.584) work which is asymptotically lower than the standard
     * O(N ^ 2) that the baseline/comba methods use.
     */
    void MultiplyKaratsubaData(const BigInt& a, const BigInt& b) {
        int a_size = a.size_;
        int b_size = b.size_;
        int size = std::min(a_size, b_size);

        if (size < KARATSUBA_MULTIPLY_THRESHOLD) {
            MultiplySimpleData(a, b);
            return;
        }

        size >>= 1;

        if (a_size * 2 < b_size || b_size * 2 < a_size) {
            const BigInt *ta = &a, *tb = &b;
            if (a_size > b_size) std::swap(ta, tb);

            BigInt t, sum = 0;
            t.data_.resize(ta->size_);
            for (int i = 0; i < tb->size_; i += ta->size_) {
                t.size_ = std::min(ta->size_, tb->size_ - i);
                for (int j = 0; j < t.size_; ++j) {
                    t.data_[j] = tb->data_[i + j];
                }
                t.DeleteLeadingZero();
                t.MultiplyData(t, *ta);
                t.ShiftLeftData(t, i * RADIX_BITS);
                sum.AddData(sum, t);
            }
            size_ = sum.size_;
            data_ = sum.data_;
            return;
        }

        BigInt x0, x1, y0, y1;
        x0.data_.resize(size); x1.data_.resize(a_size - size);
        y0.data_.resize(size); y1.data_.resize(b_size - size);
        x0.size_ = y0.size_ = size;
        x1.size_ = a_size - size;
        y1.size_ = b_size - size;

        std::copy(a.data_.begin(), a.data_.begin() + size, x0.data_.begin());
        std::copy(b.data_.begin(), b.data_.begin() + size, y0.data_.begin());
        std::copy(a.data_.begin() + size, a.data_.begin() + a_size,
                  x1.data_.begin());
        std::copy(b.data_.begin() + size, b.data_.begin() + b_size,
                  y1.data_.begin());

        x0.DeleteLeadingZero();
        y0.DeleteLeadingZero();

        BigInt x0y0, x1y1, t1;
        x0y0.MultiplyKaratsubaData(x0, y0);  // x0y0 = x0 * y0
        x1y1.MultiplyKaratsubaData(x1, y1);  // x1y1 = x1 * y1
        t1.AddData(x1, x0);                  // t1 = x1 + x0
        x0.AddData(y1, y0);                  // x0 = y1 + y0
        t1.MultiplyKaratsubaData(t1, x0);    // t1 = (x1 + x0) * (y1 + y0)
        x0.AddData(x0y0, x1y1);              // x0 = x0y0 + x1y1
        t1.SubtractData(t1, x0);             // t1 = (x1 + x0) * (y1 + y0) -
                                             //      (x1y1 + x0y0)

        t1.ShiftLeftData(t1, size * RADIX_BITS);
        x1y1.ShiftLeftData(x1y1, (size + size) * RADIX_BITS);

        t1.AddData(x0y0, t1);    // t1 = x0y0 + t1
        this->AddData(x1y1, t1); // *this = x1y1 + t1
    }

    /**
     * this.data = a.data / b.data, remainder = a.data % b.data.
     */
    void DivideData(const BigInt& a, uint32_t b,
                    uint32_t* remainder = nullptr) {
        assert(b != 0);

        // Quick outs.
        if (a.size_ == 0) {
            if (remainder != nullptr) *remainder = 0;
            SetZero();
            return;
        }
        if (b == 1) {
            if (remainder != nullptr) *remainder = 0;
            size_ = a.size_;
            data_ = a.data_;
            return;
        }

        // If b is power of 2, then quick out.
        if ((b & (b - 1)) == 0) {
            if (remainder != nullptr) *remainder = a.data_[0] & (b - 1);
            this->ShiftRightData(a, Log2(b));
            return;
        }

        uint64_t tmp = 0;
        data_.resize(a.size_);
        size_ = a.size_;

        int i = a.size_ - 1;
        for (; i - 4 >= 0; i -= 4) {
            tmp = (tmp << RADIX_BITS) | (a.data_[i]);
            data_[i] = static_cast<uint32_t>(tmp / b);
            tmp -= static_cast<uint64_t>(data_[i]) * b;

            tmp = (tmp << RADIX_BITS) | (a.data_[i - 1]);
            data_[i - 1] = static_cast<uint32_t>(tmp / b);
            tmp -= static_cast<uint64_t>(data_[i - 1]) * b;

            tmp = (tmp << RADIX_BITS) | (a.data_[i - 2]);
            data_[i - 2] = static_cast<uint32_t>(tmp / b);
            tmp -= static_cast<uint64_t>(data_[i - 2]) * b;

            tmp = (tmp << RADIX_BITS) | (a.data_[i - 3]);
            data_[i - 3] = static_cast<uint32_t>(tmp / b);
            tmp -= static_cast<uint64_t>(data_[i - 3]) * b;
        }
        for (; i >= 0; --i) {
            tmp = (tmp << RADIX_BITS) | (a.data_[i]);
            data_[i] = static_cast<uint32_t>(tmp / b);
            tmp -= static_cast<uint64_t>(data_[i]) * b;
        }
        if (remainder != nullptr) *remainder = static_cast<uint32_t>(tmp);

        DeleteLeadingZero();
    }

    /**
     * this->data = a.data / b.data, remainder = a.data % b.data.
     */
    void DivideData(const BigInt& a, const BigInt& b,
                    BigInt* remainder = nullptr) {
        assert(b.size_ != 0);

        // Quick outs.
        if (a.size_ == 0) {
            if (remainder != nullptr) remainder->size_ = 0;
            size_ = 0;
            return;
        }
        if (b.size_ == 1 && b.data_[0] == 1) {
            if (remainder != nullptr) remainder->size_ = 0;
            size_ = a.size_;
            data_ = a.data_;
            return;
        }

        // If a < b then quotient = 0, remainder = dividend.
        if (CompareData(a, b) < 0) {
            if (remainder != nullptr) {
                remainder->size_ = a.size_;
                remainder->data_ = a.data_;
            }
            size_ = 0;
            return;
        }

        if (b.size_ == 1) {
            uint32_t remainder_t;
            DivideData(a, b.data_[0], &remainder_t);
            if (remainder != nullptr) {
                if (remainder_t == 0) {
                    remainder->size_ = 0;
                } else {
                    remainder->size_ = 1;
                    remainder->data_[0] = remainder_t;
                }
            }
            return;
        }

        BigInt x = a, y = b;

        // Normalize: shift y left so that its top digit is >= RADIX_BITS / 2.
        //            Shift x left by the same amount.
        int norm = Log2(y.data_[y.size_ - 1]) + 1;
        if (norm < RADIX_BITS - 1) {
            norm = (RADIX_BITS - 1) - norm;
            x.ShiftLeft(x, norm);
            y.ShiftLeft(y, norm);
        } else {
            norm = 0;
        }

        int x_back = x.size_ - 1;
        int y_back = y.size_ - 1;
        int diff = x_back - y_back;

        y.ShiftLeftData(y, diff * RADIX_BITS);

        std::vector<uint32_t> q(a.size_ + 2, 0);
        while (CompareData(x, y) >= 0) {
            ++(q[diff]);
            x.SubtractData(x, y);
        }

        // Reset y by shifting it back down.
        y.ShiftRightData(y, diff * RADIX_BITS);

        BigInt t;
        std::vector<uint32_t> x_low(x.size_, 0);

        for (int i = x_back; i > y_back; --i) {
            if (i > x.size_) {
                continue;
            }

            int y_offset = i - y_back - 1;
            uint32_t& q_t = q[y_offset];

            if (x.data_[i] == y.data_[y_back]) {
                q_t = UINT32_MAX;
            } else {
                uint64_t tmp = static_cast<uint64_t>(x.data_[i]) << RADIX_BITS;
                tmp |= x.data_[i - 1];
                tmp /= y.data_[y_back];
                if (tmp > UINT32_MAX) tmp = UINT32_MAX;
                q_t = static_cast<uint32_t>(tmp);
            }

            BigInt t1, t2;
            t1.data_.resize(3);
            t2.data_.resize(3);
            ++q_t;
            do {
                --q_t;

                // Find left hand.
                t1.size_ = 2;
                t1.data_[0] = (y_back - 1 < 0) ? 0 : y.data_[y_back - 1];
                t1.data_[1] = y.data_[y_back];
                t1.MultiplyData(t1, q_t);

                // Find right hand.
                t2.size_ = 3;
                t2.data_[0] = (i - 2 < 0) ? 0 : x.data_[i - 2];
                t2.data_[1] = (i - 1 < 0) ? 0 : x.data_[i - 1];
                t2.data_[2] = x.data_[i];
            } while (CompareData(t1, t2) > 0);

            t.MultiplyData(y, q_t);
            for (int j = 0; j < y_offset; ++j) {
                x_low[j] = x.data_[j];
            }
            x.ShiftRightData(x, RADIX_BITS * y_offset);

            if (CompareData(x, t) < 0) {
                // This branch taken rarely.
                x.AddData(x, y);
                x.SubtractData(x, t);

                // Reset x by left shift and add the low part of x.
                x.data_.resize(x.size_ + y_offset);
                std::copy_backward(x.data_.begin(), x.data_.begin() + x.size_,
                                   x.data_.begin() + x.size_ + y_offset);
                std::copy(x_low.begin(), x_low.begin() + y_offset,
                          x.data_.begin());
                x.size_ += y_offset;
                --q_t;
            } else {
                x.SubtractData(x, t);

                // Reset x by left shift and add the low part of x.
                x.data_.resize(x.size_ + y_offset);
                std::copy_backward(x.data_.begin(), x.data_.begin() + x.size_,
                                   x.data_.begin() + x.size_ + y_offset);
                std::copy(x_low.begin(), x_low.begin() + y_offset,
                          x.data_.begin());
                x.size_ += y_offset;
            }
        }

        if (remainder != nullptr) {
            remainder->ShiftRightData(x, norm);
            remainder->DeleteLeadingZero();
        }

        size_ = a.size_ + 2;
        data_.swap(q);
        this->DeleteLeadingZero();
    }

    int sign_ = 0;               // The sign of BigInt.
    int size_ = 0;               // The used size of data.
    std::vector<uint32_t> data_; // The reverse data for BigInt.
};

/**
 * Random generator for BigInt.
 */
class BigIntRandomGenerator {
public:
    explicit BigIntRandomGenerator(int seed = 0)
        : seed_(seed) {}

    void set_seed(int seed) {
        seed_ = seed;
    }

    /**
     * Generate a random BigInt between [0, n].
     */
    BigInt operator() (const BigInt& n) const {
        return Generate(n);
    }

    /**
     * Generate a random BigInt between [0, n].
     */
    BigInt Generate(const BigInt& n) const {
        assert(n >= 0);

        if (n == 0) return BigInt(0);

        uint32_t back = n.data()[n.size() - 1];
        std::uniform_int_distribution<uint32_t> uniform(0, back);

        BigInt res;
        while (true) {
            uint32_t r = uniform(random_);
            bool flag = (r == back);
            res = r;

            int i = 0;
            for (i = n.size() - 2; i >= 0; --i) {
                uint32_t tmp = random_();
                if (flag) {
                    if (tmp > n.data()[i]) break;
                    if (tmp < n.data()[i]) flag = false;
                }

                res <<= 32;
                res += tmp;
            }
            if (i < 0) break;
        }

        return res;
    }

protected:
    uint32_t GenerateInt(uint32_t n) const {
        std::uniform_int_distribution<uint32_t> uniform(0, n);
        return uniform(random_);
    }

    int seed_;
    mutable std::mt19937 random_;
};

} // namespace cl

#endif // MATH_NUMBER_BIGINT_H_
