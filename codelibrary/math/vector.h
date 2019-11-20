//
// Copyright 2019 Yangbin Lin. All Rights Reserved.
//
// Author: yblin@jmu.edu.cn (Yangbin Lin)
//
// This file is part of the Code Library.
//

#ifndef MATH_VECTOR_H_
#define MATH_VECTOR_H_

#include <cassert>
#include <iterator>

#include "codelibrary/base/format.h"
#include "codelibrary/math/basic_linear_algebra.h"

namespace cl {

template <typename T>
class Vector : public BasicLinearAlgebra<T> {
public:
    using value_type             = T;
    using iterator               = T*;
    using const_iterator         = const T*;
    using reverse_iterator       = std::reverse_iterator<iterator>;
    using const_reverse_iterator = std::reverse_iterator<const_iterator>;

    /**
     * Default vector constructor.
     */
    Vector() : BasicLinearAlgebra<T>() {}

    /**
     * Construct a vector without initialization.
     */
    explicit Vector(int size) : BasicLinearAlgebra<T>(size) {}

    /**
     * Construct a vector with initial value.
     */
    Vector(int size, const T& v) : BasicLinearAlgebra<T>(size, v) {}

    /**
     * Construct vector from data in [first, last).
     *
     * The second template parameter is used to distinguish this function to
     * Vector(int, int).
     */
    template <typename Iter,
              typename = typename std::enable_if<std::is_convertible<
                         typename std::iterator_traits<Iter>::iterator_category,
                                  std::input_iterator_tag>::value>::type>
    Vector(Iter first, Iter last)
        : BasicLinearAlgebra<T>(first, last) {}

    /**
     * Build vector from initializer list.
     *
     * Note that this constructor can not be explicit.
     */
    Vector(std::initializer_list<T> list)
        : BasicLinearAlgebra<T>(list) {}

    virtual ~Vector() = default;

    Vector& operator *=(const T& rhs) {
        blas::Scale(this->size_, rhs, this->data_);
        return *this;
    }

    Vector& operator +=(const Vector& rhs) {
        assert(this->size_ == rhs.size_);
        blas::Add(this->size_, this->data_, rhs.data(), this->data_);
        return *this;
    }

    Vector& operator -=(const Vector& rhs) {
        assert(this->size_ == rhs.size_);

        blas::Subtract(this->size_, this->data_, rhs.data_, this->data_);
        return *this;
    }

    Vector& operator *=(const Vector& rhs) {
        assert(this->size_ == rhs.size_);

        blas::Multiply(this->size_, this->data_, rhs.data_, this->data_);
        return *this;
    }

    bool operator <(const Vector& rhs) const {
        if (this->size_ != rhs.size_) return this->size_ < rhs.size_;
        return std::lexicographical_compare(this->data_,
                                            this->data_ + this->size_,
                                            rhs.data_,
                                            rhs.data_ + this->size_);
    }

    bool operator >(const Vector& rhs) const {
        return rhs < *this;
    }

    bool operator <=(const Vector& rhs) const {
        return !(rhs < *this);
    }

    bool operator >=(const Vector& rhs) const {
        return !(*this < rhs);
    }

    /**
     * This = Vector(n, v).
     */
    void assign(int n, const T& v) {
        this->Assign(n, v);
    }

    /**
     * Assign vector from data in [first, last).
     * The second template parameter is used to distinguish this function to
     * assign(int, int).
     */
    template <typename Iter,
              typename = typename std::enable_if<std::is_convertible<
                         typename std::iterator_traits<Iter>::iterator_category,
                                  std::input_iterator_tag>::value>::type>
    void assign(Iter first, Iter last) {
        this->Assign(first, last);
    }

    /**
     * Resize the vector.
     * 
     * Note that, resize(n) is much faster than resize(n, v).
     */
    void resize(int n) {
        this->Resize(n);
    }

    /**
     * Resize the vector with the filling value.
     */
    void resize(int n, const T& v) {
        this->Resize(n, v);
    }

    /**
     * Return euclidean norm.
     */
    T norm() const {
        return blas::EuclideanNorm(this->size(), this->data());
    }

    /**
     * Return the square of the norm.
     */
    T squared_norm() const {
        return blas::SquaredEuclideanNorm(this->size(), this->data());
    }

    /**
     * Infinity norm.
     */
    T infinty_norm() const {
        return blas::InfinityNorm(this->size(), this->data());
    }

    /**
     * Compute the negative of a vector.
     */
    friend Vector operator -(const Vector& rhs) {
        Vector res = rhs;
        blas::Negate(res.size(), res.data());
        return res;
    }

    /**
     * Compute the product of a vector by a scalar.
     */
    friend Vector operator *(const Vector& lhs, const T& rhs) {
        Vector res = lhs;
        blas::Scale(res.size(), rhs, res.data());
        return res;
    }

    /**
     * Compute the product of a scale by a vector.
     */
    friend Vector operator *(const T& lhs, const Vector& rhs) {
        Vector res = rhs;
        blas::Scale(res.size(), lhs, res.data());
        return res;
    }

    /**
     * Perform element by element addition of vector a and vector b.
     */
    friend Vector operator +(const Vector& a, const Vector& b) {
        assert(a.size() == b.size());

        Vector c(a.size());
        blas::Add(c.size(), a.data(), b.data(), c.data());
        return c;
    }

    /**
     * Perform element by element subtraction of vector b from vector a.
     */
    friend Vector operator -(const Vector& a, const Vector& b) {
        assert(a.size() == b.size());

        Vector c(a.size());
        blas::Subtract(c.size(), a.data(), b.data(), c.data());
        return c;
    }

    /**
     * Perform element by element multiplication of vector a and vector b.
     */
    template <typename T1>
    friend Vector<T1> ElementwiseMultiply(const Vector<T1>& a,
                                          const Vector<T1>& b);

    friend T operator*(const Vector& a, const Vector& b) {
        assert(a.size() == b.size());
        return blas::DotProduct(a.size(), a.data(), b.data());
    }

    friend std::ostream& operator <<(std::ostream& os, const Vector& rhs) {
        os << Printer::GetInstance()->AlignFormat(rhs.begin(), rhs.end());
        return os;
    }
};

template <typename T>
Vector<T> ElementwiseMultiply(const Vector<T>& a, const Vector<T>& b) {
    assert(a.size() == b.size());

    Vector<T> c(a.size());
    blas::Multiply(c.size(), a.data(), b.data(), c.data());
    return c;
}

using FVector = Vector<float>;
using RVector = Vector<double>;

} // namespace cl

#endif // MATH_VECTOR_H_
