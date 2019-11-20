//
// Copyright 2019 Yangbin Lin. All Rights Reserved.
//
// Author: yblin@jmu.edu.cn (Yangbin Lin)
//
// This file is part of the Code Library.
//

#ifndef MATH_MATRIX_MATRIX_H_
#define MATH_MATRIX_MATRIX_H_

#include <random>

#include "codelibrary/base/array.h"
#include "codelibrary/base/array_nd.h"
#include "codelibrary/base/format.h"
#include "codelibrary/math/basic_linear_algebra.h"
#include "codelibrary/math/vector.h"

namespace cl {

/**
 * Row-major matrix.
 */
template <typename T>
class Matrix : public BasicLinearAlgebra<T> {
    static const int BLAS_N = 256;

public:
    using value_type             = T;
    using iterator               = T*;
    using const_iterator         = const T*;
    using reverse_iterator       = std::reverse_iterator<iterator>;
    using const_reverse_iterator = std::reverse_iterator<const_iterator>;

    Matrix() = default;

    Matrix(int n_rows, int n_columns, const T& v = 0)
        : BasicLinearAlgebra<T>(n_rows * n_columns, v),
          n_rows_(n_rows), n_columns_(n_columns) {
        CheckDimension(n_rows_, n_columns_);
    }

    template <typename InputIter>
    Matrix(int n_rows, int n_columns, InputIter first, InputIter last)
        : BasicLinearAlgebra<T>(first, last),
          n_rows_(n_rows), n_columns_(n_columns) {
        assert(n_rows_ * n_columns_ == this->size_);
    }

    /**
     * Assign a matrix where all elements are equal to zero, except for the
     * diagonal, whose values are equal to one.
     */
    void eye(int n) {
        assert(n >= 0);

        this->assign(n, n, 0);
        for (int i = 0; i < n; ++i) {
            this->data_[i * n + i] = 1;
        }
    }

    /**
     * Create a square diagonal matrix with the elements of vector v on the main
     * diagonal.
     */
    void diag(const Vector<T>& v) {
        int n = v.size();
        CheckDimension(n, n);

        this->assign(n, n, 0);
        int k = 0;
        for (int i = 0; i < this->size_; i += n) {
            this->data_[i] = v[k++];
        }
    }

    /**
     * Create a random matrix. The value of each element is in [0, 1]
     */
    void rand(int n, int m, int seed = 0) {
        CheckDimension(n, m);

        this->resize(n, m);
        std::uniform_real_distribution<T> random(0.0, 1.0);
        std::mt19937 engine(seed);
        for (int i = 0; i < this->size_; ++i) {
            this->data_[i] = random(engine);
        }
    }

    bool operator ==(const Matrix& rhs) const {
        return n_rows_ == rhs.n_rows_ && n_columns_ == rhs.n_columns_ &&
               std::equal(this->data_, this->data_ + this->size_, rhs.data_);
    }

    bool operator !=(const Matrix& rhs) const {
        return !(*this == rhs);
    }
    
    /**
     * Resize the matrix with the given value.
     */
    void resize(int n_rows, int n_columns) {
        CheckDimension(n_rows, n_columns);

        if (n_columns == n_columns_) {
            n_rows_ = n_rows;
            this->Resize(n_rows_ * n_columns_);
        } else {
            Array<T> tmp(this->begin(), this->end());
            this->Resize(n_rows * n_columns);

            int min_n_rows = std::min(n_rows_, n_rows);
            int min_n_columns = std::min(n_columns_, n_columns);
            for (int i = 0; i < min_n_rows; ++i) {
                std::memcpy(this->data() + i * n_columns,
                            tmp.begin() + i * n_columns_,
                            sizeof(T) * min_n_columns);
            }
            n_rows_ = n_rows;
            n_columns_ = n_columns;
        }
    }

    /**
     * Resize the matrix with the given value.
     */
    void resize(int n_rows, int n_columns, const T& value) {
        CheckDimension(n_rows, n_columns);

        if (n_columns == n_columns_) {
            n_rows_ = n_rows;
            this->Resize(n_rows_ * n_columns_, value);
        } else {
            Array<T> tmp(this->begin(), this->end());
            this->Assign(n_rows * n_columns, value);

            int min_n_rows = std::min(n_rows_, n_rows);
            int min_n_columns = std::min(n_columns_, n_columns);
            for (int i = 0; i < min_n_rows; ++i) {
                std::memcpy(this->data() + i * n_columns,
                            tmp.begin() + i * n_columns_,
                            sizeof(T) * min_n_columns);
            }
            n_rows_ = n_rows;
            n_columns_ = n_columns;
        }
    }

    /**
     * Assign the matrix.
     */
    void assign(int m, int n, const T& v) {
        resize(m, n);
        this->fill(v);
    }

    /**
     * this = [first, last).
     */
    template <typename Iterator>
    void assign(int m, int n, Iterator first, Iterator last) {
        resize(m, n);
        Iterator p = first;
        for (int i = 0; i < this->size_; ++i, ++p) {
            this->data_[i] = *p;
            assert(p != last);
        }
    }

    int n_rows() const {
        return n_rows_;
    }

    int n_columns() const {
        return n_columns_;
    }

    const T& operator() (int i, int j) const {
        return this->data_[i * n_columns_ + j];
    }

    T& operator() (int i, int j) {
        return this->data_[i * n_columns_ + j];
    }

    Matrix& operator *=(const T& rhs) {
        blas::Scale(this->size_, rhs, this->data_);
        return *this;
    }

    Matrix& operator +=(const Matrix& rhs) {
        blas::Add(this->size_, this->data_, rhs.data_, this->data_);
        return *this;
    }

    Matrix& operator -=(const Matrix& rhs) {
        blas::Subtract(this->size_, this->data_, rhs.data_, this->data_);
        return *this;
    }

    Matrix& operator *=(const Matrix& rhs) {
        Array<T> tmp(this->n_rows_ * rhs.n_columns_);
        blas::GEMM(this->n_rows_, this->n_columns_, rhs.n_columns_,
                   this->data_, rhs.data_, tmp.data_);
        this->assign(this->n_rows_, rhs.n_columns_, tmp.begin(), tmp.end());
        return *this;
    }

    Vector<T> get_column(int index) const {
        assert(index >= 0 && index < n_columns_);

        Vector<T> col(n_rows_);
        for (int i = 0; i < n_rows_; ++i, index += n_columns_) {
            col[i] = this->data_[index];
        }
        return col;
    }

    Vector<T> get_row(int index) const {
        assert(index >= 0 && index < n_rows_);

        Vector<T> row(n_columns_);
        index = n_rows_ * (index - 1);
        for (int i = 0; i < n_columns_; ++i) {
            row[i] = this->data_[index++];
        }
        return row;
    }

    /**
     * Compute the negative of a matrix.
     */
    friend Matrix operator -(const Matrix& rhs) {
        Matrix res = rhs;
        blas::Negate(res.size_, res.data_);
        return res;
    }

    /**
     * Compute the product of a matrix by a scalar.
     */
    friend Matrix operator *(const Matrix& lhs, const T& rhs) {
        Matrix res = lhs;
        blas::Scale(res.size_, rhs, res.data_);
        return res;
    }

    /**
     * Compute the product of a scale by a matrix.
     */
    friend Matrix operator *(const T& lhs, const Matrix& rhs) {
        Matrix res = rhs;
        blas::Scale(res.size_, lhs, res.data_);
        return res;
    }

    /**
     * Perform element by element addition of matrix a and matrix b.
     */
    friend Matrix operator +(const Matrix& a, const Matrix& b) {
        assert(a.n_rows_ == b.n_rows_);
        assert(a.n_columns_ == b.n_columns_);

        Matrix c(a.n_rows_, a.n_columns_);
        blas::Add(c.size_, a.data_, b.data_, c.data_);
        return c;
    }

    /**
     * Perform element by element subtraction of matrix b from matrix a.
     */
    friend Matrix operator -(const Matrix& a, const Matrix& b) {
        assert(a.n_rows_ == b.n_rows_);
        assert(a.n_columns_ == b.n_columns_);

        Matrix c(a.n_rows_, a.n_columns_);
        blas::Subtract(c.size_, a.data_, b.data_, c.data_);
        return c;
    }

    /**
     * Return Ab.
     */
    friend Vector<T> operator *(const Matrix<T>& a, const Vector<T>& b) {
        assert(a.n_columns_ == b.size());

        Vector<T> c(a.n_rows_);
        blas::GEMV(a.n_rows_, a.n_columns_, a.data_, b.data(), c.data());
        return c;
    }

    /**
     * Return A'b.
     */
    friend Vector<T> operator *(const Vector<T>& b, const Matrix<T>& a) {
        assert(a.n_rows_ == b.size());

        Vector<T> c(a.n_columns_);
        blas::GEMVTrans(a.n_rows_, a.n_columns_, a.data_, b.data(),
                        c.data());
        return c;
    }

    /**
     * Return AB.
     */
    friend Matrix operator *(const Matrix& a, const Matrix& b) {
        assert(a.n_columns_ == b.n_rows_);

        Matrix c(a.n_rows_, b.n_columns_);
        blas::GEMM(a.n_rows_, a.n_columns_, b.n_columns_, a.data_, b.data_,
                   c.data_);
        return c;
    }

    friend std::ostream& operator <<(std::ostream& os, const Matrix& rhs) {
        ArrayND<int> t(rhs.n_rows_, rhs.n_columns_);
        t.set_data(rhs.begin(), rhs.end());
        return os << t;
    }

protected:
    /**
     * Check the dimension.
     */
    static void CheckDimension(int n_rows, int n_columns) {
        assert(n_rows >= 0);
        assert(n_columns >= 0);
        assert((n_columns == 0 || n_rows <= INT_MAX / n_columns) &&
               "The given dimensions of matrix are too large.");
    }

    // Number of rows in the matrix.
    int n_rows_ = 0;
    // Number of columns in the matrix.
    int n_columns_ = 0;
};

using IMatrix = Matrix<int>;
using FMatrix = Matrix<float>;
using RMatrix = Matrix<double>;

} // namespace cl

#endif // MATH_MATRIX_MATRIX_H_
