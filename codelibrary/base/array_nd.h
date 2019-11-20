//
// Copyright 2019 Yangbin Lin. All Rights Reserved.
//
// Author: yblin@jmu.edu.cn (Yangbin Lin)
//
// This file is part of the Code Library.
//

#ifndef BASE_ARRAY_ND_H_
#define BASE_ARRAY_ND_H_

#include <cassert>
#include <ostream>

#include "codelibrary/base/array.h"
#include "codelibrary/base/format.h"

namespace cl {

template <typename T>
class ArrayND {
public:
    using value_type             = T;
    using pointer                = T*;
    using const_pointer          = const T *;
    using reference              = T&;
    using const_reference        = const T &;
    using iterator               = T*;
    using const_iterator         = const T*;
    using reverse_iterator       = std::reverse_iterator<iterator>;
    using const_reverse_iterator = std::reverse_iterator<const_iterator>;
    using size_type              = int;
    using difference_type        = int;

    ArrayND() = default;

    explicit ArrayND(int size) {
        reshape(size);
    }

    ArrayND(int size1, int size2) {
        reshape(size1, size2);
    }

    ArrayND(int size1, int size2, int size3) {
        reshape(size1, size2, size3);
    }

    explicit ArrayND(const Array<int>& shape) {
        reshape(shape);
    }

    const T& operator()(int a1) const {
        return data_[a1];
    }

    T& operator()(int a1) {
        return data_[a1];
    }

    const T& operator()(int a1, int a2) const {
        return data_[a1 * strides_[0] + a2];
    }

    T& operator()(int a1, int a2) {
        return data_[a1 * strides_[0] + a2];
    }

    const T& operator()(int a1, int a2, int a3) const {
        return data_[a1 * strides_[0] + a2 * strides_[1] + a3];
    }

    T& operator()(int a1, int a2, int a3) {
        return data_[a1 * strides_[0] + a2 * strides_[1] + a3];
    }

    const T& operator()(int a1, int a2, int a3, int a4) const {
        return data_[a1 * strides_[0] + a2 * strides_[1] +
                     a3 * strides_[2] + a4];
    }

    T& operator()(int a1, int a2, int a3, int a4) {
        return data_[a1 * strides_[0] + a2 * strides_[1] +
                     a3 * strides_[2] + a4];
    }

    const T& operator[](int index) const {
        return data_[index];
    }

    T& operator[](int index) {
        return data_[index];
    }

    /**
     * Reshape the ND array.
     */
    void reshape(int a1) {
        assert(a1 > 0);

        size_ = a1;
        shape_.resize(1);
        shape_[0] = a1;
        data_.resize(size_, T());
        strides_ = {1};
    }
    void reshape(int a1, int a2) {
        reshape({ a1, a2 }); 
    }
    void reshape(int a1, int a2, int a3) { 
        reshape({ a1, a2, a3 }); 
    }
    void reshape(int a1, int a2, int a3, int a4) { 
        reshape({ a1, a2, a3, a4 }); 
    }
    void reshape(const Array<int>& shape) {
        assert(!shape.empty());

        size_ = 1;
        for (int s : shape) {
            assert(s > 0);
            assert(size_ <= INT_MAX / s && "The given shape is too large.");

            size_ *= s;
        }
 
        shape_ = shape;
        data_.resize(size_, T());
        strides_.resize(shape_.size());
        strides_.back() = 1;
        for (int i = shape_.size() - 2; i >= 0; --i) {
            strides_[i] = strides_[i + 1] * shape_[i + 1];
        }
    }

    /**
     * Fill the array.
     */
    void fill(const T& v) {
        std::fill(data_.begin(), data_.end(), v);
    }

    iterator begin() {
        return data_.begin();
    }

    iterator end() {
        return data_.end();
    }

    const_iterator begin() const {
        return data_.begin();
    }

    const_iterator end() const {
        return data_.end();
    }

    reverse_iterator rbegin() {
        return reverse_iterator(end());
    }

    reverse_iterator rend() {
        return reverse_iterator(begin());
    }

    const_reverse_iterator rbegin() const {
        return const_reverse_iterator(end());
    }

    const_reverse_iterator rend() const {
        return const_reverse_iterator(begin());
    }

    void clear() {
        size_ = 0;
        shape_.clear();
        data_.clear();
        strides_.clear();
    }

    int n_dimension() const {
        return shape_.size();
    }

    const Array<int>& strides() const {
        return strides_;
    }

    const Array<int>& shape() const {
        return shape_;
    }

    int shape(int i) const {
        assert(i >= 0 && i < shape_.size());

        return shape_[i];
    }

    bool empty() const {
        return size_ == 0;
    }

    int size() const {
        return size_;
    }

    T* data() {
        return data_.data();
    }

    const T* data() const {
        return data_.data();
    }

    void swap(ArrayND* t) {
        assert(t);

        std::swap(size_, t->size_);
        shape_.swap(t->shape_);
        strides_.swap(t->strides_);
        data_.swap(t->data_);
    }

    void set_data(const Array<T>& data) {
        assert(data.size() == size_);

        data_ = data;
    }

    template <typename Iterator>
    void set_data(Iterator first, Iterator last) {
        auto n = std::distance(first, last);
        assert(n == size_);

        data_.assign(first, last);
    }

    /**
     * For debug.
     */
    friend std::ostream& operator <<(std::ostream& os, const ArrayND& rhs) {
        const int max_n_print_elements = 6;
        int width = ComputeWidth(rhs.data(), rhs.shape(), rhs.strides(),
                                 0, max_n_print_elements);

        Print(os, rhs.data(), rhs.shape(), rhs.strides(), 0,
              max_n_print_elements, width);
        return os;
    }

protected:
    /**
     * Recursively to print the ND array.
     */
    static void Print(std::ostream& os,
                      const T* data,
                      const Array<int>& shape,
                      const Array<int>& strides,
                      int depth,
                      int max_print,
                      int width) {
        if (depth == shape.size()) return;

        int n = shape[depth];
        if (depth + 1 == shape.size()) {
            // Print 1D sequence.
            os << Printer::GetInstance()->AlignFormat(data, data + n, width,
                                                      depth, depth);
            return;
        }

        std::string space(depth + 1, ' ');

        os << "[";

        if (n <= max_print || max_print == 0) {
            for (int i = 0; i < n; ++i) {
                if (i != 0) os << space;
                Print(os, data + strides[depth] * i, shape, strides, depth + 1,
                      max_print, width);
                if (i + 1 < n) {
                    os << ",\n";
                    if (depth + 2 < shape.size()) os << "\n";
                }
            }
            if (depth + 2 < shape.size() && depth != 0) os << "\n";
        } else {
            for (int i = 0; i < max_print / 2; ++i) {
                if (i != 0) os << space;
                Print(os, data + strides[depth] * i, shape, strides, depth + 1,
                      max_print, width);
                os << ",\n";
                if (depth + 2 < shape.size()) os << "\n";
            }

            os << space;
            os << "...,\n";
            if (depth + 2 < shape.size()) os << "\n";

            for (int i = n - (max_print + 1) / 2; i < n; ++i) {
                os << space;
                Print(os, data + strides[depth] * i, shape, strides, depth + 1,
                      max_print, width);
                if (i + 1 < n) {
                    os << ",\n";
                    if (depth + 2 < shape.size()) os << "\n";
                }
            }
            if (depth + 2 < shape.size() && depth != 0) os << "\n";
        }

        os << "]";
    }

    /**
     * Recursively to compute width for printing.
     */
    static int ComputeWidth(const T* data,
                            const Array<int>& shape,
                            const Array<int>& strides,
                            int depth,
                            int max_print) {
        if (depth == shape.size()) return 0;

        int n = shape[depth];
        if (depth + 1 == shape.size()) {
            int width = 0;
            if (n <= max_print || max_print == 0) {
                for (const T* p = data; p != data + n; ++p) {
                    std::string str = fmt::format("{}", *p);
                    assert(str.size() <= INT_MAX);
                    width = std::max(width, static_cast<int>(str.size()));
                }
            } else {
                for (int i = 0; i < (max_print + 1) / 2; ++i) {
                    std::string str = fmt::format("{}", data[i]);
                    assert(str.size() <= INT_MAX);
                    width = std::max(width, static_cast<int>(str.size()));
                }
                for (int i = (max_print + 1) / 2; i < n; ++i) {
                    std::string str = fmt::format("{}", data[i]);
                    assert(str.size() <= INT_MAX);
                    width = std::max(width, static_cast<int>(str.size()));
                }
            }
            return width;
        }

        int width = 0;
        if (n <= max_print || max_print == 0) {
            for (int i = 0; i < n; ++i) {
                width = std::max(width, ComputeWidth(data + strides[depth] * i,
                                                     shape, strides, depth + 1,
                                                     max_print));
            }
        } else {
            for (int i = 0; i < max_print / 2; ++i) {
                width = std::max(width, ComputeWidth(data + strides[depth] * i,
                                                     shape, strides, depth + 1,
                                                     max_print));
            }
            for (int i = n - (max_print + 1) / 2; i < n; ++i) {
                width = std::max(width, ComputeWidth(data + strides[depth] * i,
                                                     shape, strides, depth + 1,
                                                     max_print));
            }
        }

        return width;
    }

    // The total number of the ND array.
    int size_;

    // The number of elements in each dimension.
    Array<int> shape_;

    // Array of steps in each dimension when traversing an array.
    // strides_[i] = shape_[i + 1] * ... * shape_[nd - 1].
    Array<int> strides_;

    // The data of the ND array.
    Array<T> data_;
};

} // namespace cl

#endif // BASE_ARRAY_ND_H_