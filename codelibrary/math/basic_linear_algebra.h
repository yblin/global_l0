//
// Copyright 2019 Yangbin Lin. All Rights Reserved.
//
// Author: yblin@jmu.edu.cn (Yangbin Lin)
//
// This file is part of the Code Library.
//

#ifndef MATH_BASIC_LINEAR_ALGEBRA_H_
#define MATH_BASIC_LINEAR_ALGEBRA_H_

#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <cstring>

namespace cl {
namespace blas {

/**
 * Compute a vector-vector dot product.
 *
 * Input:
 *  n - Specify the number of elements in vectors x and y.
 *  a - Array, size at least n.
 *  b - Array, size at least n.
 *
 * Return:
 *  The result of the dot product of a and b, if n is positive.
 *  Otherwise, return 0.
 */
template <typename T>
T DotProduct(int n, const T* a, const T* b) {
    if (n <= 0) return 0;

    switch (n) {
    case 1: return a[0] * b[0];
    case 2: return a[0] * b[0] + a[1] * b[1];
    case 3: return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    case 4: return a[0] * b[0] + a[1] * b[1] + a[2] * b[2] + a[3] * b[3];
    }

    T s = 0;
    const T* p1 = a;
    const T* p2 = b;

    int i = 0;
    for (; i + 4 < n; p1 += 4, p2 += 4, i += 4) {
        s += *(p1)     * *(p2);
        s += *(p1 + 1) * *(p2 + 1);
        s += *(p1 + 2) * *(p2 + 2);
        s += *(p1 + 3) * *(p2 + 3);
    }
    for (; i < n; ++p1, ++p2, ++i)
        s += *(p1) * *(p2);
    return s;
}

/**
 * Compute the Euclidean norm of a vector.
 *
 * Input:
 *  n - Specify the number of elements in vector a.
 *  a - Array, size at least n.
 *
 * Return:
 *  The Euclidean norm of the vector a.
 */
template <typename T>
T EuclideanNorm(int n, const T* a) {
    switch (n) {
    case 0: return 0;
    case 1: return std::sqrt(a[0] * a[0]);
    case 2: return std::sqrt(a[0] * a[0] + a[1] * a[1]);
    case 3: return std::sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
    case 4: return std::sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2] +
                             a[3] * a[3]);
    }

    T sum = 0;
    const T* p = a;

    int i = 0;
    for (; i + 4 < n; p += 4, i += 4) {
        sum += *(p)     * *(p);
        sum += *(p + 1) * *(p + 1);
        sum += *(p + 2) * *(p + 2);
        sum += *(p + 3) * *(p + 3);
    }
    for (; i < n; ++p, ++i)
        sum += *(p) * *(p);

    return std::sqrt(sum);
}

/**
 * Compute the square of Euclidean norm of a vector.
 *
 * Input:
 *  n - Specify the number of elements in vector a.
 *  a - Array, size at least n.
 *
 * Return:
 *  The square Euclidean norm of the vector a.
 */
template <typename T>
T SquaredEuclideanNorm(int n, const T* a) {
    switch (n) {
    case 0: return 0;
    case 1: return a[0] * a[0];
    case 2: return a[0] * a[0] + a[1] * a[1];
    case 3: return a[0] * a[0] + a[1] * a[1] + a[2] * a[2];
    case 4: return a[0] * a[0] + a[1] * a[1] + a[2] * a[2] + a[3] * a[3];
    }

    T sum = 0;
    const T* p = a;

    int i = 0;
    for (; i + 4 < n; p += 4, i += 4) {
        sum += *(p)     * *(p);
        sum += *(p + 1) * *(p + 1);
        sum += *(p + 2) * *(p + 2);
        sum += *(p + 3) * *(p + 3);
    }
    for (; i < n; ++p, ++i)
        sum += *(p) * *(p);

    return sum;
}

/**
 * Compute the infinity norm of a vector.
 *
 * Input:
 *  n - Specify the number of elements in vector a.
 *  a - Array, size at least n.
 *
 * Return:
 *  The infinity norm of the vector a.
 */
template <typename T>
T InfinityNorm(int n, const T* a) {
    switch (n) {
    case 0: return 0;
    case 1: return std::fabs(a[0]);
    case 2: return std::max(std::fabs(a[0]), std::fabs(a[1]));
    case 3: return std::max(std::fabs(a[0]), 
                            std::max(std::fabs(a[1]), std::fabs(a[2])));
    case 4: return std::max(std::max(std::fabs(a[0]), std::fabs(a[1])),
                            std::max(std::fabs(a[2]), std::fabs(a[3])));
    }

    T max_a = 0;
    const T* p = a;

    int i = 0;
    for (; i + 4 < n; p += 4, i += 4) {
        max_a = std::max(std::fabs(*(p)),     max_a);
        max_a = std::max(std::fabs(*(p + 1)), max_a);
        max_a = std::max(std::fabs(*(p + 2)), max_a);
        max_a = std::max(std::fabs(*(p + 3)), max_a);
    }
    for (; i < n; ++p, ++i)
        max_a = std::max(std::fabs(*(p)), max_a);

    return max_a;
}

/**
 * Computes the product of a vector by a scalar, which is defined as:
 *
 *   b *= a.
 *
 * Input:
 *  n - Specify the number of elements in vector b.
 *  a - Specify the scalar a.
 *  b - Array, size at least n.
 *
 * Output:
 *  b - Update vector b.
 */
template <typename T>
void Scale(int n, T a, T* b) {
    switch (n) {
    case 0: return;
    case 1: b[0] *= a;
            return;
    case 2: b[0] *= a;
            b[1] *= a;
            return;
    case 3: b[0] *= a;
            b[1] *= a;
            b[2] *= a;
            return;
    case 4: b[0] *= a;
            b[1] *= a;
            b[2] *= a;
            b[3] *= a;
            return;
    }

    T* p = b;

    int i = 0;
    for (; i + 4 < n; p += 4, i += 4) {
        *(p)     *= a;
        *(p + 1) *= a;
        *(p + 2) *= a;
        *(p + 3) *= a;
    }
    for (; i < n; ++p, ++i)
        *(p) *= a;
}

/**
 * In-place negate.
 *
 * Input:
 *  n - Specify the number of elements in vector a.
 *  a - Pointer to array that contains the input vector a.
 *
 * Output:
 *  a - Update vector a.
 */
template <typename T>
void Negate(int n, T* a) {
    switch (n) {
    case 0: return;
    case 1: a[0] = -a[0];
            return;
    case 2: a[0] = -a[0];
            a[1] = -a[1];
            return;
    case 3: a[0] = -a[0];
            a[1] = -a[1];
            a[2] = -a[2];
            return;
    case 4: a[0] = -a[0];
            a[1] = -a[1];
            a[2] = -a[2];
            a[3] = -a[3];
            return;
    }

    T* p = a;

    int i = 0;
    for (; i + 4 < n; p += 4, i += 4) {
        *(p)     *= -*(p);
        *(p + 1) *= -*(p + 1);
        *(p + 2) *= -*(p + 2);
        *(p + 3) *= -*(p + 3);
    }
    for (; i < n; ++p, ++i)
        *(p) = -*(p);
}

/**
 * Perform element by element addition of vector a and vector b.
 *
 * Input:
 *  n   - Specify the number of elements to be calculated.
 *  a,b - Pointers to arrays that contain the input vectors a and b.
 *
 * Output:
 *  c   - Pointer to an array that contains the output vector b.
 */
template <typename T>
void Add(int n, const T* a, const T* b, T* c) {
    switch (n) {
    case 0: return;
    case 1: c[0] = a[0] + b[0];
            return;
    case 2: c[0] = a[0] + b[0];
            c[1] = a[1] + b[1];
            return;
    case 3: c[0] = a[0] + b[0];
            c[1] = a[1] + b[1];
            c[2] = a[2] + b[2];
            return;
    case 4: c[0] = a[0] + b[0];
            c[1] = a[1] + b[1];
            c[2] = a[2] + b[2];
            c[3] = a[3] + b[3];
            return;
    }

    const T* p1 = a;
    const T* p2 = b;

    T* p = c;
    int i = 0;
    
    for (; i + 4 < n; p1 += 4, p2 += 4, p += 4, i += 4) {
        *(p)     = *(p1)     + *(p2);
        *(p + 1) = *(p1 + 1) + *(p2 + 1);
        *(p + 2) = *(p1 + 2) + *(p2 + 2);
        *(p + 3) = *(p1 + 3) + *(p2 + 3);
    }
    for (; i < n; ++p1, ++p2, ++p, ++i)
        *(p) = *(p1) + *(p2);
}

/**
 * Perform element by element subtraction of vector a and vector b.
 *
 * Input:
 *  n   - Specify the number of elements to be calculated.
 *  a,b - Pointers to arrays that contain the input vectors a and b.
 *
 * Output:
 *  c   - Pointer to an array that contains the output vector b.
 */
template <typename T>
void Subtract(int n, const T* a, const T* b, T* c) {
    switch (n) {
    case 0: return;
    case 1: c[0] = a[0] - b[0];
            return;
    case 2: c[0] = a[0] - b[0];
            c[1] = a[1] - b[1];
            return;
    case 3: c[0] = a[0] - b[0];
            c[1] = a[1] - b[1];
            c[2] = a[2] - b[2];
            return;
    case 4: c[0] = a[0] - b[0];
            c[1] = a[1] - b[1];
            c[2] = a[2] - b[2];
            c[3] = a[3] - b[3];
            return;
    }
    const T* p1 = a;
    const T* p2 = b;

    T* p = c;
    int i = 0;
    for (; i + 4 < n; p1 += 4, p2 += 4, p += 4, i += 4) {
        *(p)     = *(p1)     - *(p2);
        *(p + 1) = *(p1 + 1) - *(p2 + 1);
        *(p + 2) = *(p1 + 2) - *(p2 + 2);
        *(p + 3) = *(p1 + 3) - *(p2 + 3);
    }
    for (; i < n; ++p1, ++p2, ++p, ++i)
        *(p) = *(p1) - *(p2);
}

/**
 * Perform element by element multiplication of vector a and vector b.
 *
 * Input:
 *  n   - Specify the number of elements to be calculated.
 *  a,b - Pointers to arrays that contain the input vectors a and b.
 *
 * Output:
 *  c   - Pointer to an array that contains the output vector b.
 */
template <typename T>
void Multiply(int n, const T* a, const T* b, T* c) {
    switch (n) {
    case 0: return;
    case 1: c[0] = a[0] * b[0];
            return;
    case 2: c[0] = a[0] * b[0];
            c[1] = a[1] * b[1];
            return;
    case 3: c[0] = a[0] * b[0];
            c[1] = a[1] * b[1];
            c[2] = a[2] * b[2];
            return;
    case 4: c[0] = a[0] * b[0];
            c[1] = a[1] * b[1];
            c[2] = a[2] * b[2];
            c[3] = a[3] * b[3];
            return;
    }

    const T* p1 = a;
    const T* p2 = b;

    T* p = c;
    int i = 0;
    for (; i + 4 < n; p1 += 4, p2 += 4, p += 4, i += 4) {
        *(p)     = *(p1)     * *(p2);
        *(p + 1) = *(p1 + 1) * *(p2 + 1);
        *(p + 2) = *(p1 + 2) * *(p2 + 2);
        *(p + 3) = *(p1 + 3) * *(p2 + 3);
    }
    for (; i < n; ++p1, ++p2, ++p, ++i)
        *(p) = *(p1) * *(p2);
}

/**
 * Compute a matrix-vector product using a general matrix, which is defined as:
 *
 *   c = Ab
 *
 * Input:
 *  m - The number of rows of the matrix a.
 *  n - The number of columns of the matrix a.
 *  a - Array, size m * n.
 *  b - Array, size at least n.
 *  c - Array, size at least m.
 *
 * Output:
 *  c - Update vector c.
 */
template <typename T>
void GEMV(int m, int n, const T* a, const T* b, T* c) {
    // Special treat for 1 * 1, 2 * 2, 3 * 3, and 4 * 4 matrices.
    switch (m) {
    case 0: return;
    case 1:
        if (n == 1) {
            c[0] = a[0] * b[0];
            return;
        }
        break;
    case 2:
        if (n == 2) {
            c[0] = a[0] * b[0] + a[1] * b[1];
            c[1] = a[2] * b[0] + a[3] * b[1];
            return;
        }
        break;
    case 3:
        if (n == 3) {
            c[0] = a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
            c[1] = a[3] * b[0] + a[4] * b[1] + a[5] * b[2];
            c[2] = a[6] * b[0] + a[7] * b[1] + a[8] * b[2];
            return;
        }
        break;
    case 4:
        if (n == 4) {
            c[0] = a[0]  * b[0] + a[1]  * b[1] + a[2]  * b[2] + a[3]  * b[3];
            c[1] = a[4]  * b[0] + a[5]  * b[1] + a[6]  * b[2] + a[7]  * b[3];
            c[2] = a[8]  * b[0] + a[9]  * b[1] + a[10] * b[2] + a[11] * b[3];
            c[3] = a[12] * b[0] + a[13] * b[1] + a[14] * b[2] + a[15] * b[3];
            return;
        }
        break;
    }

    std::memset(c, 0, sizeof(T) * m);
    for (int i = 0; i < m; ++i) {
        // Using cache size 4 to accelerate.
        int j = 0;
        for (; j + 3 < n; j += 4) {
            c[i] += a[i * n + j]     * b[j] +
                    a[i * n + j + 1] * b[j + 1] +
                    a[i * n + j + 2] * b[j + 2] +
                    a[i * n + j + 3] * b[j + 3];
        }
        for (; j < n; ++j) {
            c[i] += a[i * n + j] * b[j];
        }
    }
}

/**
 * Compute a matrix-vector product using a general matrix, which is defined as:
 *
 *   c = A'b
 *
 * Input:
 *  m - The number of rows of the matrix a.
 *  n - The number of columns of the matrix a.
 *  a - Array, size m * n.
 *  b - Array, size at least m.
 *  c - Array, size at least n.
 *
 * Output:
 *  c - Update vector c.
 */
template <typename T>
void GEMVTrans(int m, int n, const T* a, const T* b, T* c) {
    // Special treat for 1 * 1, 2 * 2, 3 * 3, and 4 * 4 matrices.
    switch (m) {
    case 0: return;
    case 1:
        if (n == 1) {
            c[0] = a[0] * b[0];
            return;
        }
        break;
    case 2:
        if (n == 2) {
            c[0] = a[0] * b[0] + a[2] * b[1];
            c[1] = a[1] * b[0] + a[3] * b[1];
            return;
        }
        break;
    case 3:
        if (n == 3) {
            c[0] = a[0] * b[0] + a[3] * b[1] + a[6] * b[2];
            c[1] = a[1] * b[0] + a[4] * b[1] + a[7] * b[2];
            c[2] = a[2] * b[0] + a[5] * b[1] + a[8] * b[2];
            return;
        }
        break;
    case 4:
        if (n == 4) {
            c[0] = a[0] * b[0] + a[4] * b[1] + a[8]  * b[2] + a[12] * b[3];
            c[1] = a[1] * b[0] + a[5] * b[1] + a[9]  * b[2] + a[13] * b[3];
            c[2] = a[2] * b[0] + a[6] * b[1] + a[10] * b[2] + a[14] * b[3];
            c[3] = a[3] * b[0] + a[7] * b[1] + a[11] * b[2] + a[15] * b[3];
            return;
        }
        break;
    }

    std::memset(c, 0, sizeof(T) * n);
    for (int i = 0; i < m; ++i) {
        // Using cache size 4 to accelerate.
        int j = 0;
        for (; j + 3 < n; j += 4) {
            c[j]     += a[i * n + j]     * b[i];
            c[j + 1] += a[i * n + j + 1] * b[i];
            c[j + 2] += a[i * n + j + 2] * b[i];
            c[j + 3] += a[i * n + j + 3] * b[i];
        }
        for (; j < n; ++j) {
            c[j] += a[i * n + j] * b[i];
        }
    }
}

/**
 * Computes a matrix-matrix product with general matrices, which is defined as:
 *
 *   c = a * b
 *
 * Input:
 *  m - The number of rows of the matrix a.
 *  n - The number of columns of the matrix a.
 *  k - The number of columns of the matrix b.
 *  a - Array, size m * n.
 *  b - Array, size n * k.
 *  c - Array, size m * k.
 *
 * Output:
 *  c - Update matrix c.
 */
template <typename T>
void GEMM(int m, int n, int k, const T* a, const T* b, T* c) {
    // Special treat for 1 * 1, 2 * 2, 3 * 3, and 4 * 4 matrices.
    switch (m) {
    case 0: return;
    case 1:
        if (n == 1) {
            c[0] = a[0] * b[0];
            return;
        }
        break;
    case 2:
        if (n == 2) {
            c[0] = a[0] * b[0] + a[2] * b[1];
            c[1] = a[1] * b[0] + a[3] * b[1];
            c[2] = a[0] * b[2] + a[2] * b[3];
            c[3] = a[1] * b[2] + a[3] * b[3];
            return;
        }
        break;
    case 3:
        if (n == 3) {
            c[0] = a[0] * b[0] + a[1] * b[3] + a[2] * b[6];
            c[1] = a[0] * b[1] + a[1] * b[4] + a[2] * b[7];
            c[2] = a[0] * b[2] + a[1] * b[5] + a[2] * b[8];

            c[3] = a[3] * b[0] + a[4] * b[3] + a[5] * b[6];
            c[4] = a[3] * b[1] + a[4] * b[4] + a[5] * b[7];
            c[5] = a[3] * b[2] + a[4] * b[5] + a[5] * b[8];

            c[6] = a[6] * b[0] + a[7] * b[3] + a[8] * b[6];
            c[7] = a[6] * b[1] + a[7] * b[4] + a[8] * b[7];
            c[8] = a[6] * b[2] + a[7] * b[5] + a[8] * b[8];
            return;
        }
        break;
    case 4:
        if (n == 4) {
            c[0] = a[0] * b[0] + a[1] * b[4] + a[2] * b[8]  + a[3] * b[12];
            c[1] = a[0] * b[1] + a[1] * b[5] + a[2] * b[9]  + a[3] * b[13];
            c[2] = a[0] * b[2] + a[1] * b[6] + a[2] * b[10] + a[3] * b[14];
            c[3] = a[0] * b[3] + a[1] * b[7] + a[2] * b[11] + a[3] * b[15];

            c[4] = a[4] * b[0] + a[5] * b[4] + a[6] * b[8]  + a[7] * b[12];
            c[5] = a[4] * b[1] + a[5] * b[5] + a[6] * b[9]  + a[7] * b[13];
            c[6] = a[4] * b[2] + a[5] * b[6] + a[6] * b[10] + a[7] * b[14];
            c[7] = a[4] * b[3] + a[5] * b[7] + a[6] * b[11] + a[7] * b[15];

            c[8]  = a[8] * b[0] + a[9] * b[4] + a[10] * b[8]  + a[11] * b[12];
            c[9]  = a[8] * b[1] + a[9] * b[5] + a[10] * b[9]  + a[11] * b[13];
            c[10] = a[8] * b[2] + a[9] * b[6] + a[10] * b[10] + a[11] * b[14];
            c[11] = a[8] * b[3] + a[9] * b[7] + a[10] * b[11] + a[11] * b[15];

            c[12] = a[12] * b[0] + a[13] * b[4] + a[14] * b[8]  + a[15] * b[12];
            c[13] = a[12] * b[1] + a[13] * b[5] + a[14] * b[9]  + a[15] * b[13];
            c[14] = a[12] * b[2] + a[13] * b[6] + a[14] * b[10] + a[15] * b[14];
            c[15] = a[12] * b[3] + a[13] * b[7] + a[14] * b[11] + a[15] * b[15];
            return;
        }
        break;
    }

    std::memset(c, 0, sizeof(T) * m * k);

    for (int ii = 0; ii < m; ++ii) {
        for (int kk = 0; kk < n; ++kk) {
            T* c0 = &c[ii * k];
            const T& a0 =  a[ii * n + kk];
            const T* b0 = &b[kk * k];

            int j = 0;
            for (; j + 3 < k; j += 4) {
                *c0++ += a0 * *b0++;
                *c0++ += a0 * *b0++;
                *c0++ += a0 * *b0++;
                *c0++ += a0 * *b0++;
            }
            for (; j < k; ++j) {
                *c0++ += a0 * *b0++;
            }
        }
    }
}

} // namespace blas

/**
 * Basic linear algebra object.
 */
template <typename T>
class BasicLinearAlgebra {
    static const int ALIGNMENT = 64; // Memory alignment.

public:
    using value_type             = T;
    using pointer                = T*;
    using const_pointer          = const T*;
    using reference              = T&;
    using const_reference        = const T&;
    using iterator               = T*;
    using const_iterator         = const T*;
    using reverse_iterator       = std::reverse_iterator<iterator>;
    using const_reverse_iterator = std::reverse_iterator<const_iterator>;
    using size_type              = int;
    using difference_type        = int;

    /**
     * Default BasicLinearAlgebra constructor.
     */
    BasicLinearAlgebra() = default;

    explicit BasicLinearAlgebra(int size) {
        Allocate(size);
    }

    /**
     * Construct a BasicLinearAlgebra with initialize value.
     */
    BasicLinearAlgebra(int size, const T& v) {
        Allocate(size);
        std::uninitialized_fill_n(data_, size, v);
    }

    /**
     * Construct BasicLinearAlgebra from data in [first, last).
     * The second template parameter is used to distinguish this function to
     * BasicLinearAlgebra(int, int).
     */
    template <typename Iter,
              typename = typename std::enable_if<std::is_convertible<
                         typename std::iterator_traits<Iter>::iterator_category,
                                  std::input_iterator_tag>::value>::type>
    BasicLinearAlgebra(Iter first, Iter last) {
        auto n = std::distance(first, last);
        assert(n >= 0);
        assert(n <= INT_MAX && "We only accept INT_MAX elements at most.");

        Allocate(static_cast<int>(n));
        std::uninitialized_copy(first, last, data_);
    }

    /**
     * Copy constructor.
     */
    explicit BasicLinearAlgebra(const BasicLinearAlgebra<T>& v) {
        Allocate(v.size_);
        std::uninitialized_copy(v.data_, v.data_ + v.size_, data_);
    }

    /**
     * Move constructor.
     */
    explicit BasicLinearAlgebra(BasicLinearAlgebra<T>&& v) {
        swap(&v);
    }

    /**
     * Build BasicLinearAlgebra from initializer list.
     *
     * Note that this constructor can not be explicit.
     */
    explicit BasicLinearAlgebra(std::initializer_list<T> list)
        : BasicLinearAlgebra(list.begin(), list.end()) {}

    virtual ~BasicLinearAlgebra() {
        Deallocate(data_);
    }

    /**
     * Clear the BasicLinearAlgebra.
     */
    void clear() {
        size_ = 0;
    }

    /**
     * Check if this BasicLinearAlgebra is empty.
     */
    bool empty() const {
        return size_ == 0;
    }

    /**
     * Swap the data with another BasicLinearAlgebra.
     */
    void swap(BasicLinearAlgebra<T>* rhs) {
        std::swap(data_, rhs->data_);
        std::swap(size_, rhs->size_);
        std::swap(allocated_size_, rhs->allocated_size_);
    }

    /**
     * Return the number of elements.
     */
    int size() const {
        return size_;
    }

    T* data() {
        return data_;
    }

    const T* data() const {
        return data_;
    }

    iterator begin() {
        return data_;
    }

    iterator end() {
        return data_ + size_;
    }

    const_iterator begin() const {
        return data_;
    }

    const_iterator end() const {
        return data_ + size_;
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

    T& operator[](int index) {
        return *(data_ + index);
    }

    const T& operator[](int index) const {
        return *(data_ + index);
    }

    T& at(int index) {
        assert(index >= 0 && index < size_);

        return *(data_ + index);
    }

    const T& at(int index) const {
        assert(index >= 0 && index < size_);

        return *(data_ + index);
    }

    /**
     * Move assignment.
     */
    BasicLinearAlgebra& operator =(BasicLinearAlgebra<T>&& rhs) noexcept {
        swap(&rhs);
        return *this;
    }

    /**
     * Copy assignment.
     */
    BasicLinearAlgebra& operator =(const BasicLinearAlgebra<T>& rhs) {
        if (this == &rhs) return *this; // This case happens often.

        Reallocate(rhs.size_);
        std::memcpy(data_, rhs.data(), sizeof(T) * size_);
        return *this;
    }

    bool operator ==(const BasicLinearAlgebra& rhs) const {
        return size_ == rhs.size_ &&
               std::equal(data_, data_ + size_, rhs.data_);
    }

    bool operator !=(const BasicLinearAlgebra& rhs) const {
        return !(*this == rhs);
    }

    /**
     * Fill the vector.
     */
    void fill(const T& v) {
        if (data_) std::fill_n(data_, size_, v);
    }

protected:
    /**
     * Allocate the memory for BasicLinearAlgebra.
     */
    void Allocate(int size) {
        assert(size >= 0);

        if (size > 0) {
            //  Allocates an aligned memory buffer.
            data_ = reinterpret_cast<T*>(std::malloc(size * sizeof(T)));
            assert(data_ && "Malloc error, maybe memory is not enough.");
        }
        allocated_size_ = size_ = size;
    }

    /**
     * Free the aligned memory buffer allocated by malloc.
     */
    void Deallocate(T* ptr) {
        if (ptr) std::free(ptr);
    }

    /**
     * Reallocate data.
     */
    void Reallocate(int size) {
        assert(size >= 0);

        if (size > allocated_size_) {
            Deallocate(data_);
            Allocate(size);
        } else {
            size_ = size;
        }
    }

    /**
     * Resize the BLAS without fill the value.
     */
    void Resize(int n) {
        assert(n >= 0);

        if (allocated_size_ == 0) {
            Allocate(n);
        } else {
            if (n > allocated_size_) {
                T* tmp = data_;
                int old_size = size_;

                Allocate(n);
                std::uninitialized_copy(tmp, tmp + old_size, data_);
                Deallocate(tmp);
            } else {
                size_ = n;
            }
        }
    }

    /**
     * Resize the BLAS with the filling value.
     */
    void Resize(int n, const T& v) {
        assert(n >= 0);

        if (allocated_size_ == 0) {
            Allocate(n);
            std::fill(data_, data_ + n, v);
        } else {
            if (n > allocated_size_) {
                T* tmp = data_;
                int old_size = size_;

                Allocate(n);
                std::uninitialized_copy(tmp, tmp + old_size, data_);
                std::uninitialized_fill(data_ + old_size, data_ + n, v);
                Deallocate(tmp);
            } else if (n > size_) {
                std::fill(data_ + size_, data_ + n, v);
                size_ = n;
            } else {
                size_ = n;
            }
        }
    }

    /**
     * Assign 'n' elements of value 'v' to this BasicLinearAlgebra.
     */
    void Assign(int n, const T& v) {
        assert(n >= 0);

        if (n == 0) {
            clear();
        } else if (n > allocated_size_) {
            // Need relocate the data.
            Reallocate(n);
            std::uninitialized_fill_n(data_, n, v);
        } else {
            size_ = n;
            std::fill_n(data_, n, v);
        }
    }

    /**
     * Construct BasicLinearAlgebra from data in [first, last).
     * The second template parameter is used to distinguish this function to
     * BasicLinearAlgebra(int, int).
     */
    template <typename Iter,
              typename = typename std::enable_if<std::is_convertible<
                         typename std::iterator_traits<Iter>::iterator_category,
                                  std::input_iterator_tag>::value>::type>
    void Assign(Iter first, Iter last) {
        auto n = std::distance(first, last);
        assert(n >= 0);
        assert(n <= INT_MAX && "We only accept INT_MAX elements at most.");

        Reallocate(static_cast<int>(n));
        std::copy(first, last, data_);
    }

    // Number of elements in the BasicLinearAlgebra.
    int size_ = 0;
    
    // Allocated size of BasicLinearAlgebra.
    int allocated_size_ = 0;

    // Store data.
    T* data_  = nullptr;
};

} // namespace cl

#endif // MATH_BASIC_LINEAR_ALGEBRA_H_
