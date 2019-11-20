//
// Copyright 2019 Yangbin Lin. All Rights Reserved.
//
// Author: yblin@jmu.edu.cn (Yangbin Lin)
//
// This file is part of the Code Library.
//

#ifndef MATH_MATRIX_EIGEN_SYMMETRIC_EIGEN_H_
#define MATH_MATRIX_EIGEN_SYMMETRIC_EIGEN_H_

#include "codelibrary/base/algorithm.h"
#include "codelibrary/math/matrix/matrix.h"

namespace cl {
namespace matrix {

/**
 * Compute all eigenvalues and, optionally, eigenvectors of a 2 * 2 real 
 * symmetric matrix.
 *
 * Reference:
 *   Kopp J. Efficient numerical diagonalization of hermitian 3x3 matrices[J].
 *   International Journal of Modern Physics C, 2008, 19(03): 523-548.
 *
 * Parameters:
 *  mat          - the input matrix using only the lower triangular part.
 *  eigenvalues  - the output eigenvalues will in ascending order.
 *  eigenvectors - the output eigenvectors are normalized so that ||z_i||_2 = 1.
 * 
 * Return:
 *   always return true.
 */
template <typename T>
bool SymmetricEigen2(const Matrix<T>& mat, Vector<T>* eigenvalues,
                     Matrix<T>* eigenvectors = nullptr) {
    static_assert(std::is_floating_point<T>::value,
                  "T must be floating.");

    assert(mat.n_rows() == 2 && mat.n_columns() == 2);
    assert(eigenvalues);

    // Compute the eigenvalues.
    eigenvalues->resize(2);

    Vector<T>& w = *eigenvalues;

    // [ A  B ]  =  [ cs  -sn ] [ rt1   0  ] [  cs  sn ]
    // [ B  C ]     [ sn   cs ] [  0   rt2 ] [ -sn  cs ]
    T a = mat(0, 0), b = mat(1, 0), c = mat(1, 1);
    T sm = a + c;
    T df = a - c;
    T rt = std::sqrt(df * df + b * b * 4);

    if (sm > 0) {
        w[1] = (sm + rt) / 2;
        T t = 1 / w[1];
        w[0] = (a * t) * c - (b * t) * b;
    } else if (sm < 0) {
        w[0] = (sm - rt) / 2;
        T t  = 1 / w[1];
        w[1] = (a * t) * c - (b * t) * b;
    } else {
        // This case needs to be treated separately to avoid div by 0.
        w[1] =  rt / 2;
        w[0] = -rt / 2;
    }

    if (eigenvectors) {
        // Compute eigenvectors.
        T cs = df > 0 ? df + rt : df - rt;
        T t = 0, sn = 0;
        if (std::fabs(cs) > std::fabs(b) * 2) {
            t = -b * 2 / cs;
            sn = 1 / std::sqrt(t * t + 1);
            cs = t * sn;
        } else if (std::fabs(b) == 0) {
            cs = 1;
            sn = 0;
        } else {
            t = -cs / b / 2;
            cs = 1 / std::sqrt(t * t + 1);
            sn = t * cs;
        }

        if (df > 0) {
            t = cs;
            cs = -sn;
            sn = t;
        }

        eigenvectors->resize(2, 2);
        (*eigenvectors)(0, 0) = -sn;
        (*eigenvectors)(1, 0) = cs;

        (*eigenvectors)(0, 1) = cs;
        (*eigenvectors)(1, 1) = sn;
    }
    return true;
}

/**
 * Computes all eigenvalues and, optionally, eigenvectors of a 3 * 3 real 
 * symmetric matrix.
 *
 * Reference:
 *   Kopp J. Efficient numerical diagonalization of hermitian 3x3 matrices[J].
 *   International Journal of Modern Physics C, 2008, 19(03): 523-548.
 *
 * Parameters:
 *  mat          - the input matrix using only the lower triangular part.
 *  eigenvalues  - the output eigenvalues will in ascending order.
 *  eigenvectors - the output eigenvectors are normalized so that ||z_i||_2 = 1.
 * 
 * Return:
 *   always return true.
 */
template <typename T>
bool SymmetricEigen3(const Matrix<T>& a, Vector<T>* eigenvalues,
                     Matrix<T>* eigenvectors = nullptr) {
    static_assert(std::is_floating_point<T>::value,
                  "T must be a floating point.");

    assert(a.n_rows() == 3 && a.n_columns() == 3);
    assert(eigenvalues);

    // Compute the eigenvalues.
    eigenvalues->resize(3);
    Vector<T>& w = *eigenvalues;

    T a00 = a(0, 0), a01 = a(1, 0), a02 = a(2, 0);
    T a11 = a(1, 1), a12 = a(2, 1), a22 = a(2, 2);

    T de = a01 * a12;
    T dd = a01 * a01;
    T ee = a12 * a12;
    T ff = a02 * a02;

    T m  = a00 + a11 + a22;
    T c1 = a00 * a11 + a00 * a22 + a11 * a22 - (dd + ee + ff);
    T c0 = dd * a22 + ee * a00 + ff * a11 - a00 * a11 * a22 - a02 * de * 2;

    T p = m * m - c1 * 3;
    T q = m * (p - c1 * T(1.5)) - c0 * T(13.5);
    T sqrt_p = std::sqrt(std::fabs(p));

    T phi = (c1 * c1 * (p - c1) * T(0.25) + c0 * (q + c0 * T(6.75))) * T(27);
    phi = std::atan2(std::sqrt(std::fabs(phi)), q) / 3;

    T c = sqrt_p * std::cos(phi);
    T s = sqrt_p * std::sin(phi) / std::sqrt(3);

    w[1]  = (m - c) / 3;
    w[2]  = w[1] + s;
    w[0]  = w[1] + c;
    w[1] -= s;
    if (eigenvectors) {
        // Compute eigenvectors.
        eigenvectors->resize(3, 3);
        Matrix<T>& v = *eigenvectors;

        T max_eigenvalue = std::max(std::fabs(w[0]), std::fabs(w[1]));
        max_eigenvalue = std::max(max_eigenvalue, std::fabs(w[2]));
        T epsilon = std::numeric_limits<T>::epsilon();
        T thresh = epsilon * max_eigenvalue * 8;
        thresh *= thresh;

        // Prepare calculation of eigenvectors.
        T n0tmp = a01 * a01 + a02 * a02;
        T n1tmp = a01 * a01 + a12 * a12;
        v(0, 1) = a01 * a12 - a02 * a11;
        v(1, 1) = a02 * a01 - a12 * a00;
        v(2, 1) = a01 * a01;

        // Calculate first eigenvector by the formula:
        //   v[0] = (A - w[0]).e1 x (A - w[0]).e2.
        a00 -= w[0];
        a11 -= w[0];
        v(0, 0) = v(0, 1) + a02 * w[0];
        v(1, 0) = v(1, 1) + a12 * w[0];
        v(2, 0) = a00 * a11 - v(2, 1);

        T norm = v(0, 0) * v(0, 0) + v(1, 0) * v(1, 0) + v(2, 0) * v(2, 0);
        T n0    = n0tmp + a00 * a00;
        T n1    = n1tmp + a11 * a11;
        T error = n0 * n1;

        if (n0 <= thresh) {
            v(0, 0) = 1;
            v(1, 0) = 0;
            v(2, 0) = 0;
        } else if (n1 <= thresh) {
            v(0, 0) = 0;
            v(1, 0) = 1;
            v(2, 0) = 0;
        } else if (norm < 4096 * epsilon * epsilon * error) {
            // If angle between mat_t[0] and mat_t[1] is too small, don't use
            // cross product, but calculate v ~ (1, -A0/A1, 0).
            T t = a01 * a01;
            T f = -a00 / a01;
            if (a11 * a11 > t) {
                t =  a11 * a11;
                f = -a01 / a11;
            }
            if (a12 * a12 > t)
                f = -a02 / a12;

            T norm = 1 / std::sqrt(f * f + 1);
            v(0, 0) = norm;
            v(1, 0) = f * norm;
            v(2, 0) = 0;
        } else {
            // This is the standard branch.
            T t = std::sqrt(1 / norm);
            v(0, 0) *= t;
            v(1, 0) *= t;
            v(2, 0) *= t;
        }

        // Prepare calculation of second eigenvector.
        T t = w[0] - w[1];
        if (std::fabs(t) > 8 * epsilon * max_eigenvalue) {
            // For non-degenerate eigenvalue, calculate second eigenvector by
            // the formula
            //   v[1] = (A - w[1]).e1 x (A - w[1]).e2.
            a00 += t;
            a11 += t;
            v(0, 1) += a02 * w[1];
            v(1, 1) += a12 * w[1];
            v(2, 1) = a00 * a11 - v(2, 1);

            T norm  = v(0, 1) * v(0, 1) + v(1, 1) * v(1, 1) + v(2, 1) * v(2, 1);
            T n0    = n0tmp + a00 * a00;
            T n1    = n1tmp + a11 * a11;
            T error = n0 * n1;

            if (n0 <= thresh) {
                v(0, 1) = 1;
                v(1, 1) = 0;
                v(2, 1) = 0;
            } else if (n1 <= thresh) {
                v(0, 1) = 0;
                v(1, 1) = 1;
                v(2, 1) = 0;
            } else if (norm < 4096 * epsilon * epsilon * error) {
                // If angle between mat_t[0] and mat_t[1] is too small,
                // don't use cross product, but calculate v ~ (1, -A0/A1, 0).
                T t = a01 * a01;
                T f = -a00 / a01;
                if (a11 * a11 > t) {
                    t =  a11 * a11;
                    f = -a01 / a11;
                }
                if (a12 * a12 > t) f = -a02 / a12;

                T norm = 1 / std::sqrt(f * f + 1);
                v(0, 1) = norm;
                v(1, 1) = f * norm;
                v(2, 1) = 0;
            } else {
                // This is the standard branch.
                T t = 1 / std::sqrt(1 / norm);
                v(0, 1) *= t;
                v(1, 1) *= t;
                v(2, 1) *= t;
            }
        } else {
            // For degenerate eigenvalue, calculate second eigenvector according
            // to:
            //   v[1] = v[0] x (A - w[1]).e[i].

            // Reset the mat_t to mat.
            a00 += w[0];
            a11 += w[0];
            Matrix<T> mat_t(3, 3);
            mat_t(0, 0) = a00; mat_t(0, 1) = a01; mat_t(0, 2) = a02;
            mat_t(1, 0) = a01; mat_t(1, 1) = a11; mat_t(1, 2) = a12;
            mat_t(2, 0) = a02; mat_t(2, 1) = a12; mat_t(2, 2) = a22;

            int i = 0;
            for (i = 0; i < 3; ++i) {
                mat_t(i, i) -= w[1];

                T n0 = mat_t(0, i) * mat_t(0, i) +
                       mat_t(1, i) * mat_t(1, i) +
                       mat_t(2, i) * mat_t(2, i);
                if (n0 > thresh) {
                    v(0, 1) = v(1, 0) * mat_t(2, i) - v(2, 0) * mat_t(1, i);
                    v(1, 1) = v(2, 0) * mat_t(0, i) - v(0, 0) * mat_t(2, i);
                    v(2, 1) = v(0, 0) * mat_t(1, i) - v(1, 0) * mat_t(0, i);

                    T norm = v(0, 1) * v(0, 1) + v(1, 1) * v(1, 1) +
                             v(2, 1) * v(2, 1);
                    if (norm > epsilon * epsilon * 65536 * n0) {
                        // Accept cross product only if the angle between the
                        // two vectors was not too small.
                        T t = std::sqrt(1 / norm);
                        v(0, 1) *= t;
                        v(1, 1) *= t;
                        v(2, 1) *= t;
                        break;
                    }
                }
            }

            if (i == 3) {
                // This means that any vector orthogonal to v[0] is an EV.
                for (int j = 0; j < 3; ++j) {
                    if (v(j, 0) != 0) {
                        // Find nonzero element of v[0] and swap it with the
                        // next one.
                        int k = (j + 1) % 3;
                        T t = 1 / std::sqrt(v(j, 0) * v(j, 0) + 
                                            v(k, 0) * v(k, 0));
                        v(j, 1) =  v(k, 0) * t;
                        v(k, 1) = -v(j, 0) * t;
                        v((j + 2) % 3, 1) = 0;
                        break;
                    }
                }
            }
        }

        // Calculate third eigenvector according to v[2] = v[0] x v[1].
        v(0, 2) = v(1, 0) * v(2, 1) - v(2, 0) * v(1, 1);
        v(1, 2) = v(2, 0) * v(0, 1) - v(0, 0) * v(2, 1);
        v(2, 2) = v(0, 0) * v(1, 1) - v(1, 0) * v(0, 1);

        // Sort in ascending order.
        if (w[1] < w[0]) {
            std::swap(w[0], w[1]);
            for (int i = 0; i < 3; ++i) {
                std::swap(v(i, 0), v(i, 1));
            }
        }
        if (w[2] < w[1]) {
            std::swap(w[1], w[2]);
            for (int i = 0; i < 3; ++i) {
                std::swap(v(i, 1), v(i, 2));
            }
            if (w[1] < w[0]) {
                std::swap(w[0], w[1]);
                for (int i = 0; i < 3; ++i) {
                    std::swap(v(i, 0), v(i, 1));
                }
            }
        }
    } else {
        // Sort in ascending order.
        if (w[1] < w[0]) std::swap(w[0], w[1]);
        if (w[2] < w[1]) {
            std::swap(w[1], w[2]);
            if (w[1] < w[0]) std::swap(w[0], w[1]);
        }
    }
    return true;
}

/**
 * Computes all eigenvalues and, optionally, eigenvectors of a real symmetric 
 * matrix.
 *
 * This function uses MKL to solve the eigen problem of a matrix of n > 3. For 
 * n <= 3 we uses a fast implementation proposed by [Kopp J 2008].
 *
 * Parameters:
 *  a            - the input matrix, using only the lower triangular part.
 *  eigenvalues  - the output eigenvalues will in ascending order.
 *  eigenvectors - the output eigenvectors are normalized so that ||z_i||_2 = 1.
 *                 the i-th column is the eigenvector of the i-th eigenvalue.
 *
 * Return:
 *   false if the algorithm fails to converge.
 */
template <typename T>
bool SymmetricEigen(const Matrix<T>& a, Vector<T>* eigenvalues,
                    Matrix<T>* eigenvectors = nullptr) {
    static_assert(std::is_floating_point<T>::value,
                  "T must be a floating point.");

    assert(a.n_rows() == a.n_columns());
    assert(eigenvalues);

    int n = a.n_rows();

    switch (n) {
    case 0:
        eigenvalues->clear();
        if (eigenvectors) eigenvectors->clear();
        return true;
    case 1:
        eigenvalues->assign(1, a(0, 0));
        if (eigenvectors) eigenvectors->assign(1, 1, 1);
        return true;
    case 2:
        return SymmetricEigen2(a, eigenvalues, eigenvectors);
    case 3:
        return SymmetricEigen3(a, eigenvalues, eigenvectors);
    }

    assert(false && "We currently only support n <= 3.");

    return true;
}

} // namespace matrix
} // namespace cl

#endif // MATH_MATRIX_EIGEN_SYMMETRIC_EIGEN_H_
