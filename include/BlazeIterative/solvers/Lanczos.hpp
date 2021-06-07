// Copyright (c)   2017 Tyler Olsen
//                 2018 Patrick Diehl
//                 2019 Nanmiao Wu
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BLAZE_ITERATIVE_LANCZOS_HPP
#define BLAZE_ITERATIVE_LANCZOS_HPP

#include <BlazeIterative/IterativeCommon.hpp>
#include "LanczosTag.hpp"


BLAZE_NAMESPACE_OPEN
    ITERATIVE_NAMESPACE_OPEN

        namespace detail {

            template<typename MatrixType, typename T>
            void  solve_impl(
                    DynamicVector<T> &x,
                    const  MatrixType &A,
                    const DynamicVector<T> &b,
                    LanczosTag &tag,
                    const std::size_t &n
                  )
            {

                BLAZE_INTERNAL_ASSERT(isSymmetric(A), "A must be a symmetric matrix")

                // n is dimension of Krylov subspace, n >=1;
                
                BLAZE_INTERNAL_ASSERT(n >= 1, "n must larger than or equal to 1")

                // computes a basis of the n- Krylov subspace of A:
                // the space is spanned by{b, Ab, ..., A^(n-1)b}
                // Input
                // A: m * m matrix
                // b: initial vector (length m)
                // n: number of iterations
                // Returns Q, h
                // Q: m * n matrix, the columns are an orthonormal
                // h: n * n matrix, tridiagonal real symmetric

                std::size_t m = A.columns();
                DynamicVector<T> alpha(n);
                DynamicVector<T> beta(n);
                DynamicVector<T> Av(m);
                DynamicMatrix<T> Q(m, n);
                DynamicMatrix<T> h(n, n, 0);
                DynamicVector<complex<double>> x_comp(n);


             // with an example of a diagonal matrix A
             // http://people.inf.ethz.ch/arbenz/ewp/Lnotes/chapter10.pdf

                column(Q,0) = b / norm(b);
                Av = A * column(Q,0);
                alpha[0] = trans(column(Q,0)) * Av;
                Av -= alpha[0] * column(Q,0);
                beta[0] = norm(Av);

                for(int j =1 ; j < n; ++j){
                    column(Q,j) = Av / beta[j-1];
                    Av = A * column(Q,j) - beta[j-1] * column(Q,j-1);
                    alpha[j] = trans(column(Q,j)) * Av;
                    Av -= alpha[j] * column(Q,j);
                    auto Q_1 = submatrix(Q, 0, 0, m, j+1);
                    Av -= Q_1 * (trans(Q_1) * Av);
                    beta[j] = norm(Av);
                    if(beta[j] == 0){
                        break;
                    }
                }

                auto beta_1 = subvector(beta, 0, n-1);
                band<0>(h) = alpha;
                band<1>(h) = beta_1;
                band<-1>(h) = beta_1;

                eigen(h, x_comp);
                x = real(x_comp);

            }; // end solve_imple function

        } //end namespace detail

    ITERATIVE_NAMESPACE_CLOSE
BLAZE_NAMESPACE_CLOSE

#endif //BLAZE_ITERATIVE_LANCZOS_HPP
