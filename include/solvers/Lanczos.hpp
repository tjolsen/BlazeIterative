// Copyright (c)   2017 Tyler Olsen
//                 2018 Patrick Diehl
//                 2019 Nanmiao Wu
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BLAZE_ITERATIVE_LANCZOS_HPP
#define BLAZE_ITERATIVE_LANCZOS_HPP

#include "IterativeCommon.hpp"
#include "LanczosTag.hpp"


BLAZE_NAMESPACE_OPEN
    ITERATIVE_NAMESPACE_OPEN

        namespace detail {

            template<typename MatrixType, typename T>
            void  solve_impl(
                    MatrixType &h,
                    MatrixType &Q,
                    const MatrixType &A,
                    const DynamicVector<T> &b,
                    LanczosTag &tag,
                    const std::size_t &n
                  )
            {

                BLAZE_INTERNAL_ASSERT(A.isSymmetric(), "A must be a symmetric matrix")

                // n is dimension of Krylov subspace, n >=1;
                //
                BLAZE_INTERNAL_ASSERT(n >= 1, "n must larger than or equal to 1")

                // computes a basis of the n- Krylov subspace of A:
                // the space is spanned by{b, Ab, ..., A^(n-1)b}
                // Input
                // A: m * m matrix
                // b: initial vector (length m)
                // n: dimension of Krylov subspace, n >=1
                // Returns Q, h
                // Q: m * n matrix, the columns are an orthonormal basis of the Krylov subspace
                // h: n * n matrix

                std::size_t m = A.columns();
                DynamicVector<T> alpha(n);
                DynamicVector<T> beta(n);
                DynamicVector<T> w(m);
                DynamicVector<T> Av(m);

                // Let q1 be an arbitrary vector with Euclidean norm 1
                column(Q, 0) = b / norm(b);

                // Abbreviated initial iteration step:
                beta[0] = 0;
                Av = declsym(A) * column(Q, 0);
                alpha[0] = ctrans(Av) * column(Q, 0);
                w = Av - alpha[0] * column(Q, 0);


                // for j = 1, ... n-1 do:
                for( int j = 1; j < n; ++j){
                    beta[j] = norm(w);
                    if (beta[j] != 0){
                        column(Q, j) = w / beta[j];
                    } else{
                        break;
                    }
                    Av = declsym(A) * column(Q, j);
                    alpha[j] = ctrans(Av) * column(Q, j);
                    w = Av - alpha[j] * column(Q, j) - beta[j] * column(Q, j-1);
                }

                auto sub_beta = subvector(beta, 1UL, (n-1) );
                band<0>(h) = alpha;
                band<-1>(h) = sub_beta;
                band<1>(h) = sub_beta;

            }; // end solve_imple function

        } //end namespace detail

    ITERATIVE_NAMESPACE_CLOSE
BLAZE_NAMESPACE_CLOSE

#endif //BLAZE_ITERATIVE_LANCZOS_HPP
