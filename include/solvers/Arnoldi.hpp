// Copyright (c)   2017 Tyler Olsen
//                 2018 Patrick Diehl
//                 2019 Nanmiao Wu
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BLAZE_ITERATIVE_ARNOLDI_HPP
#define BLAZE_ITERATIVE_ARNOLDI_HPP

#include "IterativeCommon.hpp"
#include "ArnoldiTag.hpp"


BLAZE_NAMESPACE_OPEN
    ITERATIVE_NAMESPACE_OPEN

        namespace detail {

            template<typename MatrixType, typename T>
            void  solve_impl(
                    MatrixType &h,
                    MatrixType &Q,
                    const MatrixType &A,
                    const DynamicVector<T> &b,
                    ArnoldiTag &tag,
                    const std::size_t &n
                   )
            {

                BLAZE_INTERNAL_ASSERT(isSymmetric(A), "A must be a symmetric matrix")

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
                // Q: m * (n+1) matrix, the columns are an orthonormal basis of the Krylov subspace
                // h: (n+1) * n matrix, A on basis Q. It is upper Hessenberg

                // A is m*m matrix
                std::size_t m = b.size();

                // Return Q and h
               // DynamicMatrix<T> Q(m, (n + 1));
               // DynamicMatrix<T> h((n + 1), n);

                // let b be an arbitrary initial vector
                column(Q, 0) = b / norm(b);

                // the next vector q_k = A* q_k_1
                DynamicVector<T> v(m);

                for (int k = 0; k < n; ++k ){
                    v = declsym(A) * column(Q, k);

                    for (int j = 0; j < k+1; ++j){
                        h(j,k) = ctrans(column(Q,j)) * v;
                        v -= h(j,k) * column(Q,j);
                    }

                    h(k+1, k) = norm(v);
                    // if h(k+1, k) = 0
                    double eps = 1e-12;
                    if (h(k+1, k) > eps){
                        column(Q,k+1) = v / h(k+1,k);
                    }
                    else{
                        break;
                    }

                }

            }; // end solve_imple function

        } //end namespace detail

ITERATIVE_NAMESPACE_CLOSE
BLAZE_NAMESPACE_CLOSE
#endif //BLAZE_ITERATIVE_ARNOLDI_HPP
