// Copyright (c)   2017 Tyler Olsen
//                 2018 Patrick Diehl
//                 2019 Nanmiao Wu
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BLAZE_ITERATIVE_LANCZOS_HPP
#define BLAZE_ITERATIVE_LANCZOS_HPP

#include <iostream>
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
                // n: number of iterations
                // Returns Q, h
                // Q: m * n matrix, the columns are an orthonormal
                // h: n * n matrix, tridiagonal real symmetric


                std::size_t m = A.columns();
                DynamicMatrix<T> V(m, (n+1));
                DynamicVector<T> alpha(n+1);
                DynamicVector<T> beta(n+1);
                DynamicVector<T> w(m);
                DynamicVector<T> Av(m);
                DynamicVector<T> q(m,0);
            //    DynamicVector<T> r(b);

                beta[0] = 0;
                alpha[0] = 0;
                column(V,0) = q;
                column(V,1) = b / norm(b);

                for(int j = 1; j <=n; ++j ){
                    Av = A * column(V,j);
                    alpha[j] = trans(column(V,j)) * Av;
                    Av -= alpha[j] * column(V,j) - beta[j-1]*column(V,j-1);
                    beta[j] = norm(Av);
                    if(beta[j] == 0){
                        break;
                    } else{
                        column(V,j+1) = Av / beta[j];
                    }
                }

                Q = submatrix(V, 0UL, 1UL, m, n);
//                int j = 0;
//                while(beta[j] !=0){
//                    column(Q,j+1) = r / beta[j];
//                    alpha[j] = trans(column(Q,j+1)) * A * column(Q,j+1);
//                    r = A * column(Q,j+1)  - alpha[j]* column(Q,j+1) - beta[j] * column(Q,j);
//                    beta[j+1] = norm(r);
//                    j += 1;
//                }



//                // Let q1 be an arbitrary vector with Euclidean norm 1
//                column(Q, 0) = b / norm(b);
//
//                // Abbreviated initial iteration step:
//                beta[0] = 0;
//                Av = A * column(Q, 0);
//                alpha[0] = ctrans(Av) * column(Q, 0);
//                w = Av - alpha[0] * column(Q, 0);
//
//
//                // for j = 1, ... n-1 do:
//                for( int j = 1; j < n; ++j){
//                    beta[j] = norm(w);
//                    if (beta[j] != 0){
//                        column(Q, j) = w / beta[j];
//                    } else{
//                        break;
//                    }
//                    Av = A * column(Q, j);
//                    alpha[j] = ctrans(Av) * column(Q, j);
//                    w = Av - alpha[j] * column(Q, j) - beta[j] * column(Q, j-1);
//                }
//
//                auto sub_beta = subvector(beta, 1UL, (n-1) );
                h = ctrans(Q) * A * Q;

                std::cout << "alpha is: " << alpha << std::endl;
                std::cout << "beta is: " << beta << std::endl;

            }; // end solve_imple function

        } //end namespace detail

    ITERATIVE_NAMESPACE_CLOSE
BLAZE_NAMESPACE_CLOSE

#endif //BLAZE_ITERATIVE_LANCZOS_HPP
