// Copyright (c)   2017 Tyler Olsen
//                 2018 Patrick Diehl
//                 2019 Nanmiao Wu
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BLAZE_ITERATIVE_GMRES_HPP
#define BLAZE_ITERATIVE_GMRES_HPP

#include <iostream>
#include "IterativeCommon.hpp"
#include "GMRESTag.hpp"
#include <utility>


BLAZE_NAMESPACE_OPEN
    ITERATIVE_NAMESPACE_OPEN

        namespace detail {

            template< typename MatrixType, typename T>
            std::pair< DynamicMatrix<T>,  DynamicVector<T> > arnoldi( const MatrixType &A, DynamicMatrix<T> &Q, int &k)
            {
                std::size_t m = A.columns();
                DynamicVector<T> q(m);
                DynamicMatrix<T> h(k+2, 1);

                q = A * column(Q,k);
                for(int i = 0; i <= k; ++i){
                    h(i,0) = ctrans(q) * column(Q,i);
                    q -= h(i,0) * column(Q,i);
                }
                h(k+1, 0) = norm(q);
                q = q / h(k+1, 0);
                return std::make_pair(h,q);
            }

            template<typename T>
            std::pair< T, T> givens_rotation( T &v1, T &v2)
            {
                T cs_k, sn_k;
                if (v1 == 0){
                    cs_k = 0;
                    sn_k = 0;

                } else {
                    auto t = sqrt(v1 * v1 + v2 * v2);
                    cs_k = abs(v1) / t;
                    sn_k = cs_k * v2 / v1;
                }

                return std::make_pair(cs_k,sn_k);
            }

            template<typename T>
            std::tuple<DynamicMatrix<T>,  T, T> apply_givens_rotation( DynamicMatrix<T> &h, DynamicVector<T> &cs, DynamicVector<T> &sn, int &k)
            {
                for(int i = 0; i < k; ++i){
                    auto temp = cs[i] * h(i,0) + sn[i] * h(i+1,0);
                    h(i+1,0) = -sn[i] * h(i,0) + cs[i] * h(i+1,0);
                    h(i,0) = temp;
                }

                T cs_k, sn_k;
                auto res = givens_rotation(h(k,0), h(k+1,0));
                cs_k = res.first;
                sn_k = res.second;

                h(k,0) = cs_k * h(k,0) + sn_k * h(k+1,0);
                h(k+1,0) = 0.0;

                return std::make_tuple(h, cs_k, sn_k);
            }

            template<typename MatrixType, typename T>
            void  solve_impl(
                    DynamicVector<T> &x,
                    const MatrixType &A,
                    const DynamicVector<T> &b,
                    const DynamicVector<T> &x0,
                    GMRESTag &tag,
                    const std::size_t &n)
            {

                BLAZE_INTERNAL_ASSERT(A.isSymmetric(), "A must be a symmetric matrix")

                BLAZE_INTERNAL_ASSERT(n >= 1, "n must larger than or equal to 1")

                // A: m * m matrix; n is max_iteration;


                std::size_t m = A.columns();
                DynamicVector<T> r(m);
                DynamicMatrix<T> H(n+1,n,0);
                DynamicMatrix<T> Q(m,n+1,0);
                DynamicVector<T> sn(n,0);
                DynamicVector<T> cs(n,0);
                DynamicVector<T> e1(m,0);
                DynamicVector<T> beta(m);
                DynamicVector<T> err_set(n+1);


                r = b - A * x0;
                auto err = norm(r) / norm(b);
                err_set[0] = err;

                e1[0] = 1;
                column(Q,0) = r / norm(r);
                beta = norm(r) * e1;

                for(int k = 0; k < n; ++k){
                    auto res_1 = arnoldi(A, Q, k);
                    column(Q,k+1) = res_1.second;

                    auto res_2 = apply_givens_rotation(res_1.first, cs, sn, k);
                    submatrix(H, 0, k, k+2, 1) = std::get<0>(res_2);
                    cs[k] = std::get<1>(res_2);
                    sn[k] = std::get<2>(res_2);

                    beta[k+1] = -sn[k] * beta[k];
                    beta[k] = cs[k] * beta[k];

                    err = abs(beta[k+1]) / norm(b);
                    err_set[k+1] = err;

                    double eps = 1e-8;
                    if (err <= eps){
                        break;
                    }

                }

                auto H_sub = submatrix(H, 0, 0, n, n);
                auto beta_sub = subvector(beta, 0, n);
                auto y = inv(H_sub) *  beta_sub;
                auto Q_sub = submatrix(Q, 0, 0, n, n);
                x += Q_sub * y;



            }; // end solve_imple function

        } //end namespace detail

    ITERATIVE_NAMESPACE_CLOSE
BLAZE_NAMESPACE_CLOSE

#endif //BLAZE_ITERATIVE_GMRES_HPP
