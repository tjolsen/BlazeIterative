// Copyright (c)   2017 Tyler Olsen
//                 2018 Patrick Diehl
//                 2019 Nanmiao Wu
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BLAZE_ITERATIVE_PRECONDITIONBICGSTABL_HPP
#define BLAZE_ITERATIVE_PRECONDITIONBICGSTABL_HPP

#include "PreconditionBiCGSTABLTag.hpp"

BLAZE_NAMESPACE_OPEN
    ITERATIVE_NAMESPACE_OPEN

        namespace detail {

/**
 *  Implementation of the decompositions provided by blaze
 */

/**
 *  Implementation of the Preconditioned BiCGSTABL method
 *  Based on the paper "BICGSTAB(l) for linear equations involving unsymmetric matrices with complex spectrum"
 */
            template<typename MatrixType, typename T>
            void solve_impl(
                    DynamicVector<T> &x,
                    const MatrixType &A,
                    const DynamicVector<T> &b,
                    const int &l,
                    PreconditionBiCGSTABLTag &tag,
                    std::string Preconditioner="")
            {


                std::size_t m = A.columns();

                DynamicVector<T> r0 = b - A * x;
                DynamicVector<T> r0_tilde(r0);

                DynamicVector<T> u_minus1(r0.size());
                DynamicVector<T> x0(x);
                DynamicVector<T> x0_hat(x);
                DynamicVector<T> sigma(l, 1);
                DynamicVector<T> gamma_dot_0(l);
                DynamicVector<T> gamma_dot_1(l);
                DynamicVector<T> gamma_dot_2(l);

                auto absolute_residual_0 = trans(r0) * r0;
                auto absolute_residual = absolute_residual_0;

                auto rho_0 = T(1);
                auto alpha = T(0);
                auto w = T(1);

                DynamicMatrix<T> r_hat(m, l + 1);
                DynamicMatrix<T> u_hat(m, l + 1);
                DynamicMatrix<T> tao(l - 1, l, 0);

                int k = -l;

                auto error = 1;
                double eps = 1e-8;

                while (error > eps ) {

                    k += l;

                    column(u_hat, 0) = u_minus1;
                    column(r_hat, 0) = r0;
                    x0_hat = x0;
                    rho_0 = - w * rho_0;

                    // Bi-CG part
                    for (int j = 0; j < l; ++j){
                        auto rho_1 = trans(column(r_hat, j)) * r0_tilde;
                        auto beta = alpha * rho_1 / rho_0;
                        rho_0 = rho_1;

                        for (int i = 0; i <=j; ++i){
                            column(u_hat, i) = column(r_hat, i) - beta * column(u_hat, i);
                        }

                        column(u_hat, j + 1) = A * column(u_hat, j);

                        auto gamma = trans(column(u_hat, j + 1)) * r0_tilde;
                        alpha = rho_0 / gamma;

                        for (int i = 0; i <= j; ++i){
                            column(r_hat, i) -= alpha * column(u_hat, i + 1);
                        }

                        column(r_hat, j + 1) = A * column(r_hat, j);
                        x0_hat += alpha * column(u_hat, 0);
                    }


                    // MR part
                    for (int j = 0; j < l; ++j){
                        for (int i = 0; i < j-1; ++i){
                            tao(i, j) = trans(column(r_hat, j)) * column(r_hat, i) / sigma[i];
                            column(r_hat, j) -= tao(i, j) * column(r_hat, i);
                        }
                        sigma[j] = trans(column(r_hat, j)) * column(r_hat, j);
                        gamma_dot_1[j] = trans(column(r_hat, 0)) * column(r_hat, j) / sigma[j];
                    }

                    gamma_dot_0[l-1] = gamma_dot_1[l-1];
                    w = gamma_dot_0[l-1];

                    for ( int j = l - 2; j >= 0; --j){
                        auto sum_1 = 0;
                        for ( int i = j+1; i < l; ++i){
                            sum_1 += tao(j, i) * gamma_dot_0[i];
                        }
                        gamma_dot_0[j] = gamma_dot_1[j] - sum_1;
                    }

                    for ( int j = 0; j <= l - 2; ++j){
                        auto sum_2 = 0;
                        for ( int i = j+1; i < l - 1; ++i){
                            sum_2 += tao(j, i) * gamma_dot_0[i + 1];
                        }
                        gamma_dot_2[j] = gamma_dot_0[j + 1] + sum_2;
                    }

                    x0_hat += gamma_dot_0[0] * column(r_hat, 0);
                    column(r_hat, 0) -= gamma_dot_1[l - 1] * column(r_hat, l);
                    column(u_hat, 0) -= gamma_dot_0[l - 1] * column(u_hat, l);

                    for ( int j = 0; j <= l-2; ++j){
                        column(u_hat, 0) -= gamma_dot_0[j] * column(u_hat, j);
                        x0_hat += gamma_dot_2[j] * column(r_hat, j);
                        column(r_hat, 0) -= gamma_dot_1[j] * column(r_hat, j);
                    }

                    u_minus1 = column(u_hat, 0);
                    r0 = column(r_hat, 0);
                    x0 = x0_hat;

                    error = norm(r0);
                }
            }

        } //end namespace detail

    ITERATIVE_NAMESPACE_CLOSE
BLAZE_NAMESPACE_CLOSE

#endif //BLAZE_ITERATIVE_PRECONDITIONBICGSTABL_HPP
