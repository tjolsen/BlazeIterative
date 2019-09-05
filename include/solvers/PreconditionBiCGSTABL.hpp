// Copyright (c)   2017 Tyler Olsen
//                 2018 Patrick Diehl
//                 2019 Nanmiao Wu
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BLAZE_ITERATIVE_PRECONDITIONBICGSTABL_HPP
#define BLAZE_ITERATIVE_PRECONDITIONBICGSTABL_HPP

#include "PreconditionBiCGSTABLTag.hpp"
#include "PreconditionCG.hpp"
#include <iostream>

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
                    const std::size_t &l,
                    PreconditionBiCGSTABLTag &tag
                    )
            {


                std::size_t m = A.columns();

//                MatrixType M;
//                preconditioner_matrix<MatrixType, T>(Preconditioner,A,M);

                DynamicVector<T> r0 = b - A * x;
                DynamicVector<T> r0_tilde(r0);

                DynamicVector<T> u0(m, 0);
                DynamicVector<T> x0(x);
                DynamicVector<T> x0_hat(x);
                DynamicVector<T> sigma(l);
                DynamicVector<T> gamma_dot_0(l);
                DynamicVector<T> gamma_dot_1(l);
                DynamicVector<T> gamma_dot_2(l-1);

                auto absolute_residual_0 = trans(r0) * r0;
                auto absolute_residual = absolute_residual_0;

                auto rho_0 = 1.0;
                auto alpha = 0.0;
                auto w = 1.0;

                DynamicMatrix<T> r_hat(m, l + 1);
                DynamicMatrix<T> u_hat(m, l + 1);
                DynamicMatrix<T> tao(l - 1, l, 0);

                int k = -l;

                std::size_t iteration{0};

                while (true) {

                    k += l;
//                    std::cout << "******************************************************************" << std::endl;
//                    std::cout << "k is: " << k << std::endl;

                    column(u_hat, 0) = u0;
                    column(r_hat, 0) = r0;
                    x0_hat = x0;
                    rho_0 = - w * rho_0;

//                    std::cout << "After increasing k, w is: " << w << std::endl;
//
//                    std::cout << "After increasing k, rho_0 is: " << rho_0 << std::endl;
//
//                    std::cout << "Before Bi-CG part, u_hat_0 is: " << column(u_hat, 0) << std::endl;
//                    std::cout << "Before Bi-CG part, r_hat_0 is: " << column(r_hat, 0) << std::endl;
//                    std::cout << "Before Bi-CG part, x0_hat is: " << x0_hat << std::endl;


                    // Bi-CG part
                    for (int j = 0; j < l; ++j){
                        auto rho_1 = trans(column(r_hat, j)) * r0_tilde;
//                        std::cout << "In Bi-CG part, when j = " << j <<", rho_0 is: " << rho_0 << std::endl;
//                        std::cout << "In Bi-CG part, when j = " << j <<", alpha is: " << alpha << std::endl;
                        if (rho_0 == 0){
                            std::cout << "break down" << std::endl;
                            break;
                        }
                        auto beta = alpha * rho_1 / rho_0;
                        rho_0 = rho_1;

//                        std::cout << "In Bi-CG part, when j = " << j <<", rho_1 is: " << rho_1 << std::endl;
//                        std::cout << "In Bi-CG part, when j = " << j <<", beta is: " << beta << std::endl;

                        for (int i = 0; i <= j; ++i){
                            column(u_hat, i) = column(r_hat, i) - beta * column(u_hat, i);
                        }


                        column(u_hat, j + 1) = A * column(u_hat, j);

//                        std::cout << "In Bi-CG part, when j = " << j <<", u_hat_[j] is: " << column(u_hat, j) << std::endl;
//                        std::cout << "In Bi-CG part, when j = " << j <<", u_hat_[j+1] is: " << column(u_hat, j + 1) << std::endl;

                        auto gamma = trans(column(u_hat, j + 1)) * r0_tilde;
                        if (gamma == 0){
                            std::cout << "break down #2" << std::endl;
                            break;
                        }
//                        std::cout << "In Bi-CG part, when j = " << j <<", gamma is: " << gamma << std::endl;
//                        std::cout << "In Bi-CG part, when j = " << j <<", rho_0 is: " << rho_0 << std::endl;
                        alpha = rho_0 / gamma;

//                        std::cout << "In Bi-CG part, when j = " << j <<", alpha is: " << alpha << std::endl;

                        for (int i = 0; i <= j; ++i){
//                            std::cout << "In Bi-CG part, when j = " << j <<" i = " << i  <<", r_hat_0 is: " << column(r_hat, i) << std::endl;
//                            std::cout << "In Bi-CG part, when j = " << j <<", alpha is: " << alpha << std::endl;
//                            std::cout << "In Bi-CG part, when j = " << j <<" i = " << i  <<", u_hat_1 is: " << column(u_hat, i + 1) << std::endl;
                            column(r_hat, i) -= alpha * column(u_hat, i + 1);
//                            std::cout << "In Bi-CG part, when j = " << j <<" i = " << i  <<", r_hat_0 is: " << column(r_hat, i) << std::endl;
                        }

                        column(r_hat, j + 1) = A * column(r_hat, j);
                        x0_hat += alpha * column(u_hat, 0);

//                        std::cout << "In Bi-CG part, when j = " << j <<", r_hat_0 is: " << column(r_hat, 0) << std::endl;
//                        std::cout << "In Bi-CG part, when j = " << j <<", x0_hat is: " << x0_hat << std::endl;
                    }
                    std::cout << "******************************************************************" << std::endl;

                    std::cout << "After Bi-CG part, u_hat_0 is: " << column(u_hat, 0) << std::endl;
                    std::cout << "After Bi-CG part, r_hat_0 is: " << column(r_hat, 0) << std::endl;
                    std::cout << "After Bi-CG part, x0_hat is: " << x0_hat << std::endl;


                    // MR part
                    for (int j = 0; j < l; ++j){
                        for (int i = 0; i < j-1; ++i){
                            tao(i, j) = trans(column(r_hat, j + 1)) * column(r_hat, i + 1) / sigma[i];
                            column(r_hat, j + 1) -= tao(i, j) * column(r_hat, i + 1);
                        }
                        sigma[j] = trans(column(r_hat, j + 1)) * column(r_hat, j + 1);
                        gamma_dot_1[j] = trans(column(r_hat, 0)) * column(r_hat, j + 1) / sigma[j];
                    }

//                    std::cout << "In MR part, gamma_dot_1[l-1] is: " << gamma_dot_1[l-1] << std::endl;


                    gamma_dot_0[l-1] = gamma_dot_1[l-1];

//                    std::cout << "In MR part, gamma_dot_0[l-1] is: " << gamma_dot_0[l-1] << std::endl;
//
//                    std::cout << "In MR part, w is: " << w << std::endl;

                    w = gamma_dot_0[l-1];

//                    std::cout << "In MR part, gamma_dot_0[l-1] is: " << gamma_dot_0[l-1] << std::endl;
//                    std::cout << "In MR part, w is: " << w << std::endl;

                    for ( int j = l - 2; j >= 0; --j){
                        auto sum_1 = 0;
                        for ( int i = j + 1; i < l; ++i){
                            sum_1 += tao(j, i) * gamma_dot_0[i];
                        }
                        gamma_dot_0[j] = gamma_dot_1[j] - sum_1;
                    }

                    for ( int j = 0; j < l - 1; ++j){
                        auto sum_2 = 0;
                        for ( int i = j+1; i < l - 1; ++i){
                            sum_2 += tao(j, i) * gamma_dot_0[i + 1];
                        }
                        gamma_dot_2[j] = gamma_dot_0[j + 1] + sum_2;
                    }

                    x0_hat += gamma_dot_0[0] * column(r_hat, 0);
                    column(r_hat, 0) -= gamma_dot_1[l - 1] * column(r_hat, l);
                    column(u_hat, 0) -= gamma_dot_0[l - 1] * column(u_hat, l);

                    for ( int j = 0; j < l - 1; ++j){
                        column(u_hat, 0) -= gamma_dot_0[j] * column(u_hat, j + 1);
                        x0_hat += gamma_dot_2[j] * column(r_hat, j + 1);
                        column(r_hat, 0) -= gamma_dot_1[j] * column(r_hat, j + 1);
                    }
//                    std::cout << "******************************************************************" << std::endl;
//
//                    std::cout << "After MR part, u_hat_0 is: " << column(u_hat, 0) << std::endl;
//                    std::cout << "After MR part, r_hat_0 is: " << column(r_hat, 0) << std::endl;
//                    std::cout << "After MR part, x0_hat is: " << x0_hat << std::endl;



                    u0 = column(u_hat, 0);
                    r0 = column(r_hat, 0);
                    x0 = x0_hat;
                    x = x0;


                    absolute_residual = norm(r0);
                    auto relative_residual = absolute_residual/absolute_residual_0;
                    if(tag.do_log()) {
                        tag.log_residual(relative_residual);
                    }

                    if(tag.terminateIteration(iteration,absolute_residual,relative_residual)) {
                        break;
                    }
                    ++iteration;

//                    std::cout << "At the end, when k is: " << k << ", rho_0 is: " << rho_0 << std::endl;
//                    std::cout << "At the end, when k is: " << k << ", w is: " << w << std::endl;

                }
            }

        } //end namespace detail

    ITERATIVE_NAMESPACE_CLOSE
BLAZE_NAMESPACE_CLOSE

#endif //BLAZE_ITERATIVE_PRECONDITIONBICGSTABL_HPP
