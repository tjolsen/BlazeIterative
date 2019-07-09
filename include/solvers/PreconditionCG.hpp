// Copyright (c)   2017 Tyler Olsen
//                 2018 Patrick Diehl
//                 2019 Nanmiao Wu
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BLAZE_ITERATIVE_PRECONDITIONCG_HPP
#define BLAZE_ITERATIVE_PRECONDITIONCG_HPP

#include "PreconditionCGTag.hpp"

BLAZE_NAMESPACE_OPEN
ITERATIVE_NAMESPACE_OPEN

    namespace detail {

        template<typename MatrixType, typename T>
        void preconditioner_matrix(std::string type, const MatrixType &A, MatrixType &M){


            if (type.compare("Jacobi preconditioning")==0){
                auto diag = diagonal( A );

                auto row_A = row(A, 0);
                std::size_t n = row_A.size();

                M =

            }

            if (type.compare("Gauss Seidel preconditioning")==0){

                M = invert<asLower>(A) + invert<asUpper>(A); // M = L + D
            }

            if (type.compare("Symmetric Gauss Seidel preconditioning")==0){
                auto D = invert<asDiagonal>(A);
                auto L = invert<asLower>(A);
                auto U = invert<asUpper>(A);

                M = (D + L) * inv(D) * (D + U);
            }

            if (type.compare("SOR preconditioning")==0){
                double omega = 1.2;  // omega is in the range of (0,2) to make sure converge
                auto D = invert<asDiagonal>(A);
                auto L = invert<asLower>(A);

                M = (D + omega * L)/omega;
            }

            if (type.compare("SSOR preconditioning")==0){
                auto D = invert<asDiagonal>(A);
                auto L = invert<asLower>(A);

                M = (D + L) * inv(D) * trans(D + L);
            }


            if (type.compare("incomplete Cholesky factorization") == 0 || type.compare("") == 0){
                auto row_A = row(A, 0);
                std::size_t n = row_A.size();
                MatrixType B = A;

                for(std::size_t k = 0; k < n; ++k){
                    B(k,k) = sqrt(B(k,k));

                    for(std::size_t i = k+1; i < n; ++i){
                        if(B(i,k) != 0)
                            B(i,k) = B(i,k)/B(k,k);
                    }

                    for(std::size_t j = k+1; j < n; ++j){
                        for(std::size_t i = j; i < n; ++i){
                            if (B(i,j) != 0)
                                B(i,j) = B(i,j) - B(i,k) * B(j,k);
                        }
                    }
                }

                for(std::size_t i = 0; i < n; ++i){
                    for(std::size_t j = i+1; j < n; ++j){
                        B(i,j) = 0;
                    }
                }

                M = B;

            }


        }

        template<typename MatrixType, typename T>
        void solve_impl(
                DynamicVector<T> &x,
                const MatrixType &A,
                const DynamicVector<T> &b,
                PreconditionCGTag &tag,
                std::string Preconditioner="")
        {

            BLAZE_INTERNAL_ASSERT(A.isSymmetric(), "A must be a symmetric matrix")

            MatrixType M;
            preconditioner_matrix<MatrixType>(Preconditioner,A,M);

            auto Minv = inv(M);

            DynamicVector<T> r = b - A * x;
            DynamicVector<T> z = Minv * r;
            DynamicVector<T> p(z);
            DynamicVector<T> Ap(p.size());

            auto absolute_residual_0 = trans(r) * r;
            auto absolute_residual = absolute_residual_0;
            //auto absolute_residual_prev = absolute_residual;

            if(tag.do_log()) {
                tag.log_residual(absolute_residual/absolute_residual_0);
            }


            std::size_t iteration{0};
            while(true) {
                //absolute_residual_prev = absolute_residual;
                Ap = declsym(A)*p;

                auto alpha = trans(r) * z/(trans(p) * Ap);
                auto precondition_residual_prev = trans(z) * r;
                x += alpha * p;
                r -= alpha * Ap;


                absolute_residual = trans(r) * r;

                if(tag.do_log()) {
                    tag.log_residual(absolute_residual/absolute_residual_0);
                }

                if(tag.terminateIteration(iteration, absolute_residual, absolute_residual/absolute_residual_0)) {
                    break;
                }

                z = Minv * r;
                auto beta = trans(z)*r/precondition_residual_prev;
                p = z + beta * p;

                ++iteration;
            }//end while


        };


    } //end namespace detail        } //end namespace detail

ITERATIVE_NAMESPACE_CLOSE
BLAZE_NAMESPACE_CLOSE

#endif //BLAZE_ITERATIVE_PRECONDITIONCG_HPP
