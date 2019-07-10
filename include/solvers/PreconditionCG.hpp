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
        // Decompose A as A = L + D + U, where L is strictly lower triangular matrix, D is diagonal matrix, and U is
        // strictly upper triangular matrix.
        template<typename MatrixType, typename T>
        void preconditioner_matrix(std::string type, const MatrixType &A, MatrixType &M){

            if (type.compare("Jacobi preconditioning")==0){
                // M = D
                std::size_t n = A.columns();
                MatrixType D(n,n);
                for(std::size_t i = 0; i < n; ++i){
                    for(std::size_t j = 0; j < n; ++j){
                        if(j == i){
                            D(i,j) = A(i,j);
                        }
                        else {
                            D(i,j) = 0;
                        }
                    }
                }
                M = D;
            }

            if (type.compare("Gauss Seidel preconditioning")==0){
                std::size_t n = A.columns();
                MatrixType L(n,n), D(n,n);
                for(std::size_t i = 0; i < n; ++i){
                    for(std::size_t j = 0; j < n; ++j){
                        if(j < i){
                            L(i,j) = A(i,j);
                        }
                        else{
                            L(i,j) = 0;
                        }
                    }
                }

                for(std::size_t i = 0; i < n; ++i){
                    for(std::size_t j = 0; j < n; ++j){
                        if(j == i){
                            D(i,j) = A(i,j);
                        }
                        else {
                            D(i,j) = 0;
                        }
                    }
                }

                M = L + D; // M = L + D
            }

            if (type.compare("Symmetric Gauss Seidel preconditioning")==0){

                std::size_t n = A.columns();
                MatrixType L(n,n), D(n,n), U(n,n);
                for(std::size_t i = 0; i < n; ++i){
                    for(std::size_t j = 0; j < n; ++j){
                        if(j < i){
                            L(i,j) = A(i,j);
                        }
                        else{
                            L(i,j) = 0;
                        }
                    }
                }

                for(std::size_t i = 0; i < n; ++i){
                    for(std::size_t j = 0; j < n; ++j){
                        if(j == i){
                            D(i,j) = A(i,j);
                        }
                        else {
                            D(i,j) = 0;
                        }
                    }
                }

                for(std::size_t i = 0; i < n; ++i){
                    for(std::size_t j = 0; j < n; ++j){
                        if(i < j){
                            U(i,j) = A(i,j);
                        }
                        else{
                            U(i,j) = 0;
                        }
                    }
                }


                M = (D + L) * inv(D) * (D + U);
            }

            if (type.compare("SOR preconditioning")==0){
                double omega = 1.2;  // omega is in the range of (0,2) to make sure converge
                std::size_t n = A.columns();
                MatrixType L(n,n), D(n,n);
                for(std::size_t i = 0; i < n; ++i){
                    for(std::size_t j = 0; j < n; ++j){
                        if(j < i){
                            L(i,j) = A(i,j);
                        }
                        else{
                            L(i,j) = 0;
                        }
                    }
                }

                for(std::size_t i = 0; i < n; ++i){
                    for(std::size_t j = 0; j < n; ++j){
                        if(j == i){
                            D(i,j) = A(i,j);
                        }
                        else {
                            D(i,j) = 0;
                        }
                    }
                }

                M = (D + omega * L)/omega;
            }

            if (type.compare("SSOR preconditioning")==0){
                std::size_t n = A.columns();
                MatrixType L(n,n), D(n,n);
                for(std::size_t i = 0; i < n; ++i){
                    for(std::size_t j = 0; j < n; ++j){
                        if(j < i){
                            L(i,j) = A(i,j);
                        }
                        else{
                            L(i,j) = 0;
                        }
                    }
                }

                for(std::size_t i = 0; i < n; ++i){
                    for(std::size_t j = 0; j < n; ++j){
                        if(j == i){
                            D(i,j) = A(i,j);
                        }
                        else {
                            D(i,j) = 0;
                        }
                    }
                }

                M = (D + L) * inv(D) * trans(D + L);
            }


            if (type.compare("incomplete Cholesky factorization") == 0 || type.compare("") == 0){
                // The Cholesky factorization of A is A = LL*, where L is a lower triangular matrix.
                // Incomplete Cholesky factorization precondition is M =KK*, where a sparse lower triangular matrix K is close to L.
                // Finding the exact Cholesky decomposition, except that any entry is set to zero
                // if the corresponding entry in A is also zero.

                MatrixType L, K;
                llh( A, L );  // L* LH decomposition of a row-major matrix
                K = L;

                std::size_t n = A.columns();
                for(std::size_t i = 0; i < n; ++i){
                    for(std::size_t j = 0; j < n; ++j){
                        if (A(i, j) == 0)
                            K(i, j) = 0;
                    }
                }

                M = K * trans(K);

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
            preconditioner_matrix<MatrixType, T>(Preconditioner,A,M);

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
