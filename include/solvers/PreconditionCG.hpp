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

        template< typename MatrixType>
        void copyStrictlyLowerPart( const MatrixType& src, MatrixType& dst )

        {

            decltype(auto) target( derestrict( ~dst ) );
            assert(dst.rows() == target.rows());
            assert(dst.columns() == target.columns());


            for(signed long int i=1UL; i<(~src).rows(); ++i ) {

                for(signed long int j=0UL; j<i; ++j ) {

                    (~target)(i,j) = (~src)(i,j);

                }

            }

        }

        // Decompose A as A = L + D + U, where L is strictly lower triangular matrix, D is diagonal matrix, and U is
        // strictly upper triangular matrix.
        template<typename MatrixType, typename T>
        void preconditioner_matrix(std::string type, const MatrixType &A, MatrixType &M){

            if (type.compare("Jacobi preconditioning")==0){

                DiagonalMatrix<MatrixType> D(A.columns());
                band<0L>(D) = band<0L>(A);
                M = D;
            }

            if (type.compare("Gauss Seidel preconditioning")==0){
                MatrixType L(A.rows(),A.columns());
                DiagonalMatrix<MatrixType> D(A.columns());
                copyStrictlyLowerPart<MatrixType>(A, L);
                band<0L>(D) = band<0L>(A);

                M = L + D;
            }

            if (type.compare("Symmetric Gauss Seidel preconditioning")==0){

                MatrixType L(A.rows(),A.columns());
                copyStrictlyLowerPart<MatrixType>(A, L);

                DiagonalMatrix<MatrixType> D(A.columns());
                band<0L>(D) = band<0L>(A);

                M = (D + L) * inv(D) * (D + trans(L));
            }

            if (type.compare("SOR preconditioning")==0){
                double omega = 1.2;  // omega is in the range of (0,2) to make sure converge
                DiagonalMatrix<MatrixType> D(A.columns());
                MatrixType L(A.rows(),A.columns());
                band<0L>(D) = band<0L>(A);
                copyStrictlyLowerPart<MatrixType>(A, L);
                M = (D + omega * L)/omega;
            }

            if (type.compare("SSOR preconditioning")==0){
                DiagonalMatrix<MatrixType> D(A.columns());
                MatrixType L(A.rows(),A.columns());
                copyStrictlyLowerPart<MatrixType>(A, L);
                band<0L>(D) = band<0L>(A);

                M = (D + L) * inv(D) * trans(D + L);
            }


            if (type.compare("incomplete Cholesky factorization") == 0 || type.compare("") == 0){
                // The Cholesky factorization of A is A = LL*, where L is a lower triangular matrix.
                // Incomplete Cholesky factorization precondition is M = KK*, where K is a sparse lower triangular matrix and is close to L.
                // Solution is: finding the exact Cholesky decomposition, except that any entry is set to zero
                // if the corresponding entry in A is also zero.

                MatrixType L(A.rows(),A.columns()), K(A.rows(),A.columns());
                llh( A, L );  // L* LH decomposition of a row-major matrix
                K = map(A, L, [](T aval, T lval) { return aval == T{0} ? T{0} : lval; });
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

            if(tag.do_log()) {
                tag.log_residual(absolute_residual/absolute_residual_0);
            }


            std::size_t iteration{0};
            while(true) {
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
