
//
// Created by tyler on 6/24/17.
//

#include "BlazeIterative.hpp"
#include <iostream>

using namespace blaze;
using namespace blaze::iterative;
// Copyright (c)   2017 Tyler Olsen
//                 2018 Patrick Diehl
//                 2019 Nanmiao Wu
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

int main() {
    // Test for ConjugateGradient
    // Test for BiCGSTAB
    // Test for PreconditionBiCGSTAB
    // Test for PreconditionCG
    // Test for GMRES

//    std::size_t N = 3;
//    DynamicMatrix<double,false> A(N,N, 0.0);
//    DynamicVector<double> b(N, 0.0);
//    for(int i=0; i<N; ++i) {
//        A(i,i) = 2.0;
//        b[i] = 1.0*(1+i);
//    }

/*
    DynamicMatrix<double,false> L(N,N,0.0);
    //copyStrictlyLowerPart<double, false>(A, L);
    std::cout << L << std::endl;
*/

    //ConjugateGradientTag tag;
    //BiCGSTABTag tag;


    //std::cout << solve(A,b,tag) << std::endl << std::endl;
/*
    PreconditionBiCGSTABTag tag;
    tag.do_log() = true;
    std::cout << solve(A,b,tag, "Cholesky") << std::endl << std::endl;
   */


//    PreconditionCGTag tag;
//    tag.do_log() = true;
//    std::cout << solve(A,b,tag, "incomplete Cholesky factorization") << std::endl << std::endl;
//
//   int iter(0);
//   for(auto r : tag.convergence_history()) {
//       std::cout << iter++ << '\t' << r << '\n';
//   }

//  Test for GMRES
    std::size_t N = 5;
    DynamicMatrix<double,false> A(N,N);
    DynamicVector<double> b(N);
    A = {{0.46, 0.60, 0.74, 0.61, 0.85}, {0.56, 0.31, 0.80, 0.94, 0.76},{0.41, 0.19, 0.15, 0.33, 0.06},{0.03, 0.92, 0.15, 0.56, 0.08},{0.09, 0.06, 0.69, 0.42, 0.96}};
    //x should be {0.1, 0.3, 0.5, 0.7, 0.9};
    b = {1.788, 1.891, 0.458, 0.818, 1.53};
    GMRESTag tag;
    tag.do_log() = true;
    DynamicVector<double> x0(N,0);
    std::size_t n = N;
    std::cout << solve(A,b,x0,tag,n) << std::endl << std::endl;


    // Test For PreconditionBiCGSTABL





//    // Test Arnoldi
//    // Test Lanczos
//    std::size_t N = 6;
//    DynamicMatrix<double,columnMajor> A(N,N,0);
//    band<0>(A) = {0, 1, 2, 3, 4, 100000};
//
//
//    DynamicVector<double> b(N, 1);
//    DynamicVector<complex<double>,columnVector> w(N); // The vector for the real eigenvalues
//    DynamicMatrix<complex<double>,rowMajor> V(N,N); // The matrix for the left eigenvectors
//    eigen(A,w,V);
//    std::cout<< "The eigenvalues of Matrix A is: " << std::endl << w << std::endl;
//
//    std::size_t n = N;
//    ArnoldiTag tag1;
//    DynamicVector<complex<double>,columnVector> w2(n);
//    DynamicMatrix<complex<double>,rowMajor> V2(n,n);
//    auto res1 = solve(A,b,tag1,n);
//    auto sub_h = submatrix( res1.second, 0UL, 0UL, (res1.second.rows()-1), res1.second.columns());
//    eigen(sub_h,w2,V2);
//    std::cout << "Arnoldi: The eigenvalues of Matrix h is: " <<std::endl << w2 << std::endl;
//
//    DynamicVector<complex<double>,columnVector> w1(n);
//    DynamicMatrix<complex<double>,rowMajor> V1(n,n);
//    LanczosTag tag;
//    auto res = solve(A,b,tag,n);
//    eigen(res.second,w1,V1);
//    std::cout << "Lanczos: The eigenvalues of Matrix h is: "  <<std::endl << w1 << std::endl;




    return 0;
}