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

   std::size_t N = 1000;
   DynamicMatrix<double,columnMajor> A(N,N, 0.0);
   DynamicVector<double> b(N, 0.0);

   for(int i=0; i<N; ++i) {
       A(i,i) = 2.0;
       b[i] = 1.0*(1+i);
   }

   //ConjugateGradientTag tag;
   //BiCGSTABTag tag;


   //std::cout << solve(A,b,tag) << std::endl << std::endl;
/*
    PreconditionBiCGSTABTag tag;
    tag.do_log() = true;
    std::cout << solve(A,b,tag, "Cholesky") << std::endl << std::endl;
   */

    PreconditionCGTag tag;
    tag.do_log() = true;
    std::cout << solve(A,b,tag, "SSOR preconditioning") << std::endl << std::endl;

   int iter(0);
   for(auto r : tag.convergence_history()) {
       std::cout << iter++ << '\t' << r << '\n';
   }


    /*
    // Test Arnoldi
    std::size_t N = 10;

    DynamicMatrix<double,columnMajor> A(N,N, 0.0);
    for(int i=1; i<N-1; ++i) {
        A(i, i-1) = -1.0;
        A(i,i) = 2.0;
        A(i,i+1) = -1.0;
    }
    A(0,1) = -1.0;
    A(0,0) = 2.0;
    A(9,8) = -1.0;
    A(9,9) = 2.0;
    std::cout<< "The Matrix A is: " << std::endl << A << std::endl;

    DynamicVector<double> b(N, 1.0);
    DynamicVector<complex<double>,columnVector> w(N); // The vector for the real eigenvalues
    DynamicMatrix<complex<double>,rowMajor> V(N,N); // The matrix for the left eigenvectors
    eigen(A,w,V);
    std::cout<< "The eigenvalues of Matrix A is: " << std::endl << w << std::endl;
    //std::cout << "The eigenvectors of Matrix A is: " << V << std::endl;

    DynamicVector<complex<double>,columnVector> w1(N); // The vector for the real eigenvalues
    DynamicMatrix<complex<double>,rowMajor> V1(N,N); // The matrix for the left eigenvectors

    std::size_t n = N-1;
    //DynamicMatrix<double> Q(N, (n + 1));
    //DynamicMatrix<double> h((n + 1), n);
    ArnoldiTag tag;
    auto res = solve(A,b,tag,n);
    //eigen(h,w1,V1);
    std::cout << "The Matrix h is: " <<std::endl << res.second << std::endl;

    */

    return 0;
}
