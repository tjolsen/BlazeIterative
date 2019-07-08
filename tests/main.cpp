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
    /*
    * // Test for ConjugateGradient
    * // Test for BiCGSTAB
    * // Test for PreconditionBiCGSTAB
   std::size_t N = 1000;
   DynamicMatrix<double,columnMajor> A(N,N, 0.0);
   DynamicVector<double> b(N, 0.0);

   for(int i=0; i<N; ++i) {
       A(i,i) = 2.0;
       b[i] = 1.0*(1+i);
   }

    *
   //ConjugateGradientTag tag;
   //BiCGSTABTag tag;


   //std::cout << solve(A,b,tag) << std::endl << std::endl;

   PreconditionBiCGSTABTag tag;
   tag.do_log() = true;
   std::cout << solve(A,b,tag, "Cholesky") << std::endl << std::endl;

   int iter(0);
   for(auto r : tag.convergence_history()) {
       std::cout << iter++ << '\t' << r << '\n';
   }
   */

    // Test Arnoldi
    std::size_t N = 3;
    DynamicMatrix<double,columnMajor> A{{-2, -4, 2},
                                        {-2, 1, 2},
                                        {4, 2, 5}};
    DynamicVector<double> b(N, 1.0);
    DynamicVector<complex<double>,columnVector> w(N); // The vector for the real eigenvalues
    DynamicMatrix<complex<double>,rowMajor> V(N,N); // The matrix for the left eigenvectors
    eigen(A,w,V);
    std::cout<< "The eigenvalues of Matrix A is: " << w << std::endl;
    std::cout << "The eigenvectors of Matrix A is: " << V << std::endl;

    DynamicVector<complex<double>,columnVector> w1(N); // The vector for the real eigenvalues
    DynamicMatrix<complex<double>,rowMajor> V1(N,N); // The matrix for the left eigenvectors

    std::size_t n = N-1;
    //DynamicMatrix<double> Q(N, (n + 1));
    //DynamicMatrix<double> h((n + 1), n);
    ArnoldiTag tag;
    auto res = solve(A,b,tag,n);
    //eigen(h,w1,V1);
    std::cout << res.second;
    std::cout<< "The eigenvalues of Matrix h is: " << w1 << std::endl;
    std::cout << "The eigenvectors of Matrix h is: " << V1 << std::endl;

    return 0;
}
