// Copyright (c)   2017 Tyler Olsen
//                 2018 Patrick Diehl
//                 2019 Nanmiao Wu
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "BlazeIterative.hpp"
#include <iostream>
#include <cstdlib>

using namespace blaze;
using namespace blaze::iterative;

int main() {

    // Test Arnoldi

    std::size_t N = 6;
    DynamicMatrix<double,columnMajor> A(N,N,0);
    band<0>(A) = {0, 1, 2, 3, 4, 100000};

    DynamicVector<double> b(N, 1);
    DynamicVector<complex<double>,columnVector> w(N);
    DynamicMatrix<complex<double>,rowMajor> V(N,N);
    DynamicVector<complex<double>,columnVector> w1(N);
    eigen(A,w,V);

    w1 = {w[5],w[0],w[4],w[1],w[2],w[3]};

    std::size_t n = N;
    ArnoldiTag tag;
    DynamicVector<complex<double>,columnVector> w2(n);
    DynamicMatrix<complex<double>,rowMajor> V2(n,n);
    w2 = solve(A,b,tag,n);


    auto error = real(norm(w1 - w2));


    if (error < EPSILON){

        std::cout << " Pass test of Arnoldi" << std::endl;
        return EXIT_SUCCESS;
    } else{
        std::cout << "Fail test of Arnoldi" << std::endl;
        return EXIT_FAILURE;
    }

}
