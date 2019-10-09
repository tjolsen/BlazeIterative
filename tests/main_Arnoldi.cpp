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
    DynamicMatrix<double,rowMajor> A(N, N, 0.0);
    band<0>(A) = {0, 1, 2, 3, 4, 100000};
    DynamicVector<double> b(N, 1);
    DynamicVector<complex<double>,columnVector> w(N);
    DynamicVector<double,columnVector> w1(N);
    DynamicVector<double,columnVector> w2(N);
    eigen(A,w);
    w1 = real(w);
    w2 = {w1[5],w1[0],w1[4],w1[1],w1[2],w1[3]};

    std::size_t n = N;
    ArnoldiTag tag;
    DynamicVector<double,columnVector> w3(n);
    w3 = solve(A,b,tag,n);
    
    auto error = norm(w2 - w3);
    
    if (error < EPSILON){
        std::cout << " Pass test of Arnoldi" << std::endl;
        return EXIT_SUCCESS;
    } else{
        std::cout << "Fail test of Arnoldi" << std::endl;
        return EXIT_FAILURE;
    }

}
