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

    // Test CG

    std::size_t N = 10;
    DynamicMatrix<double,false> A(N,N, 0.0);
    DynamicVector<double> b(N, 0.0);
    DynamicVector<double> x1(N, 0.5);
    for(int i=0; i<N; ++i) {
        A(i,i) = 2.0;
        b[i] = 1.0*(1+i);
        x1[i] += x1[i-1];
    }

    ConjugateGradientTag tag;
    tag.do_log() = true;
    auto x2 = solve(A,b,tag);

    auto error = norm(x1 - x2);




    if (error < EPSILON){
        std::cout << " Pass test of CG" << std::endl;
        return EXIT_SUCCESS;
    } else{
        std::cout << "Fail test of CG" << std::endl;
        return EXIT_FAILURE;
    }
}

