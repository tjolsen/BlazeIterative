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

    // Test GMRES

    std::size_t N = 5;
    DynamicMatrix<double,false> A(N,N);
    DynamicVector<double> b(N);
    A = {{0.46, 0.60, 0.74, 0.61, 0.85}, {0.56, 0.31, 0.80, 0.94, 0.76},{0.41, 0.19, 0.15, 0.33, 0.06},{0.03, 0.92, 0.15, 0.56, 0.08},{0.09, 0.06, 0.69, 0.42, 0.96}};
    b = {1.788, 1.891, 0.458, 0.818, 1.53};
    DynamicVector<double> x1{0.1, 0.3, 0.5, 0.7, 0.9};

    GMRESTag tag;
    tag.do_log() = true;
    std::size_t n = N;
    auto x2 = solve(A,b,tag,n);

    auto error = norm(x1 - x2);


    if (error < EPSILON){
        std::cout << " Pass test of GMRES" << std::endl;
        return EXIT_SUCCESS;
    } else{
        std::cout << "Fail test of GMRES" << std::endl;
        return EXIT_FAILURE;
    }

}

