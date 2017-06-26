//
// Created by tyler on 6/24/17.
//

#include "BlazeIterative.hpp"
#include <iostream>

using namespace blaze;
using namespace blaze::iterative;

int main() {

    std::size_t N = 10;
    DynamicMatrix<double,columnMajor> A(N,N, 0.0);
    DynamicVector<double> b(N, 0.0);

    for(int i=0; i<N; ++i) {
        A(i,i) = 2.0;
        b[i] = 1.0*(1+i);
    }

    BiCGSTABTag tag;
    tag.do_log() = true;

    std::cout << solve(A,b,tag) << std::endl << std::endl;

    int iter(0);
    for(auto r : tag.convergence_history()) {
        std::cout << iter++ << '\t' << r << '\n';
    }

    return 0;
}