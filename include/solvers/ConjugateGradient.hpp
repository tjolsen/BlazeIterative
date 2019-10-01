// Copyright (c)   2017 Tyler Olsen
//                 2018 Patrick Diehl
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BLAZE_ITERATIVE_CONJUGATEGRADIENT_HPP
#define BLAZE_ITERATIVE_CONJUGATEGRADIENT_HPP

#include "IterativeCommon.hpp"
#include "ConjugateGradientTag.hpp"


BLAZE_NAMESPACE_OPEN
ITERATIVE_NAMESPACE_OPEN

namespace detail {

template<typename MatrixType, typename T>
void solve_impl(
        DynamicVector<T> &x,
        const MatrixType &A,
        const DynamicVector<T> &b,
        ConjugateGradientTag &tag,
        std::string Preconditioner="")
{

    BLAZE_INTERNAL_ASSERT(isSymmetric(A), "A must be a symmetric matrix")

    DynamicVector<T> r = b - A*x;
    DynamicVector<T> p(r);
    DynamicVector<T> Ap(p.size());

    auto absolute_residual_0 = trans(r)*r;
    auto absolute_residual = absolute_residual_0;
    auto absolute_residual_prev = absolute_residual;

    if(tag.do_log()) {
        tag.log_residual(absolute_residual/absolute_residual_0);
    }


    std::size_t iteration{0};
    while(true) {
        absolute_residual_prev = absolute_residual;
        Ap = declsym(A)*p;

        auto alpha = absolute_residual/(trans(p)*Ap);
        x += alpha*p;
        r -= alpha*Ap;


        absolute_residual = trans(r)*r;

        if(tag.do_log()) {
            tag.log_residual(absolute_residual/absolute_residual_0);
        }

        if(tag.terminateIteration(iteration, absolute_residual, absolute_residual/absolute_residual_0)) {
            break;
        }

        auto beta = absolute_residual/absolute_residual_prev;
        p = r + beta*p;

        ++iteration;
    }//end while


};


} //end namespace detail

ITERATIVE_NAMESPACE_CLOSE
BLAZE_NAMESPACE_CLOSE

#endif //BLAZE_ITERATIVE_CONJUGATEGRADIENT_HPP
