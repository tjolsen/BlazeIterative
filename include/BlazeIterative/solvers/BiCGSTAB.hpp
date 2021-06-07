// Copyright (c)   2017 Tyler Olsen
//                 2018 Patrick Diehl
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BLAZE_ITERATIVE_BICGSTAB_HPP
#define BLAZE_ITERATIVE_BICGSTAB_HPP

#include "BiCGSTABTag.hpp"

BLAZE_NAMESPACE_OPEN
ITERATIVE_NAMESPACE_OPEN

namespace detail {

/**
 *  Implementation of the BiCGSTAB method, following the
 *  non-preconditioned version on Wikipedia
 *  (which is from Saad 2003).
 */
template<typename MatrixType, typename T>
void solve_impl(
        DynamicVector<T> &x,
        const MatrixType &A,
        const DynamicVector<T> &b,
        BiCGSTABTag &tag,
        std::string Preconditioner="")
{

    DynamicVector<T> r = b - A * x;
    DynamicVector<T> p(r);
    DynamicVector<T> v(r);
    DynamicVector<T> r0(r);
    DynamicVector<T> s(p.size());
    DynamicVector<T> t(p.size());
    DynamicVector<T> error(r);

    auto absolute_residual_0 = trans(r) * r;
    auto absolute_residual = absolute_residual_0;

    auto rho_prev = T(1);
    auto w = T(1);
    auto alpha = T(1);


    std::size_t iteration{0};
    while (true) {

        auto rho = trans(r0) * r;
        auto beta = (rho * alpha) / (rho_prev * w);

        p = r + beta * (p - w * v);
        v = A * p;
        alpha = rho / (trans(r0) * v);

        s = r - alpha * v;
        t = A*s;

        // sometimes, t will be zero, so trans(t)*t is zero.
        // This happens if the solution is exactly correct,
        // So best to set w=0, and loop will terminate below.
        auto t_dot_t = (trans(t)*t);
        if(t_dot_t == 0)
            w = 0;
        else
            w = (trans(t)*s)/t_dot_t;


        x += alpha*p + w*s;

        error = b - A*x;
        absolute_residual = trans(error)*error;
        auto relative_residual = absolute_residual/absolute_residual_0;
        if(tag.do_log()) {
            tag.log_residual(relative_residual);
        }

        if(tag.terminateIteration(iteration,absolute_residual,relative_residual)) {
            break;
        }

        r = s - w*t;
        rho_prev = rho;

        ++iteration;
    }

}

} //end namespace detail

ITERATIVE_NAMESPACE_CLOSE
BLAZE_NAMESPACE_CLOSE

#endif //BLAZE_ITERATIVE_BICGSTAB_HPP
