// Copyright (c)   2017 Tyler Olsen
//                 2018 Patrick Diehl
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BLAZE_ITERATIVE_PRECONDITIONBICGSTABTAG_HPP
#define BLAZE_ITERATIVE_PRECONDITIONBICGSTABTAG_HPP

#include "IterativeCommon.hpp"
#include "IterativeTag.hpp"

BLAZE_NAMESPACE_OPEN
ITERATIVE_NAMESPACE_OPEN

/**
 * \class PreconditionBiCGSTABTag
 * \brief Tag type to dispatch a preconditioned BiCGSTAB   solver
 *
 * BiCGSTAB is a modified version of the biconjugate gradient
 * method with improved convergence properties. It is suitable
 * for use with non-symmetric linear systems. This is a
 * non-preconditioned version of the algorithm.
 *
 */
class PreconditionBiCGSTABTag : public IterativeTag
{
public:
    PreconditionBiCGSTABTag() {
        solverName = "PreconditionBiCGSTAB";
    }
};

ITERATIVE_NAMESPACE_CLOSE
BLAZE_NAMESPACE_CLOSE


#endif //BLAZE_ITERATIVE_PRECONDITIONBICGSTABTAG_HPP
