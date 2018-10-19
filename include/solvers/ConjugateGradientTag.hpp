// Copyright (c)   2017 Tyler Olsen
//                 2018 Patrick Diehl
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BLAZE_ITERATIVE_CONJUGATEGRADIENTTAG_HPP
#define BLAZE_ITERATIVE_CONJUGATEGRADIENTTAG_HPP

#include "IterativeCommon.hpp"
#include "IterativeTag.hpp"

BLAZE_NAMESPACE_OPEN
ITERATIVE_NAMESPACE_OPEN

/**
 * \class ConjugateGradientTag
 * \brief Tag type to dispatch a Conjugate Gradient solver
 *
 * Must be used for Symmetric positive-definite (SPD) linear systems only.
 * Select a different method for non-symmetric systems.
 * If compiled without NDEBUG defined, a runtime check
 * will trigger to check that the matrix is symmetric.
 * However, no checks will be performed to ensure that
 * the matrix is positive-definite.
 */
class ConjugateGradientTag : public IterativeTag
{
public:
    ConjugateGradientTag() {
        solverName = "Conjugate Gradient";
    }
};

ITERATIVE_NAMESPACE_CLOSE
BLAZE_NAMESPACE_CLOSE

#endif //BLAZE_ITERATIVE_CONJUGATEGRADIENTTAG_HPP
