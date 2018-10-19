// Copyright (c)   2017 Tyler Olsen
//                 2018 Patrick Diehl
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BLAZE_ITERATIVE_TERMINATIONSTATUS_HPP
#define BLAZE_ITERATIVE_TERMINATIONSTATUS_HPP

#include "IterativeCommon.hpp"

BLAZE_NAMESPACE_OPEN
ITERATIVE_NAMESPACE_OPEN

enum class TerminationStatus : unsigned char {
    NOT_TERMINATED,
    CONVERGED_RELATIVE_RESIDUAL,
    CONVERGED_ABSOLUTE_RESIDUAL,
    ITERATION_LIMIT
};

ITERATIVE_NAMESPACE_CLOSE
BLAZE_NAMESPACE_CLOSE

#endif //BLAZE_ITERATIVE_TERMINATIONSTATUS_HPP
