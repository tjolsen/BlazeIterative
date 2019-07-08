// Copyright (c)   2017 Tyler Olsen
//                 2018 Patrick Diehl
//                 2019 Nanmiao Wu
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BLAZE_ITERATIVE_ARNOLDITAG_HPP
#define BLAZE_ITERATIVE_ARNOLDITAG_HPP

#include "IterativeCommon.hpp"
#include "IterativeTag.hpp"

BLAZE_NAMESPACE_OPEN
    ITERATIVE_NAMESPACE_OPEN

        class ArnoldiTag : public IterativeTag
        {
        public:
            ArnoldiTag() {
                solverName = "Arnoldi";
            }
        };

    ITERATIVE_NAMESPACE_CLOSE
BLAZE_NAMESPACE_CLOSE

#endif //BLAZE_ITERATIVE_ARNOLDITAG_HPP
