//
// Created by tyler on 6/29/17.
//

#ifndef BLAZE_ITERATIVE_BICGSTABLTAG_HPP
#define BLAZE_ITERATIVE_BICGSTABLTAG_HPP

#include "IterativeCommon.hpp"

BLAZE_NAMESPACE_OPEN
ITERATIVE_NAMESPACE_OPEN

/**
 * \class BiCGSTABTag
 * \brief Tag type to dispatch a BiCGSTAB(l) solver
 *
 * BiCGSTAB(l) is a generalized version of the BiCGSTAB
 * method with improved convergence properties. It is suitable
 * for use with non-symmetric linear systems. This is a
 * non-preconditioned version of the algorithm.
 * In short, the algorithm performs "L" BiCGSTAB steps
 * followed by a GMRES step for each iteration of the
 * BiCGSTAB(l) method.
 *
 * Based on the presentation in Sleijpen and Fokkema (1993).
 *
 */
class BiCGSTABTag : public IterativeTag
{
public:
    BiCGSTABTag() {
        solverName = "BiCGSTAB";
    }
};


ITERATIVE_NAMESPACE_CLOSE
BLAZE_NAMESPACE_CLOSE

#endif //BLAZE_ITERATIVE_BICGSTABLTAG_HPP
