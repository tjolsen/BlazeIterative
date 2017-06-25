//
// Created by tyler on 6/24/17.
//

#ifndef BLAZE_ITERATIVE_ITERATIVE_TAG_HPP
#define BLAZE_ITERATIVE_ITERATIVE_TAG_HPP

#include "TerminationStatus.hpp"


BLAZE_NAMESPACE_OPEN
ITERATIVE_NAMESPACE_OPEN

class IterativeTag
{

public:
    IterativeTag() {}

    inline bool terminateIteration(int iteration, double absolute_residual, double relative_residual)
    {
        if (isConverged(absolute_residual, relative_residual)) {
            return true;
        } else if (iteration >= maximum_iterations) {
            terminationStatus = TerminationStatus::ITERATION_LIMIT;
            return true;
        } else {
            return false;
        }
    }

    template<typename T, typename=typename std::enable_if<std::is_convertible<T, double>::value>::type>
    inline void log_residual(T residual) { convergence_history_container.push_back(residual); }

    bool &do_log() { return record_convergence_history; }

    bool do_log() const { return record_convergence_history; }

    inline TerminationStatus status() const { return terminationStatus; }

    double &relativeResidualTolerance() { return relative_residual_tolerance; }

    double relativeResidualTolerance() const { return relative_residual_tolerance; }

    double &absoluteResidualTolerance() { return absolute_residual_tolerance; }

    double absoluteResidualTolerance() const { return absolute_residual_tolerance; }

    std::size_t &maximumIterations() { return maximum_iterations; }

    std::size_t maximumIterations() const { return maximum_iterations; }

    const std::vector<double> &convergence_history() const { return convergence_history_container; }

protected:
    std::size_t maximum_iterations{20};
    double relative_residual_tolerance{1.0e-6};
    double absolute_residual_tolerance{0.0};
    std::string solverName{"Default"};
    TerminationStatus terminationStatus{TerminationStatus::NOT_TERMINATED};
    bool record_convergence_history{false};

    //container for relative residual convergence history
    std::vector<double> convergence_history_container;

    inline bool isConverged(double absolute_residual, double relative_residual)
    {
        if (std::abs(relative_residual) < relative_residual_tolerance) {
            terminationStatus = TerminationStatus::CONVERGED_RELATIVE_RESIDUAL;
            return true;
        } else if (std::abs(absolute_residual) < absolute_residual_tolerance) {
            terminationStatus = TerminationStatus::CONVERGED_ABSOLUTE_RESIDUAL;
            return true;
        } else {
            return false;
        }
    }

};


ITERATIVE_NAMESPACE_CLOSE
BLAZE_NAMESPACE_CLOSE

#endif //BLAZE_ITERATIVE_ITERATIVE_TAG_HPP
