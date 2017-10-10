cmake_minimum_required(VERSION 3.7)
project(blaze_iterative)


find_package(blaze)

set(CMAKE_CXX_STANDARD 14)

set(SOURCE_FILES
        blaze_iterative/IterativeCommon.hpp
        blaze_iterative/IterativeTag.hpp
        blaze_iterative/solve.hpp
        blaze_iterative/solvers/ConjugateGradientTag.hpp
        blaze_iterative/solvers/ConjugateGradient.hpp
        blaze_iterative/solvers/BiCGSTABTag.hpp
        blaze_iterative/solvers/BiCGSTAB.hpp
        )

add_library(BlazeIterative INTERFACE)

add_subdirectory(tests)

install(DIRECTORY blaze_iterative DESTINATION ${CMAKE_INSTALL_PREFIX} 
	PATTERN "*.hpp") 
