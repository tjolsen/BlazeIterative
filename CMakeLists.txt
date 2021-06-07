# Copyright (c) 2017 Tyler Olsen
# Copyright (c) 2018-2019 Patrick Diehl
# Copyright (c) 2019 Nanmiao Wu
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

cmake_minimum_required(VERSION 3.5)

project(BlazeIterative)

set(CMAKE_CXX_STANDARD 14)

#==========================================
# Optionally build tests (For devs)
#==========================================
option(BUILD_TESTS "Build BlazeIterative tests" OFF)

#==========================================
# Set up BlazeIterative target
#==========================================
add_library(BlazeIterative INTERFACE)

#==========================================
# Set BlazeIterative source files
#==========================================
set(SOURCE_FILES
  ${CMAKE_CURRENT_SOURCE_DIR}/include/IterativeCommon.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/include/IterativeTag.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/include/solve.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/include/solvers/ConjugateGradientTag.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/include/solvers/ConjugateGradient.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/include/solvers/BiCGSTABTag.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/include/solvers/BiCGSTAB.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/include/solvers/PreconditionBiCGSTAB.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/include/solvers/PreconditionBiCGSTABTag.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/include/solvers/Arnoldi.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/include/solvers/ArnoldiTag.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/include/solvers/PreconditionCG.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/include/solvers/PreconditionCGTag.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/include/solvers/Lanczos.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/include/solvers/LanczosTag.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/include/solvers/GMRES.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/include/solvers/GMRESTag.hpp
)


# target_sources not necessary for header-only library, but makes it more
# IDE-friendly. Need to have it trigger only in "build" mode (not install)
# in order to avoid an error triggering (cmake v 3.9)
target_sources(BlazeIterative INTERFACE $<BUILD_INTERFACE:${SOURCE_FILES}>)
target_include_directories(BlazeIterative INTERFACE
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include/>
  $<INSTALL_INTERFACE:include/>
)

#==========================================
# Import Blaze and "link" to BlazeIterative
#==========================================
find_package(blaze REQUIRED)
target_link_libraries(BlazeIterative INTERFACE blaze::blaze)
find_package(LAPACK REQUIRED)
target_link_libraries(BlazeIterative INTERFACE  ${LAPACK_LIBRARIES})

#==========================================
# Install library
#==========================================
install(TARGETS BlazeIterative EXPORT BlazeIterativeConfig)
install(
  DIRECTORY include/
  DESTINATION include/
  FILES_MATCHING PATTERN "*.hpp")
install(EXPORT BlazeIterativeConfig DESTINATION share/BlazeIterative/cmake
  NAMESPACE BlazeIterative::)

if(BUILD_TESTS)
  ENABLE_TESTING()
  add_subdirectory(tests)
endif()
