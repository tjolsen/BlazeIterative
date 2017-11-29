cmake_minimum_required(VERSION 3.5)
project(BlazeIterative)
set(CMAKE_CXX_STANDARD 14)

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
  )

# target_sources not necessary for header-only library, but makes it more
# IDE-friendly. Need to have it trigger only in "build" mode (not install)
# in order to avoid an error triggering (cmake v 3.9)
target_sources(BlazeIterative INTERFACE $<BUILD_INTERFACE:${SOURCE_FILES}>)
target_include_directories(BlazeIterative INTERFACE
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include/BlazeIterative>
  )

#==========================================
# Import Blaze and "link" to BlazeIterative
#==========================================
find_package(blaze REQUIRED)
target_link_libraries(BlazeIterative INTERFACE blaze::blaze)

#==========================================
# Install library
#==========================================
install(TARGETS BlazeIterative EXPORT BlazeIterativeConfig)
install(
  DIRECTORY include/
  DESTINATION include/BlazeIterative
  FILES_MATCHING PATTERN "*.hpp")
install(EXPORT BlazeIterativeConfig DESTINATION share/BlazeIterative/cmake
  NAMESPACE BlazeIterative::)


#==========================================
# Optionally build tests (primarily for devs)
#==========================================
option(BUILD_TESTS "Build BlazeIterative tests" OFF)
if(BUILD_TESTS)
    add_subdirectory(tests)
endif()

