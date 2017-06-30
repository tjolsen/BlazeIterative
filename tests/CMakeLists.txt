include_directories(${CMAKE_SOURCE_DIR}/blaze_iterative)

add_executable(blz_iter_test main.cpp)
target_link_libraries(blz_iter_test BlazeIterative)