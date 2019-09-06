add_executable(test_arnoldi main_Arnoldi.cpp)
target_link_libraries(test_arnoldi PRIVATE BlazeIterative)

add_executable(test_lanczos main_Lanczos.cpp)
target_link_libraries(test_lanczos PRIVATE BlazeIterative)

add_executable(test_preconditionedcg main_PreconditionedCG.cpp)
target_link_libraries(test_preconditionedcg PRIVATE BlazeIterative)

add_executable(test_gmres main_GMRES.cpp)
target_link_libraries(test_gmres PRIVATE BlazeIterative)

add_executable(test_cg main_CG.cpp)
target_link_libraries(test_cg PRIVATE BlazeIterative)


