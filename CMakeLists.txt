cmake_minimum_required(VERSION 3.13)

project(bruteForce)

set( ARGUMENTUM_BUILD_STATIC_LIBS ON INTERNAL )
find_package(OpenMP REQUIRED)
add_subdirectory(argumentum)
add_subdirectory(partsEngine)

set( CMAKE_CXX_STANDARD 17 )

include_directories("partsEngine" "argumentum/include")

add_executable(bruteForce_omp main_omp.cpp common.cpp)
target_link_libraries(bruteForce_omp PRIVATE partsEngine argumentum OpenMP::OpenMP_CXX)

add_executable(bruteForce_seq main_seq.cpp)
target_link_libraries(bruteForce_seq partsEngine)