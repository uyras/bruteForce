cmake_minimum_required(VERSION 3.0)

project(bruteForce)

find_package(OpenMP REQUIRED)

include_directories("../partsEngine")
link_directories("../partsEngine/build")

add_executable(bruteForce_omp main_omp.cpp)
target_link_libraries(bruteForce_omp partsEngine OpenMP::OpenMP_CXX)

add_executable(bruteForce_seq main_seq.cpp)
target_link_libraries(bruteForce_seq partsEngine)