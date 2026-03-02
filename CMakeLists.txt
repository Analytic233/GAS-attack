cmake_minimum_required(VERSION 3.16)
project(gas_attack LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(MSVC)
  add_compile_options(/O2)
else()
  add_compile_options(-O3)
endif()

find_package(OpenMP)

add_executable(FiLIP_verifier FiLIP_verify.cpp)
add_executable(toy_attack toy.cpp)
add_executable(output_generator output_make.cpp prng.cpp)
add_executable(attack_program attack.cpp prng.cpp verify.cpp)

if(OpenMP_CXX_FOUND)
  target_link_libraries(attack_program PRIVATE OpenMP::OpenMP_CXX)
endif()
