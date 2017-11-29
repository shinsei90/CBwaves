cmake_minimum_required(VERSION 3.0.0)

project(RK4_Test LANGUAGES CXX)

add_executable(${PROJECT_NAME} include/PDE.hpp
                               include/Corrections.hpp
                               src/Main.cpp
                               src/CorrectionsAccTerms.cpp
                               src/CorrectionsHTerms.cpp
                               src/CorrectionsHamiltonian.cpp)

target_include_directories(${PROJECT_NAME} PRIVATE include)

set_target_properties(${PROJECT_NAME} PROPERTIES CXX_STANDARD 14
                                                 CXX_STANDARD_REQUIRED ON)