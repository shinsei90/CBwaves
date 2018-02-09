cmake_minimum_required(VERSION 3.0.0)

project(cbwaves-2.0.0 LANGUAGES CXX)

add_executable( ${PROJECT_NAME} src/Main.cpp
                                src/CorrectionsAccTerms.cpp
                                src/CorrectionsHTerms.cpp
                                src/CorrectionsHamiltonian.cpp
                                include/Corrections.hpp
                                include/PDE.hpp
                                )

target_include_directories( ${PROJECT_NAME} PRIVATE include )

set_target_properties( ${PROJECT_NAME} PROPERTIES CXX_STANDARD 14
                                                  CXX_STANDARD_REQUIRED ON )