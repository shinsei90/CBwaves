cmake_minimum_required(VERSION 3.0.0)

project(RK4_Test LANGUAGES CXX)

add_executable(${PROJECT_NAME} inc/PDE.hpp
                               src/Main.cpp)

target_include_directories(${PROJECT_NAME} PRIVATE inc)

set_target_properties(${PROJECT_NAME} PROPERTIES CXX_STANDARD 14
                                                 CXX_STANDARD_REQUIRED ON)