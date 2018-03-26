cmake_minimum_required( VERSION 3.0.0 )

project( cbwaves-2.0.0 LANGUAGES CXX )

option( BUILD_TESTS "Build unit tests" ON )

add_executable( ${PROJECT_NAME} src/Main.cpp
                                src/CorrectionsAccTerms.cpp
                                src/CorrectionsHTerms.cpp
                                src/CorrectionsHamiltonian.cpp
                                src/DynamicalParams.cpp
                                include/DynamicalParams.hpp
                                include/InitParams.hpp
                                include/Config.hpp
                                include/Corrections.hpp
                                include/PDE.hpp
                                )

target_include_directories( ${PROJECT_NAME} PRIVATE include )

set_target_properties( ${PROJECT_NAME} PROPERTIES CXX_STANDARD 14
                                                  CXX_STANDARD_REQUIRED ON )

if( BUILD_TESTS )
  add_executable( PDE_Test src/PDE_Test.cpp )

  target_include_directories( PDE_Test PRIVATE include )

  set_target_properties( PDE_Test PROPERTIES CXX_STANDARD 14
                                             CXX_STANDARD_REQUIRED ON )
                                             
  find_package( GnuPlot REQUIRED )
  
  # Target that generates data file
  add_custom_command( COMMAND PDE_Test
                      WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
                      OUTPUT ${CMAKE_BINARY_DIR}/lv.dat
                      BYPRODUCTS ${CMAKE_BINARY_DIR}/lv.dat
                      DEPENDS PDE_Test
                      COMMENT "Generating data set" )
  
  add_custom_target( data DEPENDS ${CMAKE_BINARY_DIR}/lv.dat )
  
  # Target that generates data file
  add_custom_command(COMMAND ${GNUPLOT_EXECUTABLE} ${PROJECT_SOURCE_DIR}/script/lv.plt
                     WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
                     OUTPUT ${CMAKE_BINARY_DIR}/lv.svg
                     BYPRODUCTS ${CMAKE_BINARY_DIR}/lv.svg
                     DEPENDS ${PROJECT_SOURCE_DIR}/script/lv.plt data
                     COMMENT "Generating plot")
  
  add_custom_target(plot DEPENDS ${CMAKE_BINARY_DIR}/lv.svg)
endif( BUILD_TESTS )