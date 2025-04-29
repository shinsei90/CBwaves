cmake_minimum_required(VERSION 3.5)

project(CBwaves LANGUAGES CXX)

find_package(Threads REQUIRED)
find_package(Boost 1.83 COMPONENTS program_options REQUIRED)
# find_package(Boost REQUIRED)
find_package(Perl)

add_executable(${PROJECT_NAME} src/cbwaves.cxx
                               src/ConfigFile.cxx)

target_include_directories(${PROJECT_NAME} PRIVATE include
                                                   ${PROJECT_SOURCE_DIR})

target_link_libraries(${PROJECT_NAME} PRIVATE $<$<BOOL:${WIN32}>:Boost::disable_autolinking> # Boost can auto-link when used with MSVC, but relies on relative paths which CMake does not really like
                                              Threads::Threads
                                              Boost::boost)

target_compile_definitions(${PROJECT_NAME} PRIVATE $<$<BOOL:${WIN32}>:_CRT_SECURE_NO_WARNINGS>  # warning C4996: 'strcpy': This function or variable may be unsafe. Consider using strcpy_s instead.
                                                   $<$<BOOL:${WIN32}>:_SCL_SECURE_NO_WARNINGS>) # warning C4996: 'std::copy::_Unchecked_iterators::_Deprecate': Call to 'std::copy' with parameters that may be unsafe

set_target_properties(${PROJECT_NAME} PROPERTIES CXX_STANDARD 14
                                                 CXX_STANDARD_REQUIRED ON)

# if(${PERL_FOUND})

#   set(Imax 4)
#   set(Jmax 1)
#   set(script ${PROJECT_SOURCE_DIR}/examples/cbwgen.pl)
#   set(createdFiles)

#   foreach(I RANGE ${Imax})
#     foreach(J RANGE ${Jmax})
# 	  list(APPEND createdFiles ${CMAKE_CURRENT_BINARY_DIR}/test_cbwaves_${I}_${J}.ini
# 	                           ${CMAKE_CURRENT_BINARY_DIR}/test_cbwaves_${I}_${J}.des)
#     endforeach(J)
#   endforeach(I)

#   add_custom_command(COMMAND ${PERL_EXECUTABLE} ARGS ${script} ${Imax} ${Jmax}
#                      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
# 					 OUTPUT ${createdFiles}
# 					 DEPENDS ${script}
# 					 MAIN_DEPENDENCY ${script}
# 					 COMMENT "Running script examples/cbwgen.pl")

#   add_custom_target(ini ALL
#                     DEPENDS ${createdFiles}
# 					SOURCES ${script}
# 					COMMENT "Generating ini files")
# endif(${PERL_FOUND})