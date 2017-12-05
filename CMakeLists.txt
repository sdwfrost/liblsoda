cmake_minimum_required(VERSION 2.8)
message(STATUS "Running cmake version ${CMAKE_VERSION}")

# project information
project(liblsoda C CXX)

# cmake modules
include(CMakeCXXInformation)
include(CheckCXXCompilerFlag)
include(cmake/find_libsbml.cmake)

# build type
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING
    "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
    FORCE)
endif(NOT CMAKE_BUILD_TYPE)

# check C++11 availability
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
if(COMPILER_SUPPORTS_CXX11)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
else()
  message(FATAL_ERROR "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support.")
endif()

# variables
set(CMAKE_MACOSX_RPATH 1)
set(LIBLSODA_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/src)

add_subdirectory(src)
add_subdirectory(test)
add_subdirectory(test-sbml)
