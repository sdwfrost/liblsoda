cmake_minimum_required(VERSION 2.8)

# variables
set(CMAKE_MACOSX_RPATH 1)
set(LIBLSODA_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/src)

add_subdirectory(src)
add_subdirectory(test)
