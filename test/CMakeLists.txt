cmake_minimum_required(VERSION 2.8)

add_executable(testprog test.c)

# headers
include_directories(${LIBLSODA_INCLUDE_DIR})

target_link_libraries(testprog lsoda)
