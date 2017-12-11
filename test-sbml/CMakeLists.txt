cmake_minimum_required(VERSION 2.8)

add_executable(testSBML testSBML.cpp)

# headers
include_directories(${LIBLSODA_INCLUDE_DIR})
include_directories(${LIBSBML_INCLUDE_DIR})

target_link_libraries(testSBML lsoda ${LIBSBML_LIBRARIES})
