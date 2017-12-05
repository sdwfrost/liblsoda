# - Try to find libSBML
# Once done this will define
#  LIBSBML_FOUND - System has libSBML
#  LIBSBML_INCLUDE_DIRS - The libSBML include directories
#  LIBSBML_LIBRARIES - The libraries needed to use libSBML
#  LIBSBML_DEFINITIONS - Compiler switches required for using libSBML

cmake_minimum_required(VERSION 2.8)

find_package(PkgConfig)
pkg_check_modules(PC_LIBSBML QUIET libsbml)
set(LIBSBML_DEFINITIONS ${PC_LIBSBML_CFLAGS_OTHER})

find_path(LIBSBML_INCLUDE_DIR
  NAMES sbml/SBMLTypes.h
  PATHS ENV C_INCLUDE_PATH
        /usr/include /usr/local/include /opt/local/include
        ${CMAKE_SOURCE_DIR}/depend/include
  HINTS ${PC_LIBSBML_INCLUDEDIR} ${PC_LIBSBML_INCLUDE_DIRS} )

find_library(LIBSBML_LIBRARY
  NAMES libsbml.lib sbml
  PATHS ENV DYLD_LIBRARY_PATH
        ENV LD_LIBRARY_PATH
        /usr/lib64 /usr/lib /usr/local/lib /opt/local/lib
        ${CMAKE_SOURCE_DIR}/depend/lib
  HINTS ${PC_LIBSBML_LIBDIR} ${PC_LIBSBML_LIBRARY_DIRS} )

find_library(LIBBZIP2_LIBRARY
  NAMES bzip2.lib bz2 libbz2.lib
  PATHS ENV DYLD_LIBRARY_PATH
        ENV LD_LIBRARY_PATH
        /usr/lib /usr/local/lib /opt/local/lib
        ${CMAKE_SOURCE_DIR}/depend/lib
  HINTS ${PC_LIBBZIP2} ${PC_LIBBZIP2_LIBRARY_DIRS} )

set(LIBSBML_LIBRARIES ${LIBSBML_LIBRARY} ${LIBBZIP2_LIBRARY} )
set(LIBSBML_INCLUDE_DIRS ${LIBSBML_INCLUDE_DIR} )

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set LIBSBML_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(LIBSBML DEFAULT_MSG LIBSBML_LIBRARY LIBSBML_INCLUDE_DIR)
