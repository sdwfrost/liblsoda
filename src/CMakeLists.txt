cmake_minimum_required(VERSION 2.8)

set(LIBLSODA_SOURCES
  cfode.c
  cfode_static.c
  common.c
  corfailure.c
  correction.c
  daxpy.c
  ddot.c
  dgefa.c
  dgesl.c
  dscal.c
  fnorm.c
  idamax.c
  intdy.c
  lsoda.c
  methodswitch.c
  orderswitch.c
  prja.c
  scaleh.c
  solsy.c
  stoda_internal.c
  vmnorm.c
  strdup_printf.c
  )

add_library(lsoda-static STATIC ${LIBLSODA_SOURCES})
add_library(lsoda SHARED ${LIBLSODA_SOURCES})
add_library(lsoda-pic-object OBJECT ${LIBLSODA_SOURCES})
set_property(TARGET lsoda-pic-object PROPERTY POSITION_INDEPENDENT_CODE ON)

if(UNIX)
  target_link_libraries(lsoda m)
endif()
