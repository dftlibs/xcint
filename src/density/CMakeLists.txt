cmake_policy(SET CMP0074 NEW)

add_library(density STATIC "")

target_sources(
  density
  PRIVATE
    blas_interface.cpp
    blas_interface.h
    density.cpp
    compress.cpp
    compress.h
  PUBLIC
    ${CMAKE_CURRENT_LIST_DIR}/density.h
  )

find_package(BLAS REQUIRED)
target_link_libraries(
  density
  PUBLIC
    ${BLAS_LIBRARIES}
  )

target_include_directories(
  density
  PRIVATE
    ${PROJECT_BINARY_DIR}/generated
  PUBLIC
    ${CMAKE_CURRENT_LIST_DIR}
  )

set_target_properties(
  density
  PROPERTIES
    COMPILE_FLAGS "-fPIC"
  )

install(TARGETS density ARCHIVE DESTINATION lib)
