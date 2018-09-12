set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(AO_BLOCK_LENGTH "128" CACHE STRING "AO block length")
configure_file(
  ${PROJECT_SOURCE_DIR}/src/parameters.h.in
  ${PROJECT_BINARY_DIR}/generated/generated_parameters.h
  @ONLY
  )

add_subdirectory(density)

add_library(
  xcint
  STATIC
    Functional.cpp
    Functional.h
    integrator.cpp
    integrator.h
    xcint_parameters.h
  )

target_link_libraries(
  xcint
  PRIVATE
    density
    ${PROJECT_BINARY_DIR}/balboa/src/balboa-build/lib/libbalboa.so
  )

target_include_directories(
  xcint
  PRIVATE
    ${PROJECT_BINARY_DIR}/external/include
    ${PROJECT_BINARY_DIR}/generated
    ${PROJECT_BINARY_DIR}/balboa/src/balboa/balboa/
  )

add_dependencies(xcint xcfun)
add_dependencies(xcint numgrid)
add_dependencies(xcint balboa)

if(ENABLE_FC_SUPPORT)
  add_library(
    xcint_fortran
    ${PROJECT_SOURCE_DIR}/api/xcint.f90
    )
endif()

install(TARGETS xcint ARCHIVE DESTINATION lib)
