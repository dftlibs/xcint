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
    )

target_include_directories(
    xcint
    PRIVATE
       ${PROJECT_BINARY_DIR}/external/include
       ${PROJECT_BINARY_DIR}/generated
    )

set_target_properties(xcint PROPERTIES COMPILE_FLAGS "-fPIC")

add_library(
    xcint_shared
    SHARED
    empty.cpp
    )

if(APPLE)
    target_link_libraries(
        xcint_shared
        ${MATH_LIBS}
        "-Wl, -force_load"
        xcint
        density
        ${PROJECT_BINARY_DIR}/external/lib/libxcfun.a
        )
else()
    target_link_libraries(
        xcint_shared
        ${MATH_LIBS}
        "-Wl,--whole-archive"
        xcint
        density
        ${PROJECT_BINARY_DIR}/external/lib/libxcfun.a
        "-Wl,--no-whole-archive"
        )
endif()

add_dependencies(xcint xcfun)
add_dependencies(xcint numgrid)

if(ENABLE_FC_SUPPORT)
    add_library(
        xcint_fortran
        ${PROJECT_SOURCE_DIR}/api/xcint.f90
        )

    target_link_libraries(
        xcint_fortran
        xcint_shared
        )
endif()

install(TARGETS xcint ARCHIVE DESTINATION lib)
install(TARGETS xcint_shared LIBRARY DESTINATION lib)
