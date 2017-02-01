include_directories(
    ${PROJECT_SOURCE_DIR}/src
    ${PROJECT_SOURCE_DIR}/src/density
    ${PROJECT_BINARY_DIR}/external/include
    ${PROJECT_SOURCE_DIR}/api
    )

add_library(
    xcint
    src/rolex.cpp
    src/xcint_blas.cpp
    src/Functional.cpp
    src/integrator.cpp
    src/density/Basis.cpp
    src/density/AOBatch.cpp
    )

set_target_properties(xcint PROPERTIES COMPILE_FLAGS "-fPIC")

add_library(
    xcint_shared
    SHARED
    src/empty.cpp
    )

if(APPLE)
    target_link_libraries(
        xcint_shared
        ${MATH_LIBS}
        "-Wl, -force_load"
        xcint
        ${PROJECT_BINARY_DIR}/external/lib/libxcfun.a
        )
else()
    target_link_libraries(
        xcint_shared
        ${MATH_LIBS}
        "-Wl,--whole-archive"
        xcint
        ${PROJECT_BINARY_DIR}/external/lib/libxcfun.a
        "-Wl,--no-whole-archive"
        )
endif()

add_dependencies(xcint xcfun)
add_dependencies(xcint numgrid)
add_dependencies(xcint balboa)

if(ENABLE_FC_SUPPORT)
    add_library(
        xcint_fortran
        ${CMAKE_CURRENT_LIST_DIR}/../../api/xcint.f90
        )

    target_link_libraries(
        xcint_fortran
        xcint_shared
        )
endif()

install(TARGETS xcint ARCHIVE DESTINATION lib)
install(TARGETS xcint_shared LIBRARY DESTINATION lib)
