include(CTest)
enable_testing()

add_executable(
    cpp_test
    test/main.cpp
    test/energy_spherical.cpp
#   test/energy_cartesian.cpp
    )

# workaround:
# different dynamic lib suffix on mac
if(APPLE)
    set(_dyn_lib_suffix dylib)
else()
    set(_dyn_lib_suffix so)
endif()

target_link_libraries(
    cpp_test
    googletest
    xcint
    ${PROJECT_BINARY_DIR}/external/xcfun-build/libxcfun.a
    ${PROJECT_BINARY_DIR}/external/numgrid-build/lib/libnumgrid.${_dyn_lib_suffix}
    ${MATH_LIBS}
    pthread
    )

add_test(cpp_test ${PROJECT_BINARY_DIR}/bin/cpp_test ${PROJECT_SOURCE_DIR}/test)

if(ENABLE_FC_SUPPORT)
    add_executable(
        fortran_test
        test/test.F90
        )

    target_link_libraries(
        fortran_test
        xcint_fortran
        ${PROJECT_BINARY_DIR}/external/lib/libnumgrid.so
        ${PROJECT_BINARY_DIR}/external/numgrid-build/src/libnumgrid_fortran.a
        )

    add_test(fortran_test ${PROJECT_BINARY_DIR}/bin/fortran_test ${PROJECT_SOURCE_DIR}/test)
endif()
