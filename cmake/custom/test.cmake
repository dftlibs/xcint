include(CTest)
enable_testing()

add_executable(
    unit_tests
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
    unit_tests
    googletest
    xcint
    ${PROJECT_BINARY_DIR}/external/xcfun-build/libxcfun.a
    ${_numgrid_lib}
    ${PROJECT_BINARY_DIR}/external/numgrid-build/lib/libnumgrid.${_dyn_lib_suffix}
    ${MATH_LIBS}
    pthread
    )

add_test(unit_tests ${PROJECT_BINARY_DIR}/bin/unit_tests ${PROJECT_SOURCE_DIR}/test)
