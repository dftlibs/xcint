include(CTest)
enable_testing()

add_executable(
    unit_tests
    test/main.cpp
    test/energy.cpp
    )

target_link_libraries(
    unit_tests
    gtest
    xcint
    ${PROJECT_BINARY_DIR}/external/xcfun-build/libxcfun.a
    ${PROJECT_BINARY_DIR}/external/numgrid-build/lib/libnumgrid.so
    ${MATH_LIBS}
    pthread
    )

add_test(unit_tests ${PROJECT_BINARY_DIR}/bin/unit_tests)
