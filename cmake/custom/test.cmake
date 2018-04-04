option(ENABLE_UNIT_TESTS "Enable unit tests" ON)
message(STATUS "Enable testing: ${ENABLE_UNIT_TESTS}")

if(ENABLE_UNIT_TESTS)
    include(CTest)
    enable_testing()

    include(ExternalProject)

    ExternalProject_Add(
        gtest
        PREFIX "${PROJECT_BINARY_DIR}/gtest"
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG release-1.8.0
        INSTALL_COMMAND true  # currently no install command
        )

    include_directories(${PROJECT_BINARY_DIR}/gtest/src/gtest/googletest/include)
    include_directories(${PROJECT_SOURCE_DIR}/src)

    link_directories(${PROJECT_BINARY_DIR}/gtest/src/gtest-build/googlemock/gtest/)

    add_executable(
        cpp_test
        test/main.cpp
        test/energy_spherical.cpp
    #   test/energy_cartesian.cpp
        )

    add_dependencies(cpp_test gtest)

    # workaround:
    # different dynamic lib suffix on mac
    if(APPLE)
        set(_dyn_lib_suffix dylib)
    else()
        set(_dyn_lib_suffix so)
    endif()

    target_link_libraries(
        cpp_test
        libgtest.a
        xcint
        density
        ${PROJECT_BINARY_DIR}/external/xcfun-build/libxcfun.a
        ${NUMGRID_LIBS}
        ${PROJECT_BINARY_DIR}/balboa/src/balboa-build/lib/libbalboa.so
        ${MATH_LIBS}
        pthread
        )

    add_test(cpp_test ${PROJECT_BINARY_DIR}/bin/cpp_test ${PROJECT_SOURCE_DIR}/test)

    if(ENABLE_FC_SUPPORT)
        add_executable(
            fortran_test
            test/test.f90
            )

        target_link_libraries(
            fortran_test
            xcint_fortran
            ${PROJECT_BINARY_DIR}/balboa/src/balboa-build/lib/libbalboa.so
            ${NUMGRID_LIBS}
            )

        set_target_properties(fortran_test PROPERTIES LINKER_LANGUAGE Fortran)

        add_test(fortran_test ${PROJECT_BINARY_DIR}/bin/fortran_test ${PROJECT_SOURCE_DIR}/test)
    endif()
endif()
