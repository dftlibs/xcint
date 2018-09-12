# we will use the network to fetch Google Test sources
# make it possible to disable unit tests when not on network
option(ENABLE_UNIT_TESTS "Enable unit tests" ON)
message(STATUS "Enable testing: ${ENABLE_UNIT_TESTS}")

if(ENABLE_UNIT_TESTS)
  # the following code to fetch googletest
  # is inspired by and adapted after:
  #   - https://cmake.org/cmake/help/v3.11/module/FetchContent.html
  include(FetchContent)

  FetchContent_Declare(
    googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG        release-1.8.1
  )

  FetchContent_GetProperties(googletest)

  if(NOT googletest_POPULATED)
    FetchContent_Populate(googletest)

    # Prevent GoogleTest from overriding our compiler/linker options
    # when building with Visual Studio
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
    # Prevent GoogleTest from using PThreads
    set(gtest_disable_pthreads ON CACHE BOOL "" FORCE)

    # adds the targers: gtest, gtest_main, gmock, gmock_main
    add_subdirectory(
      ${googletest_SOURCE_DIR}
      ${googletest_BINARY_DIR}
      )

    # Silence std::tr1 warning on MSVC
    if(MSVC)
      foreach(_tgt gtest gtest_main gmock gmock_main)
        target_compile_definitions(${_tgt}
          PRIVATE
            "_SILENCE_TR1_NAMESPACE_DEPRECATION_WARNING"
          )
      endforeach()
    endif()
  endif()

  add_executable(
    cpp_test
    main.cpp
    energy_spherical.cpp
    )

  find_package(BLAS REQUIRED)

  target_link_libraries(
    cpp_test
    gtest_main
    xcint
    density
    ${PROJECT_BINARY_DIR}/external/xcfun-build/libxcfun.a
    ${NUMGRID_LIBS}
    ${PROJECT_BINARY_DIR}/balboa/src/balboa-build/lib/libbalboa.so
    ${BLAS_LIBRARIES}
    )

  add_test(
    NAME cpp_test
    COMMAND $<TARGET_FILE:cpp_test> ${PROJECT_SOURCE_DIR}/test
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
    )

  if(ENABLE_FC_SUPPORT)
    add_executable(
      fortran_test
      test.f90
      )

    target_link_libraries(
      fortran_test
      xcint_fortran
      xcint
      density
      ${PROJECT_BINARY_DIR}/external/xcfun-build/libxcfun.a
      ${NUMGRID_LIBS}
      ${PROJECT_BINARY_DIR}/balboa/src/balboa-build/lib/libbalboa.so
      ${BLAS_LIBRARIES}
      )

    set_target_properties(fortran_test PROPERTIES LINKER_LANGUAGE Fortran)

    add_test(
      NAME fortran_test
      COMMAND $<TARGET_FILE:fortran_test> ${PROJECT_SOURCE_DIR}/test
      WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
      )
  endif()
endif()
