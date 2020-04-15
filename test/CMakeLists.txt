# we will use the network to fetch Google Test sources
# make it possible to disable unit tests when not on network
option(ENABLE_XCINT_UNIT_TESTS "Enable unit tests" ON)
message(STATUS "Enable testing: ${ENABLE_XCINT_UNIT_TESTS}")

if(ENABLE_XCINT_UNIT_TESTS)
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

  FetchContent_Declare(numgrid_sources
    QUIET
    URL
      https://github.com/dftlibs/numgrid/archive/v1.1.1.tar.gz
    )

  FetchContent_GetProperties(numgrid_sources)

  set(ENABLE_UNIT_TESTS FALSE CACHE BOOL "")
  set(ENABLE_FC_SUPPORT FALSE CACHE BOOL "")

  if(NOT numgrid_sources_POPULATED)
    FetchContent_Populate(numgrid_sources)

    add_subdirectory(
      ${numgrid_sources_SOURCE_DIR}
      ${numgrid_sources_BINARY_DIR}
      )
  endif()

  add_executable(
    cpp_test
    main.cpp
    )

  find_package(BLAS REQUIRED)

  add_library(
    energy_spherical
    energy_spherical.cpp
    )

  target_link_libraries(
    energy_spherical
    gtest_main
    numgrid-objects
    )

  target_link_libraries(
    cpp_test
    xcint
    density
    balboa
    energy_spherical
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
      balboa
      numgrid-shared
      numgrid_fortran_bindings
      ${BLAS_LIBRARIES}
      )

    target_include_directories(
      fortran_test
      PRIVATE
        ${PROJECT_BINARY_DIR}/_deps/numgrid_sources-build/modules
      )

    add_dependencies(fortran_test numgrid_fortran_bindings)

    set_target_properties(fortran_test PROPERTIES LINKER_LANGUAGE Fortran)

    add_test(
      NAME fortran_test
      COMMAND $<TARGET_FILE:fortran_test> ${PROJECT_SOURCE_DIR}/test
      WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
      )
  endif()
endif()
