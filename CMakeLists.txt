# set minimum cmake version
cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

# project name
project(xcint)

# do not rebuild if rules (compiler flags) change
set(CMAKE_SKIP_RULE_DEPENDENCY TRUE)

# if CMAKE_BUILD_TYPE undefined, we set it to Debug
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Debug")
endif()

# directory which holds enabled cmake modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
    ${PROJECT_SOURCE_DIR}/cmake/modules)

# included cmake modules
include(autocmake_python)
include(autocmake_cxx)
include(autocmake_GNU.CXX)
include(autocmake_Intel.CXX)
include(autocmake_code_coverage)
include(autocmake_omp)
include(autocmake_math_libs)
include(autocmake_safeguards)
include(autocmake_default_build_paths)
include(autocmake_ConfigExternal)
include(autocmake_api)
include(autocmake_custom)
include(autocmake_gtest)
include(autocmake_test)
