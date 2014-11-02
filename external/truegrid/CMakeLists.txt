cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

project(truegrid)

enable_language(C CXX)

set(CMAKE_MODULE_PATH
    ${CMAKE_MODULE_PATH}
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/compilers
    )

include(OptimizeForArchitecture)
OptimizeForArchitecture()
if(AVX_FOUND)
    add_definitions(-DCPU_MEM_ALIGN=32)
    message(STATUS "Setting CPU_MEM_ALIGN=32")
elseif(SSE3_FOUND)
    add_definitions(-DCPU_MEM_ALIGN=16)
    message(STATUS "Setting CPU_MEM_ALIGN=16")
else()
    message(FATAL_ERROR "neither SSE3 nor AVX detected; adapt CMakeLists.txt")
endif()

include(ConfigSafeGuards)
include(ConfigCompilerFlags)

include_directories(
    ${PROJECT_SOURCE_DIR}/src
    ${PROJECT_SOURCE_DIR}/interfaces
    ${PROJECT_SOURCE_DIR}/external/lebedev
    )

add_library(
    truegrid
    src/becke_partitioning.cpp
    src/grid_radial.cpp
    src/Grid.cpp
    src/io.cpp
    src/rolex.cpp
    src/MemAllocator.cpp
    interfaces/truegrid_c_interface.cpp
    external/lebedev/sphere_lebedev_rule.cpp
    )

install(TARGETS truegrid ARCHIVE DESTINATION lib)
