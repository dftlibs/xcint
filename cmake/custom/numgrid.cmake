include(ExternalProject)

set(numgrid_args
    -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
    -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
    -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
    -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
    -DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR}/external
    -DENABLE_FC_SUPPORT=${ENABLE_FC_SUPPORT}
    )

ExternalProject_Add(numgrid
    SOURCE_DIR ${PROJECT_SOURCE_DIR}/external/numgrid
    BINARY_DIR ${PROJECT_BINARY_DIR}/external/numgrid-build
    STAMP_DIR ${PROJECT_BINARY_DIR}/external/numgrid-stamp
    TMP_DIR ${PROJECT_BINARY_DIR}/external/numgrid-tmp
    INSTALL_DIR ${PROJECT_BINARY_DIR}/external
    CMAKE_ARGS ${numgrid_args}
    )

include_directories(${PROJECT_SOURCE_DIR}/external/numgrid/api)
include_directories(${PROJECT_BINARY_DIR}/external/numgrid-build/modules)
