include(ExternalProject)

set(ExternalProjectCMakeArgs
    -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
    -DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR}/external
    -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
    -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
    -DENABLE_FORTRAN_INTERFACE=OFF
    -DENABLE_TESTALL=OFF
    -DENABLE_STATIC_LINKING=ON # we need the -fPIC
    -DXC_MAX_ORDER=6
    )

ExternalProject_Add(xcfun
    SOURCE_DIR ${PROJECT_SOURCE_DIR}/external/xcfun
    BINARY_DIR ${PROJECT_BINARY_DIR}/external/xcfun-build
    STAMP_DIR ${PROJECT_BINARY_DIR}/external/xcfun-stamp
    TMP_DIR ${PROJECT_BINARY_DIR}/external/xcfun-tmp
    INSTALL_DIR ${PROJECT_BINARY_DIR}/external
    CMAKE_ARGS ${ExternalProjectCMakeArgs}
    )
