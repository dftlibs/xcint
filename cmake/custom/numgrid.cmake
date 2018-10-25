include(ExternalProject)

ExternalProject_Add(
  numgrid
  PREFIX "${PROJECT_BINARY_DIR}/numgrid"
  GIT_REPOSITORY https://github.com/dftlibs/numgrid.git
  GIT_TAG v1.0.2
  INSTALL_COMMAND true  # currently no install command
  CMAKE_ARGS
    -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
    -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
    -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
    -DENABLE_UNIT_TESTS=OFF
  )

include(GNUInstallDirs)
link_directories(${PROJECT_BINARY_DIR}/numgrid/src/numgrid-build/${CMAKE_INSTALL_LIBDIR})

# workaround:
# different dynamic lib suffix on mac
if(APPLE)
    set(_dyn_lib_suffix dylib)
else()
    set(_dyn_lib_suffix so)
endif()

set(NUMGRID_LIBS
    libnumgrid.${_dyn_lib_suffix}
    )

include_directories(${PROJECT_BINARY_DIR}/numgrid/src/numgrid/numgrid)
include_directories(${PROJECT_BINARY_DIR}/numgrid/src/numgrid-build/modules)
include_directories(${PROJECT_BINARY_DIR}/numgrid/src/numgrid-build/numgrid)
include_directories(${PROJECT_BINARY_DIR}/numgrid/src/numgrid-build/include)
