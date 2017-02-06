include(ExternalProject)

ExternalProject_Add(
    numgrid
    PREFIX "${PROJECT_BINARY_DIR}/numgrid"
    GIT_REPOSITORY https://github.com/dftlibs/numgrid.git
    GIT_TAG 19c27ea532b6cbe1c2059afde3d7b225881af00e
    INSTALL_COMMAND true  # currently no install command
    )

link_directories(${PROJECT_BINARY_DIR}/numgrid/src/numgrid-build/src)
link_directories(${PROJECT_BINARY_DIR}/numgrid/src/numgrid-build/lib)

# workaround:
# different dynamic lib suffix on mac
if(APPLE)
    set(_dyn_lib_suffix dylib)
else()
    set(_dyn_lib_suffix so)
endif()

set(NUMGRID_LIBS
    libnumgrid.${_dyn_lib_suffix}
    liblebedev.${_dyn_lib_suffix}
    libnumgrid_fortran.a
    )

include_directories(${PROJECT_BINARY_DIR}/numgrid/src/numgrid/api)
include_directories(${PROJECT_BINARY_DIR}/numgrid/src/numgrid-build/include)
include_directories(${PROJECT_BINARY_DIR}/numgrid/src/numgrid-build/modules)
