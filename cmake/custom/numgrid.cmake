include(ExternalProject)

ExternalProject_Add(
    numgrid
    PREFIX "${PROJECT_BINARY_DIR}/numgrid"
    GIT_REPOSITORY https://github.com/dftlibs/numgrid.git
    GIT_TAG v0.5.0
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

# consider renaming to include instead of modules
include_directories(${PROJECT_BINARY_DIR}/numgrid/src/numgrid-build/modules)
