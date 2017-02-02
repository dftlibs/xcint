include(ExternalProject)

ExternalProject_Add(
    balboa
    PREFIX "${PROJECT_BINARY_DIR}/balboa"
    GIT_REPOSITORY https://github.com/bast/balboa.git
    GIT_TAG 1b27595a32c8bbbd22b5a380056cce578ea3299a
    INSTALL_COMMAND true  # currently no install command
    )

link_directories(${PROJECT_BINARY_DIR}/balboa/src/balboa-build/lib)

# workaround:
# different dynamic lib suffix on mac
if(APPLE)
    set(_dyn_lib_suffix dylib)
else()
    set(_dyn_lib_suffix so)
endif()

set(BALBOA_LIBS
    libbalboa.${_dyn_lib_suffix}
    )

include_directories(${PROJECT_BINARY_DIR}/balboa/src/balboa/api)
