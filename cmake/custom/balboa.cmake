include(ExternalProject)

ExternalProject_Add(
    balboa
    PREFIX "${PROJECT_BINARY_DIR}/balboa"
    GIT_REPOSITORY https://github.com/bast/balboa.git
    GIT_TAG a26be6df03c3736c96e401b1a3b0c2d9f4f8272c
    INSTALL_COMMAND true  # currently no install command
    )

link_directories(${PROJECT_BINARY_DIR}/balboa/src/balboa-build/src)
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

# consider renaming to include instead of modules
include_directories(${PROJECT_BINARY_DIR}/balboa/src/balboa-build/modules)
