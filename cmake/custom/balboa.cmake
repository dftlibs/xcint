include(ExternalProject)

ExternalProject_Add(
    balboa
    PREFIX "${PROJECT_BINARY_DIR}/balboa"
    GIT_REPOSITORY https://github.com/bast/balboa.git
    GIT_TAG c095e8ccf8fa797a55114723b1705e11751158b8
    INSTALL_COMMAND true  # currently no install command
    CMAKE_ARGS "-DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}"
    )

# workaround:
# different dynamic lib suffix on mac
if(APPLE)
    set(_dyn_lib_suffix dylib)
else()
    set(_dyn_lib_suffix so)
endif()

include_directories(${PROJECT_BINARY_DIR}/balboa/src/balboa/balboa)
include_directories(${PROJECT_BINARY_DIR}/balboa/src/balboa-build/include)
