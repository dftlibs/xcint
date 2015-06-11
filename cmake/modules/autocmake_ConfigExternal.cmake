
include(ExternalProject)

macro(add_external _project)

    ExternalProject_Add(${_project}
        DOWNLOAD_COMMAND ${UPDATE_COMMAND}
        DOWNLOAD_DIR ${PROJECT_SOURCE_DIR}
        SOURCE_DIR ${PROJECT_SOURCE_DIR}/external/${_project}
        BINARY_DIR ${PROJECT_BINARY_DIR}/external/${_project}-build
        STAMP_DIR ${PROJECT_BINARY_DIR}/external/${_project}-stamp
        TMP_DIR ${PROJECT_BINARY_DIR}/external/${_project}-tmp
        INSTALL_DIR ${PROJECT_BINARY_DIR}/external
        CMAKE_ARGS ${ExternalProjectCMakeArgs}
        )
    include_directories(${PROJECT_BINARY_DIR}/external/${_project}-build)
    include_directories(${PROJECT_BINARY_DIR}/external/${_project}-build/modules)
    link_directories(${PROJECT_BINARY_DIR}/external/lib)

    if(ALWAYS_RESET_EXTERNAL_BUILDS)
        # remove stamps for external builds so that they are rebuilt every time
        add_custom_command(
            TARGET ${_project}
            PRE_BUILD
            COMMAND rm -rf ${PROJECT_BINARY_DIR}/external/${_project}-stamp
            WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
            )
    endif()
endmacro()
