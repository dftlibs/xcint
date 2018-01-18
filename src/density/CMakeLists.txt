set(AO_BLOCK_LENGTH "128" CACHE STRING "AO block length")
set(AO_CHUNK_LENGTH "32" CACHE STRING "AO chunk length")
set(MAX_GEO_DIFF_ORDER "5" CACHE STRING "Maximum geometric differentiation order")
set(MAX_L_VALUE "5" CACHE STRING "Maximum L value")

configure_file(
    parameters.h.in
    ${PROJECT_BINARY_DIR}/generated/density_parameters.h
    @ONLY
    )

add_library(
    density
    SHARED
    blas_interface.cpp
    blas_interface.h
    density.cpp
    density.h
    compress.cpp
    compress.h
    )

target_include_directories(
    density
    PRIVATE
    ${PROJECT_BINARY_DIR}/generated
    )

set_target_properties(density PROPERTIES COMPILE_FLAGS "-fPIC")
