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
