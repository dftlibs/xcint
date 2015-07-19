option(ENABLE_CBLAS "Enable CBLAS interface instead of the regular Fortran BLAS" OFF)

if(ENABLE_CBLAS)
    add_definitions(-DENABLE_CBLAS)
endif()
