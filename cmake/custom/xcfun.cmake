find_package(XCFun CONFIG QUIET)

if(TARGET XCFun::xcfun)
  message(STATUS "Found XCFun version ${XCFun_VERSION}")
else()
  include(FetchContent)

  set(XCFun_XC_MAX_ORDER 6)

  FetchContent_Populate(xcfun_sources
    QUIET
    GIT_REPOSITORY
      https://github.com/dftlibs/xcfun.git
    GIT_TAG
      v2.0.1
    )

  add_subdirectory(
    ${xcfun_sources_SOURCE_DIR}
    ${xcfun_sources_BINARY_DIR}
    )

  unset(XCFun_XC_MAX_ORDER)
endif()
