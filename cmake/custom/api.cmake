include_directories(${PROJECT_SOURCE_DIR}/api)

file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/include)

file(COPY ${PROJECT_SOURCE_DIR}/api/xcint.h DESTINATION ${PROJECT_BINARY_DIR}/include)
