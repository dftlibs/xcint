include_directories(${PROJECT_BINARY_DIR}/generated)

file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/generated)

# generate ave_contributions.h
add_custom_command(
    OUTPUT
        ${PROJECT_BINARY_DIR}/generated/ave_contributions.h
    COMMAND
        ${PYTHON_EXECUTABLE} ${PROJECT_SOURCE_DIR}/src/generate_ave_contributions.py > ${PROJECT_BINARY_DIR}/generated/ave_contributions.h
    DEPENDS
        src/generate_ave_contributions.py
    )
add_custom_target(
    generate_ave
    ALL
    DEPENDS
        ${PROJECT_BINARY_DIR}/generated/ave_contributions.h
    )

add_dependencies(xcint generate_ave)
