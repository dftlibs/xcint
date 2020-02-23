enable_testing()

add_test(
    NAME
        test_main
    COMMAND
        py.test -vv -s ${PROJECT_SOURCE_DIR}/test/test.py
    )

set_property(
    TEST
        test_main
    PROPERTY
        ENVIRONMENT
            BALBOA_LIBRARY_DIR=${PROJECT_BINARY_DIR}/lib
            BALBOA_INCLUDE_DIR=${PROJECT_SOURCE_DIR}/balboa
            PYTHONPATH=${PROJECT_SOURCE_DIR}
    )

add_test(
    NAME
        test_generate
    COMMAND
        py.test -vv -s ${PROJECT_SOURCE_DIR}/balboa/generate.py
    )

set_property(
    TEST
        test_generate
    PROPERTY
        ENVIRONMENT
            BALBOA_LIBRARY_DIR=${PROJECT_BINARY_DIR}/lib
            BALBOA_INCLUDE_DIR=${PROJECT_SOURCE_DIR}/balboa
            PYTHONPATH=${PROJECT_SOURCE_DIR}
    )
