#include <cstdlib>

#include "gtest/gtest.h"

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    if (argc == 1)
    {
        fputs("ERROR: call this binary with XCINT_TEST_DIRECTORY as first argument!\n", stderr);
        abort();
    }
    setenv("XCINT_TEST_DIRECTORY", argv[1], true);
    return RUN_ALL_TESTS();
}
