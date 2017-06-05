#include "octree.hpp"

#include "gtest/gtest.h"

TEST(TestCase1, SqrTest)
{
    EXPECT_EQ(4.0, sqr(2.0));
    EXPECT_EQ(9.0, sqr(3.0));
}
