#include "octree.hpp"

#include "gtest/gtest.h"

TEST(Octree, Instantiation)
{
	ASSERT_NO_THROW(Octree());
	Position c1(0.0, 0.0, 0.0);
	Position c2(10.0, 10.0, 10.0);
	ASSERT_NO_THROW(Octree(c1, c2));
}

TEST(Octree, Adding)
{
	Octree oct(
		Position(-1.0, -1.0, -1.0),
		Position(1.0, 1.0, 1.0)
	);
	ASSERT_EQ(oct.count(), 0);
	OctreeElement e1(0.0, 0.0, 0.0, 1.0);
	ASSERT_NO_THROW(oct.add(e1));
	ASSERT_EQ(oct.count(), 1);
}
