#include "octree.hpp"

#include "test-utils.hpp"

#include "gtest/gtest.h"

#include <iostream>
#include <fstream>

TEST(OctreeNode, DistToNode)
{
	OctreeNode n(Position(10.0, 20.0, 30.0), 2.0);
	auto d1 = n.getDistsToNode(Position(10.0, 20.0, 30.0));
	EXPECT_NEAR(d1.farest, sqrt(3.0), 1e-6);
	EXPECT_NEAR(d1.nearest, sqrt(3.0), 1e-6);

	auto d2 = n.getDistsToNode(Position(12.0, 22.0, 32.0));
	EXPECT_NEAR(d2.nearest, sqrt(3.0), 1e-6);
	EXPECT_NEAR(d2.farest, 3*sqrt(3.0), 1e-6);
}

TEST(OctreeBase, Instantiation)
{
	ASSERT_NO_THROW(Octree());
	Position c(28.0, -43.2212, 1.23e50);
	ASSERT_NO_THROW(Octree(c, 17.8));
	ASSERT_NO_THROW(Octree(22.0));
	ASSERT_NO_THROW(Octree());
}

TEST(OctreeBase, InitWithAdding)
{
}

TEST(OctreeBase, Adding)
{
	Octree oct(
		Position(0.0, 0.0, 0.0),
		2
	);
	ASSERT_EQ(oct.count(), 0);
	// Level 0
	OctreeElement e1(-0.01, -0.01, -0.01, 1.0);
	ASSERT_NO_THROW(oct.add(e1));
	ASSERT_EQ(oct.count(), 1);

	// Level 1
	OctreeElement e2(0.01, 0.01, 0.01, 1.0);
	ASSERT_NO_THROW(oct.add(e2));
	EXPECT_TRUE(oct.root().subnodes[0] != nullptr);
	EXPECT_TRUE(oct.root().subnodes[1] == nullptr);
	EXPECT_TRUE(oct.root().subnodes[2] == nullptr);
	EXPECT_TRUE(oct.root().subnodes[3] == nullptr);
	EXPECT_TRUE(oct.root().subnodes[4] == nullptr);
	EXPECT_TRUE(oct.root().subnodes[5] == nullptr);
	EXPECT_TRUE(oct.root().subnodes[6] == nullptr);
	EXPECT_TRUE(oct.root().subnodes[7] != nullptr);

	// Level 2
	OctreeElement e3(-0.011, -0.011, -0.011, 1.0);
	ASSERT_NO_THROW(oct.add(e3));
	ASSERT_EQ(oct.count(), 3);
	EXPECT_TRUE(oct.root().subnodes[0]->subnodes[0] == nullptr);
	EXPECT_TRUE(oct.root().subnodes[0]->subnodes[1] == nullptr);
	EXPECT_TRUE(oct.root().subnodes[0]->subnodes[2] == nullptr);
	EXPECT_TRUE(oct.root().subnodes[0]->subnodes[3] == nullptr);
	EXPECT_TRUE(oct.root().subnodes[0]->subnodes[4] == nullptr);
	EXPECT_TRUE(oct.root().subnodes[0]->subnodes[5] == nullptr);
	EXPECT_TRUE(oct.root().subnodes[0]->subnodes[6] == nullptr);
	EXPECT_TRUE(oct.root().subnodes[0]->subnodes[7] != nullptr);

	OctreeElement e4(1.0, 1.0, 1.0, 1.0);
	ASSERT_NO_THROW(oct.add(e4));
	OctreeElement e5(-1.0, 1.0, -1.0, 1.0);
	ASSERT_NO_THROW(oct.add(e5));
}

TEST(OctreeBase, DbgOutput)
{
	Octree oct(
		Position(0.0, 0.0, 0.0),
		2
	);
	OctreeElement e1(-0.01, -0.01, -0.01, 1.0);
	oct.add(e1);
	OctreeElement e2(0.01, 0.01, 0.01, 1.0);
	oct.add(e2);
	OctreeElement e3(-0.011, -0.011, -0.011, 1.0);
	oct.add(e3);
	std::ofstream file("dbg-out-test.txt", std::ios::out);
	ASSERT_NO_THROW(oct.dbgOutCoords(file));
}

class OctreeAccess : public ::testing::Test
{
public:
	void addGrid(int n)
	{
		for (int i=0; i<n; i++)
			for (int j=0; j<n; j++)
				for (int k=0; k<n; k++)
				{
					Position p(
						-size/2.0 + size / (n-1) * i,
						-size/2.0 + size / (n-1) * j,
						-size/2.0 + size / (n-1) * k
					);
					addElement(p);
				}
	}

	Position& findNearestBruteForce(const Position& pos)
	{
		double distMin = (positions.front() - pos).len();
		Position* p = &positions.front();
		for (auto it = positions.begin(); it != positions.end(); ++it)
		{
			double d = (*it - pos).len();
			if (d < distMin)
			{
				distMin = d;
				p = &(*it);
			}
		}
		return *p;
	}

	void addElement(const Position& pos)
	{
		oct.add(OctreeElement(pos, 1.0));
		//std::cout << p.str() << std::endl;
		positions.push_back(pos);
	}
	double size = 2.0;
	Octree oct{Position(0.0, 0.0, 0.0), size};
	std::vector<Position> positions;
};

TEST_F(OctreeAccess, Initialization)
{
	ASSERT_NO_THROW(addGrid(10));
}

TEST_F(OctreeAccess, FindNearest0)
{
	Position target(0.1, -23, 876);
	ASSERT_ANY_THROW(oct.getNearest(target) );
}

TEST_F(OctreeAccess, FindNearest1)
{
	addElement(Position(0.2, -0.8, 1.0));
	Position target(0.1, -23, 876);
	OctreeElement *ans = nullptr;
	Position *brute = &(findNearestBruteForce(target));
	ASSERT_NO_THROW(ans = & oct.getNearest(target) );
	ASSERT_EQ(*brute, ans->pos);
}

TEST_F(OctreeAccess, FindNearest4)
{
	addGrid(4);
	Position target(0.1, -0.8, 0.5);
	OctreeElement *ans = nullptr;
	Position *brute = &(findNearestBruteForce(target));
	ASSERT_NO_THROW(ans = & oct.getNearest(target) );
	ASSERT_EQ(*brute, ans->pos);
}

TEST_F(OctreeAccess, FindNearest20)
{
	addGrid(2);
	{
		Position target(0.1, -0.8, 0.5);
		OctreeElement *ans = nullptr;
		Position *brute = &(findNearestBruteForce(target));
		ASSERT_NO_THROW(ans = & oct.getNearest(target) );
		ASSERT_EQ(*brute, ans->pos);
	}
	{
		Position target(10, -678, -0.0001);
		OctreeElement *ans = nullptr;
		Position *brute = &(findNearestBruteForce(target));
		ASSERT_NO_THROW(ans = & oct.getNearest(target) );
		ASSERT_EQ(*brute, ans->pos);
	}
	{
		Position target(1.0, -0.8, 0.5);
		OctreeElement *ans = nullptr;
		Position *brute = &(findNearestBruteForce(target));
		ASSERT_NO_THROW(ans = & oct.getNearest(target) );
		ASSERT_EQ(*brute, ans->pos);
	}
	{
		Position target(positions.front());
		OctreeElement *ans = nullptr;
		Position *brute = &(findNearestBruteForce(target));
		ASSERT_NO_THROW(ans = & oct.getNearest(target) );
		ASSERT_EQ(*brute, ans->pos);
	}
	// @todo: empty, only one
}

TEST_F(OctreeAccess, Performance)
{
	addGrid(100);
	Position target(0.1, -0.8, 0.5);
	auto tb = TimeMesurment::execution(
		[this, &target](){ findNearestBruteForce(target); }
	);
	auto t  = TimeMesurment::execution(
		[this, &target](){ oct.getNearest(target); }
	);
	double part = double(t) / double(tb);
	ASSERT_TRUE(part < 0.1);
	//std::cout << "part: " << part<< std::endl;
}
