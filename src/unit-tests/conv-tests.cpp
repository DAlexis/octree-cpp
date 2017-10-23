#include "octree.hpp"

#include "test-utils.hpp"

#include "gtest/gtest.h"

#include <iostream>

using namespace std;

TEST(ConvolutionInternals, FindScales)
{
    Convolution c;
    c.addScale(20, 2);
    c.addScale(10, 1);
    c.addScale(50, 5);
    c.addScale(40, 4);
    c.addScale(30, 3);
    ASSERT_NO_THROW(c.sortDistsScales());
    ASSERT_EQ(c.findScale(15), 1);
    ASSERT_EQ(c.findScale(90), 5);
    ASSERT_EQ(c.findScale(31), 3);
}

class ConvolutionTests : public ::testing::Test
{
public:
    void addSomePoints()
    {
        oct.add(OctreeElement(Position(2.0, 3.0, -8.0), 3.0));
        oct.add(OctreeElement(Position(0.0, 0.0, 0.0), 1.0));
        oct.add(OctreeElement(Position(8.0, 9.0, 9.0), 1.0));
        oct.add(OctreeElement(Position(-3.0, -9.0, -4.0), 2.0));
        oct.add(OctreeElement(Position(-7.0, -9.0, -4.0), 1.0));
        oct.add(OctreeElement(Position(-1.0, -4.0, -2.0), 1.0));
    }

    int callsCounter = 0;
    double somePointsMass = 9;
    Octree oct{Position(0.0, 0.0, 0.0), 20};
    Convolution conv;
    Convolution::Visitor massSumVisitor =
        [this](const Position& target, const Position& object, double mass)
        {
            callsCounter++;
            return mass;
        };
};

TEST_F(ConvolutionTests, ConvoluteNoScale)
{
    addSomePoints();
    ASSERT_EQ(oct.mass(), 9.0);
    ASSERT_NO_THROW(conv.convolute(oct, Position(0.0, 0.0, 0.0), massSumVisitor));
    ASSERT_EQ(callsCounter, 6);
    ASSERT_EQ(conv.convolute(oct, Position(0.0, 0.0, 0.0), massSumVisitor), somePointsMass);
}

TEST_F(ConvolutionTests, ConvoluteOneScalingZone)
{
    addSomePoints();
    callsCounter = 0;
    ASSERT_EQ(conv.convolute(oct, Position(0.0, 0.0, 0.0), massSumVisitor), somePointsMass);
    int cc = callsCounter;
    conv.addScale(0.1, 1000);
    callsCounter = 0;
    ASSERT_EQ(conv.convolute(oct, Position(15.0, 15.0, 15.0), massSumVisitor), somePointsMass);
    ASSERT_EQ(callsCounter, 1);
    callsCounter = 0;
    ASSERT_EQ(conv.convolute(oct, Position(9.0, 9.0, 9.0), massSumVisitor), somePointsMass);
    ASSERT_LT(callsCounter, cc);
}
