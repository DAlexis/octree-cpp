#include "octree.hpp"

#include "test-utils.hpp"

#include "gtest/gtest.h"

#include <iostream>

using namespace std;
using namespace octree;

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
        oct.add(make_shared<ElementValue>(Position(2.0, 3.0, -8.0), 3.0));
        oct.add(make_shared<ElementValue>(Position(0.0, 0.0, 0.0), 1.0));
        oct.add(make_shared<ElementValue>(Position(8.0, 9.0, 9.0), 1.0));
        oct.add(make_shared<ElementValue>(Position(-3.0, -9.0, -4.0), 2.0));
        oct.add(make_shared<ElementValue>(Position(-7.0, -9.0, -4.0), 1.0));
        oct.add(make_shared<ElementValue>(Position(-1.0, -4.0, -2.0), 1.0));
    }

    void addManyPoints()
    {
        PointsGenerator::addGrid(10, 10, oct, &positions);
    }

    double getCoulombFieldBruteForce(const Position& target)
    {
        double result = 0;
        for (auto &it: positions)
        {
            result += coulomb(target, it, 1.0);
        }
        return result;
    }

    int callsCounter = 0;
    double somePointsMass = 9;

    Octree oct{Position(0.0, 0.0, 0.0), 20};

    // Vector for Coulomb
    std::vector<Position> positions;
    Convolution conv;
    Convolution::Visitor massSumVisitor =
        [this](const Position& target, const Position& object, double mass)
        {
            callsCounter++;
            return mass;
        };
    Convolution::Visitor coulomb =
        [this](const Position& target, const Position& object, double mass)
        {
            return mass/(target-object).len();
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

TEST_F(ConvolutionTests, ConvoluteCoulomb1)
{
    addManyPoints();
    Position p1 = {1.123, 2.345, 3.456};
    double realField = getCoulombFieldBruteForce(p1);
    double convField = conv.convolute(oct, p1, coulomb);
    ASSERT_NEAR(realField, convField, 1e-8) << "Convolution gives bad ansver with ideal precision";

    conv.addScale(5, 3);
    convField = conv.convolute(oct, p1, coulomb);
    ASSERT_NEAR(realField, convField, 1e-3*realField) << "Convolution result has bad precision";
    //cout << realField << " " << convField << endl;

    conv.addScale(7, 10);
    convField = conv.convolute(oct, p1, coulomb);
    ASSERT_NEAR(realField, convField, 3e-3*realField) << "Convolution result has bad precision";
    //cout << realField << " " << convField << endl;
}

TEST_F(ConvolutionTests, ConvoluteCoulomb2)
{
    addManyPoints();
    Position p1 = {11.23, -23.45, -34.56};
    double realField = getCoulombFieldBruteForce(p1);
    double convField = conv.convolute(oct, p1, coulomb);
    ASSERT_NEAR(realField, convField, 1e-8) << "Convolution gives bad ansver with ideal precision";

    conv.addScale(5, 3);
    convField = conv.convolute(oct, p1, coulomb);
    ASSERT_NEAR(realField, convField, 1e-3*realField) << "Convolution result has bad precision";
    //cout << realField << " " << convField << endl;

    conv.addScale(7, 10);
    convField = conv.convolute(oct, p1, coulomb);
    ASSERT_NEAR(realField, convField, 3e-3*realField) << "Convolution result has bad precision";
    //cout << realField << " " << convField << endl;

    conv.addScale(10, 1000);
    convField = conv.convolute(oct, p1, coulomb);
    ASSERT_NEAR(realField, convField, 3e-3*realField) << "Convolution result has bad precision";
    //cout << realField << " " << convField << endl;
}

