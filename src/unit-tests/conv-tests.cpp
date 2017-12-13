#include "octree.hpp"

#include "test-utils.hpp"

#include "gtest/gtest.h"

#include <iostream>

using namespace std;
using namespace octree;

TEST(DiscreteScales, FindScales)
{
    DiscreteScales c;
    c.addScale(20, 2);
    c.addScale(10, 1);
    c.addScale(50, 5);
    c.addScale(40, 4);
    c.addScale(30, 3);
    ASSERT_EQ(c.findScale(15), 1);
    ASSERT_EQ(c.findScale(90), 5);
    ASSERT_EQ(c.findScale(31), 3);
}

class ConvolutionTests : public ::testing::Test
{
public:
    shared_ptr<ElementValue> addSomePoints()
    {
        positions.push_back(Position(2.0, 3.0, -8.1));
        masses.push_back(1.0);
        oct.add(make_shared<ElementValue>(positions.back(), masses.back()));

        positions.push_back(Position(0.0, 0.0, 0.2));
        masses.push_back(2.0);
        oct.add(make_shared<ElementValue>(positions.back(), masses.back()));

        positions.push_back(Position(8.0, 9.0, 9.3));
        masses.push_back(3.0);
        oct.add(make_shared<ElementValue>(positions.back(), masses.back()));

        positions.push_back(Position(-3.0, -9.0, -4.4));
        masses.push_back(4.0);
        oct.add(make_shared<ElementValue>(positions.back(), masses.back()));

        positions.push_back(Position(-7.0, -9.0, -4.5));
        masses.push_back(5.0);
        oct.add(make_shared<ElementValue>(positions.back(), masses.back()));

        positions.push_back(Position(-1.0, -4.0, -2.6));
        masses.push_back(6.0);
        auto last = make_shared<ElementValue>(positions.back(), masses.back());
        oct.add(last);
        return last;
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
    double somePointsMass = 21;

    Octree oct{Position(0.0, 0.0, 0.0), 20};

    // Vector for Coulomb
    std::vector<Position> positions;
    std::vector<double> masses;

    DiscreteScales scales;
    Convolution<double> conv{scales};
    Convolution<double>::Visitor massSumVisitor =
        [this](const Position& target, const Position& object, double mass)
        {
            callsCounter++;
            return mass;
        };
    Convolution<double>::Visitor coulomb =
        [this](const Position& target, const Position& object, double mass)
        {
            double d = (target-object).len();
            if (d == 0.0)
                return 0.0;
            return mass/d;
        };
};

TEST_F(ConvolutionTests, ConvoluteEmptyTree)
{
    double result = 0.0;
    ASSERT_NO_THROW(result = conv.convolute(oct, Position(0.0, 0.0, 0.0), massSumVisitor));
    ASSERT_EQ(callsCounter, 0);
    ASSERT_EQ(result, 0.0);
}

TEST_F(ConvolutionTests, ConvoluteNoScale)
{
    addSomePoints();
    ASSERT_EQ(oct.mass(), somePointsMass);
    ASSERT_NO_THROW(conv.convolute(oct, Position(0.0, 0.0, 0.0), massSumVisitor));
    ASSERT_EQ(callsCounter, 6);
    ASSERT_EQ(conv.convolute(oct, Position(0.0, 0.0, 0.0), massSumVisitor), somePointsMass);
}


TEST_F(ConvolutionTests, ConvoluteWithScale)
{
    auto last = addSomePoints();
    double result = 0;
    double trueResult = 0;
    auto target = last->pos;
    // Convolution in position of last without last
    ASSERT_NO_THROW(result = conv.convolute(oct, last->pos, coulomb));

    // Manually calculated convolution without last
    for (size_t i=0; i!=positions.size()-1; i++)
    {
        trueResult += masses[i]/(positions[i]-target).len();
    }
    ASSERT_NEAR_RELATIVE(result, trueResult, 1e-6);

    //ASSERT_EQ(conv.convolute(oct, Position(0.0, 0.0, 0.0), massSumVisitor), somePointsMass);
}

TEST_F(ConvolutionTests, ConvoluteOneScalingZone)
{
    addSomePoints();
    callsCounter = 0;
    ASSERT_EQ(conv.convolute(oct, Position(0.0, 0.0, 0.0), massSumVisitor), somePointsMass);
    int cc = callsCounter;
    scales.addScale(0.1, 1000);
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

    scales.addScale(5, 3);
    convField = conv.convolute(oct, p1, coulomb);
    ASSERT_NEAR(realField, convField, 1e-3*realField) << "Convolution result has bad precision";
    //cout << realField << " " << convField << endl;

    scales.addScale(7, 10);
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

    scales.addScale(5, 3);
    convField = conv.convolute(oct, p1, coulomb);
    ASSERT_NEAR(realField, convField, 1e-3*realField) << "Convolution result has bad precision";
    //cout << realField << " " << convField << endl;

    scales.addScale(7, 10);
    convField = conv.convolute(oct, p1, coulomb);
    ASSERT_NEAR(realField, convField, 3e-3*realField) << "Convolution result has bad precision";
    //cout << realField << " " << convField << endl;

    scales.addScale(10, 1000);
    convField = conv.convolute(oct, p1, coulomb);
    ASSERT_NEAR(realField, convField, 3e-3*realField) << "Convolution result has bad precision";
    //cout << realField << " " << convField << endl;
}


class ConvolutionTestsTempated : public ::testing::Test
{
public:
    void addManyPoints()
    {
        PointsGenerator::addGrid(10, 10, oct, &positions);
    }

    FullEField getCoulombFieldBruteForce(const Position& target)
    {
        FullEField result;
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
    DiscreteScales scales;
    Convolution<FullEField> conv{scales};

    Convolution<FullEField>::Visitor coulomb =
        [this](const Position& target, const Position& object, double mass)
        {
            FullEField result;
            double dist = (target-object).len();
            double dist3 = dist*dist*dist;
            double tmp = mass/dist3;
            result.potential = mass/dist;
            result.E[0] = (target - object).x[0]*tmp;
            result.E[1] = (target - object).x[1]*tmp;
            result.E[2] = (target - object).x[2]*tmp;
            return result;
        };
};


TEST_F(ConvolutionTestsTempated, ConvoluteCoulomb1)
{
    addManyPoints();
    Position p1 = {4.23, -3.45, -1.56};
    FullEField realField = getCoulombFieldBruteForce(p1);
    FullEField convField = conv.convolute(oct, p1, coulomb);

    ASSERT_NEAR_RELATIVE(realField.E[0], convField.E[0], 1e-3);
    ASSERT_NEAR_RELATIVE(realField.E[1], convField.E[1], 1e-3);
    ASSERT_NEAR_RELATIVE(realField.E[2], convField.E[2], 1e-3);
    ASSERT_NEAR_RELATIVE(realField.potential, convField.potential, 1e-3);

    scales.addScale(4, 3);
    convField = conv.convolute(oct, p1, coulomb);

    ASSERT_NEAR_RELATIVE(realField.E[0], convField.E[0], 5e-3);
    ASSERT_NEAR_RELATIVE(realField.E[1], convField.E[1], 5e-3);
    ASSERT_NEAR_RELATIVE(realField.E[2], convField.E[2], 5e-3);
    ASSERT_NEAR_RELATIVE(realField.potential, convField.potential, 5e-3);

    scales.addScale(7, 11);
    convField = conv.convolute(oct, p1, coulomb);

    ASSERT_NEAR_RELATIVE(realField.E[0], convField.E[0], 1e-2);
    ASSERT_NEAR_RELATIVE(realField.E[1], convField.E[1], 1e-2);
    ASSERT_NEAR_RELATIVE(realField.E[2], convField.E[2], 1e-2);
    ASSERT_NEAR_RELATIVE(realField.potential, convField.potential, 1e-2);

    p1 = {14.23, -23.45, -1.56};
    realField = getCoulombFieldBruteForce(p1);
    scales.addScale(10, 30);
    convField = conv.convolute(oct, p1, coulomb);
    ASSERT_NEAR_RELATIVE(realField.E[0], convField.E[0], 1e-1);
    ASSERT_NEAR_RELATIVE(realField.E[1], convField.E[1], 1e-1);
    ASSERT_NEAR_RELATIVE(realField.E[2], convField.E[2], 1e-1);
    ASSERT_NEAR_RELATIVE(realField.potential, convField.potential, 1e-1);
}
