/*
 * test-utils.hpp
 *
 *  Created on: 9 июн. 2017 г.
 *      Author: dalexies
 */

#ifndef UNIT_TESTS_TEST_UTILS_HPP_
#define UNIT_TESTS_TEST_UTILS_HPP_

#include "octree.hpp"
#include <chrono>
#include <vector>
#include <cmath>

using namespace octree;

struct TimeMesurment
{
    template<typename F, typename ...Args>
    static std::chrono::microseconds::rep execution(F&& func, Args&&... args)
    {
        auto start = std::chrono::steady_clock::now();
        std::forward<decltype(func)>(func)(std::forward<Args>(args)...);
        auto duration = std::chrono::duration_cast< std::chrono::microseconds >
                            (std::chrono::steady_clock::now() - start);
        return duration.count();
    }
};

class PointsGenerator
{
public:
    static void addGrid(int n, double size, Octree& oct, std::vector<Position>* positions);
    static Position& findNearestBruteForce(const Position& pos, std::vector<Position>& positions);
};

class FullEField
{
public:
    double E[3] = {0.0, 0.0, 0.0};
    double potential = 0.0;

    FullEField& operator+=(const FullEField& right)
    {
        E[0] += right.E[0];
        E[1] += right.E[1];
        E[2] += right.E[2];
        potential += right.potential;
        return *this;
    }
};

#define ASSERT_NEAR_RELATIVE(a, b, err)    ASSERT_NEAR((a), (b), fabs(err*a))

#endif /* UNIT_TESTS_TEST_UTILS_HPP_ */
