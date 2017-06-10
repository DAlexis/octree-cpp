/*
 * test-utils.hpp
 *
 *  Created on: 9 июн. 2017 г.
 *      Author: dalexies
 */

#ifndef UNIT_TESTS_TEST_UTILS_HPP_
#define UNIT_TESTS_TEST_UTILS_HPP_

#include <chrono>

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

#endif /* UNIT_TESTS_TEST_UTILS_HPP_ */
