#include "octree.hpp"
#include <iostream>

void func()
{
    auto f = []() {
        std::cout << "Hello, cmake fan" << std::endl;
    };
    
    f();
}

double sqr(double x)
{
    return x*x;
}
