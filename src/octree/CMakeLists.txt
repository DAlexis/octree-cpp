cmake_minimum_required(VERSION 2.8)

project(octree)

set(LIB_SOURCE
    octree.cpp
    octree.hpp
    geom-vector.hpp
)

add_library(${PROJECT_NAME} STATIC ${LIB_SOURCE})

target_include_directories(${PROJECT_NAME} PUBLIC ${PROJECT_SOURCE_DIR})

set_property(TARGET ${PROJECT_NAME} PROPERTY CXX_STANDARD 11)
set_property(TARGET ${PROJECT_NAME} PROPERTY CXX_STANDARD_REQUIRED ON)
