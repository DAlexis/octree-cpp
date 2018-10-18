cmake_minimum_required(VERSION 2.8)

project(octree-tests)

include_directories (
    ${GTEST_INCLUDE_DIR}
    ${octree_INCLUDE_DIRS}
)

set(EXE_SOURCES
    octree-tests.cpp
    conv-tests.cpp
    test-utils.cpp
    test-utils.hpp
)

add_executable(${PROJECT_NAME} ${EXE_SOURCES})

target_link_libraries (${PROJECT_NAME}
    gtest
    gtest_main
    ${octree_LIBRARY}
)

add_test(NAME octree_tests
         COMMAND ${PROJECT_NAME})
