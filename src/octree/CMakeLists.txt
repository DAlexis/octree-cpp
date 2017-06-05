cmake_minimum_required(VERSION 2.8)

project(octree)

set(LIB_SOURCE
    octree.cpp
    octree.hpp
)

set(${PROJECT_NAME}_USED_INCDIRS
    # includes here
)

include_directories(${${PROJECT_NAME}_USED_INCDIRS})

add_library(${PROJECT_NAME} STATIC ${LIB_SOURCE})


set_property(TARGET ${PROJECT_NAME} PROPERTY CXX_STANDARD 11)
set_property(TARGET ${PROJECT_NAME} PROPERTY CXX_STANDARD_REQUIRED ON)

export_include_dirs(
    ${PROJECT_SOURCE_DIR}
    ${${PROJECT_NAME}_USED_INCDIRS}
)

export_libs(
    ${PROJECT_NAME}
)
