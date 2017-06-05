#!/bin/bash

set -e

echo "Usage: ./generate-eclipse-project.sh <debug|release>. Debug is default"

workspace="eclipse-workspace"
project="octree"
dir=$workspace/$project

mkdir -p $dir

cfg="Debug"
if [ "$1" == "release" ];
then
    cfg="Release"
    build_dir=$build_dir_prefix"/release"
fi

(
cd $dir
cmake -G"Eclipse CDT4 - Unix Makefiles" -D CMAKE_BUILD_TYPE=Debug ../../src

# Patching definitions for C++11 support
sed -i s/199711L/201103L/ .cproject
)
