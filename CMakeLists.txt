cmake_minimum_required (VERSION 2.6)
project (3DMCGlb)

#for MacOSX
cmake_policy(SET CMP0042 NEW)

set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_HOME_DIRECTORY}")

option (test "Build Unit tests" OFF)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O3 -Wall -std=c++11 -DDOCTEST_CONFIG_DISABLE -DNDEBUG")

if (test)
    set(CMAKE_CXX_FLAGS "-g -Wall -std=c++11 -DDOCTEST_CONFIG_IMPLEMENT_WITH_MAIN")
endif()

add_subdirectory (src)
