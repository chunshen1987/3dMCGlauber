cmake_minimum_required (VERSION 2.6)
project (3DMCGlb)

#for MacOSX
cmake_policy(SET CMP0042 NEW)

set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_HOME_DIRECTORY}")
set (CMAKE_INSTALL_PREFIX "${CMAKE_HOME_DIRECTORY}")
set (CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
message("prefix = " ${CMAKE_INSTALL_PREFIX})
message("rpath = " ${CMAKE_INSTALL_RPATH})

option (test "Build Unit tests" OFF)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O3 -Wall -std=c++11 -fsanitize=address")


if (test)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DDOCTEST_CONFIG_IMPLEMENT_WITH_MAIN")
else (test)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DDOCTEST_CONFIG_DISABLE -DNDEBUG")
endif (test)

add_subdirectory (src)
