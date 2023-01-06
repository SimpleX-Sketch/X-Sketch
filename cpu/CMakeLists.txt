cmake_minimum_required(VERSION 3.14)
project(cpu)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "-O2 -g")

include_directories(.)
include_directories(Common)

if(CMAKE_COMPILER_IS_GNUCXX)
    message(STATUS "GCC detected, adding compile flags")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pedantic -Wall -Wextra -Wredundant-decls")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse2 -mssse3 -msse4.1 -msse4.2 -mavx -mbmi -march=native")
    set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -g")
    ADD_DEFINITIONS(-DUNIX)
endif(CMAKE_COMPILER_IS_GNUCXX)

if (UNIX)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
else()
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++11")
endif(UNIX)


add_executable(cpu main.cc Common/BOBHash.h)
