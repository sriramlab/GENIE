cmake_minimum_required(VERSION 3.12...3.30)
project(FLEX LANGUAGES C CXX)


include_directories(include/)

set(CMAKE_CXX_FLAGS "-std=c++0x -O3 -lrt -lpthread ")

IF (NOT DEFINED SSE_SUPPORT)
    SET(SSE_SUPPORT 0)
ENDIF()
ADD_DEFINITIONS(-DSSE_SUPPORT=${SSE_SUPPORT})

IF (NOT DEFINED DEBUG)
    SET(DEBUG 0)
ENDIF()
ADD_DEFINITIONS(-DDEBUG=${DEBUG})


add_executable(FLEX_h2 src/flex_h2.cpp src/genotype.cpp src/matmult.cpp src/storage.cpp)
