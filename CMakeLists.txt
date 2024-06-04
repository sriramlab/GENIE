cmake_minimum_required(VERSION 3.2)

project(EM_C++ LANGUAGES CXX)

include_directories(include/)

# Set the CMake standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Set a default build type to Release (-O3 -DNDEBUG for GCC)
if (NOT CMAKE_BUILD_TYPE)
    message(STATUS "Setting build type to 'Release' as none was specified.")
    message(STATUS "  You can specify with the flag -DCMAKE_BUILD_TYPE=<Debug|Release|MinSizeRel|RelWithDebInfo>")
    message(STATUS "  'Release' will build optimised binaries, but 'Debug' may be better while developing.")
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
endif ()

message(STATUS "CMake build type is set to ${CMAKE_BUILD_TYPE}")

set(CMAKE_CXX_FLAGS "-lrt -lpthread ")

# Find required dependencies
find_package(Boost REQUIRED)
find_package(GSL REQUIRED)

IF (NOT DEFINED SSE_SUPPORT)
    SET(SSE_SUPPORT 0)
ENDIF()
ADD_DEFINITIONS(-DSSE_SUPPORT=${SSE_SUPPORT})

IF (NOT DEFINED DEBUG)
    SET(DEBUG 0)
ENDIF()
ADD_DEFINITIONS(-DDEBUG=${DEBUG})


add_executable(GENIE_mem src/ge_mem_flexible.cpp src/genotype.cpp src/matmult.cpp src/storage.cpp)
add_executable(GENIE src/ge_hetro_flexible.cpp src/genotype.cpp src/matmult.cpp src/storage.cpp)
add_executable(GENIE_multi_pheno src/ge_hetro_flexible_multi_pheno.cpp src/genotype.cpp src/matmult.cpp src/storage.cpp)

# Link executables against required dependencies
target_link_libraries(GENIE PRIVATE Boost::headers GSL::gsl)
target_link_libraries(GENIE_mem PRIVATE Boost::headers GSL::gsl)
target_link_libraries(GENIE_multi_pheno PRIVATE Boost::headers GSL::gsl)
