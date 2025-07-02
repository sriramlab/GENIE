cmake_minimum_required(VERSION 3.5)

project(GENIE LANGUAGES CXX)

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

# Find required dependencies
find_package(Threads REQUIRED)

# Determine system architecture
if(CMAKE_SYSTEM_PROCESSOR MATCHES "^(x86_64|AMD64)")
    set(ARCHITECTURE "amd64")
elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "^(aarch64|arm64)")
    set(ARCHITECTURE "arm64")
else()
    set(ARCHITECTURE "unknown")
endif()

# Print the detected architecture
message(STATUS "Detected architecture: ${ARCHITECTURE}")

# Handle user preference for SSE support (default: OFF)
option(ENABLE_SSE "Enable SSE support" OFF)

# Set SSE_SUPPORT based on the ENABLE_SSE option and architecture
if (ENABLE_SSE)
    if (ARCHITECTURE STREQUAL "amd64")
        set(SSE_SUPPORT 1)
    else()
        message(WARNING "SSE support currently only available on amd64 platforms. Turning OFF SSE support.")
        set(SSE_SUPPORT 0)
    endif()
else()
    set(SSE_SUPPORT 0)
endif()

ADD_DEFINITIONS(-DSSE_SUPPORT=${SSE_SUPPORT})

# Check whether we need to link against the rt library
include(CheckCXXSourceCompiles)
check_cxx_source_compiles("
#include <time.h>
int main() {
    struct timespec ts;
    return clock_gettime(CLOCK_REALTIME, &ts);
}
" HAVE_CLOCK_GETTIME)

set(genie_lib_header
        include/auxillary.h
        include/arguments.h
        include/genotype.h
        include/genomult.h
        include/helper.h
        include/io.h
        include/mailman.h
        include/matmult.h
        include/printable.h
        include/std.h
        include/storage.h
		include/functions.h
		include/vectorfn.h
		include/statsfn.h
)

set(genie_lib_source
		src/auxillary.cpp
        src/genotype.cpp
        src/genomult.cpp
        src/matmult.cpp
        src/storage.cpp
		src/functions.cpp
)

add_library(genie_lib STATIC ${genie_lib_header} ${genie_lib_source})
target_include_directories(genie_lib INTERFACE include)
target_link_libraries(genie_lib PRIVATE Threads::Threads)
IF (NOT HAVE_CLOCK_GETTIME)
    target_link_libraries(genie_lib PRIVATE rt)
ENDIF ()

add_executable(GENIE_mem src/ge_mem_flexible.cpp)
add_executable(GENIE src/ge_flexible.cpp)

# Link executables against required dependencies
target_link_libraries(GENIE_mem PRIVATE genie_lib)
target_link_libraries(GENIE PRIVATE genie_lib)

# Handle user preference for SSE support (default: OFF)
option(BUILD_PYTHON_MODULE "Whether the python module is being built from the scikit-build-core context" OFF)
if (BUILD_PYTHON_MODULE)
    install(TARGETS GENIE_mem GENIE DESTINATION rhe)
endif ()
