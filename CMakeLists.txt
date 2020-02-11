cmake_minimum_required(VERSION 3.15)
project(Tinker LANGUAGES Fortran CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS ON)

set(CMAKE_VERBOSE_MAKEFILE ON)

option(ENABLE_UNITTEST "Enable unit test cases" ON)
option(ENABLE_NATIVE "Enable native cpu instruction support" ON)

file(STRINGS "AUTHOR" CACANA_AUTHOR)
file(READ "ScriptSyntaxSummary" SCRIPT_SYNTAX_SUMMARY)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/externals/sanitizers-cmake/cmake" ${CMAKE_MODULE_PATH})
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/externals/findFFTW" ${CMAKE_MODULE_PATH})
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

# gromacs
find_package(PkgConfig REQUIRED)
pkg_search_module(GROMACS REQUIRED libgromacs>=5.1.4 libgromacs_mpi>=5.1.4)
list(GET GROMACS_INCLUDE_DIRS 0 GROMACS_INCLUDE)

# NetCDF
find_package(NetCDF REQUIRED)

# Google Sanitizer
find_package(Sanitizers)

configure_file(config.h.in config.h)

add_subdirectory(src)
if (ENABLE_UNITTEST)
    include(CTest)
    add_subdirectory(unittest)
endif ()

