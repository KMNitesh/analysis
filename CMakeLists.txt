cmake_minimum_required(VERSION 3.15)
project(Tinker LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS ON)

set(CMAKE_SKIP_BUILD_RPATH OFF)
set(CMAKE_BUILD_WITH_INSTALL_RPATH ON)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH ON)

option(ENABLE_UNITTEST "Enable unit test cases" ON)
option(ENABLE_NATIVE "Enable native cpu instruction support" ON)

file(STRINGS "AUTHOR" CACANA_AUTHOR)
file(READ "ScriptSyntaxSummary" SCRIPT_SYNTAX_SUMMARY)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/externals/sanitizers-cmake/cmake" ${CMAKE_MODULE_PATH})
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/externals/findFFTW" ${CMAKE_MODULE_PATH})
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

if (CMAKE_CXX_COMPILER_ID MATCHES GNU)
    string(APPEND CMAKE_CXX_FLAGS_DEBUG " -Wall -Wextra")
    string(APPEND CMAKE_CXX_FLAGS_RELEASE " -ffast-math -flto")
    SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -flto")
    SET(CMAKE_AR  "gcc-ar")
    SET(CMAKE_RANLIB "gcc-ranlib")
    SET(CMAKE_CXX_ARCHIVE_CREATE "<CMAKE_AR> qcs <TARGET> <LINK_FLAGS> <OBJECTS>")
    SET(CMAKE_CXX_ARCHIVE_FINISH "<CMAKE_RANLIB> <TARGET>")
    if (ENABLE_NATIVE)
        string(APPEND CMAKE_CXX_FLAGS_RELEASE " -march=native")
    endif()
endif()

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

