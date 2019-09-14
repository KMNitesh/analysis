cmake_minimum_required(VERSION 3.13)
project(Tinker)

enable_language(Fortran C CXX)

set(CMAKE_CXX_STANDARD 17)

#set(CMAKE_CXX_COMPILER_LAUNCHER ccache)

option(ENABLE_UNITTEST "Enable unit test cases" ON)
option(ENABLE_NATIVE "Enable native cpu instruction support" ON)

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    set(CMAKE_CXX_FLAGS_DEBUG "-O0 -pipe -ggdb3 -Wall -Wextra")
    set(CMAKE_CXX_FLAGS_RELEASE "-DNDEBUG -DBOOST_DISABLE_ASSERTS -pipe -O3 -ggdb3 -ffast-math  -fopenmp")
    if (ENABLE_NATIVE)
        set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -march=native")
    endif ()
endif ()

if ("${CMAKE_C_COMPILER_ID}" STREQUAL "GNU")
    set(CMAKE_C_FLAGS_DEBUG "-O0 -pipe -ggdb3 -Wall -Wextra")
    set(CMAKE_C_FLAGS_RELEASE "-DNDEBUG -pipe -O3 -ggdb3 -ffast-math -fopenmp")
    if (ENABLE_NATIVE)
        set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -march=native")
    endif ()
endif ()

if ("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU")
    set(CMAKE_Fotran_FLAGS_DEBUG "-O0 -pipe -ggdb3 -Wall -Wextra")
    set(CMAKE_Fortran_FLAGS_RELEASE "-DNDEBUG -pipe -O3 -ffast-math -fopenmp")
    if (ENABLE_NATIVE)
        set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -march=native")
    endif ()
endif ()

set(CMAKE_VERBOSE_MAKEFILE ON)

# Set global property (all targets are impacted)
#set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE "${CMAKE_COMMAND} -E time")


if (NOT GROMACS_SRC_ROOT)
    set(GROMACS_SRC_ROOT "/public/home/xiamr/prog/gromacs-5.1.4")
endif ()

if (NOT EXISTS ${GROMACS_SRC_ROOT}/src)
    message(FATAL_ERROR "gromacs 5.1.4 source root not found, please set GROMACS_SRC_ROOT")
endif ()

message("gromacs 5.1.4 source root : " ${GROMACS_SRC_ROOT})

include_directories(${GROMACS_SRC_ROOT}/src)
include_directories(lib)


if (DEFINED ENV{LIBRARY_PATH})
    string(REPLACE ":" ";" LIBRARY_DIRS $ENV{LIBRARY_PATH})
endif ()
if (DEFINED ENV{CPATH})
    string(REPLACE ":" ";" CPATH_DIRS $ENV{CPATH})
endif ()
if (DEFINED ENV{C_INCLUDE_PATH})
    string(REPLACE ":" ";" C_INCLUDE_DIRS $ENV{C_INCLUDE_PATH})
endif ()
if (DEFINED ENV{CPLUS_INCLUDE_PATH})
    string(REPLACE ":" ";" CPLUS_INCLUDE_DIRS $ENV{CPLUS_INCLUDE_PATH})
endif ()
set(DEFAULT_INCLUDE_DIRS "/usr/include;/usr/local/include")

set(INCLUES_PATH "${CPATH_DIRS} ${C_INCLUDE_DIRS} ${CPLUS_INCLUDE_DIRS} ${DEFAULT_INCLUDE_DIRS}")

# netcdf
find_path(NETCDF_INCLUDES netcdf.h HITS ${INCLUES_PATH})
if (NOT NETCDF_INCLUDES)
    message(FATAL_ERROR "netcdf.h not found")
endif ()
include_directories(${NETCDF_INCLUDES})

find_library(NETCDF_LIBS netcdf PATHS ${LIBRARY_DIRS})
if (NOT NETCDF_LIBS)
    message(FATAL_ERROR "netcdf library not found")
endif ()
link_libraries(${NETCDF_LIBS})


# fftw3
find_path(FFTW3_INCLUDES fftw3.h HITS ${INCLUES_PATH})
if (NOT FFTW3_INCLUDES)
    message(FATAL_ERROR "fftw3.h not found")
endif ()
include_directories(${FFTW3_INCLUDES})

find_library(FFTW3_threads_LIBS fftw3_threads HITS ${LIBRARY_DIRS})
if (NOT FFTW3_threads_LIBS)
    message(FATAL_ERROR "fftw3_threads library not found")
endif ()
link_libraries(fftw3 ${FFTW3_threads_LIBS})


# intel tbb
find_path(TBB_INCLUDES "tbb/tbb.h" HITS ${INCLUES_PATH})
if (NOT TBB_INCLUDES)
    message(FATAL_ERROR "tbb header not found")
endif ()
include_directories(${TBB_INCLUDES})

if (${CMAKE_BUILD_TYPE} STREQUAL "Debug")
    set(TBB_LIB_NAME tbb_debug)
    set(TBBMALLOC_LIB_NAME tbbmalloc_debug)
else ()
    set(TBB_LIB_NAME tbb)
    set(TBBMALLOC_LIB_NAME tbbmalloc)
endif ()

find_library(TBB_LIBS ${TBB_LIB_NAME} PATHS ${LIBRARY_DIRS})
if (NOT TBB_LIBS)
    message(FATAL_ERROR "tbb library not found")
endif ()
link_libraries(${TBB_LIBS})

find_library(TBBMALLOC_LIBS ${TBBMALLOC_LIB_NAME} PATHS ${LIBRARY_DIRS})
if (NOT TBBMALLOC_LIBS)
    message(FATAL_ERROR "tbbmalloc library not found")
endif ()
link_libraries(${TBBMALLOC_LIBS})

# Boost
set(BOOST_INCLUDEDIR ${INCLUES_PATH})
set(BOOST_LIBRARYDIR ${LIBRARY_DIRS})
find_package(Boost 1.67 REQUIRED COMPONENTS regex program_options filesystem container)
include_directories(${Boost_INCLUDE_DIRS})
link_libraries(${Boost_LIBRARIES})

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/externals/sanitizers-cmake/cmake" ${CMAKE_MODULE_PATH})
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

# Google Sanitizer
find_package(Sanitizers)

# gromacs
find_library(GROMACS_LIBS NAMES gromacs gromacs_mpi PATHS ${LIBRARY_DIRS})
if (NOT GROMACS_LIBS)
    message(FATAL_ERROR "gromacs library not found")
endif ()

# Eigen3
find_package(Eigen3 3.3 REQUIRED)
link_libraries(Eigen3::Eigen)
message(STATUS "FOUND Eigen3: ${EIGEN3_INCLUDE_DIR}")

enable_testing()

file(REMOVE_RECURSE ${PROJECT_BINARY_DIR}/*.gcda)

include_directories(include)

add_subdirectory(lib)
add_subdirectory(src)
if (ENABLE_UNITTEST)
    add_subdirectory(unittest)
endif ()

