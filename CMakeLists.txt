cmake_minimum_required(VERSION 3.13)
project(Tinker)
enable_language(Fortran C CXX)

set(CMAKE_CXX_STANDARD 14)

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
    set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g -Wall -Wextra")
    set(CMAKE_CXX_FLAGS_RELEASE "-DNDEBUG -DINTEL -O3 -no-prec-div -xHost -axCORE-AVX2 -qopenmp")
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    set(CMAKE_CXX_FLAGS_DEBUG "-O0 -ggdb3 -Wall -Wextra")
    set(CMAKE_CXX_FLAGS_RELEASE "-DNDEBUG -O3 -ffast-math -march=native -fopenmp")
endif()

if ("${CMAKE_C_COMPILER_ID}" STREQUAL "Intel")
    set(CMAKE_C_FLAGS_DEBUG "-O0 -g -Wall -Wextra")
    set(CMAKE_C_FLAGS_RELEASE "-DNDEBUG -DINTEL -O3 -no-prec-div -xHost -axCORE-AVX2 -qopenmp")
elseif("${CMAKE_C_COMPILER_ID}" STREQUAL "GNU")
    set(CMAKE_C_FLAGS_DEBUG "-O0 -ggdb3 -Wall -Wextra")
    set(CMAKE_C_FLAGS_RELEASE "-DNDEBUG -O3 -ffast-math -march=native -fopenmp")
endif()

if ("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Intel")
    set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -g -Wall -Wextra")
    set(CMAKE_Fortran_FLAGS_RELEASE "-DNDEBUG -DINTEL -O3 -no-prec-div -xHost -axCORE-AVX2 -qopenmp")
elseif("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU")
    set(CMAKE_Fotran_FLAGS_DEBUG "-O0 -ggdb3 -Wall -Wextra")
    set(CMAKE_Fortran_FLAGS_RELEASE "-DNDEBUG -O3 -ffast-math -march=native -fopenmp")
endif()


set(CMAKE_VERBOSE_MAKEFILE ON)

if (NOT GROMACS_SRC_ROOT)
    set(GROMACS_SRC_ROOT "/public/home/xiamr/prog/gromacs-5.1.4")
endif()

if (NOT EXISTS ${GROMACS_SRC_ROOT}/src)
    message(FATAL_ERROR "gromacs 5.1.4 source root not found, please set GROMACS_SRC_ROOT" )
endif()

message("gromacs 5.1.4 source root : " ${GROMACS_SRC_ROOT})

include_directories(${GROMACS_SRC_ROOT}/src)
include_directories(lib)

file(GLOB LIBTINKER_SRC  lib/*.f lib/*.f90 lib/*.cpp)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

string(REPLACE ":" ";" LIBRARY_DIRS $ENV{LIBRARY_PATH})
string(REPLACE ":" ";" CPATH_DIRS $ENV{CPATH})
if ($ENV{C_INCLUDE_PATH})
    string(REPLACE ":" ";" C_INCLUDE_DIRS $ENV{C_INCLUDE_PATH})
endif()
if ($ENV{CPLUS_INCLUDE_PATH})
    string(REPLACE ":" ";" CPLUS_INCLUDE_DIRS $ENV{CPLUS_INCLUDE_PATH})
endif ()
set(DEFAULT_INCLUDE_DIRS "/usr/include;/usr/local/include")

set(INCLUES_PATH "${CPATH_DIRS} ${C_INCLUDE_DIRS} ${CPLUS_INCLUDE_DIRS} ${DEFAULT_INCLUDE_DIRS}")

# netcdf
find_path(NETCDF_INCLUDES netcdf.h PATHS ${INCLUES_PATH})
if(NOT NETCDF_INCLUDES)
    message(FATAL_ERROR "netcdf.h not found")
endif()
include_directories(${NETCDF_INCLUDES})

find_library(NETCDF_LIBS netcdf PATHS ${LIBRARY_DIRS})
if(NOT NETCDF_LIBS)
    message(FATAL_ERROR "netcdf library not found")
endif()
link_libraries(${NETCDF_LIBS})


# fmt
find_path(FMT_INCLUDES "fmt/printf.h" PATHS ${INCLUES_PATH})
if(NOT FMT_INCLUDES)
    message(FATAL_ERROR "fmt header not found")
endif()
include_directories(${FMT_INCLUDES})

find_library(FMT_LIBS fmt PATHS ${LIBRARY_DIRS})
if(NOT FMT_LIBS)
    message(FATAL_ERROR "fmt library not found")
endif()
link_libraries(${FMT_LIBS})

# fftw3
find_path(FFTW3_INCLUDES fftw3.h PATHS ${INCLUES_PATH})
if(NOT FFTW3_INCLUDES)
    message(FATAL_ERROR "fftw3.h not found")
endif()
include_directories(${FFTW3_INCLUDES})

find_library(FFTW3_threads_LIBS fftw3_threads PATHS ${LIBRARY_DIRS})
if(NOT FFTW3_threads_LIBS)
    message(FATAL_ERROR "fftw3_threads library not found")
endif()
link_libraries(fftw3 ${FFTW3_threads_LIBS})

# readline
find_path(READLINE_INCLUDES "readline/readline.h" PATHS ${INCLUES_PATH})
if(NOT READLINE_INCLUDES)
    message(FATAL_ERROR "readline header not found")
endif()
include_directories(${READLINE_INCLUDES})

find_library(READLINE_LIBS readline PATHS ${LIBRARY_DIRS})
if(NOT READLINE_LIBS)
    message(FATAL_ERROR "readline library not found")
endif()
link_libraries(${READLINE_LIBS})


# intel tbb
find_path(TBB_INCLUDES "tbb/tbb.h" PATHS ${INCLUES_PATH})
if(NOT TBB_INCLUDES)
    message(FATAL_ERROR "tbb header not found")
endif()
include_directories(${TBB_INCLUDES})

find_library(TBB_LIBS tbb PATHS ${LIBRARY_DIRS})
if(NOT TBB_LIBS)
    message(FATAL_ERROR "tbb library not found")
endif()
link_libraries(${TBB_LIBS})

find_library(TBBMALLOC_LIBS tbbmalloc PATHS ${LIBRARY_DIRS})
if(NOT TBBMALLOC_LIBS)
    message(FATAL_ERROR "tbbmalloc library not found")
endif()
link_libraries(${TBBMALLOC_LIBS})

# Boost
set(BOOST_INCLUDEDIR ${INCLUES_PATH})
set(BOOST_LIBRARYDIR ${LIBRARY_DIRS})
find_package(Boost 1.67 REQUIRED COMPONENTS regex program_options)
include_directories(${Boost_INCLUDE_DIRS})
link_libraries(${Boost_LIBRARIES})

# gromacs
find_library(GROMACS_LIBS NAMES gromacs gromacs_mpi PATHS ${LIBRARY_DIRS})
if(NOT GROMACS_LIBS)
    message(FATAL_ERROR "gromacs library not found")
endif()

add_library(libtinker_shared SHARED ${LIBTINKER_SRC} )
set_target_properties(libtinker_shared PROPERTIES OUTPUT_NAME "tinker_shared")
target_link_libraries(libtinker_shared ${GROMACS_LIBS})

link_libraries(libtinker_shared)
add_executable(dynamic dynamic.f)
add_executable(analyze analyze.f)
add_executable(archive archive.f)
add_executable(radial radial.f)
add_executable(diffuse diffuse.f)
add_executable(anneal anneal.f)
add_executable(bar bar.f)
add_executable(xyzpdb xyzpdb.f)
add_executable(pdbxyz pdbxyz.f)
add_executable(minimize minimize.f)
add_executable(analysis analysis.cpp)

