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
include_directories(
   "/public/home/xiamr/prog/gromacs-5.1.4/src"
   "/public/soft/software/netcdf-4.3.0/include"
    "lib"
)

link_directories(
   "/public/soft/software/gromacs-5.1.4/lib64"
)

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


# Boost
set(BOOST_INCLUDEDIR ${CPATH_DIRS} ${C_INCLUDE_DIRS} ${CPLUS_INCLUDE_DIRS} ${DEFAULT_INCLUDE_DIRS})
set(BOOST_LIBRARYDIR ${LIBRARY_DIRS})
find_package(Boost 1.67 REQUIRED COMPONENTS regex program_options)
include_directories(${Boost_INCLUDE_DIRS})
link_libraries(${Boost_LIBRARIES})


add_library(libtinker_shared SHARED ${LIBTINKER_SRC} )
set_target_properties(libtinker_shared PROPERTIES OUTPUT_NAME "tinker_shared")
target_link_libraries(libtinker_shared gromacs_mpi)

link_libraries(fftw3_threads fftw3 netcdf libtinker_shared fmt tbb readline tbbmalloc ncurses)
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


