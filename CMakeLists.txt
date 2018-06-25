cmake_minimum_required(VERSION 3.1)
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
   "/public/soft/software/boost-1.67.0/include"
    "lib"
)
link_directories(
   "/public/soft/software/gromacs-5.1.4/lib64"
   "/public/soft/software/boost-1.67.0/lib/"
)

file(GLOB LIBTINKER_SRC  lib/*.f lib/*.f90 lib/*.cpp)

ADD_LIBRARY(libtinker_shared SHARED ${LIBTINKER_SRC} )
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
target_link_libraries(analysis boost_program_options)
