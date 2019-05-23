include_directories(${PROJECT_SOURCE_DIR}/lib ${PROJECT_BINARY_DIR}/lib)
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

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    set(CMAKE_C_FLAGS_DEBUG "-O0 -g -Wall -Wextra")
    set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g -Wall -Wextra")
endif()


add_library(base memory.cpp atom.cpp common.cpp molecule.cpp frame.cpp 
                   trajectoryreader.cpp forcefield.cpp gro_writer.cpp
                   xtc_writer.cpp trr_writer.cpp netcdf_writer.cpp)

add_executable(analysis analysis.cpp)
target_link_libraries(analysis base)


