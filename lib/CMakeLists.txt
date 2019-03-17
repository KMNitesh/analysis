file(GLOB LIBTINKER_SRC  *.f *.f90 *.cpp)
add_library(libtinker_shared SHARED ${LIBTINKER_SRC} )
set_target_properties(libtinker_shared PROPERTIES OUTPUT_NAME "tinker_shared")
target_link_libraries(libtinker_shared ${GROMACS_LIBS})