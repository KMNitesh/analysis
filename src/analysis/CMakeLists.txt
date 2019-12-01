
file(GLOB_RECURSE ANALYSISSRCS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.hpp *.cpp)
list(REMOVE_ITEM ANALYSISSRCS program/analysis.cpp)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
add_library(base ${ANALYSISSRCS})
add_sanitizers(base)

target_link_libraries(base ${GROMACS_LIBS})

add_executable(analysis program/analysis.cpp)
target_link_libraries(analysis base)
add_sanitizers(analysis)



