
file(GLOB ANALYSISSRCS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.hpp *.cpp)
list(REMOVE_ITEM ANALYSISSRCS analysis.cpp)

add_library(base ${ANALYSISSRCS})
add_sanitizers(base)

add_executable(analysis analysis.cpp)
target_link_libraries(analysis base)
add_sanitizers(analysis)



