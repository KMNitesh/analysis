include_directories(${PROJECT_SOURCE_DIR}/lib ${PROJECT_BINARY_DIR}/lib)
link_libraries(libtinker_shared)

add_subdirectory(tinker)
add_subdirectory(analysis)