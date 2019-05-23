include_directories(${PROJECT_SOURCE_DIR}/src ${PROJECT_BINARY_DIR}/src)

add_executable(test_grammar test_grammar.cpp)
target_link_libraries(test_grammar base)

add_executable(test_trajectoryreader_tpr test_trajectoryreader_tpr.cpp)
target_link_libraries(test_trajectoryreader_tpr base)

add_test(NAME "TestGrammar" COMMAND test_grammar)
add_test(NAME "TestTrajectoryReaderTPR" COMMAND test_trajectoryreader_tpr WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})