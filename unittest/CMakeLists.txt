include_directories(${PROJECT_SOURCE_DIR}/src ${PROJECT_BINARY_DIR}/src)

add_executable(test_grammar test_grammar.cpp)
target_link_libraries(test_grammar base)

add_executable(test_trajectoryreader_tpr test_trajectoryreader_tpr.cpp)
target_link_libraries(test_trajectoryreader_tpr base)

add_executable(test_center_selection_grammar test_center_selection_grammar.cpp)
target_link_libraries(test_center_selection_grammar base)

add_test(NAME "TestGrammar" COMMAND test_grammar)
add_test(NAME "TestCenterRuleGrammer" COMMAND test_center_selection_grammar)
add_test(NAME "TestTrajectoryReaderTPR" COMMAND test_trajectoryreader_tpr WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})


find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})
include(GoogleTest)
add_executable(test_make_program_options test_make_program_options.cpp)
target_link_libraries(test_make_program_options base ${GTEST_BOTH_LIBRARIES})
gtest_discover_tests(test_make_program_options)
