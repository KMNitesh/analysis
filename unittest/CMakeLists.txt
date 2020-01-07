include_directories(${PROJECT_SOURCE_DIR}/src/analysis ${PROJECT_BINARY_DIR}/src)

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
include_directories(SYSTEM ${GTEST_INCLUDE_DIRS})
include(GoogleTest)

file(GLOB GTEST_SRCS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cpp)
list(REMOVE_ITEM GTEST_SRCS test_grammar.cpp test_trajectoryreader_tpr.cpp test_center_selection_grammar.cpp)

foreach (GTEST_TEST_CASE_SRC ${GTEST_SRCS})
    get_filename_component(BASE_NAME ${GTEST_TEST_CASE_SRC} NAME_WE)
    add_executable(${BASE_NAME} ${GTEST_TEST_CASE_SRC})
    target_link_libraries(${BASE_NAME} base ${GTEST_BOTH_LIBRARIES} gmock gmock_main pthread)
    gtest_discover_tests(${BASE_NAME})
endforeach ()
