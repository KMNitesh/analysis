include_directories(${PROJECT_SOURCE_DIR}/src/analysis ${PROJECT_BINARY_DIR}/src)

find_package(GTest REQUIRED)
include_directories(SYSTEM ${GTEST_INCLUDE_DIRS})
include(GoogleTest)

file(GLOB GTEST_SRCS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cpp)

foreach (GTEST_TEST_CASE_SRC ${GTEST_SRCS})
    get_filename_component(BASE_NAME ${GTEST_TEST_CASE_SRC} NAME_WE)
    add_executable(${BASE_NAME} ${GTEST_TEST_CASE_SRC})
    target_link_libraries(${BASE_NAME} base ${GTEST_BOTH_LIBRARIES} gmock gmock_main pthread)
    gtest_discover_tests(${BASE_NAME} WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/test_files)
endforeach ()
