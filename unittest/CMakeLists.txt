
set(CMAKE_FIND_PACKAGE_PREFER_CONFIG TRUE)
find_package(GTest 1.10.0 REQUIRED)

file(GLOB GTEST_SRCS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cpp)

include(GoogleTest)
foreach (GTEST_TEST_CASE_SRC ${GTEST_SRCS})
    get_filename_component(BASE_NAME ${GTEST_TEST_CASE_SRC} NAME_WE)
    add_executable(${BASE_NAME} ${GTEST_TEST_CASE_SRC})
    target_link_libraries(${BASE_NAME} PRIVATE base GTest::gtest GTest::gtest_main GTest::gmock GTest::gmock_main)
    gtest_discover_tests(${BASE_NAME} WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/test_files)
    add_sanitizers(${BASE_NAME})
endforeach ()
