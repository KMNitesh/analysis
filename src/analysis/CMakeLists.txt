
find_package(Threads REQUIRED)

find_package(Boost 1.72 REQUIRED CONFIG COMPONENTS regex program_options filesystem iostreams)

find_package(TBB 2020.1 REQUIRED COMPONENTS tbb tbbmalloc)
target_compile_definitions(TBB::tbb INTERFACE TBB_SUPPRESS_DEPRECATED_MESSAGES)

find_package(pugixml REQUIRED)

find_package(Eigen3 3.3 REQUIRED)

find_package(Readline REQUIRED)

find_package(nlohmann_json 3.7.3 REQUIRED)

file(GLOB_RECURSE ANALYSISSRCS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} CONFIGURE_DEPENDS *.hpp *.cpp)
list(REMOVE_ITEM ANALYSISSRCS program/analysis.cpp)

add_library(base STATIC ${ANALYSISSRCS})

target_include_directories(base
    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}
    PUBLIC ${PROJECT_BINARY_DIR}
    )

target_include_directories(base SYSTEM
    PUBLIC ${GROMACS_INCLUDE}
    )

target_link_libraries(base
    PUBLIC Boost::regex Boost::program_options Boost::filesystem Boost::iostreams
    PUBLIC TBB::tbb TBB::tbbmalloc
    PUBLIC pugixml-shared
    PUBLIC Eigen3::Eigen
    PUBLIC nlohmann_json::nlohmann_json
    PUBLIC Readline::Readline
    PRIVATE NetCDF::NetCDF
    PRIVATE ${GROMACS_LIBRARIES}
    PRIVATE ${CMAKE_DL_LIBS}
    PRIVATE Threads::Threads
    )

target_link_directories(base
    PUBLIC ${GROMACS_LIBRARY_DIRS}
    )

target_compile_definitions(base
    PRIVATE $<$<CONFIG:Release>:BOOST_DISABLE_ASSERTS>
    )

add_executable(analysis program/analysis.cpp)
target_link_libraries(analysis PRIVATE base)
add_sanitizers(analysis)



