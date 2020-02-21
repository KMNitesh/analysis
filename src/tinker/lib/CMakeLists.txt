# OpenMP
find_package(OpenMP REQUIRED)
# fftw3
find_package(FFTW REQUIRED COMPONENTS DOUBLE_THREADS_LIB)

file(GLOB LIBTINKER_SRC CONFIGURE_DEPENDS *.f *.f90 *.cpp)

add_library(libtinker_shared SHARED ${LIBTINKER_SRC})
set_property(TARGET libtinker_shared PROPERTY OUTPUT_NAME "tinker-shared")

target_include_directories(libtinker_shared
    PUBLIC ${CMAKE_CURRENT_BINARY_DIR}
)

target_include_directories(libtinker_shared SYSTEM
    PRIVATE ${GROMACS_INCLUDE}
)

target_link_libraries(libtinker_shared
    PRIVATE FFTW::DoubleThreads
    PUBLIC OpenMP::OpenMP_Fortran
    PRIVATE NetCDF::NetCDF
    PRIVATE ${GROMACS_LIBRARIES}
)

