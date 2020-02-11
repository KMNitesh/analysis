# OpenMP
find_package(OpenMP REQUIRED)
# fftw3
find_package(FFTW REQUIRED COMPONENTS DOUBLE_THREADS_LIB)

file(GLOB LIBTINKER_SRC *.f *.f90 *.cpp)

add_library(libtinker_shared SHARED ${LIBTINKER_SRC})
set_property(TARGET libtinker_shared PROPERTY OUTPUTNAME "libtinker_shared")

target_include_directories(libtinker_shared
        PUBLIC ${CMAKE_CURRENT_BINARY_DIR}
        )

target_include_directories(libtinker_shared SYSTEM
        PUBLIC ${GROMACS_INCLUDE}
        )

target_link_libraries(libtinker_shared
        PUBLIC FFTW::DoubleThreads
        PUBLIC OpenMP::OpenMP_Fortran
        PRIVATE NetCDF::NetCDF
        PRIVATE ${GROMACS_LIBRARIES}
        )

target_compile_options(libtinker_shared
        PUBLIC $<$<AND:$<CONFIG:Debug>,$<OR:$<COMPILE_LANG_AND_ID:Fortran,GNU>,$<COMPILE_LANG_AND_ID:CXX,GNU>>>:-Wall -Wextra>
        PUBLIC $<$<AND:$<CONFIG:Release>,$<OR:$<COMPILE_LANG_AND_ID:Fortran,GNU>,$<COMPILE_LANG_AND_ID:CXX,GNU>>>:-ffast-math>
        PUBLIC $<$<AND:$<CONFIG:Release>,$<OR:$<COMPILE_LANG_AND_ID:Fortran,GNU>,$<COMPILE_LANG_AND_ID:CXX,GNU>>,$<BOOL:${ENABLE_NATIVE}>>:-march=native>
        )
