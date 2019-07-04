include_directories(${PROJECT_SOURCE_DIR}/lib ${PROJECT_BINARY_DIR}/lib)
link_libraries(libtinker_shared)

add_executable(dynamic dynamic.f)
add_executable(analyze analyze.f)
add_executable(archive archive.f)
add_executable(radial radial.f)
add_executable(diffuse diffuse.f)
add_executable(anneal anneal.f)
add_executable(bar bar.f)
add_executable(xyzpdb xyzpdb.f)
add_executable(pdbxyz pdbxyz.f)
add_executable(minimize minimize.f)

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    set(CMAKE_C_FLAGS_DEBUG "-O0 -g -Wall -Wextra")
    set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g -Wall -Wextra")
endif()


add_library(base memory.cpp atom.cpp common.cpp molecule.cpp frame.cpp 
                   trajectoryreader.cpp forcefield.cpp gro_writer.cpp
        xtc_writer.cpp trr_writer.cpp netcdf_writer.cpp
        Trajconv.cpp Distance.cpp CoordinateNumPerFrame.cpp FirstCoordExchangeSearch.cpp
        RadicalDistribtuionFunction.cpp RMSDCal.cpp RMSFCal.cpp HBond.cpp taskMenu.cpp
        ResidenceTime.cpp GreenKubo.cpp Cluster.cpp AminoTop.cpp NMRRange.cpp
        RotAcfCutoff.cpp DiffuseCutoff.cpp Diffuse.cpp DipoleAngle.cpp
        DipoleAngleSingleDistanceNormal.cpp DipoleAngleVolumeNormal.cpp
        DipoleAngle2Gibbs.cpp ShellDensity.cpp SearchInteractionResidue.cpp
        FindMinBetweenTwoGroups.cpp DemixIndexOfTwoGroup.cpp PrintTopolgy.cpp mainUtils.cpp
        TrajectoryWriterFactoryImpl.cpp PBCUtils.cpp DipoleAngleWithDistanceRange.cpp
        DipoleAxisDistribution.cpp EquatorialAngle.cpp PlaneAngle.cpp AngleWat.cpp
        RotAcf.cpp DistanceAngle.cpp Histrogram.cpp Histrogram2D.cpp)

add_executable(analysis analysis.cpp)
target_link_libraries(analysis base)


