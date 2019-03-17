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
add_executable(analysis analysis.cpp memory.cpp atom.cpp common.cpp
        molecule.cpp frame.cpp trajectoryreader.cpp forcefield.cpp)


add_executable(test_grammar test_grammar.cpp atom.cpp common.cpp)
add_test(NAME "TestGrammar" COMMAND test_grammar)
