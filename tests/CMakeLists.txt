# test-level CMakeLists.txt

add_executable(StuffedMolecule_test StuffedMolecule_test.cpp)

# For now, I still have include everything from LibInt2 and eigen3, but in later versions I will use a more OOP-approach so that this isn't necessary.
target_include_directories(StuffedMolecule_test PRIVATE /usr/local/libint/2.2.0/include)
target_include_directories(StuffedMolecule_test PRIVATE /usr/local/libint/2.2.0/include/libint2)
target_include_directories(StuffedMolecule_test PRIVATE /opt/local/include/eigen3)

# Include the Boost headers
target_include_directories(StuffedMolecule_test PRIVATE ${Boost_INCLUDE_DIRS})

# Include this project's headers
target_include_directories(StuffedMolecule_test PRIVATE ../include)

# Link with this project's library
target_link_libraries(StuffedMolecule_test PRIVATE libint-eigen-wrapper)

add_test(NAME StuffedMolecule_test COMMAND StuffedMolecule_test)