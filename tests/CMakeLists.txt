# test-level CMakeLists.txt

add_executable(one_two_test one_two_test.cpp)

# For now, I still have include everything from LibInt2 and eigen3, but in later versions I will use a more OOP-approach so that this isn't necessary.
target_include_directories(one_two_test PRIVATE /usr/local/libint/2.2.0/include)
target_include_directories(one_two_test PRIVATE /usr/local/libint/2.2.0/include/libint2)
target_include_directories(one_two_test PRIVATE /opt/local/include/eigen3)

# Include the Boost headers
target_include_directories(one_two_test PRIVATE ${Boost_INCLUDE_DIRS})

# Include this project's headers
target_include_directories(one_two_test PRIVATE ../include)

# Link with this project's library
target_link_libraries(one_two_test PRIVATE libint-eigen-wrapper)

add_test(NAME one_two_test COMMAND one_two_test)