# test-level CMakeLists.txt

# Set-up boost testing
find_package(Boost COMPONENTS unit_test_framework REQUIRED)
include_directories(${Boost_INCLUDE_DIRS}
        ../src)


# Include LibInt2 and the eigen3 headers
include_directories(/usr/local/libint/2.2.0/include)
include_directories(/usr/local/libint/2.2.0/include/libint2)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I/opt/local/include/eigen3")


# Add source and header files to compile
add_executable(some_test
        some_test.cpp)


# Link the LibInt2 library to the object files
target_link_libraries(some_test /usr/local/libint/2.2.0/lib/libint2.a)


# Specifying test options
add_test(NAME test_name COMMAND test_command)