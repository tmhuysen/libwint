# test-level CMakeLists.txt


# For every *.cpp-file (every file corresponding to one set of tests) ...
file(GLOB test_sources ./*.cpp)

foreach(test_source ${test_sources})
    # ... extract the filename without extension as a name for our executable ...
    get_filename_component(test_name ${test_source} NAME_WE)

    # ... add an executable based on the test source ...
    add_executable(${test_name} ${test_source})

    # ... add the boost headers (for testing) ...
    target_include_directories(${test_name} PRIVATE ${Boost_INCLUDE_DIRS})

    # ... add this project's headers ...
    target_include_directories(${test_name} PRIVATE ../include)

    # ... add the eigen3 headers ...
    target_include_directories(${test_name} PRIVATE /opt/local/include/eigen3)

    # ... add the libint2 headers ...
    target_include_directories(${test_name} PRIVATE /usr/local/libint/2.2.0/include)
    target_include_directories(${test_name} PRIVATE /usr/local/libint/2.2.0/include/libint2)

    # ... link to this project's library ...
    target_link_libraries(${test_name} PRIVATE libint-eigen)

    # ... and finally make a test ...
    add_test(NAME ${test_name} COMMAND ${test_name})
endforeach()
