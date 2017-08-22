# src-level CMakeLists.txt

# Include LibInt2 and the eigen3 headers
include_directories(/usr/local/libint/2.2.0/include)
include_directories(/usr/local/libint/2.2.0/include/libint2)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I/opt/local/include/eigen3")


# Add source files and header files to compile
add_executable(libint-eigen
        main.cpp
        libint-wrapper.cpp
        libint-wrapper.hpp)


# Link the LibInt2 library to the object files
target_link_libraries(libint-eigen /usr/local/libint/2.2.0/lib/libint2.a)