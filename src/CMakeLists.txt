# src-level CMakeLists.txt


add_library(libint-eigen-wrapper libint-wrapper.cpp)
target_include_directories(libint-eigen-wrapper PRIVATE ../include)  # for the target 'libint-eigen'

target_include_directories(libint-eigen-wrapper PRIVATE /usr/local/libint/2.2.0/include)
target_include_directories(libint-eigen-wrapper PRIVATE /usr/local/libint/2.2.0/include/libint2)        # FIXME: test if we can leave this open
target_include_directories(libint-eigen-wrapper PRIVATE /opt/local/include/eigen3)
target_link_libraries(libint-eigen-wrapper /usr/local/libint/2.2.0/lib/libint2.a)