# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.10

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /opt/local/bin/cmake

# The command to remove a file.
RM = /opt/local/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/Forks/libwint

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/Forks/libwint/bbb

# Include any dependencies generated for this target.
include tests/CMakeFiles/AOBasis_test.dir/depend.make

# Include the progress variables for this target.
include tests/CMakeFiles/AOBasis_test.dir/progress.make

# Include the compile flags for this target's objects.
include tests/CMakeFiles/AOBasis_test.dir/flags.make

tests/CMakeFiles/AOBasis_test.dir/AOBasis_test.cpp.o: tests/CMakeFiles/AOBasis_test.dir/flags.make
tests/CMakeFiles/AOBasis_test.dir/AOBasis_test.cpp.o: ../tests/AOBasis_test.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/Forks/libwint/bbb/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object tests/CMakeFiles/AOBasis_test.dir/AOBasis_test.cpp.o"
	cd /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/Forks/libwint/bbb/tests && /opt/local/bin/clang++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/AOBasis_test.dir/AOBasis_test.cpp.o -c /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/Forks/libwint/tests/AOBasis_test.cpp

tests/CMakeFiles/AOBasis_test.dir/AOBasis_test.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/AOBasis_test.dir/AOBasis_test.cpp.i"
	cd /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/Forks/libwint/bbb/tests && /opt/local/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/Forks/libwint/tests/AOBasis_test.cpp > CMakeFiles/AOBasis_test.dir/AOBasis_test.cpp.i

tests/CMakeFiles/AOBasis_test.dir/AOBasis_test.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/AOBasis_test.dir/AOBasis_test.cpp.s"
	cd /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/Forks/libwint/bbb/tests && /opt/local/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/Forks/libwint/tests/AOBasis_test.cpp -o CMakeFiles/AOBasis_test.dir/AOBasis_test.cpp.s

tests/CMakeFiles/AOBasis_test.dir/AOBasis_test.cpp.o.requires:

.PHONY : tests/CMakeFiles/AOBasis_test.dir/AOBasis_test.cpp.o.requires

tests/CMakeFiles/AOBasis_test.dir/AOBasis_test.cpp.o.provides: tests/CMakeFiles/AOBasis_test.dir/AOBasis_test.cpp.o.requires
	$(MAKE) -f tests/CMakeFiles/AOBasis_test.dir/build.make tests/CMakeFiles/AOBasis_test.dir/AOBasis_test.cpp.o.provides.build
.PHONY : tests/CMakeFiles/AOBasis_test.dir/AOBasis_test.cpp.o.provides

tests/CMakeFiles/AOBasis_test.dir/AOBasis_test.cpp.o.provides.build: tests/CMakeFiles/AOBasis_test.dir/AOBasis_test.cpp.o


# Object files for target AOBasis_test
AOBasis_test_OBJECTS = \
"CMakeFiles/AOBasis_test.dir/AOBasis_test.cpp.o"

# External object files for target AOBasis_test
AOBasis_test_EXTERNAL_OBJECTS =

tests/AOBasis_test: tests/CMakeFiles/AOBasis_test.dir/AOBasis_test.cpp.o
tests/AOBasis_test: tests/CMakeFiles/AOBasis_test.dir/build.make
tests/AOBasis_test: /usr/local/cpputil/lib/libcpputil.a
tests/AOBasis_test: src/libwint.a
tests/AOBasis_test: /usr/local/libint/2.3.1/lib/libint2.a
tests/AOBasis_test: tests/CMakeFiles/AOBasis_test.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/Forks/libwint/bbb/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable AOBasis_test"
	cd /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/Forks/libwint/bbb/tests && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/AOBasis_test.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
tests/CMakeFiles/AOBasis_test.dir/build: tests/AOBasis_test

.PHONY : tests/CMakeFiles/AOBasis_test.dir/build

tests/CMakeFiles/AOBasis_test.dir/requires: tests/CMakeFiles/AOBasis_test.dir/AOBasis_test.cpp.o.requires

.PHONY : tests/CMakeFiles/AOBasis_test.dir/requires

tests/CMakeFiles/AOBasis_test.dir/clean:
	cd /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/Forks/libwint/bbb/tests && $(CMAKE_COMMAND) -P CMakeFiles/AOBasis_test.dir/cmake_clean.cmake
.PHONY : tests/CMakeFiles/AOBasis_test.dir/clean

tests/CMakeFiles/AOBasis_test.dir/depend:
	cd /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/Forks/libwint/bbb && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/Forks/libwint /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/Forks/libwint/tests /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/Forks/libwint/bbb /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/Forks/libwint/bbb/tests /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/Forks/libwint/bbb/tests/CMakeFiles/AOBasis_test.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : tests/CMakeFiles/AOBasis_test.dir/depend
