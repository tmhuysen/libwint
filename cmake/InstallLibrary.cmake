# In this CMake file, we will provide the necessary commands to install the library

# Installation of the target library and its associated header files
#   The library (liblibint-wrapper.a) is installed in ${INSTALLATION_PREFIX}/lib
#   The header files are installed in ${INSTALLATION_PREFIX}/include/libint-wrapper

# 1. Install the library
#   The target of this project is a library called libint-wrapper. We'll also have to export this target, to be able to use it in other CMake-based projects.
#   To specify that this target should also be exported, we add the EXPORT option. This is used in conjuction with the install(EXPORT) command below
#   The ARCHIVE option specifies that we're working with a static library
install(TARGETS libint-wrapper
        EXPORT libint-wrapper ARCHIVE
        DESTINATION ${INSTALLATION_PREFIX}/lib)

# 2. Install the header files
install(DIRECTORY ../include/ DESTINATION ${INSTALLATION_PREFIX}/include/libint-wrapper)

# 3. Export the target library: this creates a file libint-wrapper.cmake, to be loaded by an outside project
install(EXPORT libint-wrapper DESTINATION ${INSTALLATION_PREFIX}/include/libint-wrapper)