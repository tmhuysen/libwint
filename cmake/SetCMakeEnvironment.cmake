# In this CMake file, all CMake variables will be set


# Parse the project name into uppercase and lowercase.
string(TOUPPER ${PROJECT_NAME} PROJECT_NAME_UPPERCASE)  # Uppercase is needed in version.hpp.in
string(TOLOWER ${PROJECT_NAME} PROJECT_NAME_LOWERCASE)



# The name of the library should be equal to the project name
if(NOT LIBRARY_NAME)
    set(LIBRARY_NAME ${PROJECT_NAME})
endif()

# We want to make a static library
set(LIBRARY_TYPE STATIC)
set(EXPORT_TYPE ARCHIVE)



# Find the source folder
set(PROJECT_SOURCE_FOLDER ${CMAKE_SOURCE_DIR}/src)

# Find the source files
file(GLOB PROJECT_SOURCE_FILES ${PROJECT_SOURCE_FOLDER}/*.cpp)

# Find the header folder
set(PROJECT_INCLUDE_FOLDER ${CMAKE_SOURCE_DIR}/include)

# Find the tests folder
set(PROJECT_TESTS_FOLDER ${CMAKE_SOURCE_DIR}/tests)

# Find the source files for the tests
file(GLOB PROJECT_TEST_SOURCE_FILES ${PROJECT_TESTS_FOLDER}/*.cpp)



#   LIBINT_PREFIX
# FIXME: this should be changed into find_package(libint)
if(NOT LIBINT_PREFIX)
    if(DEFINED ENV{LIBINT_PREFIX})
        set(LIBINT_PREFIX $ENV{LIBINT_PREFIX})
    else()
        message(WARNING "You did not specify -DLIBINT_PREFIX, nor is LIBINT_PREFIX set in your environment. Proceeding with default value for LIBINT_PREFIX.")
        set(LIBINT_PREFIX /usr/local/libint/2.3.1)
    endif()
endif()



# Give the user the option to specify an installation prefix. If not given as -DINSTALLATION_PREFIX, defaults to /usr/local.
if(NOT INSTALLATION_PREFIX)
    set(INSTALLATION_PREFIX ${CMAKE_INSTALL_PREFIX})
endif()
