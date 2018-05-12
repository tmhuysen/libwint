#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "libwint" for configuration "Release"
set_property(TARGET libwint APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(libwint PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "CXX"
  IMPORTED_LOCATION_RELEASE "/usr/local/libwint/lib/libwint.a"
  )

list(APPEND _IMPORT_CHECK_TARGETS libwint )
list(APPEND _IMPORT_CHECK_FILES_FOR_libwint "/usr/local/libwint/lib/libwint.a" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
