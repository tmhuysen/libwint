# Install script for directory: /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/Forks/libwint/src

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/usr/local/libwint/lib/libwint.a")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/usr/local/libwint/lib" TYPE STATIC_LIBRARY FILES "/Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/Forks/libwint/bb/src/libwint.a")
  if(EXISTS "$ENV{DESTDIR}/usr/local/libwint/lib/libwint.a" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/usr/local/libwint/lib/libwint.a")
    execute_process(COMMAND "/opt/local/bin/ranlib" "$ENV{DESTDIR}/usr/local/libwint/lib/libwint.a")
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/usr/local/libwint/include/AOBasis.hpp;/usr/local/libwint/include/LibintCommunicator.hpp;/usr/local/libwint/include/Molecule.hpp;/usr/local/libwint/include/SOBasis.hpp;/usr/local/libwint/include/SOMullikenBasis.hpp;/usr/local/libwint/include/libwint.hpp;/usr/local/libwint/include/transformations.hpp;/usr/local/libwint/include/version.hpp")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/usr/local/libwint/include" TYPE FILE FILES
    "/Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/Forks/libwint/include/AOBasis.hpp"
    "/Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/Forks/libwint/include/LibintCommunicator.hpp"
    "/Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/Forks/libwint/include/Molecule.hpp"
    "/Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/Forks/libwint/include/SOBasis.hpp"
    "/Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/Forks/libwint/include/SOMullikenBasis.hpp"
    "/Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/Forks/libwint/include/libwint.hpp"
    "/Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/Forks/libwint/include/transformations.hpp"
    "/Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/Forks/libwint/include/version.hpp"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}/usr/local/libwint/cmake/libwintTargets.cmake")
    file(DIFFERENT EXPORT_FILE_CHANGED FILES
         "$ENV{DESTDIR}/usr/local/libwint/cmake/libwintTargets.cmake"
         "/Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/Forks/libwint/bb/src/CMakeFiles/Export/_usr/local/libwint/cmake/libwintTargets.cmake")
    if(EXPORT_FILE_CHANGED)
      file(GLOB OLD_CONFIG_FILES "$ENV{DESTDIR}/usr/local/libwint/cmake/libwintTargets-*.cmake")
      if(OLD_CONFIG_FILES)
        message(STATUS "Old export file \"$ENV{DESTDIR}/usr/local/libwint/cmake/libwintTargets.cmake\" will be replaced.  Removing files [${OLD_CONFIG_FILES}].")
        file(REMOVE ${OLD_CONFIG_FILES})
      endif()
    endif()
  endif()
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/usr/local/libwint/cmake/libwintTargets.cmake")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/usr/local/libwint/cmake" TYPE FILE FILES "/Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/Forks/libwint/bb/src/CMakeFiles/Export/_usr/local/libwint/cmake/libwintTargets.cmake")
  if("${CMAKE_INSTALL_CONFIG_NAME}" MATCHES "^([Rr][Ee][Ll][Ee][Aa][Ss][Ee])$")
    list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
     "/usr/local/libwint/cmake/libwintTargets-release.cmake")
    if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
        message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
    endif()
    if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
        message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
    endif()
file(INSTALL DESTINATION "/usr/local/libwint/cmake" TYPE FILE FILES "/Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/Forks/libwint/bb/src/CMakeFiles/Export/_usr/local/libwint/cmake/libwintTargets-release.cmake")
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/usr/local/libwint/cmake/libwintConfig.cmake;/usr/local/libwint/cmake/libwintConfigVersion.cmake")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/usr/local/libwint/cmake" TYPE FILE FILES
    "/Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/Forks/libwint/cmake/libwintConfig.cmake"
    "/Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/Forks/libwint/cmake/libwintConfigVersion.cmake"
    )
endif()

