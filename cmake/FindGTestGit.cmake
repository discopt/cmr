# FindGTestGit
# ------------
#
# Download the sources of the Google C++ Testing Framework, build them and add them in a subdirectory.
#
# Sets variables from the subproject.

message(STATUS "Using GTest sources from git repository.")

# Create a cmake subproject for downloading in gtestgit-control/.
file(WRITE ${CMAKE_BINARY_DIR}/gtestgit-control/CMakeLists.txt "\
cmake_minimum_required(VERSION 3.5.0)
project(gtestgit-control NONE)
include(ExternalProject)
ExternalProject_Add(googletest
   GIT_REPOSITORY https://github.com/google/googletest.git
   GIT_TAG main
   SOURCE_DIR \"${CMAKE_BINARY_DIR}/gtestgit-src\"
   BINARY_DIR \"${CMAKE_BINARY_DIR}/gtestgit-build\"
   CONFIGURE_COMMAND \"\"
   BUILD_COMMAND \"\"
   INSTALL_COMMAND \"\"
   TEST_COMMAND \"\"
)")

# Run cmake in gtestgit-control.
execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
   RESULT_VARIABLE result
   WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/gtestgit-control )
if(result)
   message(FATAL_ERROR "CMake step for FindGTestGit failed: ${result}")
endif()

# Build gtestgit-control.
execute_process(COMMAND ${CMAKE_COMMAND} --build .
   RESULT_VARIABLE result
   WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/gtestgit-control )
if(result)
   message(FATAL_ERROR "Build step for FindGTestGit failed: ${result}")
endif()

# Add subdirectory with sources in gtest-src and build in gtest-build.
add_subdirectory(${CMAKE_BINARY_DIR}/gtestgit-src
   ${CMAKE_BINARY_DIR}/gtestgit-build
   EXCLUDE_FROM_ALL)
