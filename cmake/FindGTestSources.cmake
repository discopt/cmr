# FindGTestSources
# ----------------
#
# Locate the sources of the Google C++ Testing Framework, build them and add them in a subdirectory.
#
# Result variables
# ^^^^^^^^^^^^^^^^
#
# This module will set the following variables in your project:
#
# ``GTESTSOURCES_FOUND``
#   Found the sources of the Google Testing framework
# ``GTESTSOURCES_BASE``
#   the directory containing the sources of the Google Test headers
#
# In addition also sets variables from the subproject.

find_path(GTESTSOURCES_BASE
   googletest/src/gtest-all.cc
   PATHS
      $ENV{GTESTSOURCES_ROOT}
      ${GTESTSOURCES_ROOT}
   HINTS
      /usr/src/googletest
      /usr/src/gtest
)
if(GTESTSOURCES_BASE)
   message(STATUS "Found GTest sources in ${GTESTSOURCES_BASE}.")
   set(GTESTSOURCES_FOUND TRUE)
   set(OLD_CMAKE_POLICY_DEFAULT_CMP0048 ${CMAKE_POLICY_DEFAULT_CMP0048})
   set(CMAKE_POLICY_DEFAULT_CMP0048 NEW)
   
   set(OLD_CMAKE_POLICY_DEFAULT_CMP0054 ${CMAKE_POLICY_DEFAULT_CMP0054})
   set(CMAKE_POLICY_DEFAULT_CMP0054 NEW)

   file(WRITE ${CMAKE_BINARY_DIR}/gtestsources-control/CMakeLists.txt "\
cmake_minimum_required(VERSION 3.5.0)
project(gtestsources-control NONE)
include(ExternalProject)
ExternalProject_Add(googletest
   SOURCE_DIR \"${GTESTSOURCES_BASE}\"
   BINARY_DIR \"${CMAKE_BINARY_DIR}/gtestsources-build\"
   CONFIGURE_COMMAND \"\"
   BUILD_COMMAND \"\"
   INSTALL_COMMAND \"\"
   TEST_COMMAND \"\"
)")
   add_subdirectory(${GTESTSOURCES_BASE}
      ${CMAKE_BINARY_DIR}/gtestsources-build
      EXCLUDE_FROM_ALL)
   set(CMAKE_POLICY_DEFAULT_CMP0054 ${OLD_CMAKE_POLICY_DEFAULT_CMP0054})
   set(CMAKE_POLICY_DEFAULT_CMP0048 ${OLD_CMAKE_POLICY_DEFAULT_CMP0048})
else()
   message(STATUS "Could not find GTestSources.")
   set(GTESTSOURCES_FOUND FALSE)
endif()
