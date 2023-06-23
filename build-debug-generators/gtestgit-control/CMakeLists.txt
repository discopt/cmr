cmake_minimum_required(VERSION 3.5.0)
project(gtestgit-control NONE)
include(ExternalProject)
ExternalProject_Add(googletest
   GIT_REPOSITORY https://github.com/google/googletest.git
   GIT_TAG main
   SOURCE_DIR "/home/runner/work/cmr/cmr/build-debug-generators/gtestgit-src"
   BINARY_DIR "/home/runner/work/cmr/cmr/build-debug-generators/gtestgit-build"
   CONFIGURE_COMMAND ""
   BUILD_COMMAND ""
   INSTALL_COMMAND ""
   TEST_COMMAND ""
)