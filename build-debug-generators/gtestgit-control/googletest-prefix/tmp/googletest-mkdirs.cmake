# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION 3.5)

file(MAKE_DIRECTORY
  "/home/runner/work/cmr/cmr/build-debug-generators/gtestgit-src"
  "/home/runner/work/cmr/cmr/build-debug-generators/gtestgit-build"
  "/home/runner/work/cmr/cmr/build-debug-generators/gtestgit-control/googletest-prefix"
  "/home/runner/work/cmr/cmr/build-debug-generators/gtestgit-control/googletest-prefix/tmp"
  "/home/runner/work/cmr/cmr/build-debug-generators/gtestgit-control/googletest-prefix/src/googletest-stamp"
  "/home/runner/work/cmr/cmr/build-debug-generators/gtestgit-control/googletest-prefix/src"
  "/home/runner/work/cmr/cmr/build-debug-generators/gtestgit-control/googletest-prefix/src/googletest-stamp"
)

set(configSubDirs )
foreach(subDir IN LISTS configSubDirs)
    file(MAKE_DIRECTORY "/home/runner/work/cmr/cmr/build-debug-generators/gtestgit-control/googletest-prefix/src/googletest-stamp/${subDir}")
endforeach()
if(cfgdir)
  file(MAKE_DIRECTORY "/home/runner/work/cmr/cmr/build-debug-generators/gtestgit-control/googletest-prefix/src/googletest-stamp${cfgdir}") # cfgdir has leading slash
endif()
