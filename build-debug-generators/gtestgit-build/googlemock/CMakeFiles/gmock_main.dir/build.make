# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.26

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/local/bin/cmake

# The command to remove a file.
RM = /usr/local/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/runner/work/cmr/cmr

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/runner/work/cmr/cmr/build-debug-generators

# Include any dependencies generated for this target.
include gtestgit-build/googlemock/CMakeFiles/gmock_main.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include gtestgit-build/googlemock/CMakeFiles/gmock_main.dir/compiler_depend.make

# Include the progress variables for this target.
include gtestgit-build/googlemock/CMakeFiles/gmock_main.dir/progress.make

# Include the compile flags for this target's objects.
include gtestgit-build/googlemock/CMakeFiles/gmock_main.dir/flags.make

gtestgit-build/googlemock/CMakeFiles/gmock_main.dir/src/gmock_main.cc.o: gtestgit-build/googlemock/CMakeFiles/gmock_main.dir/flags.make
gtestgit-build/googlemock/CMakeFiles/gmock_main.dir/src/gmock_main.cc.o: gtestgit-src/googlemock/src/gmock_main.cc
gtestgit-build/googlemock/CMakeFiles/gmock_main.dir/src/gmock_main.cc.o: gtestgit-build/googlemock/CMakeFiles/gmock_main.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/runner/work/cmr/cmr/build-debug-generators/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object gtestgit-build/googlemock/CMakeFiles/gmock_main.dir/src/gmock_main.cc.o"
	cd /home/runner/work/cmr/cmr/build-debug-generators/gtestgit-build/googlemock && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT gtestgit-build/googlemock/CMakeFiles/gmock_main.dir/src/gmock_main.cc.o -MF CMakeFiles/gmock_main.dir/src/gmock_main.cc.o.d -o CMakeFiles/gmock_main.dir/src/gmock_main.cc.o -c /home/runner/work/cmr/cmr/build-debug-generators/gtestgit-src/googlemock/src/gmock_main.cc

gtestgit-build/googlemock/CMakeFiles/gmock_main.dir/src/gmock_main.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/gmock_main.dir/src/gmock_main.cc.i"
	cd /home/runner/work/cmr/cmr/build-debug-generators/gtestgit-build/googlemock && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/runner/work/cmr/cmr/build-debug-generators/gtestgit-src/googlemock/src/gmock_main.cc > CMakeFiles/gmock_main.dir/src/gmock_main.cc.i

gtestgit-build/googlemock/CMakeFiles/gmock_main.dir/src/gmock_main.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/gmock_main.dir/src/gmock_main.cc.s"
	cd /home/runner/work/cmr/cmr/build-debug-generators/gtestgit-build/googlemock && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/runner/work/cmr/cmr/build-debug-generators/gtestgit-src/googlemock/src/gmock_main.cc -o CMakeFiles/gmock_main.dir/src/gmock_main.cc.s

# Object files for target gmock_main
gmock_main_OBJECTS = \
"CMakeFiles/gmock_main.dir/src/gmock_main.cc.o"

# External object files for target gmock_main
gmock_main_EXTERNAL_OBJECTS =

lib/libgmock_main.so.1.13.0: gtestgit-build/googlemock/CMakeFiles/gmock_main.dir/src/gmock_main.cc.o
lib/libgmock_main.so.1.13.0: gtestgit-build/googlemock/CMakeFiles/gmock_main.dir/build.make
lib/libgmock_main.so.1.13.0: lib/libgmock.so.1.13.0
lib/libgmock_main.so.1.13.0: lib/libgtest.so.1.13.0
lib/libgmock_main.so.1.13.0: gtestgit-build/googlemock/CMakeFiles/gmock_main.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/runner/work/cmr/cmr/build-debug-generators/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX shared library ../../lib/libgmock_main.so"
	cd /home/runner/work/cmr/cmr/build-debug-generators/gtestgit-build/googlemock && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/gmock_main.dir/link.txt --verbose=$(VERBOSE)
	cd /home/runner/work/cmr/cmr/build-debug-generators/gtestgit-build/googlemock && $(CMAKE_COMMAND) -E cmake_symlink_library ../../lib/libgmock_main.so.1.13.0 ../../lib/libgmock_main.so.1.13.0 ../../lib/libgmock_main.so

lib/libgmock_main.so: lib/libgmock_main.so.1.13.0
	@$(CMAKE_COMMAND) -E touch_nocreate lib/libgmock_main.so

# Rule to build all files generated by this target.
gtestgit-build/googlemock/CMakeFiles/gmock_main.dir/build: lib/libgmock_main.so
.PHONY : gtestgit-build/googlemock/CMakeFiles/gmock_main.dir/build

gtestgit-build/googlemock/CMakeFiles/gmock_main.dir/clean:
	cd /home/runner/work/cmr/cmr/build-debug-generators/gtestgit-build/googlemock && $(CMAKE_COMMAND) -P CMakeFiles/gmock_main.dir/cmake_clean.cmake
.PHONY : gtestgit-build/googlemock/CMakeFiles/gmock_main.dir/clean

gtestgit-build/googlemock/CMakeFiles/gmock_main.dir/depend:
	cd /home/runner/work/cmr/cmr/build-debug-generators && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/runner/work/cmr/cmr /home/runner/work/cmr/cmr/build-debug-generators/gtestgit-src/googlemock /home/runner/work/cmr/cmr/build-debug-generators /home/runner/work/cmr/cmr/build-debug-generators/gtestgit-build/googlemock /home/runner/work/cmr/cmr/build-debug-generators/gtestgit-build/googlemock/CMakeFiles/gmock_main.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : gtestgit-build/googlemock/CMakeFiles/gmock_main.dir/depend

