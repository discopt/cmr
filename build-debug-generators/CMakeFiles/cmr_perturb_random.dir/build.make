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
include CMakeFiles/cmr_perturb_random.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/cmr_perturb_random.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/cmr_perturb_random.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/cmr_perturb_random.dir/flags.make

CMakeFiles/cmr_perturb_random.dir/src/gen/random_perturb.c.o: CMakeFiles/cmr_perturb_random.dir/flags.make
CMakeFiles/cmr_perturb_random.dir/src/gen/random_perturb.c.o: /home/runner/work/cmr/cmr/src/gen/random_perturb.c
CMakeFiles/cmr_perturb_random.dir/src/gen/random_perturb.c.o: CMakeFiles/cmr_perturb_random.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/runner/work/cmr/cmr/build-debug-generators/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object CMakeFiles/cmr_perturb_random.dir/src/gen/random_perturb.c.o"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT CMakeFiles/cmr_perturb_random.dir/src/gen/random_perturb.c.o -MF CMakeFiles/cmr_perturb_random.dir/src/gen/random_perturb.c.o.d -o CMakeFiles/cmr_perturb_random.dir/src/gen/random_perturb.c.o -c /home/runner/work/cmr/cmr/src/gen/random_perturb.c

CMakeFiles/cmr_perturb_random.dir/src/gen/random_perturb.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/cmr_perturb_random.dir/src/gen/random_perturb.c.i"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/runner/work/cmr/cmr/src/gen/random_perturb.c > CMakeFiles/cmr_perturb_random.dir/src/gen/random_perturb.c.i

CMakeFiles/cmr_perturb_random.dir/src/gen/random_perturb.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/cmr_perturb_random.dir/src/gen/random_perturb.c.s"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/runner/work/cmr/cmr/src/gen/random_perturb.c -o CMakeFiles/cmr_perturb_random.dir/src/gen/random_perturb.c.s

# Object files for target cmr_perturb_random
cmr_perturb_random_OBJECTS = \
"CMakeFiles/cmr_perturb_random.dir/src/gen/random_perturb.c.o"

# External object files for target cmr_perturb_random
cmr_perturb_random_EXTERNAL_OBJECTS =

cmr-perturb-random: CMakeFiles/cmr_perturb_random.dir/src/gen/random_perturb.c.o
cmr-perturb-random: CMakeFiles/cmr_perturb_random.dir/build.make
cmr-perturb-random: libcmr.so
cmr-perturb-random: CMakeFiles/cmr_perturb_random.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/runner/work/cmr/cmr/build-debug-generators/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking C executable cmr-perturb-random"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/cmr_perturb_random.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/cmr_perturb_random.dir/build: cmr-perturb-random
.PHONY : CMakeFiles/cmr_perturb_random.dir/build

CMakeFiles/cmr_perturb_random.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/cmr_perturb_random.dir/cmake_clean.cmake
.PHONY : CMakeFiles/cmr_perturb_random.dir/clean

CMakeFiles/cmr_perturb_random.dir/depend:
	cd /home/runner/work/cmr/cmr/build-debug-generators && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/runner/work/cmr/cmr /home/runner/work/cmr/cmr /home/runner/work/cmr/cmr/build-debug-generators /home/runner/work/cmr/cmr/build-debug-generators /home/runner/work/cmr/cmr/build-debug-generators/CMakeFiles/cmr_perturb_random.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/cmr_perturb_random.dir/depend

