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
include CMakeFiles/cmr_matrix.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/cmr_matrix.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/cmr_matrix.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/cmr_matrix.dir/flags.make

CMakeFiles/cmr_matrix.dir/src/main/matrix_main.c.o: CMakeFiles/cmr_matrix.dir/flags.make
CMakeFiles/cmr_matrix.dir/src/main/matrix_main.c.o: /home/runner/work/cmr/cmr/src/main/matrix_main.c
CMakeFiles/cmr_matrix.dir/src/main/matrix_main.c.o: CMakeFiles/cmr_matrix.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/runner/work/cmr/cmr/build-debug-generators/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object CMakeFiles/cmr_matrix.dir/src/main/matrix_main.c.o"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT CMakeFiles/cmr_matrix.dir/src/main/matrix_main.c.o -MF CMakeFiles/cmr_matrix.dir/src/main/matrix_main.c.o.d -o CMakeFiles/cmr_matrix.dir/src/main/matrix_main.c.o -c /home/runner/work/cmr/cmr/src/main/matrix_main.c

CMakeFiles/cmr_matrix.dir/src/main/matrix_main.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/cmr_matrix.dir/src/main/matrix_main.c.i"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/runner/work/cmr/cmr/src/main/matrix_main.c > CMakeFiles/cmr_matrix.dir/src/main/matrix_main.c.i

CMakeFiles/cmr_matrix.dir/src/main/matrix_main.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/cmr_matrix.dir/src/main/matrix_main.c.s"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/runner/work/cmr/cmr/src/main/matrix_main.c -o CMakeFiles/cmr_matrix.dir/src/main/matrix_main.c.s

# Object files for target cmr_matrix
cmr_matrix_OBJECTS = \
"CMakeFiles/cmr_matrix.dir/src/main/matrix_main.c.o"

# External object files for target cmr_matrix
cmr_matrix_EXTERNAL_OBJECTS =

cmr-matrix: CMakeFiles/cmr_matrix.dir/src/main/matrix_main.c.o
cmr-matrix: CMakeFiles/cmr_matrix.dir/build.make
cmr-matrix: libcmr.so
cmr-matrix: CMakeFiles/cmr_matrix.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/runner/work/cmr/cmr/build-debug-generators/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking C executable cmr-matrix"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/cmr_matrix.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/cmr_matrix.dir/build: cmr-matrix
.PHONY : CMakeFiles/cmr_matrix.dir/build

CMakeFiles/cmr_matrix.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/cmr_matrix.dir/cmake_clean.cmake
.PHONY : CMakeFiles/cmr_matrix.dir/clean

CMakeFiles/cmr_matrix.dir/depend:
	cd /home/runner/work/cmr/cmr/build-debug-generators && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/runner/work/cmr/cmr /home/runner/work/cmr/cmr /home/runner/work/cmr/cmr/build-debug-generators /home/runner/work/cmr/cmr/build-debug-generators /home/runner/work/cmr/cmr/build-debug-generators/CMakeFiles/cmr_matrix.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/cmr_matrix.dir/depend

