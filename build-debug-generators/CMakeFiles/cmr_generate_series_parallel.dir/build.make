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
include CMakeFiles/cmr_generate_series_parallel.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/cmr_generate_series_parallel.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/cmr_generate_series_parallel.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/cmr_generate_series_parallel.dir/flags.make

CMakeFiles/cmr_generate_series_parallel.dir/src/gen/series_parallel_gen.c.o: CMakeFiles/cmr_generate_series_parallel.dir/flags.make
CMakeFiles/cmr_generate_series_parallel.dir/src/gen/series_parallel_gen.c.o: /home/runner/work/cmr/cmr/src/gen/series_parallel_gen.c
CMakeFiles/cmr_generate_series_parallel.dir/src/gen/series_parallel_gen.c.o: CMakeFiles/cmr_generate_series_parallel.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/runner/work/cmr/cmr/build-debug-generators/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object CMakeFiles/cmr_generate_series_parallel.dir/src/gen/series_parallel_gen.c.o"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT CMakeFiles/cmr_generate_series_parallel.dir/src/gen/series_parallel_gen.c.o -MF CMakeFiles/cmr_generate_series_parallel.dir/src/gen/series_parallel_gen.c.o.d -o CMakeFiles/cmr_generate_series_parallel.dir/src/gen/series_parallel_gen.c.o -c /home/runner/work/cmr/cmr/src/gen/series_parallel_gen.c

CMakeFiles/cmr_generate_series_parallel.dir/src/gen/series_parallel_gen.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/cmr_generate_series_parallel.dir/src/gen/series_parallel_gen.c.i"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/runner/work/cmr/cmr/src/gen/series_parallel_gen.c > CMakeFiles/cmr_generate_series_parallel.dir/src/gen/series_parallel_gen.c.i

CMakeFiles/cmr_generate_series_parallel.dir/src/gen/series_parallel_gen.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/cmr_generate_series_parallel.dir/src/gen/series_parallel_gen.c.s"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/runner/work/cmr/cmr/src/gen/series_parallel_gen.c -o CMakeFiles/cmr_generate_series_parallel.dir/src/gen/series_parallel_gen.c.s

# Object files for target cmr_generate_series_parallel
cmr_generate_series_parallel_OBJECTS = \
"CMakeFiles/cmr_generate_series_parallel.dir/src/gen/series_parallel_gen.c.o"

# External object files for target cmr_generate_series_parallel
cmr_generate_series_parallel_EXTERNAL_OBJECTS =

cmr-generate-series-parallel: CMakeFiles/cmr_generate_series_parallel.dir/src/gen/series_parallel_gen.c.o
cmr-generate-series-parallel: CMakeFiles/cmr_generate_series_parallel.dir/build.make
cmr-generate-series-parallel: libcmr.so
cmr-generate-series-parallel: CMakeFiles/cmr_generate_series_parallel.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/runner/work/cmr/cmr/build-debug-generators/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking C executable cmr-generate-series-parallel"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/cmr_generate_series_parallel.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/cmr_generate_series_parallel.dir/build: cmr-generate-series-parallel
.PHONY : CMakeFiles/cmr_generate_series_parallel.dir/build

CMakeFiles/cmr_generate_series_parallel.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/cmr_generate_series_parallel.dir/cmake_clean.cmake
.PHONY : CMakeFiles/cmr_generate_series_parallel.dir/clean

CMakeFiles/cmr_generate_series_parallel.dir/depend:
	cd /home/runner/work/cmr/cmr/build-debug-generators && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/runner/work/cmr/cmr /home/runner/work/cmr/cmr /home/runner/work/cmr/cmr/build-debug-generators /home/runner/work/cmr/cmr/build-debug-generators /home/runner/work/cmr/cmr/build-debug-generators/CMakeFiles/cmr_generate_series_parallel.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/cmr_generate_series_parallel.dir/depend

