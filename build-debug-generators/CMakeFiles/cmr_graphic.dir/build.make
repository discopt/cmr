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
include CMakeFiles/cmr_graphic.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/cmr_graphic.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/cmr_graphic.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/cmr_graphic.dir/flags.make

CMakeFiles/cmr_graphic.dir/src/main/graphic_main.c.o: CMakeFiles/cmr_graphic.dir/flags.make
CMakeFiles/cmr_graphic.dir/src/main/graphic_main.c.o: /home/runner/work/cmr/cmr/src/main/graphic_main.c
CMakeFiles/cmr_graphic.dir/src/main/graphic_main.c.o: CMakeFiles/cmr_graphic.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/runner/work/cmr/cmr/build-debug-generators/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object CMakeFiles/cmr_graphic.dir/src/main/graphic_main.c.o"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT CMakeFiles/cmr_graphic.dir/src/main/graphic_main.c.o -MF CMakeFiles/cmr_graphic.dir/src/main/graphic_main.c.o.d -o CMakeFiles/cmr_graphic.dir/src/main/graphic_main.c.o -c /home/runner/work/cmr/cmr/src/main/graphic_main.c

CMakeFiles/cmr_graphic.dir/src/main/graphic_main.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/cmr_graphic.dir/src/main/graphic_main.c.i"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/runner/work/cmr/cmr/src/main/graphic_main.c > CMakeFiles/cmr_graphic.dir/src/main/graphic_main.c.i

CMakeFiles/cmr_graphic.dir/src/main/graphic_main.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/cmr_graphic.dir/src/main/graphic_main.c.s"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/runner/work/cmr/cmr/src/main/graphic_main.c -o CMakeFiles/cmr_graphic.dir/src/main/graphic_main.c.s

# Object files for target cmr_graphic
cmr_graphic_OBJECTS = \
"CMakeFiles/cmr_graphic.dir/src/main/graphic_main.c.o"

# External object files for target cmr_graphic
cmr_graphic_EXTERNAL_OBJECTS =

cmr-graphic: CMakeFiles/cmr_graphic.dir/src/main/graphic_main.c.o
cmr-graphic: CMakeFiles/cmr_graphic.dir/build.make
cmr-graphic: libcmr.so
cmr-graphic: CMakeFiles/cmr_graphic.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/runner/work/cmr/cmr/build-debug-generators/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking C executable cmr-graphic"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/cmr_graphic.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/cmr_graphic.dir/build: cmr-graphic
.PHONY : CMakeFiles/cmr_graphic.dir/build

CMakeFiles/cmr_graphic.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/cmr_graphic.dir/cmake_clean.cmake
.PHONY : CMakeFiles/cmr_graphic.dir/clean

CMakeFiles/cmr_graphic.dir/depend:
	cd /home/runner/work/cmr/cmr/build-debug-generators && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/runner/work/cmr/cmr /home/runner/work/cmr/cmr /home/runner/work/cmr/cmr/build-debug-generators /home/runner/work/cmr/cmr/build-debug-generators /home/runner/work/cmr/cmr/build-debug-generators/CMakeFiles/cmr_graphic.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/cmr_graphic.dir/depend

