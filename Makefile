# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.5

# Default target executed when no arguments are given to make.
default_target: all

.PHONY : default_target

# Allow only one "make -f Makefile2" at a time, but pass parallelism.
.NOTPARALLEL:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/daniel/ownCloud/fem/polyFEM

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/daniel/ownCloud/fem/polyFEM

#=============================================================================
# Targets provided globally by CMake.

# Special rule for the target rebuild_cache
rebuild_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running CMake to regenerate build system..."
	/usr/bin/cmake -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR)
.PHONY : rebuild_cache

# Special rule for the target rebuild_cache
rebuild_cache/fast: rebuild_cache

.PHONY : rebuild_cache/fast

# Special rule for the target edit_cache
edit_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "No interactive CMake dialog available..."
	/usr/bin/cmake -E echo No\ interactive\ CMake\ dialog\ available.
.PHONY : edit_cache

# Special rule for the target edit_cache
edit_cache/fast: edit_cache

.PHONY : edit_cache/fast

# The main all target
all: cmake_check_build_system
	$(CMAKE_COMMAND) -E cmake_progress_start /home/daniel/ownCloud/fem/polyFEM/CMakeFiles /home/daniel/ownCloud/fem/polyFEM/CMakeFiles/progress.marks
	$(MAKE) -f CMakeFiles/Makefile2 all
	$(CMAKE_COMMAND) -E cmake_progress_start /home/daniel/ownCloud/fem/polyFEM/CMakeFiles 0
.PHONY : all

# The main clean target
clean:
	$(MAKE) -f CMakeFiles/Makefile2 clean
.PHONY : clean

# The main clean target
clean/fast: clean

.PHONY : clean/fast

# Prepare targets for installation.
preinstall: all
	$(MAKE) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall

# Prepare targets for installation.
preinstall/fast:
	$(MAKE) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall/fast

# clear depends
depend:
	$(CMAKE_COMMAND) -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 1
.PHONY : depend

#=============================================================================
# Target rules for targets named main

# Build rule for target.
main: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 main
.PHONY : main

# fast build rule for target.
main/fast:
	$(MAKE) -f CMakeFiles/main.dir/build.make CMakeFiles/main.dir/build
.PHONY : main/fast

Delta.o: Delta.cpp.o

.PHONY : Delta.o

# target to build an object file
Delta.cpp.o:
	$(MAKE) -f CMakeFiles/main.dir/build.make CMakeFiles/main.dir/Delta.cpp.o
.PHONY : Delta.cpp.o

Delta.i: Delta.cpp.i

.PHONY : Delta.i

# target to preprocess a source file
Delta.cpp.i:
	$(MAKE) -f CMakeFiles/main.dir/build.make CMakeFiles/main.dir/Delta.cpp.i
.PHONY : Delta.cpp.i

Delta.s: Delta.cpp.s

.PHONY : Delta.s

# target to generate assembly for a file
Delta.cpp.s:
	$(MAKE) -f CMakeFiles/main.dir/build.make CMakeFiles/main.dir/Delta.cpp.s
.PHONY : Delta.cpp.s

draw.o: draw.cpp.o

.PHONY : draw.o

# target to build an object file
draw.cpp.o:
	$(MAKE) -f CMakeFiles/main.dir/build.make CMakeFiles/main.dir/draw.cpp.o
.PHONY : draw.cpp.o

draw.i: draw.cpp.i

.PHONY : draw.i

# target to preprocess a source file
draw.cpp.i:
	$(MAKE) -f CMakeFiles/main.dir/build.make CMakeFiles/main.dir/draw.cpp.i
.PHONY : draw.cpp.i

draw.s: draw.cpp.s

.PHONY : draw.s

# target to generate assembly for a file
draw.cpp.s:
	$(MAKE) -f CMakeFiles/main.dir/build.make CMakeFiles/main.dir/draw.cpp.s
.PHONY : draw.cpp.s

fields.o: fields.cpp.o

.PHONY : fields.o

# target to build an object file
fields.cpp.o:
	$(MAKE) -f CMakeFiles/main.dir/build.make CMakeFiles/main.dir/fields.cpp.o
.PHONY : fields.cpp.o

fields.i: fields.cpp.i

.PHONY : fields.i

# target to preprocess a source file
fields.cpp.i:
	$(MAKE) -f CMakeFiles/main.dir/build.make CMakeFiles/main.dir/fields.cpp.i
.PHONY : fields.cpp.i

fields.s: fields.cpp.s

.PHONY : fields.s

# target to generate assembly for a file
fields.cpp.s:
	$(MAKE) -f CMakeFiles/main.dir/build.make CMakeFiles/main.dir/fields.cpp.s
.PHONY : fields.cpp.s

gradient.o: gradient.cpp.o

.PHONY : gradient.o

# target to build an object file
gradient.cpp.o:
	$(MAKE) -f CMakeFiles/main.dir/build.make CMakeFiles/main.dir/gradient.cpp.o
.PHONY : gradient.cpp.o

gradient.i: gradient.cpp.i

.PHONY : gradient.i

# target to preprocess a source file
gradient.cpp.i:
	$(MAKE) -f CMakeFiles/main.dir/build.make CMakeFiles/main.dir/gradient.cpp.i
.PHONY : gradient.cpp.i

gradient.s: gradient.cpp.s

.PHONY : gradient.s

# target to generate assembly for a file
gradient.cpp.s:
	$(MAKE) -f CMakeFiles/main.dir/build.make CMakeFiles/main.dir/gradient.cpp.s
.PHONY : gradient.cpp.s

linear.o: linear.cpp.o

.PHONY : linear.o

# target to build an object file
linear.cpp.o:
	$(MAKE) -f CMakeFiles/main.dir/build.make CMakeFiles/main.dir/linear.cpp.o
.PHONY : linear.cpp.o

linear.i: linear.cpp.i

.PHONY : linear.i

# target to preprocess a source file
linear.cpp.i:
	$(MAKE) -f CMakeFiles/main.dir/build.make CMakeFiles/main.dir/linear.cpp.i
.PHONY : linear.cpp.i

linear.s: linear.cpp.s

.PHONY : linear.s

# target to generate assembly for a file
linear.cpp.s:
	$(MAKE) -f CMakeFiles/main.dir/build.make CMakeFiles/main.dir/linear.cpp.s
.PHONY : linear.cpp.s

main.o: main.cpp.o

.PHONY : main.o

# target to build an object file
main.cpp.o:
	$(MAKE) -f CMakeFiles/main.dir/build.make CMakeFiles/main.dir/main.cpp.o
.PHONY : main.cpp.o

main.i: main.cpp.i

.PHONY : main.i

# target to preprocess a source file
main.cpp.i:
	$(MAKE) -f CMakeFiles/main.dir/build.make CMakeFiles/main.dir/main.cpp.i
.PHONY : main.cpp.i

main.s: main.cpp.s

.PHONY : main.s

# target to generate assembly for a file
main.cpp.s:
	$(MAKE) -f CMakeFiles/main.dir/build.make CMakeFiles/main.dir/main.cpp.s
.PHONY : main.cpp.s

move.o: move.cpp.o

.PHONY : move.o

# target to build an object file
move.cpp.o:
	$(MAKE) -f CMakeFiles/main.dir/build.make CMakeFiles/main.dir/move.cpp.o
.PHONY : move.cpp.o

move.i: move.cpp.i

.PHONY : move.i

# target to preprocess a source file
move.cpp.i:
	$(MAKE) -f CMakeFiles/main.dir/build.make CMakeFiles/main.dir/move.cpp.i
.PHONY : move.cpp.i

move.s: move.cpp.s

.PHONY : move.s

# target to generate assembly for a file
move.cpp.s:
	$(MAKE) -f CMakeFiles/main.dir/build.make CMakeFiles/main.dir/move.cpp.s
.PHONY : move.cpp.s

nabla.o: nabla.cpp.o

.PHONY : nabla.o

# target to build an object file
nabla.cpp.o:
	$(MAKE) -f CMakeFiles/main.dir/build.make CMakeFiles/main.dir/nabla.cpp.o
.PHONY : nabla.cpp.o

nabla.i: nabla.cpp.i

.PHONY : nabla.i

# target to preprocess a source file
nabla.cpp.i:
	$(MAKE) -f CMakeFiles/main.dir/build.make CMakeFiles/main.dir/nabla.cpp.i
.PHONY : nabla.cpp.i

nabla.s: nabla.cpp.s

.PHONY : nabla.s

# target to generate assembly for a file
nabla.cpp.s:
	$(MAKE) -f CMakeFiles/main.dir/build.make CMakeFiles/main.dir/nabla.cpp.s
.PHONY : nabla.cpp.s

onto_from_mesh.o: onto_from_mesh.cpp.o

.PHONY : onto_from_mesh.o

# target to build an object file
onto_from_mesh.cpp.o:
	$(MAKE) -f CMakeFiles/main.dir/build.make CMakeFiles/main.dir/onto_from_mesh.cpp.o
.PHONY : onto_from_mesh.cpp.o

onto_from_mesh.i: onto_from_mesh.cpp.i

.PHONY : onto_from_mesh.i

# target to preprocess a source file
onto_from_mesh.cpp.i:
	$(MAKE) -f CMakeFiles/main.dir/build.make CMakeFiles/main.dir/onto_from_mesh.cpp.i
.PHONY : onto_from_mesh.cpp.i

onto_from_mesh.s: onto_from_mesh.cpp.s

.PHONY : onto_from_mesh.s

# target to generate assembly for a file
onto_from_mesh.cpp.s:
	$(MAKE) -f CMakeFiles/main.dir/build.make CMakeFiles/main.dir/onto_from_mesh.cpp.s
.PHONY : onto_from_mesh.cpp.s

periodic.o: periodic.cpp.o

.PHONY : periodic.o

# target to build an object file
periodic.cpp.o:
	$(MAKE) -f CMakeFiles/main.dir/build.make CMakeFiles/main.dir/periodic.cpp.o
.PHONY : periodic.cpp.o

periodic.i: periodic.cpp.i

.PHONY : periodic.i

# target to preprocess a source file
periodic.cpp.i:
	$(MAKE) -f CMakeFiles/main.dir/build.make CMakeFiles/main.dir/periodic.cpp.i
.PHONY : periodic.cpp.i

periodic.s: periodic.cpp.s

.PHONY : periodic.s

# target to generate assembly for a file
periodic.cpp.s:
	$(MAKE) -f CMakeFiles/main.dir/build.make CMakeFiles/main.dir/periodic.cpp.s
.PHONY : periodic.cpp.s

quad_coeffs.o: quad_coeffs.cpp.o

.PHONY : quad_coeffs.o

# target to build an object file
quad_coeffs.cpp.o:
	$(MAKE) -f CMakeFiles/main.dir/build.make CMakeFiles/main.dir/quad_coeffs.cpp.o
.PHONY : quad_coeffs.cpp.o

quad_coeffs.i: quad_coeffs.cpp.i

.PHONY : quad_coeffs.i

# target to preprocess a source file
quad_coeffs.cpp.i:
	$(MAKE) -f CMakeFiles/main.dir/build.make CMakeFiles/main.dir/quad_coeffs.cpp.i
.PHONY : quad_coeffs.cpp.i

quad_coeffs.s: quad_coeffs.cpp.s

.PHONY : quad_coeffs.s

# target to generate assembly for a file
quad_coeffs.cpp.s:
	$(MAKE) -f CMakeFiles/main.dir/build.make CMakeFiles/main.dir/quad_coeffs.cpp.s
.PHONY : quad_coeffs.cpp.s

volumes.o: volumes.cpp.o

.PHONY : volumes.o

# target to build an object file
volumes.cpp.o:
	$(MAKE) -f CMakeFiles/main.dir/build.make CMakeFiles/main.dir/volumes.cpp.o
.PHONY : volumes.cpp.o

volumes.i: volumes.cpp.i

.PHONY : volumes.i

# target to preprocess a source file
volumes.cpp.i:
	$(MAKE) -f CMakeFiles/main.dir/build.make CMakeFiles/main.dir/volumes.cpp.i
.PHONY : volumes.cpp.i

volumes.s: volumes.cpp.s

.PHONY : volumes.s

# target to generate assembly for a file
volumes.cpp.s:
	$(MAKE) -f CMakeFiles/main.dir/build.make CMakeFiles/main.dir/volumes.cpp.s
.PHONY : volumes.cpp.s

# Help Target
help:
	@echo "The following are some of the valid targets for this Makefile:"
	@echo "... all (the default if no target is provided)"
	@echo "... clean"
	@echo "... depend"
	@echo "... rebuild_cache"
	@echo "... edit_cache"
	@echo "... main"
	@echo "... Delta.o"
	@echo "... Delta.i"
	@echo "... Delta.s"
	@echo "... draw.o"
	@echo "... draw.i"
	@echo "... draw.s"
	@echo "... fields.o"
	@echo "... fields.i"
	@echo "... fields.s"
	@echo "... gradient.o"
	@echo "... gradient.i"
	@echo "... gradient.s"
	@echo "... linear.o"
	@echo "... linear.i"
	@echo "... linear.s"
	@echo "... main.o"
	@echo "... main.i"
	@echo "... main.s"
	@echo "... move.o"
	@echo "... move.i"
	@echo "... move.s"
	@echo "... nabla.o"
	@echo "... nabla.i"
	@echo "... nabla.s"
	@echo "... onto_from_mesh.o"
	@echo "... onto_from_mesh.i"
	@echo "... onto_from_mesh.s"
	@echo "... periodic.o"
	@echo "... periodic.i"
	@echo "... periodic.s"
	@echo "... quad_coeffs.o"
	@echo "... quad_coeffs.i"
	@echo "... quad_coeffs.s"
	@echo "... volumes.o"
	@echo "... volumes.i"
	@echo "... volumes.s"
.PHONY : help



#=============================================================================
# Special targets to cleanup operation of make.

# Special rule to run CMake to check the build system integrity.
# No rule that depends on this can have commands that come from listfiles
# because they might be regenerated.
cmake_check_build_system:
	$(CMAKE_COMMAND) -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 0
.PHONY : cmake_check_build_system

