# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.16

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


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
CMAKE_SOURCE_DIR = /home/oem/PRG/CPP/Paulo/MEFPoro/v6_static_NewHessiana_07_06_2022

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/oem/PRG/CPP/Paulo/MEFPoro/v6_static_NewHessiana_07_06_2022/build

# Include any dependencies generated for this target.
include Examples/CMakeFiles/BasicFEM.dir/depend.make

# Include the progress variables for this target.
include Examples/CMakeFiles/BasicFEM.dir/progress.make

# Include the compile flags for this target's objects.
include Examples/CMakeFiles/BasicFEM.dir/flags.make

Examples/CMakeFiles/BasicFEM.dir/main.cpp.o: Examples/CMakeFiles/BasicFEM.dir/flags.make
Examples/CMakeFiles/BasicFEM.dir/main.cpp.o: ../Examples/main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/oem/PRG/CPP/Paulo/MEFPoro/v6_static_NewHessiana_07_06_2022/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object Examples/CMakeFiles/BasicFEM.dir/main.cpp.o"
	cd /home/oem/PRG/CPP/Paulo/MEFPoro/v6_static_NewHessiana_07_06_2022/build/Examples && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/BasicFEM.dir/main.cpp.o -c /home/oem/PRG/CPP/Paulo/MEFPoro/v6_static_NewHessiana_07_06_2022/Examples/main.cpp

Examples/CMakeFiles/BasicFEM.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/BasicFEM.dir/main.cpp.i"
	cd /home/oem/PRG/CPP/Paulo/MEFPoro/v6_static_NewHessiana_07_06_2022/build/Examples && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/oem/PRG/CPP/Paulo/MEFPoro/v6_static_NewHessiana_07_06_2022/Examples/main.cpp > CMakeFiles/BasicFEM.dir/main.cpp.i

Examples/CMakeFiles/BasicFEM.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/BasicFEM.dir/main.cpp.s"
	cd /home/oem/PRG/CPP/Paulo/MEFPoro/v6_static_NewHessiana_07_06_2022/build/Examples && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/oem/PRG/CPP/Paulo/MEFPoro/v6_static_NewHessiana_07_06_2022/Examples/main.cpp -o CMakeFiles/BasicFEM.dir/main.cpp.s

# Object files for target BasicFEM
BasicFEM_OBJECTS = \
"CMakeFiles/BasicFEM.dir/main.cpp.o"

# External object files for target BasicFEM
BasicFEM_EXTERNAL_OBJECTS =

Examples/BasicFEM: Examples/CMakeFiles/BasicFEM.dir/main.cpp.o
Examples/BasicFEM: Examples/CMakeFiles/BasicFEM.dir/build.make
Examples/BasicFEM: libNewFEMlib.a
Examples/BasicFEM: /usr/lib/x86_64-linux-gnu/openmpi/lib/libmpi_cxx.so
Examples/BasicFEM: /usr/lib/x86_64-linux-gnu/openmpi/lib/libmpi.so
Examples/BasicFEM: /home/oem/petsc/arch-linux2-c-debug/lib/libpetsc.so
Examples/BasicFEM: /home/oem/petsc/arch-linux2-c-debug/lib/libcmumps.a
Examples/BasicFEM: /home/oem/petsc/arch-linux2-c-debug/lib/libdmumps.a
Examples/BasicFEM: /home/oem/petsc/arch-linux2-c-debug/lib/libsmumps.a
Examples/BasicFEM: /home/oem/petsc/arch-linux2-c-debug/lib/libzmumps.a
Examples/BasicFEM: /home/oem/petsc/arch-linux2-c-debug/lib/libmumps_common.a
Examples/BasicFEM: /home/oem/petsc/arch-linux2-c-debug/lib/libpord.a
Examples/BasicFEM: /home/oem/petsc/arch-linux2-c-debug/lib/libscalapack.a
Examples/BasicFEM: /home/oem/petsc/arch-linux2-c-debug/lib/libflapack.a
Examples/BasicFEM: /home/oem/petsc/arch-linux2-c-debug/lib/libfblas.a
Examples/BasicFEM: /home/oem/petsc/arch-linux2-c-debug/lib/libptesmumps.a
Examples/BasicFEM: /home/oem/petsc/arch-linux2-c-debug/lib/libptscotchparmetis.a
Examples/BasicFEM: /home/oem/petsc/arch-linux2-c-debug/lib/libptscotch.a
Examples/BasicFEM: /home/oem/petsc/arch-linux2-c-debug/lib/libptscotcherr.a
Examples/BasicFEM: /home/oem/petsc/arch-linux2-c-debug/lib/libesmumps.a
Examples/BasicFEM: /home/oem/petsc/arch-linux2-c-debug/lib/libscotch.a
Examples/BasicFEM: /home/oem/petsc/arch-linux2-c-debug/lib/libscotcherr.a
Examples/BasicFEM: /home/oem/petsc/arch-linux2-c-debug/lib/libhdf5_hl.so
Examples/BasicFEM: /home/oem/petsc/arch-linux2-c-debug/lib/libhdf5.so
Examples/BasicFEM: /home/oem/petsc/arch-linux2-c-debug/lib/libparmetis.so
Examples/BasicFEM: /home/oem/petsc/arch-linux2-c-debug/lib/libmetis.so
Examples/BasicFEM: /usr/lib/x86_64-linux-gnu/openmpi/lib/libmpi_usempif08.so
Examples/BasicFEM: /usr/lib/x86_64-linux-gnu/openmpi/lib/libmpi_usempi_ignore_tkr.so
Examples/BasicFEM: /usr/lib/x86_64-linux-gnu/openmpi/lib/libmpi_mpifh.so
Examples/BasicFEM: /usr/lib/x86_64-linux-gnu/openmpi/lib/libmpi.so
Examples/BasicFEM: /usr/lib/gcc/x86_64-linux-gnu/9/libgfortran.so
Examples/BasicFEM: /usr/lib/x86_64-linux-gnu/libm.so
Examples/BasicFEM: /usr/lib/gcc/x86_64-linux-gnu/9/libgcc_s.so
Examples/BasicFEM: /usr/lib/x86_64-linux-gnu/libpthread.so
Examples/BasicFEM: /usr/lib/x86_64-linux-gnu/librt.so
Examples/BasicFEM: /usr/lib/gcc/x86_64-linux-gnu/9/libquadmath.so
Examples/BasicFEM: /usr/lib/gcc/x86_64-linux-gnu/9/libstdc++.so
Examples/BasicFEM: /usr/lib/x86_64-linux-gnu/libdl.so
Examples/BasicFEM: /home/oem/petsc/arch-linux2-c-debug/lib/libpetsc.so
Examples/BasicFEM: /home/oem/petsc/arch-linux2-c-debug/lib/libcmumps.a
Examples/BasicFEM: /home/oem/petsc/arch-linux2-c-debug/lib/libdmumps.a
Examples/BasicFEM: /home/oem/petsc/arch-linux2-c-debug/lib/libsmumps.a
Examples/BasicFEM: /home/oem/petsc/arch-linux2-c-debug/lib/libzmumps.a
Examples/BasicFEM: /home/oem/petsc/arch-linux2-c-debug/lib/libmumps_common.a
Examples/BasicFEM: /home/oem/petsc/arch-linux2-c-debug/lib/libpord.a
Examples/BasicFEM: /home/oem/petsc/arch-linux2-c-debug/lib/libscalapack.a
Examples/BasicFEM: /home/oem/petsc/arch-linux2-c-debug/lib/libflapack.a
Examples/BasicFEM: /home/oem/petsc/arch-linux2-c-debug/lib/libfblas.a
Examples/BasicFEM: /home/oem/petsc/arch-linux2-c-debug/lib/libptesmumps.a
Examples/BasicFEM: /home/oem/petsc/arch-linux2-c-debug/lib/libptscotchparmetis.a
Examples/BasicFEM: /home/oem/petsc/arch-linux2-c-debug/lib/libptscotch.a
Examples/BasicFEM: /home/oem/petsc/arch-linux2-c-debug/lib/libptscotcherr.a
Examples/BasicFEM: /home/oem/petsc/arch-linux2-c-debug/lib/libesmumps.a
Examples/BasicFEM: /home/oem/petsc/arch-linux2-c-debug/lib/libscotch.a
Examples/BasicFEM: /home/oem/petsc/arch-linux2-c-debug/lib/libscotcherr.a
Examples/BasicFEM: /home/oem/petsc/arch-linux2-c-debug/lib/libhdf5_hl.so
Examples/BasicFEM: /home/oem/petsc/arch-linux2-c-debug/lib/libhdf5.so
Examples/BasicFEM: /home/oem/petsc/arch-linux2-c-debug/lib/libparmetis.so
Examples/BasicFEM: /home/oem/petsc/arch-linux2-c-debug/lib/libmetis.so
Examples/BasicFEM: /usr/lib/x86_64-linux-gnu/openmpi/lib/libmpi_usempif08.so
Examples/BasicFEM: /usr/lib/x86_64-linux-gnu/openmpi/lib/libmpi_usempi_ignore_tkr.so
Examples/BasicFEM: /usr/lib/x86_64-linux-gnu/openmpi/lib/libmpi_mpifh.so
Examples/BasicFEM: /usr/lib/gcc/x86_64-linux-gnu/9/libgfortran.so
Examples/BasicFEM: /usr/lib/x86_64-linux-gnu/libm.so
Examples/BasicFEM: /usr/lib/gcc/x86_64-linux-gnu/9/libgcc_s.so
Examples/BasicFEM: /usr/lib/x86_64-linux-gnu/libpthread.so
Examples/BasicFEM: /usr/lib/x86_64-linux-gnu/librt.so
Examples/BasicFEM: /usr/lib/gcc/x86_64-linux-gnu/9/libquadmath.so
Examples/BasicFEM: /usr/lib/gcc/x86_64-linux-gnu/9/libstdc++.so
Examples/BasicFEM: /usr/lib/x86_64-linux-gnu/libdl.so
Examples/BasicFEM: Examples/CMakeFiles/BasicFEM.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/oem/PRG/CPP/Paulo/MEFPoro/v6_static_NewHessiana_07_06_2022/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable BasicFEM"
	cd /home/oem/PRG/CPP/Paulo/MEFPoro/v6_static_NewHessiana_07_06_2022/build/Examples && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/BasicFEM.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
Examples/CMakeFiles/BasicFEM.dir/build: Examples/BasicFEM

.PHONY : Examples/CMakeFiles/BasicFEM.dir/build

Examples/CMakeFiles/BasicFEM.dir/clean:
	cd /home/oem/PRG/CPP/Paulo/MEFPoro/v6_static_NewHessiana_07_06_2022/build/Examples && $(CMAKE_COMMAND) -P CMakeFiles/BasicFEM.dir/cmake_clean.cmake
.PHONY : Examples/CMakeFiles/BasicFEM.dir/clean

Examples/CMakeFiles/BasicFEM.dir/depend:
	cd /home/oem/PRG/CPP/Paulo/MEFPoro/v6_static_NewHessiana_07_06_2022/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/oem/PRG/CPP/Paulo/MEFPoro/v6_static_NewHessiana_07_06_2022 /home/oem/PRG/CPP/Paulo/MEFPoro/v6_static_NewHessiana_07_06_2022/Examples /home/oem/PRG/CPP/Paulo/MEFPoro/v6_static_NewHessiana_07_06_2022/build /home/oem/PRG/CPP/Paulo/MEFPoro/v6_static_NewHessiana_07_06_2022/build/Examples /home/oem/PRG/CPP/Paulo/MEFPoro/v6_static_NewHessiana_07_06_2022/build/Examples/CMakeFiles/BasicFEM.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : Examples/CMakeFiles/BasicFEM.dir/depend
