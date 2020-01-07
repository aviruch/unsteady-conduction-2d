# Header
cmake_minimum_required(VERSION 2.8)
project(unsteady-conduction-2d)
enable_language(Fortran)
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_VERBOSE_MAKEFILE TRUE)

# Append .f90 files
file(GLOB_RECURSE SRCS RELATIVE ${PROJECT_SOURCE_DIR} *.f90)

# Compile
add_executable(amain amain.f90 ${SRCS})