# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

project( pFEM )

cmake_minimum_required(VERSION 2.6.2)

find_package(CGAL REQUIRED QUIET OPTIONAL_COMPONENTS Core )

find_package(Eigen3)
#include_directories(${EIGEN3_INCLUDE_DIR})
include( ${EIGEN3_USE_FILE} )

#find_library(LAPACK_LIB NAMES lapack)

#find_library(BLAS_LIB NAMES blas)

find_library(FFTW_LIB NAMES fftw3)
#find_library(FFTW_OMP_LIB NAMES fftw3_omp)


add_executable( pFEM
	main.cpp linear.cpp gradient.cpp nabla.cpp
	quad_coeffs.cpp
	periodic.cpp fields.cpp Delta.cpp move.cpp
	draw.cpp onto_from_mesh.cpp volumes.cpp CH_FFT.cpp
	fftwpp/fftw++.cc )

target_link_libraries( pFEM   ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES}   ${FFTW_LIB} ${FFTW_OMP_LIB})

#create_single_source_cgal_program( "main.cpp" )

