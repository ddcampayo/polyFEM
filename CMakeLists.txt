# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

project( embedded_ )

cmake_minimum_required(VERSION 2.6.2)
if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}" VERSION_GREATER 2.6)
  if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION}" VERSION_GREATER 2.8.3)
    cmake_policy(VERSION 2.8.4)
  else()
    cmake_policy(VERSION 2.6)
  endif()
endif()

# These don't work ...
#set(CGAL_DIR "${HOME}/CGAL-4.6.3")
#set(CMAKE_BUILD_TYPE "Release")

find_package(CGAL QUIET COMPONENTS Core )

include( ${CGAL_USE_FILE} )

#set(EIGEN3_INCLUDE_DIR "/usr/local/include/eigen")
set(EIGEN3_INCLUDE_DIR "/usr/local/include/eigen3")

find_package(Eigen3) 

if(EIGEN3_FOUND) 
   message(STATUS "NOTICE: Eigen library found.") 
   include( ${EIGEN3_USE_FILE} ) 
else() 
   message(STATUS "NOTICE: Eigen library is not found.") 
endif() 


set(CMAKE_MODULE_PATH "/usr/local/include/eigen3/cmake/;${CMAKE_MODULE_PATH}")
#set(CMAKE_MODULE_PATH "/usr/include/eigen3/cmake/;${CMAKE_MODULE_PATH}")
#set(CHOLMOD_LIBRARIES "/usr/include/suitesparse/")

#find_package( Cholmod )

find_library(LAPACK_LIB NAMES lapack)

find_library(BLAS_LIB NAMES blas)

find_library(FFTW_LIB NAMES fftw3)
find_library(FFTW_OMP_LIB NAMES fftw3_omp)

# -lfftw3 -lfftw3_omp


#find_library(SS_LIB NAMES suitesparseconfig)

include_directories( ${CHOLMOD_INCLUDES})

#ADD_LIBRARY(${LAPACK_LIB} ${BLAS_LIB} ${CHOLMOD_LIBRARIES} )

include( CGAL_CreateSingleSourceCGALProgram )

create_single_source_cgal_program("main.cpp" "linear.cpp"   "gradient.cpp"  "nabla.cpp" "quad_coeffs.cpp" "periodic.cpp" "fields.cpp" "Delta.cpp" "move.cpp" "draw.cpp" "onto_from_mesh.cpp" "volumes.cpp" "CH_FFT.cpp" "fftwpp/fftw++.cc" )

target_link_libraries(main ${FFTW_LIB} ${FFTW_OMP_LIB})

#target_link_libraries(main ${CHOLMOD_LIBRARIES} ${LAPACK_LIB} ${BLAS_LIB} ${SS_LIB}) 
#target_link_libraries(main ${CHOLMOD_LIBRARIES} ${LAPACK_LIB} ${BLAS_LIB} )
#target_link_libraries(main ${CHOLMOD_LIBRARIES} ${LAPACK_LIB} ${BLAS_LIB} /usr/local/lib/libsuitesparseconfig.a )
