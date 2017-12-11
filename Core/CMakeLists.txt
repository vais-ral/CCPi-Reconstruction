#   Copyright 2017 Edoardo Pasca
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.

cmake_minimum_required (VERSION 3.0)

project(IterativeReconstructionLibray)
#https://stackoverflow.com/questions/13298504/using-cmake-with-setup-py

# The version number.

set (CIL_VERSION $ENV{CIL_VERSION} CACHE INTERNAL "Core Imaging Library version" FORCE)

## Build the regularizers package as a library
message("Creating Iterative Reconstruction CGLS as shared library")
include(GenerateExportHeader)

#include_directories(${Boost_INCLUDE_DIRS}) 
add_library(cgls_lib SHARED
	${CMAKE_CURRENT_SOURCE_DIR}/src/mpi.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/utils.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/instruments.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/results.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/voxels.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/total_v.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/parallel.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/cone.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/diamond.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/tv_core.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/p2D.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/c2D.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/timer.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/tikhonov.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/ui_calls.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/Algorithms/cgls.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/Algorithms/mlem.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/Algorithms/sirt.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/Algorithms/tv_reg.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/Readers/tiff.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/Readers/xtek.cpp
	#${CMAKE_CURRENT_SOURCE_DIR}/module/diamond_wrapper.cpp
	)

target_link_libraries(cgls_lib ${BOOST_ROOT})
include_directories(cgls_lib PUBLIC $ENV{LIBRARY_INC}/include ${CMAKE_CURRENT_SOURCE_DIR}/src/)
	

GENERATE_EXPORT_HEADER(cgls_lib)
#install (TARGETS cgls_lib DESTINATION $ENV{PREFIX}/lib)

## test
#add_executable(regularizer_test ${CMAKE_CURRENT_SOURCE_DIR}/test/test_regularizer.cpp)
#target_link_libraries (regularizer_test LINK_PUBLIC regularizers_lib)