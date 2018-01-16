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

project(IterativeReconstructionLibrary)
#https://stackoverflow.com/questions/13298504/using-cmake-with-setup-py

# The version number.

set (CIL_VERSION $ENV{CIL_VERSION} CACHE INTERNAL "Core Imaging Library version" FORCE)
#set (LIBRARY_LIB C:/Apps/Miniconda2/envs/cil/Library/lib)
#set (LIBRARY_INC C:/Apps/Miniconda2/envs/cil/Library/)

# conda orchestrated build

message("CIL_VERSION ${CIL_VERSION}")
#include(find_package)
if (WIN32)
	#set (BOOST_ROOT ${CONDA_PREFIX}\Library\lib)
endif()
message("BOOST_ROOT " ${BOOST_ROOT})
#set(Boost_NO_SYSTEM_PATHS ON)
find_package(Boost  REQUIRED
	COMPONENTS python numpy)

## Build the regularizers package as a library
message("Creating Iterative Reconstruction CGLS as shared library")
message("CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS}")

set(CMAKE_BUILD_TYPE "Release")

if(WIN32)
 set (FLAGS "/DWIN32 /EHsc /DBOOST_ALL_NO_LIB /openmp /DHAS_TIFF /DCCPiReconstructionIterative_EXPORTS")
 
 #set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /VERBOSE:LIB")
# [100%] Linking CXX shared library cgls_lib.dll
#   Creating library cgls_lib.lib and object cgls_lib.exp
#LINK : warning LNK4098: defaultlib 'MSVCRT' conflicts with use of other libs; use #/NODEFAULTLIB:library
#   Creating library cgls_lib.lib and object cgls_lib.exp
#LINK : warning LNK4098: defaultlib 'MSVCRT' conflicts with use of other libs; use #/NODEFAULTLIB:library
#[100%] Built target cgls_lib
  set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /NODEFAULTLIB:MSVCRT.lib")
  
  set (EXTRA_LIBRARIES 
      ${Boost_LIBRARIES}
		#"boost_python3-vc140-mt-1_64" 
		#"boost_numpy3-vc140-mt-1_64"
		"tiff"
		)
		
  message("library lib: ${LIBRARY_LIB}")
  
elseif(UNIX)
   set (FLAGS "-fopenmp -O2 -funsigned-char -Wall  -Wl,--no-undefined -DHAS_TIFF -DCCPiReconstructionIterative_EXPORTS")  
   set (EXTRA_LIBRARIES 
		${Boost_LIBRARIES}
		#"boost_python3" 
		#"boost_numpy3"
		"tiff"
		"gomp"
		)
   
endif()

link_directories(${LIBRARY_LIB})
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${FLAGS}")
message("CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS}")

add_library(cilrec SHARED
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

target_link_libraries(cilrec ${EXTRA_LIBRARIES} )

include_directories(cilrec PUBLIC 
                      ${LIBRARY_INC}/include 
		      ${CMAKE_CURRENT_SOURCE_DIR}/src/  )
#
if (UNIX)
message ("I'd install into ${CMAKE_INSTALL_PREFIX}/lib")
install(TARGETS cilrec 
	LIBRARY DESTINATION lib
	CONFIGURATIONS ${CMAKE_BUILD_TYPE} 
	)
elseif(WIN32)
message ("I'd install into ${CMAKE_INSTALL_PREFIX} lib bin")
  install(TARGETS cilrec 
	RUNTIME DESTINATION bin
	ARCHIVE DESTINATION lib
	CONFIGURATIONS ${CMAKE_BUILD_TYPE} 
	)
endif()
## test
#add_executable(regularizer_test ${CMAKE_CURRENT_SOURCE_DIR}/test/test_regularizer.cpp)
