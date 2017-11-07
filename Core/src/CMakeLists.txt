cmake_minimum_required(VERSION 2.8.8 FATAL_ERROR)
PROJECT (CCPiReconstructionIterative)


# Allow the developer to select if Dynamic or Static libraries are built
OPTION (BUILD_SHARED_LIBS "Build Shared Libraries" OFF)
# Set the LIB_TYPE variable to STATIC
SET (LIB_TYPE STATIC)
IF (BUILD_SHARED_LIBS)
  # User wants to build Dynamic Libraries, so change the LIB_TYPE variable to CMake keyword 'SHARED'
  SET (LIB_TYPE SHARED)
ENDIF (BUILD_SHARED_LIBS)

ADD_LIBRARY(CCPiReconstructionIterative ${LIB_TYPE} 
										"mpi.cpp"
										"utils.cpp"
										"instruments.cpp"
										"results.cpp"
										"voxels.cpp"
										"Algorithms/cgls.cpp"
										"Algorithms/mlem.cpp"
										"Algorithms/sirt.cpp"
										"total_v.cpp"
										"parallel.cpp"
										"cone.cpp"
										"diamond.cpp"
										"Readers/xtek.cpp"
										"Algorithms/tv_reg.cpp"
										"tv_core.cpp"
										"p2D.cpp"
										"c2D.cpp"
										"Readers/tiff.cpp"
										"timer.cpp"
										"tikhonov.cpp"
										"ui_calls.cpp"
			)
			
#Boost Libraries
SET(Boost_USE_STATIC_LIBS        OFF) # only find static libs
SET(Boost_USE_MULTITHREADED      ON)
find_package(Boost 1.63.0 COMPONENTS date_time filesystem system)
if(Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIRS})
  target_link_libraries(CCPiReconstructionIterative ${Boost_LIBRARIES})
endif()
message(${Boost_INCLUDE_DIRS})
#OpenMP Settings
find_package(OpenMP)
if (OPENMP_FOUND)
    SET (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()
			
#Include iterative directories
include_directories("Readers" "Algorithms" ".")
IF( CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64" )
  SET_TARGET_PROPERTIES(CCPiReconstructionIterative PROPERTIES COMPILE_FLAGS "-fPIC")
ENDIF( CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64" )
			
							
SET_TARGET_PROPERTIES(CCPiReconstructionIterative  PROPERTIES PUBLIC_HEADER "base_types.hpp;CCPiUserApplicationInterface.h;")							


INSTALL(TARGETS CCPiReconstructionIterative EXPORT CCPiReconstructionTargets
        LIBRARY  DESTINATION "${INSTALL_LIB_DIR}" COMPONENT shlib
	PUBLIC_HEADER DESTINATION "${INSTALL_INCLUDE_DIR}" COMPONENT dev
	RUNTIME  DESTINATION "${INSTALL_BIN_DIR}" COMPONENT bin
	ARCHIVE  DESTINATION "${INSTALL_LIB_DIR}" COMPONENT lib
       )


############Expose include files
set(${PROJECT_NAME}_INCLUDE_DIRS ${PROJECT_SOURCE_DIR}
                                 ${C_INCLUDE_DIRS}
    CACHE INTERNAL "${PROJECT_NAME}: Include Directories" FORCE)


