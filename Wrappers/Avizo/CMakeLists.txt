cmake_minimum_required(VERSION 2.8)
project(CCPiReconstructionAvizo)

# Create the source files
configure_file(${CMAKE_SOURCE_DIR}/CMake/version.h.in ${CMAKE_CURRENT_SOURCE_DIR}/src/version.h)
configure_file(${CMAKE_SOURCE_DIR}/CMake/api.h.in ${CMAKE_CURRENT_SOURCE_DIR}/src/CCPiReconstruction/api.h)
	#List of function to wrap in Avizo
	set(Avizo_Function_List XTek_recon
							Parallel_Beam_recon
							NeXus_normalise
							Cone_Beam_recon
							CCPiXtekNoShiftReconstruction
							CCPiXtekAvizoPrepareFilter
							CCPiXradiaAvizoPrepareFilter		
	)

	set(Avizo_Wrapper_Exit_Function_String "")
        set(Avizo_Wrapper_Extern_Function_String "")
        set(Avizo_Wrapper_Init_Function_String "")
	foreach(fname ${Avizo_Function_List})
		set(Avizo_Wrapper_Init_Function_String "${Avizo_Wrapper_Init_Function_String}    mcInitClass_${fname}();")
		set(Avizo_Wrapper_Exit_Function_String "${Avizo_Wrapper_Exit_Function_String}    mcExitClass_${fname}();")
		set(Avizo_Wrapper_Extern_Function_String "${Avizo_Wrapper_Extern_Function_String}void mcExitClass_${fname}();")
		if(WIN32)
			set(Avizo_Wrapper_Init_Function_String "${Avizo_Wrapper_Init_Function_String}\r\n")
			set(Avizo_Wrapper_Exit_Function_String "${Avizo_Wrapper_Exit_Function_String}\r\n")
			set(Avizo_Wrapper_Extern_Function_String "${Avizo_Wrapper_Extern_Function_String}\r\n")
		else()
			set(Avizo_Wrapper_Init_Function_String "${Avizo_Wrapper_Init_Function_String}\n")
			set(Avizo_Wrapper_Exit_Function_String "${Avizo_Wrapper_Exit_Function_String}\n")
			set(Avizo_Wrapper_Extern_Function_String "${Avizo_Wrapper_Extern_Function_String}\n")
		endif()	
		set(Avizo_Wrapper_Extern_Function_String "${Avizo_Wrapper_Extern_Function_String}void mcInitClass_${fname}();")
		if(WIN32)
			set(Avizo_Wrapper_Extern_Function_String "${Avizo_Wrapper_Extern_Function_String}\r\n")
		else()
			set(Avizo_Wrapper_Extern_Function_String "${Avizo_Wrapper_Extern_Function_String}\n")
		endif()			
	endforeach()
configure_file(${CMAKE_SOURCE_DIR}/CMake/init.cpp.in ${CMAKE_CURRENT_SOURCE_DIR}/src/CCPiReconstruction/init.cpp)

# Set all the source files in the src directory
FILE(GLOB CCPiReconstructionAvizo_SRC ${PROJECT_SOURCE_DIR}/src/CCPiReconstruction/*.cpp ${PROJECT_SOURCE_DIR}/src/CCPiReconstruction/*.h)


#MESSAGE( "Avizo Source Files:" ${CCPiAvizo_SRC})
add_library(CCPiReconstructionAvizo SHARED ${CCPiReconstructionAvizo_SRC})

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/CMake)



find_package (Avizo)
IF (Avizo_FOUND)
  include_directories( ${Avizo_INCLUDE_DIRS} )
  target_link_libraries (CCPiReconstructionAvizo  ${Avizo_LIBRARIES})
  add_definitions(${Avizo_DEFINITIONS})
ENDIF (Avizo_FOUND)

#Custom definitions
add_definitions(-DPACKAGEID=amirapackage_CCPiReconstructionAvizo)
find_package(OpenMP)
if (OPENMP_FOUND)
    SET (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()
IF (WIN32)

ELSE()
set(CMAKE_CXX_FLAGS "-shared -fPIC -Wall -Wno-unknown-pragmas -mmmx -msse -msse2 -mfpmath=sse -fmessage-length=0 -fexceptions -frounding-math")
ENDIF()

include_directories(${CCPiReconstructionIterative_INCLUDE_DIRS})
target_link_libraries (CCPiReconstructionAvizo CCPiReconstructionIterative)

