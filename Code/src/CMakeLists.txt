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
				
			)
			
IF( CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64" )
  SET_TARGET_PROPERTIES(CCPiReconstructionIterative PROPERTIES COMPILE_FLAGS "-fPIC")
ENDIF( CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64" )
			


TARGET_LINK_LIBRARIES(CCPiReconstructionIterative  ${CCPI_ITK_LIBS} ${CCPI_VTK_LIBS})							
SET_TARGET_PROPERTIES(CCPiReconstructionIterative  PROPERTIES PUBLIC_HEADER "CCPiAccessibleVolumeInputImages.h;CCPiConsoleUserInterface.h;CCPiImageData.h;CCPiLabelQuantificationResult.h;CCPiParticleTracker.h;CCPiTrack.h;CCPiAccessibleVolumeITKImpl.h;CCPiDefines.h;CCPiLabelQuantificationITKImpl.h;CCPiParticle.h;CCPiSimpleHistogramThresholdingITKImpl.h;CCPiUserApplicationInterface.h;")							


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


