#############################################################################
# CMakeLists.txt file for Building CCPiReconstruction iterative 
############################################################################
cmake_minimum_required(VERSION 2.8)
project(CCPiReconstruction)

set(CCPIRECONSTRUCTION_MAJOR_VERSION 0)
set(CCPIRECONSTRUCTION_MINOR_VERSION 9)
set(CCPIRECONSTRUCTION_PATCH_VERSION 1)
set(CCPIRECONSTRUCTION_VERSION 
	${CCPIRECONSTRUCTION_MAJOR_VERSION}.${CCPIRECONSTRUCTION_MINOR_VERSION}.${CCPIRECONSTRUCTION_PATCH_VERSION})

set(INSTALL_LIB_DIR lib CACHE PATH "Installation directory for libraries")
set(INSTALL_BIN_DIR bin CACHE PATH "Installation directory for executables")
set(INSTALL_INCLUDE_DIR include CACHE PATH "Installation directory for header files")

if(WIN32 AND NOT CYGWIN)
	SET(DEF_INSTALL_CMAKE_DIR cmake)
else()
	SET(DEF_INSTALL_CMAKE_DIR lib/cmake/CCPiReconstruction)
endif()
set(INSTALL_CMAKE_DIR ${DEF_INSTALL_CMAKE_DIR} CACHE PATH
	"Installation directory for CMake files")

#Make relative paths absolute
foreach(p LIB BIN INCLUDE CMAKE)
	set(var INSTALL_{p}_DIR)
	if(NOT IS_ABSOLUTE "${${var}}")
		set(${var} "${CMAKE_INSTALL_PREFIX}/${${var}}")
	endif()
endforeach()


#setup include directories
include_directories(
	"${PROJECT_SOURCE_DIR}"
	"${PROJECT_BINARY_DIR}")

# Add sub-directories
add_subdirectory(Code/src)

OPTION(BUILD_AVIZO_PLUGIN "Build Avizo plugin" OFF)
IF(BUILD_AVIZO_PLUGIN)
	add_subdirectory(Wrappers/Avizo)
ENDIF()


SET(CCPIRECONSTRUCTION_BUILD_LIBRARIES "CCPiReconstructionIterative")
export(TARGETS CCPiReconstructionIterative
	FILE "${PROJECT_BINARY_DIR}/CCPiReconstructionTargets.cmake")

export(PACKAGE CCPiReconstruction)
file(RELATIVE_PATH REL_INCLUDE_DIR "${INSTALL_CMAKE_DIR}"
	"${INSTALL_INCLUDE_DIR}")
# for build tree
set(CONF_INCLUDE_DIRS "${PROJECT_SOURCE_DIR}" "${PROJECT_BINARY_DIR}")
configure_file(${CMAKE_SOURCE_DIR}/CMake/CCPiReconstructionConfig.cmake.in  
		"${PROJECT_BINARY_DIR}/CCPiReconstructionConfig.cmake" @ONLY)
# for install tree
set(CONF_INCLUDE_DIRS "\${CCPIRECONSTRUCTION_CMAKE_DIR}/${REL_INCLUDE_DIR}")
configure_file(${CMAKE_SOURCE_DIR}/CMake/CCPiReconstructionConfig.cmake.in
		"${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CCPiReconstructionConfig.cmake" @ONLY)
# for both
configure_file(${CMAKE_SOURCE_DIR}/CMake/CCPiReconstructionConfigVersion.cmake.in
		"${PROJECT_BINARY_DIR}/CCPiReconstructionConfigVersion.cmake" @ONLY)

#Install
install(FILES
	"${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CCPiReconstructionConfig.cmake"
	"${PROJECT_BINARY_DIR}/CCPiReconstructionConfigVersion.cmake"
	DESTINATION "${INSTALL_CMAKE_DIR}" COMPONENT dev)

# Install the export set
install(EXPORT CCPiReconstructionTargets DESTINATION
	"${INSTALL_CMAKE_DIR}" COMPONENT dev)		

