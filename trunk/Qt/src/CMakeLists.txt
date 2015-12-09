cmake_minimum_required(VERSION 2.8.8 FATAL_ERROR)

PROJECT(CCPiReconstructionQt)

FIND_PACKAGE(Qt4 REQUIRED)
FIND_PACKAGE( Boost 1.55 COMPONENTS system filesystem REQUIRED )
FIND_PACKAGE(OpenMP)
FIND_PACKAGE(TIFF REQUIRED)

IF (OPENMP_FOUND)
	SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")	
ENDIF()	

IF(TIFF_FOUND)
	ADD_DEFINITIONS(-DHAS_TIFF)
ENDIF(TIFF_FOUND)

SET(CCPiReconstructionQt_SOURCES ReconstructionWizardImpl.cpp UserInterfaceQT.cpp ProgressWizardPage.cpp SavingWizardPage.cpp)
SET(CCPiReconstructionQt_HEADERS ReconstructionWizardImpl.h UserInterfaceQT.h ProgressWizardPage.h SavingWizardPage.h)
SET(CCPiReconstruction_FORMS ReconstructionWizard.ui)

QT4_WRAP_CPP(CCPiReconstructionQt_HEADERS_MOC ${CCPiReconstructionQt_HEADERS})
QT4_WRAP_UI(CCPiReconstructionQt_FORMS_HEADERS ${CCPiReconstruction_FORMS})

INCLUDE(${QT_USE_FILE})
ADD_DEFINITIONS(${QT_DEFINITIONS})
include_directories( ${Boost_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR} ../../src ${CMAKE_CURRENT_BINARY_DIR} ${TIFF_INCLUDE_DIR})
link_directories(${Boost_LIBRARY_DIRS})
if(UNIX AND NOT APPLE)
	string(TOLOWER ${PROJECT_NAME} PROJECT_NAME_LOWERCASE)
	set(BIN_INSTALL_DIR "bin")
	set(DOC_INSTALL_DIR "share/doc/${PROJECT_NAME_LOWERCASE}/")
else()
	set(BIN_INSTALL_DIR ".")
	set(DOC_INSTALL_DIR ".")
endif()
SET(CCPiReconstruction_SOURCES 	../../src/cgls.cpp
								../../src/c2D.cpp
								../../src/cone.cpp
								../../src/diamond.cpp
								../../src/instruments.cpp
								../../src/mlem.cpp
								../../src/p2D.cpp
								../../src/parallel.cpp
								../../src/mpi.cpp
								../../src/tikhonov.cpp
								../../src/sirt.cpp
								../../src/tv_core.cpp
								../../src/tv_reg.cpp
								../../src/voxels.cpp
								../../src/xtek.cpp
								../../src/timer.cpp
								../../src/utils.cpp
								../../src/total_v.cpp
								../../src/tiff.cpp
								../../src/results.cpp
								)
								
ADD_EXECUTABLE(CCPiReconstructionQt ${CCPiReconstructionQt_SOURCES} 
				${CCPiReconstructionQt_FORMS_HEADERS} 
				${CCPiReconstructionQt_HEADERS_MOC} ${CCPiReconstruction_SOURCES} ReconstructionImpl.cpp ui_calls.cpp main.cpp)
IF(WIN32)
   SET_TARGET_PROPERTIES(CCPiReconstructionQt PROPERTIES LINK_FLAGS_RELEASE "/SUBSYSTEM:WINDOWS /ENTRY:mainCRTStartup")
ENDIF(WIN32)				
TARGET_LINK_LIBRARIES(CCPiReconstructionQt ${QT_LIBRARIES} ${Boost_FILESYSTEM_LIBRARY} ${Boost_SYSTEM_LIBRARY} ${TIFF_LIBRARIES})


#--------------------------------------------------------------------------------
# Now the installation stuff below
#--------------------------------------------------------------------------------
SET(plugin_dest_dir bin)
SET(qtconf_dest_dir bin)
SET(APPS "\${CMAKE_INSTALL_PREFIX}/bin/CCPiReconstructionQt")
IF(APPLE)
  SET(plugin_dest_dir CCPiReconstructionQt.app/Contents/MacOS)
  SET(qtconf_dest_dir CCPiReconstructionQt.app/Contents/Resources)
  SET(APPS "\${CMAKE_INSTALL_PREFIX}/CCPiReconstructionQt.app")
ENDIF(APPLE)
IF(WIN32)
  SET(APPS "\${CMAKE_INSTALL_PREFIX}/bin/CCPiReconstructionQt.exe")
ENDIF(WIN32)



INSTALL(TARGETS CCPiReconstructionQt
    BUNDLE DESTINATION . COMPONENT Runtime
    RUNTIME DESTINATION bin COMPONENT Runtime
)
   
SET(CPACK_PACKAGE_NAME "CCPiReconstructionQt")
SET(CPACK_PACKAGE_VENDOR "CCPi")
SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "CCPi Reconstruction GUI using QT")
SET(CPACK_PACKAGE_VERSION "1.0.0")
SET(CPACK_PACKAGE_MAJOR "1")
SET(CPACK_VERSION_MINOR "0")
SET(CPACK_VERSION_PATCH "0")
SET(CPACK_PACKAGE_INSTALL_DIRECTORY "CCPi")
SET(CPACK_COMPONENT_APPLICATIONS_DISPLAY_NAME "CCPi Reconstruction Qt")

SET(CPACK_COMPONENTS_ALL Runtime)

INCLUDE(InstallRequiredSystemLibraries)
#--------------------------------------------------------------------------------
# Use BundleUtilities to get all other dependencies for the application to work.
# It takes a bundle or executable along with possible plugins and inspects it
# for dependencies.  If they are not system dependencies, they are copied.

# directories to look for dependencies
SET(DIRS ${QT_LIBRARY_DIRS} ${Boost_LIBRARY_DIRS} ${TIFF_LIBRARY})

# Now the work of copying dependencies into the bundle/package
# The quotes are escaped and variables to use at install time have their $ escaped
# An alternative is the do a configure_file() on a script and use install(SCRIPT  ...).
# Note that the image plugins depend on QtSvg and QtXml, and it got those copied
# over.
INSTALL(CODE "
    include(BundleUtilities)
    fixup_bundle(\"${APPS}\"  \"\" \"${DIRS}\")
    " COMPONENT Runtime)

IF(WIN32)
	INSTALL(FILES ${CMAKE_INSTALL_SYSTEM_RUNTIME_LIBS} DESTINATION bin COMPONENT Runtime)
	IF (OPENMP_FOUND AND (MSVC_VERSION EQUAL 1600))
		INSTALL(FILES "C:/Program Files (x86)/Microsoft Visual Studio 10.0/VC/redist/x64/Microsoft.VC100.OPENMP/vcomp100.dll" DESTINATION bin COMPONENT Runtime)
	ENDIF()
ENDIF(WIN32)

IF(APPLE)
# To Create a package, one can run "cpack -G DragNDrop CPackConfig.cmake" on Mac OS X
# where CPackConfig.cmake is created by including CPack
# And then there's ways to customize this as well
set(CPACK_BINARY_DRAGNDROP ON)
ENDIF()
INCLUDE(CPack)
