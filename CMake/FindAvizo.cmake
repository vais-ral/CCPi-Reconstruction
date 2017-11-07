# - Try to find Avizo
# Once done, this will define
#
# Avizo_FOUND - system has Avizo
# Avizo_INCLUDE_DIRS - Avizo include directories
# Avizo_LIBRARIES - link these to build for Avizo plugin
# Avizo_DEFINITIONS - definition for Avizio
	# Root Avizo Dir
set(Avizo_DIR CACHE PATH "Avizo Installed Root directory")
set(AVIZO_MAJOR_VERSION CACHE STRING "Avizo Major Version Number")
set(AVIZO_MINOR_VERSION CACHE STRING "Avizo Minor Version Number")
IF((NOT DEFINED Avizo_DIR_old OR NOT (Avizo_DIR STREQUAL Avizo_DIR_old)) AND EXISTS ${Avizo_DIR})
	# Include dir
	find_path(Avizo_INCLUDE_DIR NAMES hxcore/HxBase.h PATHS ${Avizo_DIR} PATH_SUFFIXES include)
	# Find Qt header in Avizo
	find_path(Avizo_INCLUDE_Qt_DIR NAMES QtGui/QtGui PATHS ${Avizo_INCLUDE_DIR} PATH_SUFFIXES arch-LinuxAMD64/qt NO_DEFAULT_PATH)	
#	find_path(Avizo_INCLUDE_qt_commercial_charts_DIR NAMES qt-commercial-charts PATHS ${Avizo_INCLUDE_DIR}/arch-LinuxAMD64 PATH_SUFFIXES qt-commercial-charts)
	
	#Update the Qt headers when Qt Dir is changed
	IF(NOT DEFINED Avizo_INCLUDE_Qt_DIR_old OR NOT (Avizo_INCLUDE_Qt_DIR STREQUAL Avizo_INCLUDE_Qt_DIR_old))
		
		UNSET(Avizo_INCLUDE_Qt_QtDBus_DIR CACHE)
		UNSET(Avizo_INCLUDE_Qt_QtGui_DIR CACHE)
		UNSET(Avizo_INCLUDE_Qt_QtCore_DIR CACHE)
		UNSET(Avizo_INCLUDE_Qt_QtDesigner_DIR CACHE)
		UNSET(Avizo_INCLUDE_Qt_QtNetwork_DIR CACHE)
		UNSET(Avizo_INCLUDE_Qt_QtOpenGL_DIR CACHE)
		UNSET(Avizo_INCLUDE_Qt_QtScript_DIR CACHE)
		UNSET(Avizo_INCLUDE_Qt_QtSql_DIR CACHE)
		UNSET(Avizo_INCLUDE_Qt_QtSvg_DIR CACHE)
		UNSET(Avizo_INCLUDE_Qt_QtUiTools_DIR CACHE)
		UNSET(Avizo_INCLUDE_Qt_QtXml_DIR CACHE)
		UNSET(Avizo_INCLUDE_Qt_QtTest_DIR CACHE)
		UNSET(Avizo_INCLUDE_Qt_QtXmlPatterns_DIR CACHE)
		UNSET(Avizo_INCLUDE_Qt_QtDeclarative_DIR CACHE)
		UNSET(Avizo_INCLUDE_Qt_QtWebKit_DIR CACHE)
		UNSET(Avizo_INCLUDE_Qt_qt_commercial_charts_DIR CACHE)
		set(Avizo_INCLUDE_Qt_Qt_DIR ${Avizo_INCLUDE_Qt_DIR}/Qt CACHE PATH "Qt include path")
		set(Avizo_INCLUDE_Qt_QtDBus_DIR ${Avizo_INCLUDE_Qt_DIR}/QtDBus CACHE PATH "QtDBus include path")
		set(Avizo_INCLUDE_Qt_QtGui_DIR ${Avizo_INCLUDE_Qt_DIR}/QtGui CACHE PATH "QtGui include path")
		set(Avizo_INCLUDE_Qt_QtCore_DIR ${Avizo_INCLUDE_Qt_DIR}/QtCore CACHE PATH "QtCore include path")
		set(Avizo_INCLUDE_Qt_QtDesigner_DIR ${Avizo_INCLUDE_Qt_DIR}/QtDesigner CACHE PATH "QtDesigner include path")
		set(Avizo_INCLUDE_Qt_QtNetwork_DIR ${Avizo_INCLUDE_Qt_DIR}/QtNetwork CACHE PATH "QtNetwork include path")
		set(Avizo_INCLUDE_Qt_QtOpenGL_DIR ${Avizo_INCLUDE_Qt_DIR}/QtOpenGL CACHE PATH "QtOpenGL include path")
		set(Avizo_INCLUDE_Qt_QtScript_DIR ${Avizo_INCLUDE_Qt_DIR}/QtScript CACHE PATH "QtScript include path")
		set(Avizo_INCLUDE_Qt_QtSql_DIR ${Avizo_INCLUDE_Qt_DIR}/QtSql CACHE PATH "QtSql include path")
		set(Avizo_INCLUDE_Qt_QtSvg_DIR ${Avizo_INCLUDE_Qt_DIR}/QtSvg CACHE PATH "QtSvg include path")
		set(Avizo_INCLUDE_Qt_QtTest_DIR ${Avizo_INCLUDE_Qt_DIR}/QtTest CACHE PATH "QtTest include path")
		set(Avizo_INCLUDE_Qt_QtUiTools_DIR ${Avizo_INCLUDE_Qt_DIR}/QtUiTools CACHE PATH "QtUiTools path")
		set(Avizo_INCLUDE_Qt_QtXml_DIR ${Avizo_INCLUDE_Qt_DIR}/QtXml CACHE PATH "QtXml include path")
		set(Avizo_INCLUDE_Qt_QtXmlPatterns_DIR ${Avizo_INCLUDE_Qt_DIR}/QtXmlPatterns CACHE PATH "QtXmlPatterns include path")
		set(Avizo_INCLUDE_Qt_QtDeclarative_DIR ${Avizo_INCLUDE_Qt_DIR}/QtDeclarative CACHE PATH "QtDeclarative include path")
		set(Avizo_INCLUDE_Qt_QtWebKit_DIR ${Avizo_INCLUDE_Qt_DIR}/QtWebKit CACHE PATH "QtWebKit include path")
	
	ENDIF()
	set(Avizo_INCLUDE_Qt_DIR_old ${Avizo_INCLUDE_Qt_DIR} CACHE INTERNAL "Copy of Qt Dir")

	# Find OpenInventor in Avizo
	find_path(Avizo_INCLUDE_OIV_DIR NAMES Inventor/SoType.h PATHS ${Avizo_INCLUDE_DIR} PATH_SUFFIXES arch-LinuxAMD64/oiv)
	# Find Boost header in Avizo
	find_path(Avizo_INCLUDE_Boost_DIR NAMES boost/type.hpp PATHS ${Avizo_INCLUDE_DIR}  PATH_SUFFIXES arch-LinuxAMD64)

	unset(Avizo_LIBRARY_HX_X11)
	unset(Avizo_LIBRARY_HX_Tcl)
	unset(Avizo_LIBRARY_QT_Core)
	unset(Avizo_LIBRARY_QT_DesignerComponents)
  	unset(Avizo_LIBRARY_QT_Designer)
	unset(Avizo_LIBRARY_QT_Gui)
	unset(Avizo_LIBRARY_QT_Network)
	unset(Avizo_LIBRARY_QT_OpenGL)
	unset(Avizo_LIBRARY_QT_Script)
	unset(Avizo_LIBRARY_QT_Sql)
	unset(Avizo_LIBRARY_QT_Svg)
	unset(Avizo_LIBRARY_QT_Test)
	unset(Avizo_LIBRARY_QT_Xml)
	unset(Avizo_LIBRARY_QT_Declarative)
	unset(Avizo_LIBRARY_QT_XmlPatterns)
	unset(Avizo_LIBRARY_QT_WebKit)

	# Finally Avizo libraries
	find_library(Avizo_LIBRARY NAMES hxcore PATHS ${Avizo_DIR}/lib/arch-LinuxAMD64-Optimize)
	find_library(Avizo_LIBRARY_HX_Plot NAMES hxplot PATHS ${Avizo_DIR}/lib/arch-LinuxAMD64-Optimize)
	find_library(Avizo_LIBRARY_HX_Time NAMES hxtime PATHS ${Avizo_DIR}/lib/arch-LinuxAMD64-Optimize)
	find_library(Avizo_LIBRARY_HX_Surface NAMES hxsurface PATHS ${Avizo_DIR}/lib/arch-LinuxAMD64-Optimize)
	find_library(Avizo_LIBRARY_HX_Color NAMES hxcolor PATHS ${Avizo_DIR}/lib/arch-LinuxAMD64-Optimize)
	find_library(Avizo_LIBRARY_HX_Field NAMES hxfield PATHS ${Avizo_DIR}/lib/arch-LinuxAMD64-Optimize)
	find_library(Avizo_LIBRARY_HX_Spreadsheet NAMES hxspreadsheet PATHS ${Avizo_DIR}/lib/arch-LinuxAMD64-Optimize)
	find_library(Avizo_LIBRARY_HX_AmiraMesh NAMES amiramesh PATHS ${Avizo_DIR}/lib/arch-LinuxAMD64-Optimize)
	find_library(Avizo_LIBRARY_HX_Mclib NAMES mclib PATHS ${Avizo_DIR}/lib/arch-LinuxAMD64-Optimize)
	find_library(Avizo_LIBRARY_HX_InventorBase NAMES InventorBase PATHS ${Avizo_DIR}/lib/arch-LinuxAMD64-Optimize)
	find_library(Avizo_LIBRARY_HX_InventorGL NAMES InventorGL PATHS ${Avizo_DIR}/lib/arch-LinuxAMD64-Optimize)
	find_library(Avizo_LIBRARY_HX_Inventor NAMES Inventor PATHS ${Avizo_DIR}/lib/arch-LinuxAMD64-Optimize)
	find_library(Avizo_LIBRARY_HX_InventorQt4 NAMES InventorQt4 PATHS ${Avizo_DIR}/lib/arch-LinuxAMD64-Optimize)
	find_library(Avizo_LIBRARY_HX_HardCopy NAMES HardCopy PATHS ${Avizo_DIR}/lib/arch-LinuxAMD64-Optimize)
	find_library(Avizo_LIBRARY_HX_MeshViz NAMES MeshViz PATHS ${Avizo_DIR}/lib/arch-LinuxAMD64-Optimize)
	find_library(Avizo_LIBRARY_HX_MeshVizExtractor NAMES MeshVizExtractor PATHS ${Avizo_DIR}/lib/arch-LinuxAMD64-Optimize)
	find_library(Avizo_LIBRARY_HX_MeshVizImpl NAMES MeshVizImpl PATHS ${Avizo_DIR}/lib/arch-LinuxAMD64-Optimize)
	find_library(Avizo_LIBRARY_HX_MeshVizDataMapping NAMES MeshVizDataMapping PATHS ${Avizo_DIR}/lib/arch-LinuxAMD64-Optimize)
	find_library(Avizo_LIBRARY_HX_VolumeViz NAMES VolumeViz PATHS ${Avizo_DIR}/lib/arch-LinuxAMD64-Optimize)
	find_library(Avizo_LIBRARY_HX_LDM NAMES LDM PATHS ${Avizo_DIR}/lib/arch-LinuxAMD64-Optimize)
	find_library(Avizo_LIBRARY_HX_X11 NAMES X11 PATHS ${Avizo_DIR}/lib/arch-LinuxAMD64-Optimize)
	find_library(Avizo_LIBRARY_HX_Tcl NAMES tcl8.5 PATHS ${Avizo_DIR}/lib/arch-LinuxAMD64-Optimize NO_DEFAULT_PATH)

	find_library(Avizo_LIBRARY_QT_Core NAMES QtCore PATHS ${Avizo_DIR}/lib/arch-LinuxAMD64-Optimize NO_DEFAULT_PATH) 
	find_library(Avizo_LIBRARY_QT_DesignerComponents NAMES QtDesignerComponents PATHS ${Avizo_DIR}/lib/arch-LinuxAMD64-Optimize NO_DEFAULT_PATH)
	find_library(Avizo_LIBRARY_QT_Designer NAMES QtDesigner PATHS ${Avizo_DIR}/lib/arch-LinuxAMD64-Optimize NO_DEFAULT_PATH)
	find_library(Avizo_LIBRARY_QT_Gui NAMES QtGui PATHS ${Avizo_DIR}/lib/arch-LinuxAMD64-Optimize NO_DEFAULT_PATH)
	find_library(Avizo_LIBRARY_QT_Network NAMES QtNetwork PATHS ${Avizo_DIR}/lib/arch-LinuxAMD64-Optimize NO_DEFAULT_PATH)
	find_library(Avizo_LIBRARY_QT_OpenGL NAMES QtOpenGL PATHS ${Avizo_DIR}/lib/arch-LinuxAMD64-Optimize NO_DEFAULT_PATH)
	find_library(Avizo_LIBRARY_QT_Script NAMES QtScript PATHS ${Avizo_DIR}/lib/arch-LinuxAMD64-Optimize NO_DEFAULT_PATH)
	find_library(Avizo_LIBRARY_QT_Sql NAMES QtSql PATHS ${Avizo_DIR}/lib/arch-LinuxAMD64-Optimize NO_DEFAULT_PATH)
	find_library(Avizo_LIBRARY_QT_Svg NAMES QtSvg PATHS ${Avizo_DIR}/lib/arch-LinuxAMD64-Optimize NO_DEFAULT_PATH)
	find_library(Avizo_LIBRARY_QT_Test NAMES QtTest PATHS ${Avizo_DIR}/lib/arch-LinuxAMD64-Optimize NO_DEFAULT_PATH)
	find_library(Avizo_LIBRARY_QT_Xml NAMES QtXml PATHS ${Avizo_DIR}/lib/arch-LinuxAMD64-Optimize NO_DEFAULT_PATH)
	find_library(Avizo_LIBRARY_QT_Declarative NAMES QtDeclarative PATHS ${Avizo_DIR}/lib/arch-LinuxAMD64-Optimize NO_DEFAULT_PATH)
	find_library(Avizo_LIBRARY_QT_XmlPatterns NAMES QtXmlPatterns PATHS ${Avizo_DIR}/lib/arch-LinuxAMD64-Optimize NO_DEFAULT_PATH)
	find_library(Avizo_LIBRARY_QT_WebKit NAMES QtWebKit PATHS ${Avizo_DIR}/lib/arch-LinuxAMD64-Optimize NO_DEFAULT_PATH)

ENDIF()

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(Avizo DEFAULT_MSG Avizo_LIBRARY Avizo_INCLUDE_DIR)
mark_as_advanced(	Avizo_LIBRARY_QT_Core
		  	Avizo_LIBRARY_QT_DesignerComponents
		  	Avizo_LIBRARY_QT_Designer
		  	Avizo_LIBRARY_QT_Gui
		  	Avizo_LIBRARY_QT_Network
		  	Avizo_LIBRARY_QT_OpenGL
		  	Avizo_LIBRARY_QT_Script
		  	Avizo_LIBRARY_QT_Sql
		  	Avizo_LIBRARY_QT_Svg
		  	Avizo_LIBRARY_QT_Test
		  	Avizo_LIBRARY_QT_Xml
		  	Avizo_LIBRARY_QT_Declarative
		  	Avizo_LIBRARY_QT_XmlPatterns
		  	Avizo_LIBRARY_QT_WebKit

			Avizo_LIBRARY_HX_Plot
			Avizo_LIBRARY_HX_Time
			Avizo_LIBRARY_HX_Surface
			Avizo_LIBRARY_HX_Color
			Avizo_LIBRARY_HX_Field
			Avizo_LIBRARY_HX_Spreadsheet
			Avizo_LIBRARY_HX_AmiraMesh
			Avizo_LIBRARY_HX_Mclib
			Avizo_LIBRARY_HX_InventorBase
			Avizo_LIBRARY_HX_InventorGL
			Avizo_LIBRARY_HX_Inventor
			Avizo_LIBRARY_HX_InventorQt4
			Avizo_LIBRARY_HX_HardCopy
			Avizo_LIBRARY_HX_MeshViz
			Avizo_LIBRARY_HX_MeshVizExtractor
			Avizo_LIBRARY_HX_MeshVizImpl
			Avizo_LIBRARY_HX_MeshVizDataMapping
			Avizo_LIBRARY_HX_VolumeViz
			Avizo_LIBRARY_HX_LDM
			Avizo_LIBRARY_HX_X11
			Avizo_LIBRARY_HX_Tcl
Avizo_INCLUDE_Boost_DIR Avizo_INCLUDE_OIV_DIR Avizo_INCLUDE_Qt_DIR Avizo_INCLUDE_Qt_Qt_DIR Avizo_INCLUDE_Qt_QtCore_DIR Avizo_INCLUDE_Qt_QtDBus_DIR Avizo_INCLUDE_Qt_QtDesigner_DIR Avizo_INCLUDE_Qt_QtGui_DIR  Avizo_INCLUDE_Qt_QtNetwork_DIR Avizo_INCLUDE_Qt_QtOpenGL_DIR Avizo_INCLUDE_Qt_QtScript_DIR Avizo_INCLUDE_Qt_QtSql_DIR Avizo_INCLUDE_Qt_QtSvg_DIR Avizo_INCLUDE_Qt_QtTest_DIR Avizo_INCLUDE_Qt_QtUiTools_DIR Avizo_INCLUDE_Qt_QtXml_DIR Avizo_INCLUDE_Qt_QtXmlPatterns_DIR Avizo_INCLUDE_Qt_QtDeclarative_DIR Avizo_INCLUDE_Qt_QtWebKit_DIR
)


set(Avizo_INCLUDE_DIRS ${Avizo_INCLUDE_DIR} ${Avizo_INCLUDE_Boost_DIR} ${Avizo_INCLUDE_OIV_DIR} ${Avizo_INCLUDE_Qt_DIR} ${Avizo_INCLUDE_Qt_QtCore_DIR} ${Avizo_INCLUDE_Qt_QtDBus_DIR} ${Avizo_INCLUDE_Qt_QtDesigner_DIR} ${Avizo_INCLUDE_Qt_QtGui_DIR}  ${Avizo_INCLUDE_Qt_QtNetwork_DIR} ${Avizo_INCLUDE_Qt_QtOpenGL_DIR} ${Avizo_INCLUDE_Qt_QtScript_DIR} ${Avizo_INCLUDE_Qt_QtSql_DIR} ${Avizo_INCLUDE_Qt_QtSvg_DIR} ${Avizo_INCLUDE_Qt_QtTest_DIR} ${Avizo_INCLUDE_Qt_QtUiTools_DIR} ${Avizo_INCLUDE_Qt_QtXml_DIR} ${Avizo_INCLUDE_Qt_QtXmlPatterns_DIR} ${Avizo_INCLUDE_Qt_QtDeclarative_DIR} ${Avizo_INCLUDE_Qt_QtWebKit_DIR} )#${Avizo_INCLUDE_qt_commercial_charts_DIR} )
set(Avizo_LIBRARIES    	${Avizo_LIBRARY}
		  	${Avizo_LIBRARY_QT_Core}
		  	${Avizo_LIBRARY_QT_DesignerComponents}
		  	${Avizo_LIBRARY_QT_Designer}
		  	${Avizo_LIBRARY_QT_Gui}
		  	${Avizo_LIBRARY_QT_Network}
		  	${Avizo_LIBRARY_QT_OpenGL}
		  	${Avizo_LIBRARY_QT_Script}
		  	${Avizo_LIBRARY_QT_Sql}
		  	${Avizo_LIBRARY_QT_Svg}
		  	${Avizo_LIBRARY_QT_Test}
		  	${Avizo_LIBRARY_QT_Xml}
		  	${Avizo_LIBRARY_QT_Declarative}
		  	${Avizo_LIBRARY_QT_XmlPatterns}
		  	${Avizo_LIBRARY_QT_WebKit}

			${Avizo_LIBRARY_HX_Plot}
			${Avizo_LIBRARY_HX_Time}
			${Avizo_LIBRARY_HX_Surface}
			${Avizo_LIBRARY_HX_Color}
			${Avizo_LIBRARY_HX_Field}
			${Avizo_LIBRARY_HX_Spreadsheet}
			${Avizo_LIBRARY_HX_AmiraMesh}
			${Avizo_LIBRARY_HX_Mclib}
			${Avizo_LIBRARY_HX_InventorBase}
			${Avizo_LIBRARY_HX_InventorGL}
			${Avizo_LIBRARY_HX_Inventor}
			${Avizo_LIBRARY_HX_InventorQt4}
			${Avizo_LIBRARY_HX_HardCopy}
			${Avizo_LIBRARY_HX_MeshViz}
			${Avizo_LIBRARY_HX_MeshVizExtractor}
			${Avizo_LIBRARY_HX_MeshVizImpl}
			${Avizo_LIBRARY_HX_MeshVizDataMapping}
			${Avizo_LIBRARY_HX_VolumeViz}
			${Avizo_LIBRARY_HX_LDM}
			${Avizo_LIBRARY_HX_X11}
			${Avizo_LIBRARY_HX_Tcl}		
)
set(Avizo_DEFINITIONS  	-DHX_CORE_GIT -DQT_RESTRICTED_CAST_FROM_ASCII -DHX_INTERVAL_WARN -DHX_INTERNAL_CLASS_WARN
		       	-DBOOST_FILESYSTEM_VERSION=2 -DQT_THREAD_SUPPORT -DQT_CLEAN_NAMESPACE -DAMIRA_RELEASE
		       	-DFORTRAN_UNDERLINE -D_REENTRANT -DAMIRA64 -DHX_OS_LINUX -DHX_ARCH_LINUXAMD64 -DHX_LITTLE_ENDIAN
			-DHX_NO_MATHF -DHX_HAS_X11 -DUNICODE_CHECK -DOIV_MULTI_THREADS -DUSE_NON_CONST)
        
set(Avizo_DIR_old ${Avizo_DIR} CACHE INTERNAL "Copy of Avizo Root Dir")
