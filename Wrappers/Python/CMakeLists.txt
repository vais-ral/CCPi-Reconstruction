#   Copyright 2018 Edoardo Pasca
cmake_minimum_required (VERSION 3.0)

project(ReconstructionPython)
#https://stackoverflow.com/questions/13298504/using-cmake-with-setup-py

# The version number.
set (CIL_VERSION $ENV{CIL_VERSION} CACHE INTERNAL "Core Imaging Library version" FORCE)

message("Building python wrappers")
# conda orchestrated build
message("CIL_VERSION ${CIL_VERSION}")
#include (GenerateExportHeader)

find_package(PythonInterp REQUIRED)

find_package(Boost REQUIRED COMPONENTS python36 numpy36)

if (Boost_FOUND)
  message("Boost version " ${Boost_VERSION})
  message("Boost include dir " ${Boost_INCLUDE_DIRS})
  message("Boost library dir " ${Boost_LIBRARY_DIRS})
  message("Boost libraries " ${Boost_LIBRARIES})
endif()

message ("found " ${BOOST_PYTHON_LIB})
message ("found " ${BOOST_NUMPY_LIB})

find_package(OpenMP)
	
## Build the regularizers package as a library
message("Creating Reconstruction Wrapper")

message("CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS}")

set(CMAKE_BUILD_TYPE "Release")

if (WIN32)
  list(GET Boost_LIBRARIES 0 place )
  get_filename_component(BOOST_PYTHON_LIB ${place} NAME_WE )
  list(GET Boost_LIBRARIES 1 place )
  get_filename_component(BOOST_NUMPY_LIB ${place} NAME_WE )
else()
  # on linux the library looks like libboost_numpy3.so:
  # we need to get rid of the lib at the beginning
  list(GET Boost_LIBRARIES 0 place )
  get_filename_component(place2 ${place} NAME_WE )
  string(REGEX REPLACE "^lib(.*)" "\\1" BOOST_PYTHON_LIB "${place2}")
  
  list(GET Boost_LIBRARIES 1 place )
  get_filename_component(place2 ${place} NAME_WE )
  string(REGEX REPLACE "^lib(.*)" "\\1" BOOST_NUMPY_LIB "${place2}")
endif()



set(LIBRARY_INC ${CMAKE_SOURCE_DIR}/src/include)
set(LIBRARY_LIB ${CMAKE_BINARY_DIR}/)

set(SETUP_PY_IN "${CMAKE_CURRENT_SOURCE_DIR}/setup.py.in")
set(SETUP_PY    "${CMAKE_CURRENT_BINARY_DIR}/setup.py")
set(OUTPUT      "${CMAKE_CURRENT_BINARY_DIR}/timestamp")

configure_file(${SETUP_PY_IN} ${SETUP_PY})

file(GLOB_RECURSE DEPS ${CMAKE_CURRENT_SOURCE_DIR}/ccpi/*.py )

if (CONDA_BUILD)
  add_custom_target(pythonsetup ALL
                    COMMAND ${CMAKE_COMMAND} -E env CIL_VERSION=${CIL_VERSION}
                    ${PYTHON_EXECUTABLE} ${SETUP_PY} -vv install
                    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                    COMMAND ${CMAKE_COMMAND} -E touch ${OUTPUT}
                    DEPENDS cilrec
                    )

else()
  if (WIN32)
    add_custom_target(pythonsetup ALL
                      COMMAND ${CMAKE_COMMAND} -E env CIL_VERSION=${CIL_VERSION}
            
                      PREFIX=${CMAKE_SOURCE_DIR}/src/
                      LIBRARY_INC=${CMAKE_SOURCE_DIR}/src/include
                      LIBRARY_LIB=${CMAKE_BINARY_DIR}/${CMAKE_BUILD_TYPE}
                      ${PYTHON_EXECUTABLE} ${SETUP_PY} build_ext
                      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                      COMMAND ${CMAKE_COMMAND} -E touch ${OUTPUT}
                      DEPENDS cilrec
                      )
  else()
    add_custom_target(pythonsetup ALL  
                      COMMAND ${CMAKE_COMMAND} -E env CIL_VERSION=${CIL_VERSION}
                      PREFIX=${CMAKE_SOURCE_DIR}/src/
                      LIBRARY_INC=${CMAKE_SOURCE_DIR}/src/include
                      LIBRARY_LIB=${CMAKE_BINARY_DIR}/
                      ${PYTHON_EXECUTABLE} ${SETUP_PY} build_ext --verbose --build-lib=${CMAKE_CURRENT_BINARY_DIR}/build/lib
                      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                      COMMAND ${CMAKE_COMMAND} -E touch ${OUTPUT}
                      DEPENDS cilrec 
                      )
  endif()

  message("install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/build/lib/ccpi DESTINATION ${PYTHON_DEST} )")
  install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/build/lib/ccpi DESTINATION ${PYTHON_DEST} )
endif()

add_custom_target(PythonWrapper ALL DEPENDS pythonsetup)