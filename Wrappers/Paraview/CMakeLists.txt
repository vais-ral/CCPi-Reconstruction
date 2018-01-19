cmake_minimum_required(VERSION 2.8)
project(CCPiReconstructionParaview)

if (NOT ParaView_BINARY_DIR)
  find_package(ParaView REQUIRED)
  include(${PARAVIEW_USE_FILE})
endif()

include(ParaViewPlugins)

# create a paraview plugin containing server manager xml and the server
# manager classes to build
# this plugin can be loaded on the server side

ADD_PARAVIEW_PLUGIN(CCPiReconstruction(CCPi) "1.0"
  SERVER_MANAGER_XML CCPiReconstruction.xml 
  SERVER_MANAGER_SOURCES CCPiReconstructionParaviewImpl.cxx
  SERVER_SOURCES CCPiParaviewUserInterface.cpp
)

#Boost Libraries
SET(Boost_USE_STATIC_LIBS        OFF) # only find dynamic libs
SET(Boost_USE_MULTITHREADED  ON)
find_package(Boost)
if(Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIRS})
  target_link_libraries(CCPiReconstruction ${Boost_LIBRARIES})
endif()

INCLUDE_DIRECTORIES(${CCPiReconstructionIterative_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS})

TARGET_LINK_LIBRARIES(CCPiReconstruction CCPiReconstructionIterative ${Boost_LIBRARIES})
