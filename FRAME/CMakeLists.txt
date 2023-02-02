cmake_minimum_required(VERSION 2.8.10)

project(FRAME)

#==============================================================================
#Versionnage du projet

set(FRAME_VERSION_MAJOR 1)

set(FRAME_VERSION_MINOR 0)

set(FRAME_VERSION_PATCH 0)

#==============================================================================
# DEBUG LEVEL

option(DEBUG_FRAME  "Build with extra debug messages.")

if (DEBUG_FRAME)
   add_definitions (-DDEBUG_FRAME)
endif(DEBUG_FRAME)

enable_testing()

#==============================================================================

set(DEFAULT_PLATFORM_NAME ${CMAKE_SYSTEM_NAME}_${CMAKE_SYSTEM_PROCESSOR})


set(PLATFORM_NAME ${DEFAULT_PLATFORM_NAME} CACHE STRING "Platform named used in
the installation directory name")


#==============================================================================
# LIBRARIES DEPLOYMENT LOCATION
#==============================================================================

# Headers files to be installed
set (INSTALL_FRAME_DIR
 ${CMAKE_INSTALL_PREFIX}/${FRAME_VERSION_MAJOR}.${FRAME_VERSION_MINOR}.${FRAME_VERSION_PATCH}/${PLATFORM_NAME}/${CMAKE_BUILD_TYPE}
    CACHE STRING "FRAME Install directory"   
    )

set (INSTALL_FRAME_LIB ${INSTALL_FRAME_DIR}/lib)
set (INSTALL_FRAME_INC ${INSTALL_FRAME_DIR}/include)

if(CMAKE_SYSTEM MATCHES Windows)
  add_definitions(-DWIN32)
endif()



if(CMAKE_SYSTEM MATCHES Linux)
  add_definitions(-std=c++0x )
endif()


#=============== CORE UTILITIES ======================

set(FRAME_CORE ON CACHE BOOL "Frame Field Core utility module")

if (FRAME_CORE)
  add_subdirectory(Core)
  include_directories(${CMAKE_CURRENT_SOURCE_DIR}/Core/inc)
endif (FRAME_CORE)

#=============== FRAME FIELD GENERATION ======================

set(FRAME_FIELD_GEN ON CACHE BOOL "Frame Field Generation (2D and 3D) module")

if (FRAME_FIELD_GEN)
  add_subdirectory(frameFieldGen)
  include_directories(${CMAKE_CURRENT_SOURCE_DIR}/frameFieldGen/inc)
endif (FRAME_FIELD_GEN)


#=============== SINGULARITY GRAPH ======================

set(SingGraphBuilder ON CACHE BOOL "Singularity graph module")
if (SingGraphBuilder)
  add_subdirectory(SingGraphBuilder)
  include_directories(${CMAKE_CURRENT_SOURCE_DIR}/SingGraphBuilder/inc)
endif (SingGraphBuilder)
 
#=============== HEX DOM MESHER ======================

set(FHEDO ON CACHE BOOL "Frame-based hex-dominant mesher module")
if (FHEDO)
  add_subdirectory(FHeDo)
  include_directories(${CMAKE_CURRENT_SOURCE_DIR}/FHeDo/inc)
endif (FHEDO)
 
#=============== MESH GENERATION VIA PATCHING ===================
set(PATCHING OFF CACHE BOOL "Patching mesher (2D quad and 3D hex-dominant) module")

if (PATCHING)
  add_subdirectory(Patching)
  include_directories(${CMAKE_CURRENT_SOURCE_DIR}/Patching/inc)
endif (PATCHING)


 
