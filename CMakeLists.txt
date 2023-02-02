cmake_minimum_required(VERSION 3.14.0)

project(GMDSSuite)

#==============================================================================
#Versionnage du projet

set(GMDS_SUITE_VERSION_MAJOR 0)
set(GMDS_SUITE_VERSION_MINOR 7)
set(GMDS_SUITE_VERSION_PATCH 2)

set(GMDS_VERSION
        ${GMDS_SUITE_VERSION_MAJOR}.${GMDS_SUITE_VERSION_MINOR}.${GMDS_SUITE_VERSION_PATCH})

#==============================================================================

# DEBUG LEVEL

#option(DEBUG_GMDS  "Build with extra debug messages.")

#option(DEBUG_GMDS2 "Build with extra vtk ouptut file when required.")

#==============================================================================
# CAMAL
#==============================================================================
#IF(NOT CAMAL_LIB_PATH)
#  FIND_PATH(CAMAL_LIB_PATH NAMES closed
#    HINTS $ENV{CAMAL_LIB_PATH}
#    DOC "Path to CAMAL libraries")
#ENDIF(NOT CAMAL_LIB_PATH)
#
#IF(NOT CAMAL_INC_PATH)
#  FIND_PATH(CAMAL_INC_PATH NAMES closed
#    HINTS $ENV{CAMAL_INC_PATH}
#    DOC "Path to CAMAL include files")
#ENDIF(NOT CAMAL_INC_PATH)
#
#
#if (DEBUG_GMDS)
#   add_definitions (-DDEBUG_GMDS)
#endif(DEBUG_GMDS)
#
#if (DEBUG_GMDS2)
#   add_definitions (-DDEBUG_GMDS2)
#endif (DEBUG_GMDS2)

#enable_testing()

#==============================================================================

set(DEFAULT_PLATFORM_NAME ${CMAKE_SYSTEM_NAME}_${CMAKE_SYSTEM_PROCESSOR})


set(PLATFORM_NAME ${DEFAULT_PLATFORM_NAME} CACHE STRING "Platform named used in
the installation directory name")


#==============================================================================
# LIBRARIES DEPLOYMENT LOCATION
#==============================================================================

# Headers files to be installed
#set (INSTALL_GMDS_SUITE_DIR
# ${CMAKE_INSTALL_PREFIX}/${GMDS_SUITE_VERSION_MAJOR}.${GMDS_SUITE_VERSION_MINOR}.${GMDS_SUITE_VERSION_PATCH}/${PLATFORM_NAME}/${CMAKE_BUILD_TYPE}
#    CACHE STRING "GMDS Suite Install directory"
#    )
#
#set (INSTALL_GMDS_SUITE_LIB ${INSTALL_GMDS_SUITE_DIR}/lib)
#set (INSTALL_GMDS_SUITE_INC ${INSTALL_GMDS_SUITE_DIR}/include)

if(CMAKE_SYSTEM MATCHES Windows)
  add_definitions(-DWIN32)
endif()

#==============================================================================
# THIRD PARTY LIBRARIES
#==============================================================================

#add_subdirectory(3rd-party)
#include_directories(${gmock_SOURCE_DIR}/include)
#include_directories(${gmock_SOURCE_DIR}/gtest/include)
#==============================================================================
#FOR EIGEN
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/3rd-party)
#==============================================================================
##FOR ALGLIB
#include_directories(${CMAKE_CURRENT_SOURCE_DIR}/3rd-party)
#link_directories(${CMAKE_BINARY_DIR}/3rd-party/alglib)
#
##==============================================================================
##FOR HLBFGS
#include_directories(${CMAKE_CURRENT_SOURCE_DIR}/3rd-party/HLBFGS)
#link_directories(${CMAKE_BINARY_DIR}/3rd-party/HLBFGS)
#
#
##==============================================================================
##FOR OpenNL
#include_directories(${CMAKE_CURRENT_SOURCE_DIR}/3rd-party/OpenNL)
#link_directories(${CMAKE_BINARY_DIR}/3rd-party/OpenNL)
#
#
##==============================================================================
##FOR inih
#include_directories(${CMAKE_CURRENT_SOURCE_DIR}/3rd-party/inih/inc)
#link_directories(${CMAKE_BINARY_DIR}/3rd-party/inih)
#
##==============================================================================
##FOR Predicates_psm
#include_directories(${CMAKE_CURRENT_SOURCE_DIR}/3rd-party/Predicates_psm)
#link_directories(${CMAKE_BINARY_DIR}/3rd-party/Predicates_psm)
#
##==============================================================================
##FOR Tetgen
#include_directories(${CMAKE_CURRENT_SOURCE_DIR}/3rd-party/tetgen)
#link_directories(${CMAKE_BINARY_DIR}/3rd-party/tetgen)
#
#
##==============================================================================
##FOR ANN
#include_directories(${CMAKE_CURRENT_SOURCE_DIR}/3rd-party/ANN_1.1.2/include)
#link_directories(${CMAKE_BINARY_DIR}/3rd-party/ANN_1.1.2/lib)

#==============================================================================
# GMDS CORE
#==============================================================================
#if(CMAKE_SYSTEM MATCHES Linux)
#  add_definitions(-std=c++0x )
#endif()
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/gen/)

add_subdirectory(gmds)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/gmds/inc)
link_directories(${CMAKE_BINARY_DIR}/gmds/Debug)

#==============================================================================
# FRAME COMPONENT
#==============================================================================
#set(FRAME OFF CACHE BOOL "FRAME Component (quad/hex meshing via Frame Fields)")
#
#if (FRAME)
#  add_subdirectory(FRAME)
#endif (FRAME)
#
#
#
#
##==============================================================================
## GIIM COMPONENT (Geometry Insertion In Mesh)
##==============================================================================
#set(GIIM_FIELD_GEN OFF CACHE BOOL "Geometry Insertion module")
#
#if (GIIM_FIELD_GEN)
#  add_subdirectory(GeomInsertion)
#  include_directories(${CMAKE_CURRENT_SOURCE_DIR}/GeomInsertion/inc)
#endif (GIIM_FIELD_GEN)
#
#
##==============================================================================
## MESQUITE COMPONENT (Mesquite smoother GMDS interface)
##==============================================================================
#set(MESQUITE_INTERFACE OFF CACHE BOOL "Mesquite Interface")
#
#if (MESQUITE_INTERFACE)
#  add_subdirectory(MesquiteInterface)
#  include_directories(${CMAKE_CURRENT_SOURCE_DIR}/MesquiteInterface/inc)
#endif (MESQUITE_INTERFACE)
#
##==============================================================================
## SUBMAPPING COMPONENT (SubMapping)
##==============================================================================
#set(SUBMAPPING OFF CACHE BOOL "SubMapping")
#
#if (SUBMAPPING)
#  add_subdirectory(SubMapping)
#  include_directories(${CMAKE_CURRENT_SOURCE_DIR}/SubMapping/inc)
#endif (SUBMAPPING)
