cmake_minimum_required(VERSION 2.8.10)

project(GeomInsertion)

#==============================================================================
#Versionnage du projet

set(FRAME_VERSION_MAJOR 1)

set(FRAME_VERSION_MINOR 0)

set(FRAME_VERSION_PATCH 0)

#==============================================================================
# EXTERNAL DEPENDENCIES

# GLPK LIBRARY
set(GEOMINSERTION_GTS_LIB_DIR "" CACHE PATH "path to gts library")
set(GEOMINSERTION_GTS_INCLUDE_DIR "" CACHE PATH "path to gts includes")

#==============================================================================
# DEBUG LEVEL

option(DEBUG_GIIM  "Build with extra debug messages.")

if (DEBUG_GIIM)
   add_definitions (-DDEBUG_GIIM)
endif(DEBUG_GIIM)

enable_testing()

#==============================================================================
# CODE STRUCTURE
#==============================================================================

include_directories("${CMAKE_CURRENT_SOURCE_DIR}/inc")

#file(GLOB_RECURSE GIIM_SRC_FILES src/*.cpp)
set(GIIM_SRC_FILES
	inc/GMDS/DummyClass.h
	inc/GMDS/FacetedMeshIntersectionService.h
	inc/GMDS/GeomMeshIntersectionService.h
	inc/GMDS/LaplacianSmoothing.h
	inc/GMDS/MeshGraph.h
	inc/GMDS/MeshInsertDetail.h
	inc/GMDS/MeshInsertDetailInOut.h
	inc/GMDS/MeshModelAlgo.h
	inc/GMDS/MeshSplitter3Refinement.h
	inc/GMDS/MeshSplitter.h
	inc/GMDS/NormalLaplacianSmoothing.h
	inc/GMDS/OrderedSmartLaplacianSmoothing.h
	inc/GMDS/SmartDampenedLaplacianSmoothing.h
	inc/GMDS/SmartLaplacianSmoothing.h
	inc/GMDS/Smoothing.h
	src/DummyClass.cpp
        src/FacetedMeshIntersectionService.cpp
        src/GeomMeshIntersectionService.cpp
        src/LaplacianSmoothing.cpp
        src/MeshGraph.cpp
        src/MeshInsertDetail.cpp
        src/MeshInsertDetailInOut.cpp
        src/MeshModelAlgo.cpp
        src/MeshSplitter3Refinement.cpp
        src/MeshSplitter.cpp
	src/NormalLaplacianSmoothing.cpp
        src/OrderedSmartLaplacianSmoothing.cpp
	src/SmartDampenedLaplacianSmoothing.cpp
        src/SmartLaplacianSmoothing.cpp
)

include_directories(
	/usr/include/glib-2.0
	/usr/lib64/glib-2.0/include
	${GEOMINSERTION_GTS_INCLUDE_DIR}
)
link_directories(${GEOMINSERTION_GTS_LIB_DIR} /usr/lib64)


#==============================================================================

set(DEFAULT_PLATFORM_NAME ${CMAKE_SYSTEM_NAME}_${CMAKE_SYSTEM_PROCESSOR})


set(PLATFORM_NAME ${DEFAULT_PLATFORM_NAME} CACHE STRING "Platform named used in
the installation directory name")

#==============================================================================
# TEST MODULE
#==============================================================================

set(GIIM_TEST OFF CACHE BOOL "GIIM Algo test module")

if (GIIM_TEST)
  add_subdirectory(test)
endif ()

#==============================================================================
# LIBRARIES DEPLOYMENT LOCATION
#==============================================================================

# Headers files to be installed
set (INSTALL_GIIM_DIR
 ${CMAKE_INSTALL_PREFIX}/${GIIM_VERSION_MAJOR}.${GIIM_VERSION_MINOR}.${GIIM_VERSION_PATCH}/${PLATFORM_NAME}/${CMAKE_BUILD_TYPE}
    CACHE STRING "GIIM Install directory"   
    )

set (INSTALL_GIIM_LIB ${INSTALL_GIIM_DIR}/lib)
set (INSTALL_GIIM_INC ${INSTALL_GIIM_DIR}/include)

if(CMAKE_SYSTEM MATCHES Windows)
  add_definitions(-DWIN32)
endif()



if(CMAKE_SYSTEM MATCHES Linux)
  add_definitions(-std=c++0x )
endif()

#==============================================================================
# LIBRARY BUILDING
#==============================================================================


add_library(GIIM SHARED ${GIIM_SRC_FILES})


target_link_libraries(GIIM gmds)


add_custom_command(
   TARGET GIIM 
   POST_BUILD
   COMMAND ${CMAKE_COMMAND} ARGS -E copy $<TARGET_FILE:GIIM> ${CMAKE_BINARY_DIR}/gen/
   COMMENT "Copy GIIM executable"
   VERBATIM
   )

#==============================================================================
# LIBRARY INSTALLATION
#==============================================================================

install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/inc/"
	DESTINATION  ${INSTALL_GMDS_SUITE_INC}
	DIRECTORY_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ
	GROUP_EXECUTE WORLD_READ WORLD_EXECUTE SETGID
	FILE_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ
	GROUP_EXECUTE WORLD_READ WORLD_EXECUTE SETGID
	)
