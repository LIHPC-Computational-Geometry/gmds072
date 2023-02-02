cmake_minimum_required(VERSION 2.8.10)

project(SubMapping)

#==============================================================================
#Versionnage du projet

set(FRAME_VERSION_MAJOR 1)

set(FRAME_VERSION_MINOR 0)

set(FRAME_VERSION_PATCH 0)

#==============================================================================
# EXTERNAL DEPENDENCIES

# GLPK LIBRARY
set(SUBMAPPING_GLPK_LIB_DIR "" CACHE PATH "path to glpk library")
set(SUBMAPPING_GLPK_INCLUDE_DIR "" CACHE PATH "path to glpk includes")

#==============================================================================
# DEBUG LEVEL

option(DEBUG_SUBMAPPING  "Build with extra debug messages.")

if (DEBUG_SUBMAPPING)
   add_definitions (-DDEBUG_SUBMAPPING)
endif(DEBUG_SUBMAPPING)

enable_testing()

#==============================================================================
# CODE STRUCTURE
#==============================================================================

include_directories("${CMAKE_CURRENT_SOURCE_DIR}/inc")

set(SUBMAPPING_SRC_FILES
	inc/GMDS/SubMapping/DummyClassSubMapping.h
	inc/GMDS/SubMapping/SubMapping.h
	inc/GMDS/SubMapping/SubMappingCommon.h
	src/DummyClassSubMapping.cpp
	src/SubMapping.cpp
)

include_directories(${SUBMAPPING_GLPK_INCLUDE_DIR})
link_directories(${SUBMAPPING_GLPK_LIB_DIR})

#==============================================================================

set(DEFAULT_PLATFORM_NAME ${CMAKE_SYSTEM_NAME}_${CMAKE_SYSTEM_PROCESSOR})


set(PLATFORM_NAME ${DEFAULT_PLATFORM_NAME} CACHE STRING "Platform named used in
the installation directory name")

#==============================================================================
# TEST MODULE
#==============================================================================

set(SUBMAPPING_TEST OFF CACHE BOOL "SUBMAPPING Algo test module")

if (SUBMAPPING_TEST)
  add_subdirectory(test)
endif ()

#==============================================================================
# LIBRARIES DEPLOYMENT LOCATION
#==============================================================================

# Headers files to be installed
set (INSTALL_SUBMAPPING_DIR
 ${CMAKE_INSTALL_PREFIX}/${SUBMAPPING_VERSION_MAJOR}.${SUBMAPPING_VERSION_MINOR}.${SUBMAPPING_VERSION_PATCH}/${PLATFORM_NAME}/${CMAKE_BUILD_TYPE}
    CACHE STRING "SUBMAPPING Install directory"   
    )

set (INSTALL_SUBMAPPING_LIB ${INSTALL_SUBMAPPING_DIR}/lib)
set (INSTALL_SUBMAPPING_INC ${INSTALL_SUBMAPPING_DIR}/include)

if(CMAKE_SYSTEM MATCHES Windows)
  add_definitions(-DWIN32)
endif()



if(CMAKE_SYSTEM MATCHES Linux)
  add_definitions(-std=c++0x )
endif()

#==============================================================================
# LIBRARY BUILDING
#==============================================================================


add_library(SUBMAPPING SHARED ${SUBMAPPING_SRC_FILES})


target_link_libraries(SUBMAPPING gmds glpk)


add_custom_command(
   TARGET SUBMAPPING 
   POST_BUILD
   COMMAND ${CMAKE_COMMAND} ARGS -E copy $<TARGET_FILE:SUBMAPPING> ${CMAKE_BINARY_DIR}/gen/
   COMMENT "Copy SUBMAPPING executable"
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
