cmake_minimum_required(VERSION 2.8.10)

project(MesquiteInterface)

#==============================================================================
#Versionnage du projet

set(FRAME_VERSION_MAJOR 1)

set(FRAME_VERSION_MINOR 0)

set(FRAME_VERSION_PATCH 0)

#==============================================================================
# DEBUG LEVEL

option(DEBUG_MESQUITEINTERFACE  "Build with extra debug messages.")

if (DEBUG_MESQUITEINTERFACE)
   add_definitions (-DDEBUG_MESQUITEINTERFACE)
endif(DEBUG_MESQUITEINTERFACE)

enable_testing()

#==============================================================================
# CODE STRUCTURE
#==============================================================================

include_directories("${CMAKE_CURRENT_SOURCE_DIR}/inc")

set(MESQUITEINTERFACE_SRC_DIR  src)
file(GLOB_RECURSE SRC_FILES src/*.cpp)

include_directories(/homePOYOP/travail/utils/mesquite-2.3.0_build/include
/usr/include/glib-2.0
/usr/lib64/glib-2.0/include
/homePOYOP/travail/utils/gts-0.7.6_build/include)
link_directories(/homePOYOP/travail/utils/mesquite-2.3.0_build/lib64
/homePOYOP/travail/utils/gts-0.7.6_build/lib
/homePOYOP/travail/utils/gts-0.7.6_build/lib64)


#==============================================================================

set(DEFAULT_PLATFORM_NAME ${CMAKE_SYSTEM_NAME}_${CMAKE_SYSTEM_PROCESSOR})


set(PLATFORM_NAME ${DEFAULT_PLATFORM_NAME} CACHE STRING "Platform named used in
the installation directory name")

#==============================================================================
# TEST MODULE
#==============================================================================

set(MESQUITEINTERFACE_TEST OFF CACHE BOOL "MESQUITEINTERFACE Algo test module")

if (MESQUITEINTERFACE_TEST)
  add_subdirectory(test)
endif ()

#==============================================================================
# LIBRARIES DEPLOYMENT LOCATION
#==============================================================================

# Headers files to be installed
set (INSTALL_MESQUITEINTERFACE_DIR
 ${CMAKE_INSTALL_PREFIX}/${MESQUITEINTERFACE_VERSION_MAJOR}.${MESQUITEINTERFACE_VERSION_MINOR}.${MESQUITEINTERFACE_VERSION_PATCH}/${PLATFORM_NAME}/${CMAKE_BUILD_TYPE}
    CACHE STRING "MESQUITEINTERFACE Install directory"   
    )

set (INSTALL_MESQUITEINTERFACE_LIB ${INSTALL_MESQUITEINTERFACE_DIR}/lib)
set (INSTALL_MESQUITEINTERFACE_INC ${INSTALL_MESQUITEINTERFACE_DIR}/include)

if(CMAKE_SYSTEM MATCHES Windows)
  add_definitions(-DWIN32)
endif()



if(CMAKE_SYSTEM MATCHES Linux)
  add_definitions(-std=c++0x )
endif()

#==============================================================================
# LIBRARY BUILDING
#==============================================================================


add_library(MESQUITEINTERFACE SHARED ${SRC_FILES})


target_link_libraries(MESQUITEINTERFACE gmds GIIM gts)


add_custom_command(
   TARGET MESQUITEINTERFACE 
   POST_BUILD
   COMMAND ${CMAKE_COMMAND} ARGS -E copy $<TARGET_FILE:MESQUITEINTERFACE> ${CMAKE_BINARY_DIR}/gen/
   COMMENT "Copy MESQUITEINTERFACE executable"
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
