#==============================================================================
# GMDS FRAME FIELD GEN COMPONENT
#==============================================================================
cmake_minimum_required(VERSION 2.8.10)

project(frameFieldGen)

#==============================================================================
# CAMAL
#==============================================================================
include_directories(${CAMAL_INC_PATH})
link_directories(${CAMAL_LIB_PATH})

#==============================================================================
# DEBUG LEVEL
#==============================================================================

option(DEBUG_FFG "Build with extra debug messages.")
option(DEBUG_FFG "Build with extra vtk ouptut file when required.")

if (DEBUG_FFG)
  add_definitions (-DDEBUG_FFG)
endif(DEBUG_FFG)

if (DEBUG_FFG2)
  add_definitions (-DDEBUG_FFG2)
endif (DEBUG_FFG2)

#==============================================================================
# CODE STRUCTURE
#==============================================================================

include_directories("${CMAKE_CURRENT_SOURCE_DIR}/inc")

set(FRAME_FIELD_GEN_SRC_DIR  src)
file(GLOB_RECURSE SRC_FILES inc/*.h src/*.cpp)

#==============================================================================
# TEST MODULE
#==============================================================================

set(FRAME_FIELD_GEN_TEST OFF CACHE BOOL "FFG Algo test module")

if (FRAME_FIELD_GEN_TEST)
  add_subdirectory(test)
endif ()
 
#==============================================================================
# LIBRARY BUILDING
#==============================================================================

add_definitions(-DTETLIBRARY)

add_executable(frameFieldGenExe ${SRC_FILES})
add_library(frameFieldGen SHARED ${SRC_FILES})


target_link_libraries(frameFieldGenExe gmds inih alglib hlbfgs frame_core)
target_link_libraries(frameFieldGen inih gmds alglib hlbfgs frame_core)


add_custom_command(
   TARGET frameFieldGen 
   POST_BUILD
   COMMAND ${CMAKE_COMMAND} ARGS -E copy $<TARGET_FILE:frameFieldGen> ${CMAKE_BINARY_DIR}/gen/
   COMMENT "Copy frameFieldGen executable"
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


install(TARGETS frameFieldGenExe
	RUNTIME
	DESTINATION  ${INSTALL_GMDS_SUITE_LIB}
	PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ
	GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
	) 
	#replace LIBRARY by RUNTIME on Win

install(TARGETS frameFieldGenExe 
	DESTINATION  ${INSTALL_GMDS_SUITE_LIB}
	EXPORT frameFieldGen-targets
	) 
	
install(EXPORT frameFieldGen-targets DESTINATION  ${INSTALL_GMDS_SUITE_DIR})
