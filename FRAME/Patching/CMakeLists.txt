#==============================================================================
# GMDS PATCHING COMPONENT
#==============================================================================
cmake_minimum_required(VERSION 2.8.10)

project(Patching)

#==============================================================================
# DEBUG LEVEL
#==============================================================================
option(DEBUG_PATCHING  "Build with extra debug messages.")
option(DEBUG_PATCHING2 "Build with extra vtk ouptut file when required.")

if (DEBUG_PATCHING)
  add_definitions (-DDEBUG_PATCHING)
endif(DEBUG_PATCHING)

if (DEBUG_PATCHING2)
  add_definitions (-DDEBUG_PATCHING2)
endif (DEBUG_PATCHING2)

#==============================================================================
# CODE STRUCTURE
#==============================================================================
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/inc")

set(PATCHING_SRC_DIR  src)
file(GLOB_RECURSE SRC_FILES inc/*.h src/*.cpp)

#==============================================================================
# TEST MODULE 
#==============================================================================
set(PATCHING_TEST OFF CACHE BOOL "Patching test module")
if (PATCHING_TEST)
  add_subdirectory(test)
endif ()

#==============================================================================
# LIBRARY BUILDING
#==============================================================================

add_executable(PatchingExe ${SRC_FILES})
add_library   (Patching    ${SRC_FILES})

target_link_libraries(PatchingExe gmds frame_core)
target_link_libraries(Patching    gmds frame_core)



add_custom_command(
   TARGET Patching
   POST_BUILD
   COMMAND ${CMAKE_COMMAND} ARGS -E copy $<TARGET_FILE:Patching> ${CMAKE_BINARY_DIR}/gen/
   COMMENT "Copy Patching executable"
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

install(TARGETS PatchingExe
	RUNTIME
	DESTINATION  ${INSTALL_GMDS_SUITE_LIB}
	PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ
	GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
	) 
	#replace LIBRARY by RUNTIME on Win

install(TARGETS PatchingExe
	DESTINATION  ${INSTALL_GMDS_SUITE_LIB}
	EXPORT Patching-targets
	) 

install(EXPORT Patching-targets DESTINATION  ${INSTALL_GMDS_SUITE_DIR})
