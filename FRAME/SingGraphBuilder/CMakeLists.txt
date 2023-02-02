#==============================================================================
# GMDS SINGULARITY GRAPH BUILDER COMPONENT
#==============================================================================
cmake_minimum_required(VERSION 2.8.10)

project(SingGraphBuilder)

#==============================================================================
# DEBUG LEVEL
#==============================================================================
option(DEBUG_SGB "Build with extra debug messages.")
option(DEBUG_SGB "Build with extra vtk ouptut file when required.")

if (DEBUG_SGB)
  add_definitions (-DDEBUG_SGB)
endif(DEBUG_SGB)

if (DEBUG_SGB2)
  add_definitions (-DDEBUG_SGB2)
endif (DEBUG_SGB2)

#==============================================================================
# CODE STRUCTURE
#==============================================================================
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/inc")

set(FRAME_FIELD_GEN_SRC_DIR  src)
file(GLOB_RECURSE SRC_FILES inc/*.h src/*.cpp)

#==============================================================================
# TEST MODULE 
#==============================================================================
set(SG_BUILDER_TEST OFF CACHE BOOL "SGB Algo test module")
if (SG_BUILDER_TEST)
  add_subdirectory(test)
endif ()

#==============================================================================
# LIBRARY BUILDING
#==============================================================================

add_executable(singularityGraphBuilderExe ${SRC_FILES})
add_library   (singularityGraphBuilder    ${SRC_FILES})

target_link_libraries(singularityGraphBuilderExe gmds frame_core)
target_link_libraries(singularityGraphBuilder    gmds frame_core)



add_custom_command(
   TARGET singularityGraphBuilder
   POST_BUILD
   COMMAND ${CMAKE_COMMAND} ARGS -E copy $<TARGET_FILE:singularityGraphBuilder> ${CMAKE_BINARY_DIR}/gen/
   COMMENT "Copy singularityGraphBuilder executable"
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

install(TARGETS singularityGraphBuilderExe
	RUNTIME
	DESTINATION  ${INSTALL_GMDS_SUITE_LIB}
	PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ
	GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
	) 
	#replace LIBRARY by RUNTIME on Win

install(TARGETS singularityGraphBuilderExe
	DESTINATION  ${INSTALL_GMDS_SUITE_LIB}
	EXPORT singularityGraphBuilder-targets
	) 

install(EXPORT singularityGraphBuilder-targets DESTINATION  ${INSTALL_GMDS_SUITE_DIR})
