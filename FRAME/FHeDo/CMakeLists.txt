#==============================================================================
# FHeDo COMPONENT
#==============================================================================
cmake_minimum_required(VERSION 2.8.10)

project(FHeDo)

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

if (DEBUG_FHEDO)
  add_definitions (-DDEBUG_FHEDO)
endif(DEBUG_FHEDO)

#==============================================================================
# CODE STRUCTURE
#==============================================================================

include_directories("${CMAKE_CURRENT_SOURCE_DIR}/inc")

set(FHEDO_SRC_DIR  src)
file(GLOB_RECURSE SRC_FILES inc/*.h src/*.cpp)

#==============================================================================
# TEST MODULE
#==============================================================================

set(FHEDO_TEST OFF CACHE BOOL "FHEDO test module")

if (FHEDO_TEST)
  add_subdirectory(test)
endif ()
 
#==============================================================================
# LIBRARY BUILDING
#==============================================================================

add_definitions(-DTETLIBRARY)

add_executable(FHeDo ${SRC_FILES})
add_library(FHeDoLib SHARED ${SRC_FILES})

#### FRANCK, for Lin64, we had to add the dl library in the 2 link lines below.  If this messes you up on MAC, let us know.
#### - Matt Staten
target_link_libraries(FHeDo gmds dl inih alglib hlbfgs frame_core singularityGraphBuilder OpenNL Predicates_psm tet ANN CMLPave CMLShr CMLUtil)
target_link_libraries(FHeDoLib inih gmds dl alglib hlbfgs frame_core singularityGraphBuilder OpenNL Predicates_psm  tet ANN CMLPave CMLShr CMLUtil)


add_custom_command(
   TARGET FHeDoLib 
   POST_BUILD
   COMMAND ${CMAKE_COMMAND} ARGS -E copy $<TARGET_FILE:FHeDoLib> ${CMAKE_BINARY_DIR}/gen/
   COMMENT "Copy FHeDo executable"
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


install(TARGETS FHeDo
	RUNTIME
	DESTINATION  ${INSTALL_GMDS_SUITE_LIB}
	PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ
	GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
	) 
	#replace LIBRARY by RUNTIME on Win

install(TARGETS FHeDo 
	DESTINATION  ${INSTALL_GMDS_SUITE_LIB}
	EXPORT FHeDoLib-targets
	) 
	
install(EXPORT FHeDoLib-targets DESTINATION  ${INSTALL_GMDS_SUITE_DIR})
