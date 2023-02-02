cmake_minimum_required(VERSION 2.8.11)

project(gmdsVtk)

#==============================================================================
# CODE STRUCTURE
#==============================================================================
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/inc")
set(GMDSVTK_SRC_DIR  src/IOVTK)
file(GLOB_RECURSE SRC_FILES src/*.cpp)

#==============================================================================
# TEST MODULE 
#==============================================================================
set(GMDS_VTK_TEST OFF CACHE BOOL "GMDS VTK Read and Write test module")
if (GMDS_VTK_TEST)
add_subdirectory(test)
endif ()

#==============================================================================
# LIBRARY BUILDING
#==============================================================================
add_library(gmdsVtk SHARED ${SRC_FILES})
target_link_libraries(gmdsVtk gmds vtkIO vtkCommon vtksys vtkFiltering vtktiff Xt pthread)


#==============================================================================
# LIBRARY INSTALLATION
#==============================================================================
install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/inc/GMDS"
	DESTINATION ${INSTALL_GMDS_SUITE_INC}
	DIRECTORY_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ
	GROUP_EXECUTE WORLD_READ WORLD_EXECUTE SETGID
	FILE_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ
	GROUP_EXECUTE WORLD_READ WORLD_EXECUTE SETGID
	)

install(TARGETS gmdsVtk
	LIBRARY 
	DESTINATION ${INSTALL_GMDS_SUITE_LIB}
	PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ
	GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
	) #replace LIBRARY by RUNTIME on Win

install(TARGETS gmdsVtk
	DESTINATION ${INSTALL_GMDS_SUITE_LIB}
	EXPORT gmdsVtk-targets
	) #replace LIBRARY by RUNTIME on Win

install(EXPORT gmdsVtk-targets DESTINATION  ${INSTALL_GMDS_SUITE_DIR})
