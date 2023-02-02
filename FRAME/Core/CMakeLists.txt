#==============================================================================
# FRAME CORE COMPONENT
#==============================================================================
cmake_minimum_required(VERSION 2.8.10)

project(frame_core)

#==============================================================================
# DEBUG LEVEL
#==============================================================================
option(DEBUG_FRAME_CORE  "Build with extra debug messages.")

if (DEBUG_FRAME_CORE)
  add_definitions (-DDEBUG_FRAME_CORE)
endif(DEBUG_FRAME_CORE)
#==============================================================================
# CODE STRUCTURE
#==============================================================================
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/inc")

set(FRAME_CORE_SRC_DIR  src)
file(GLOB_RECURSE SRC_FILES inc/*.h src/*.cpp)

#==============================================================================
# LIBRARY BUILDING
#==============================================================================

add_library   (frame_core    ${SRC_FILES})

target_link_libraries(frame_core gmds Predicates_psm)

