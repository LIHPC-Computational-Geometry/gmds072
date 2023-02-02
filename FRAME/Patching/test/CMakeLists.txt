cmake_minimum_required(VERSION 2.8.10)

project(Patching_testSuite)

#==============================================================================
# CODE STRUCTURE
#==============================================================================
file(GLOB_RECURSE SRC_FILES *.h *.cpp)

#==============================================================================
# LIBRARY BUILDING
#==============================================================================

add_executable(UT_Patching main.cpp)

target_link_libraries(UT_Patching Patching gmds gmock)

add_custom_command(
   TARGET UT_Patching
   POST_BUILD
   COMMAND ${CMAKE_COMMAND} ARGS -E copy $<TARGET_FILE:UT_Patching> ../../Testing/UT_Patching
   COMMENT "Copy Patching unit tests"
   )

add_test (NAME    UT_Patching
	  COMMAND $<TARGET_FILE:UT_Patching>
	   WORKING_DIRECTORY ${CMAKE_BINARY_DIR})

add_subdirectory(Samples)
