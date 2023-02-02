cmake_minimum_required(VERSION 2.8.11)

project(gmdsVtk_testSuite)

add_executable(UT_gmdsVTK main.cpp)

target_link_libraries(UT_gmdsVTK gmock gmds gmdsVtk vtkIO vtkCommon vtksys vtkFiltering vtktiff imf Xt pthread)

add_custom_command(
   TARGET UT_gmdsVTK
   POST_BUILD
   COMMAND ${CMAKE_COMMAND} ARGS -E copy $<TARGET_FILE:UT_gmdsVTK> ../../Testing/UT_gmdsVTK
   COMMENT "Copy gmdsVTK unit tests"
   )

add_test( NAME    UT_gmdsVTK 
	  COMMAND  $<TARGET_FILE:UT_gmdsVTK>
	  WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
 
