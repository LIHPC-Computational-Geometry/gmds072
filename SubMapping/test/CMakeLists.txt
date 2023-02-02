cmake_minimum_required(VERSION 2.8.10)


add_executable(UT_SUBMAPPING main.cpp)

target_link_libraries(UT_SUBMAPPING SUBMAPPING gmock gmds glpk)

add_custom_command(
   TARGET UT_SUBMAPPING
   POST_BUILD
   COMMAND ${CMAKE_COMMAND} ARGS -E copy $<TARGET_FILE:UT_SUBMAPPING>
${CMAKE_BINARY_DIR}/Testing/UT_SUBMAPPING
   COMMENT "Copy SUBMAPPING unit tests"
   VERBATIM
   )


add_test (NAME    UT_SUBMAPPING 
	  COMMAND  $<TARGET_FILE:UT_SUBMAPPING>
	  WORKING_DIRECTORY ${CMAKE_BINARY_DIR})

