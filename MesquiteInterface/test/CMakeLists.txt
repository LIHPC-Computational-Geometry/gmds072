cmake_minimum_required(VERSION 2.8.10)


add_executable(UT_MI main.cpp)

target_link_libraries(UT_MI MESQUITEINTERFACE gmock gmds mesquite)

add_custom_command(
   TARGET UT_MI
   POST_BUILD
   COMMAND ${CMAKE_COMMAND} ARGS -E copy $<TARGET_FILE:UT_MI> ${CMAKE_BINARY_DIR}/Testing/UT_MI
   COMMENT "Copy MI unit tests"
   VERBATIM
   )


add_test (NAME    UT_MI 
	  COMMAND  $<TARGET_FILE:UT_MI>
	  WORKING_DIRECTORY ${CMAKE_BINARY_DIR})

