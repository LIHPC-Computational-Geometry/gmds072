cmake_minimum_required(VERSION 2.8.10)


add_executable(UT_GIIM main.cpp)

target_link_libraries(UT_GIIM GIIM gmock gmds gts glib)

add_custom_command(
   TARGET UT_GIIM
   POST_BUILD
   COMMAND ${CMAKE_COMMAND} ARGS -E copy $<TARGET_FILE:UT_GIIM> ${CMAKE_BINARY_DIR}/Testing/UT_GIIM   
   COMMENT "Copy GIIM unit tests"
   VERBATIM
   )


add_test (NAME    UT_GIIM 
	  COMMAND  $<TARGET_FILE:UT_GIIM>
	  WORKING_DIRECTORY ${CMAKE_BINARY_DIR})

