cmake_minimum_required(VERSION 2.8.10)


add_executable(UT_gmds main.cpp)

target_link_libraries(UT_gmds gmock gmds)

add_custom_command(
   TARGET UT_gmds
   POST_BUILD
   COMMAND ${CMAKE_COMMAND} ARGS -E copy $<TARGET_FILE:UT_gmds> ${CMAKE_BINARY_DIR}/Testing/UT_gmds   
   COMMENT "Copy gmds unit tests"
   VERBATIM
   )


add_test (NAME    UT_gmds 
	  COMMAND  $<TARGET_FILE:UT_gmds>
	  WORKING_DIRECTORY ${CMAKE_BINARY_DIR})


add_subdirectory(Samples)
