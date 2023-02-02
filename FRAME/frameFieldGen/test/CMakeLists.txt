cmake_minimum_required(VERSION 2.8.10)

project(frameFieldGen_testSuite)

add_executable(UT_FrameField main.cpp)

target_link_libraries(UT_FrameField frameFieldGen gmds gmock)

add_custom_command(
   TARGET UT_FrameField
   POST_BUILD
   COMMAND ${CMAKE_COMMAND} ARGS -E copy $<TARGET_FILE:UT_FrameField> ../../Testing/UT_FrameField
   COMMENT "Copy FrameField unit tests"
   )

add_test (NAME    UT_FrameField
	  COMMAND $<TARGET_FILE:UT_FrameField>
	   WORKING_DIRECTORY ${CMAKE_BINARY_DIR})

add_subdirectory(Samples)
