cmake_minimum_required(VERSION 2.8.10)

project(SingularityGraph_testSuite)

add_executable(UT_SingularityGraph main.cpp)

target_link_libraries(UT_SingularityGraph singularityGraphBuilder gmds gmock)

add_custom_command(
   TARGET UT_SingularityGraph
   POST_BUILD
   COMMAND ${CMAKE_COMMAND} ARGS -E copy $<TARGET_FILE:UT_SingularityGraph> ../../Testing/UT_SingularityGraph
   COMMENT "Copy SingularityGraph unit tests"
   )

add_test (NAME    UT_SingularityGraph
	  COMMAND $<TARGET_FILE:UT_SingularityGraph>
	   WORKING_DIRECTORY ${CMAKE_BINARY_DIR})

add_subdirectory(Samples)
