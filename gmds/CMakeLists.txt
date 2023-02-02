#==============================================================================
# GMDS CORE MODULE
#==============================================================================
cmake_minimum_required(VERSION 3.14.0)

project(gmds CXX)

#add_definitions(-DDLLEXPORT)

#==============================================================================
# CODE STRUCTURE
#==============================================================================

include_directories("${CMAKE_CURRENT_SOURCE_DIR}/inc")

set(GMDS_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src/IG  ${CMAKE_CURRENT_SOURCE_DIR}/src/Utils ${CMAKE_CURRENT_SOURCE_DIR}/src/Math ${CMAKE_CURRENT_SOURCE_DIR}/src/Algo ${CMAKE_CURRENT_SOURCE_DIR}/src/CAD ${CMAKE_CURRENT_SOURCE_DIR}/src/IO)

set(GMDS_HEADER_FILES_Algo
	./inc/GMDS/Algo/BoundaryOperator.h
	./inc/GMDS/Algo/DistanceFieldBuilder2D.h
	./inc/GMDS/Algo/DistanceFieldBuilder3D.h
	./inc/GMDS/Algo/GETMe.h
	./inc/GMDS/Algo/Pillowing.h
	./inc/GMDS/Algo/SheetOperator.h
	./inc/GMDS/Algo/SubMeshExtractor.h
)
set(GMDS_HEADER_FILES_CAD
	./inc/GMDS/CAD/FacetedCurve.h
	./inc/GMDS/CAD/FacetedGeomManager.h
	./inc/GMDS/CAD/FacetedPoint.h
	./inc/GMDS/CAD/FacetedSurface.h
	./inc/GMDS/CAD/FacetedTriangulationService.h
	./inc/GMDS/CAD/FacetedVolume.h
	./inc/GMDS/CAD/FACFacetedGeomReadAndWrite.h
	./inc/GMDS/CAD/GeomCurve.h
	./inc/GMDS/CAD/GeomEntity.h
	./inc/GMDS/CAD/GeomManager.h
	./inc/GMDS/CAD/GeomPoint.h
	./inc/GMDS/CAD/GeomServiceAbstractFactory.h
	./inc/GMDS/CAD/GeomService.h
	./inc/GMDS/CAD/GeomSurface.h
	./inc/GMDS/CAD/GeomTriangulationService.h
	./inc/GMDS/CAD/GeomVolume.h
)
set(GMDS_HEADER_FILES_IG
	./inc/GMDS/IG/CellGroup_def.h
	./inc/GMDS/IG/CellGroup.h
	./inc/GMDS/IG/Cell.h
	./inc/GMDS/IG/EdgeContainer.h
	./inc/GMDS/IG/Edge.h
	./inc/GMDS/IG/FaceContainer.h
	./inc/GMDS/IG/Face.h
	./inc/GMDS/IG/IG.h
	./inc/GMDS/IG/IGMeshDoctor.h
	./inc/GMDS/IG/IGMesh.h
	./inc/GMDS/IG/IGMeshQualityEvaluation.h
	./inc/GMDS/IG/NodeContainer.h
	./inc/GMDS/IG/Node.h
	./inc/GMDS/IG/RegionContainer.h
	./inc/GMDS/IG/Region.h
)
set(GMDS_HEADER_FILES_IO
	./inc/GMDS/IO/CubitFacReader_def.h
	./inc/GMDS/IO/CubitFacReader.h
	./inc/GMDS/IO/IReader_def.h
	./inc/GMDS/IO/IReader.h
	./inc/GMDS/IO/MatrixMarketWriter_def.h
	./inc/GMDS/IO/MatrixMarketWriter.h
	./inc/GMDS/IO/MeditReader_def.h
	./inc/GMDS/IO/MeditReader.h
	./inc/GMDS/IO/VTKFacetedGeomReadAndWrite.h
	./inc/GMDS/IO/VTKReader_def.h
	./inc/GMDS/IO/VTKReader.h
	./inc/GMDS/IO/VTKWriter_def.h
	./inc/GMDS/IO/VTKWriter.h
)
set(GMDS_HEADER_FILES_Math
	./inc/GMDS/Math/AxisAngleRotation.h
	./inc/GMDS/Math/BezierCurve.h
	./inc/GMDS/Math/Chart.h
	./inc/GMDS/Math/Constants.h
	./inc/GMDS/Math/Cross2D.h
	./inc/GMDS/Math/Cross.h
	./inc/GMDS/Math/Hexahedron.h
	./inc/GMDS/Math/Line.h
	./inc/GMDS/Math/Matrix_def.h
	./inc/GMDS/Math/Matrix.h
	./inc/GMDS/Math/Numerics.h
	./inc/GMDS/Math/P1Elements.h
	./inc/GMDS/Math/Plane.h
	./inc/GMDS/Math/Point.h
	./inc/GMDS/Math/Prism3.h
	./inc/GMDS/Math/Pyramid.h
	./inc/GMDS/Math/Quadrilateral.h
	./inc/GMDS/Math/Quaternion.h
	./inc/GMDS/Math/Ray.h
	./inc/GMDS/Math/Segment.h
	./inc/GMDS/Math/SHarmonicL4.h
	./inc/GMDS/Math/SHFrame.h
	./inc/GMDS/Math/Tetrahedron.h
	./inc/GMDS/Math/Triangle.h
	./inc/GMDS/Math/VectorDyn.h
	./inc/GMDS/Math/Vector.h
	./inc/GMDS/Math/VectorND.h
)
set(GMDS_HEADER_FILES_Parallel
	./inc/GMDS/Parallel/DistributedCellData.h
	./inc/GMDS/Parallel/DistributedMemoryManager.h
	./inc/GMDS/Parallel/Message.h
	./inc/GMDS/Parallel/MigrationTool.h
)
set(GMDS_HEADER_FILES_Utils
	./inc/GMDS/Utils/Assert.h
	./inc/GMDS/Utils/ChunkMemoryAllocator.h
	./inc/GMDS/Utils/CommonFlags.h
	./inc/GMDS/Utils/CommonTypes.h
	./inc/GMDS/Utils/Exception.h
	./inc/GMDS/Utils/IndexedVector_def.h
	./inc/GMDS/Utils/IndexedVector.h
	./inc/GMDS/Utils/Log.h
	./inc/GMDS/Utils/LogStream.h
	./inc/GMDS/Utils/Marks32.h
	./inc/GMDS/Utils/OrientedGraph.h
	./inc/GMDS/Utils/Parameters.h
	./inc/GMDS/Utils/PlatformUtils.h
	./inc/GMDS/Utils/RandomGenerator.h
	./inc/GMDS/Utils/SmartBitVector.h
	./inc/GMDS/Utils/SmartVector_def.h
	./inc/GMDS/Utils/SmartVector.h
	./inc/GMDS/Utils/Timer.h
	./inc/GMDS/Utils/Variable_def.h
	./inc/GMDS/Utils/Variable.h
	./inc/GMDS/Utils/VariableManager.h
)

#set(GMDS_HEADER_FILES
list(APPEND GMDS_HEADER_FILES
		${GMDS_HEADER_FILES_Algo}
		${GMDS_HEADER_FILES_CAD}
		${GMDS_HEADER_FILES_IG}
		${GMDS_HEADER_FILES_IO}
		${GMDS_HEADER_FILES_Math}
		${GMDS_HEADER_FILES_Parallel}
		${GMDS_HEADER_FILELS_Utils}
)

set(GMDS_SRC_FILES
	./src/Algo/BoundaryOperator.cpp
	./src/Algo/DistanceFieldBuilder2D.cpp
	./src/Algo/DistanceFieldBuilder3D.cpp
	./src/Algo/GETMe.cpp
	./src/Algo/Pillowing.cpp
	./src/Algo/SheetOperator.cpp
	./src/Algo/SubMeshExtractor.cpp
	./src/CAD/FacetedCurve.cpp
	./src/CAD/FacetedGeomManager.cpp
	./src/CAD/FacetedPoint.cpp
	./src/CAD/FacetedSurface.cpp
	./src/CAD/FacetedTriangulationService.cpp
	./src/CAD/FacetedVolume.cpp
	./src/CAD/FACFacetedGeomReadAndWrite.cpp
	./src/CAD/GeomCurve.cpp
	./src/CAD/GeomPoint.cpp
	./src/CAD/GeomSurface.cpp
	./src/IG/Cell.cpp
	./src/IG/CellGroupImpl.cpp
	./src/IG/EdgeContainer.cpp
	./src/IG/Edge.cpp
	./src/IG/FaceContainer.cpp
	./src/IG/Face.cpp
	./src/IG/IGMesh.cpp
	./src/IG/IGMeshDoctor.cpp
	./src/IG/IGMeshQualityEvaluation.cpp
	./src/IG/NodeContainer.cpp
	./src/IG/Node.cpp
	./src/IG/RegionContainer.cpp
	./src/IG/Region.cpp
	./src/IO/IGMeshWriterImpl.cpp
	./src/IO/VTKFacetedGeomReadAndWrite.cpp
	./src/Math/AxisAngleRotation.cpp
	./src/Math/BezierCurve.cpp
	./src/Math/Chart.cpp
	./src/Math/Constants.cpp
	./src/Math/Cross2D.cpp
	./src/Math/Cross.cpp
	./src/Math/Hexahedron.cpp
	./src/Math/Line.cpp
	./src/Math/Numerics.cpp
	./src/Math/P1Elements.cpp
	./src/Math/Plane.cpp
	./src/Math/Point.cpp
	./src/Math/Prism3.cpp
	./src/Math/Pyramid.cpp
	./src/Math/Quadrilateral.cpp
	./src/Math/Quaternion.cpp
	./src/Math/Ray.cpp
	./src/Math/Segment.cpp
	./src/Math/SHarmonicL4.cpp
	./src/Math/Tetrahedron.cpp
	./src/Math/Triangle.cpp
	./src/Math/Vector.cpp
	./src/Math/VectorDyn.cpp
	./src/Parallel/DistributedCellData.cpp
	./src/Parallel/DistributedMemoryManager.cpp
	./src/Parallel/Message.cpp
	./src/Parallel/MigrationTool.cpp
	./src/Utils/Assert.cpp
	./src/Utils/CommonTypes.cpp
	./src/Utils/Log.cpp
	./src/Utils/OrientedGraph.cpp
	#./src/Utils/Parameters.cpp
	./src/Utils/RandomGenerator.cpp
	./src/Utils/SmartBitVector.cpp
	./src/Utils/SmartVector.cpp
	./src/Utils/Timer.cpp
	./src/Utils/Variable.cpp
	./src/Utils/VariableManager.cpp
)

#==============================================================================
# TEST MODULE
#==============================================================================

#set(GMDS_TEST ON CACHE BOOL "GMDS Kernel test module")
#if (GMDS_TEST)
#  add_subdirectory(test)
#endif ()

#==============================================================================
# LIBRARY BUILDING
#==============================================================================

add_library(gmds SHARED ${GMDS_HEADER_FILES} ${GMDS_SRC_FILES})

#set_target_properties(gmds PROPERTIES PUBLIC_HEADER "${GMDS_HEADER_FILES}")


target_compile_features(gmds PUBLIC cxx_std_11)


target_link_libraries(gmds)

#add_custom_command(
#   TARGET gmds
#   POST_BUILD
#   COMMAND ${CMAKE_COMMAND} ARGS -E copy_directory  $<TARGET_FILE:gmds> ${CMAKE_BINARY_DIR}/gen/
#   COMMENT "Copy gmds lib"
#   VERBATIM
#   )

#target_include_directories(gmds INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>/generated>)
target_include_directories(gmds PUBLIC "$<INSTALL_INTERFACE:include>")

install(TARGETS gmds
		EXPORT gmds DESTINATION
		LIBRARY DESTINATION lib
		PUBLIC_HEADER DESTINATION include
		)
#		PUBLIC_HEADER DESTINATION include)

install(FILES
		${GMDS_HEADER_FILES_Algo} DESTINATION include/GMDS/Algo)
install(FILES
		${GMDS_HEADER_FILES_CAD} DESTINATION include/GMDS/CAD)
install(FILES
		${GMDS_HEADER_FILES_IG} DESTINATION include/GMDS/IG)
install(FILES
		${GMDS_HEADER_FILES_IO} DESTINATION include/GMDS/IO)
install(FILES
		${GMDS_HEADER_FILES_Math} DESTINATION include/GMDS/Math)
install(FILES
		${GMDS_HEADER_FILES_Parallel} DESTINATION include/GMDS/Parallel)
install(FILES
		${GMDS_HEADER_FILES_Utils} DESTINATION include/GMDS/Utils)

#install(FILES ${GMDS_HEADER_FILES} DESTINATION include)


set(INCLUDE_INSTALL_DIR include/ )
set(LIB_INSTALL_DIR lib/ )
set(SYSCONFIG_INSTALL_DIR etc/gmds/ )

include(CMakePackageConfigHelpers)
configure_package_config_file(
		gmdsConfig.cmake.in ${CMAKE_CURRENT_BINARY_DIR}/gmdsConfig.cmake
		INSTALL_DESTINATION ${LIB_INSTALL_DIR}/gmds/cmake
		PATH_VARS INCLUDE_INSTALL_DIR SYSCONFIG_INSTALL_DIR
		)

write_basic_package_version_file(
		${CMAKE_CURRENT_BINARY_DIR}/gmdsConfigVersion.cmake
		VERSION 0.7.2
		COMPATIBILITY SameMajorVersion
)

install(FILES
		${CMAKE_CURRENT_BINARY_DIR}/gmdsConfig.cmake
		${CMAKE_CURRENT_BINARY_DIR}/gmdsConfigVersion.cmake
		DESTINATION ${LIB_INSTALL_DIR}/gmds/cmake
		)


install(EXPORT gmds DESTINATION lib/gmds)
