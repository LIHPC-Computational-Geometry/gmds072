/*----------------------------------------------------------------------------*/
#include <string>
#include <iostream>
#include <fstream>
/*----------------------------------------------------------------------------*/
#include "SingularityGraphBuilder.h"
#include "SingularityGraph.h"
/*----------------------------------------------------------------------------*/
#include <GMDS/IG/IGMesh.h>
#include <GMDS/IG/IGMeshDoctor.h>
#include <GMDS/IO/MeditReader.h>
#include <GMDS/IO/VTKReader.h>
#include <GMDS/IO/VTKWriter.h>
#include <GMDS/Algo/BoundaryOperator.h>
#include <GMDS/Algo/DistanceFieldBuilder3D.h>
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/
int main(int argc, char* argv[])
{
	std::cout << "==== Singularity Graph Builder ====" << std::endl;
	MeshModel model(DIM3 | R | F | E | N | R2N | F2N | E2N | R2F | F2R | F2E | E2F | N2R | N2F | N2E);
	IGMesh mesh(model);

	//==================================================================
	// MESH READING
	//==================================================================
	std::cout << "Reading " << std::endl;
	try
	{
		std::string fIn, fOut;
		fOut = fIn;

		if (argc < 2)
            throw gmds::GMDSException("Wrong parameters");

		if (argc >= 2)
		{
			fIn = std::string(argv[1]);
			std::cout << "INPUT FILE: " << fIn << std::endl;
			fOut = fIn;
		}
		if (argc == 3){
			fOut = std::string(argv[2]);
			std::cout << "OUTPUT DIR: " << fOut << std::endl;
		}

		std::cout << "Start reading file " << fIn << std::endl;
		VTKReader<IGMesh> reader(mesh);
		reader.read(fIn);
		std::cout << "    DONE" << std::endl;
		std::cout << "IN - (R,F,E,N) = (" << mesh.getNbRegions()
			<< ", " << mesh.getNbFaces()
			<< ", " << mesh.getNbEdges()
			<< ", " << mesh.getNbNodes() << ")" << std::endl;

		//==================================================================
		// MESH TOPOLOGY PREPARATION
		//==================================================================
		std::cout << "Start mesh correction" << std::endl;
		IGMeshDoctor doc(&mesh);
		std::cout << "Build Faces" << std::endl;
		doc.buildFacesAndR2F();
		std::cout << "    DONE" << std::endl;
		std::cout << "Build Edges" << std::endl;
		doc.buildEdgesAndX2E();
		std::cout << "    DONE" << std::endl;
		std::cout << "Update mesh connectivity" << std::endl;
		doc.updateUpwardConnectivity();
		std::cout << "    DONE" << std::endl;
		std::cout << "IN - (R,F,E,N) = (" << mesh.getNbRegions()
			<< ", " << mesh.getNbFaces()
			<< ", " << mesh.getNbEdges()
			<< ", " << mesh.getNbNodes() << ")" << std::endl;

		//==================================================================
		// SINGULARITY GRAPH EXTRACTION
		//==================================================================
		std::cout << "Start singularity graph extraction" << std::endl;
		SingularityGraphBuilder builder(&mesh);
		builder.setDebugDirectory(fOut);
		builder.execute();
	}
	catch (exception & e)
	{
		std::cout << e.what();
		std::cout << std::endl;
	}
	return 0;
}
/*----------------------------------------------------------------------------*/
