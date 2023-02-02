/*----------------------------------------------------------------------------*/
/*
 * VTKFacetedGeomReader.cpp
 *
 *  Created on: 29 août 2014
 *      Author: ledouxf
 */
/*----------------------------------------------------------------------------*/
#include "GMDS/IOVTK/VTKFacetedGeomReadAndWrite.h"
#include "GMDS/IOVTK/VTKWriter.h"
#include "GMDS/IOVTK/VTKReader.h"
/*----------------------------------------------------------------------------*/
//#include <>
/*----------------------------------------------------------------------------*/
namespace gmds {
/*----------------------------------------------------------------------------*/
VTKFacetedGeomReadAndWrite::VTKFacetedGeomReadAndWrite()
{}
/*----------------------------------------------------------------------------*/
VTKFacetedGeomReadAndWrite::~VTKFacetedGeomReadAndWrite()
{}
/*----------------------------------------------------------------------------*/
void VTKFacetedGeomReadAndWrite::import(
			 geom::FacetedGeomManager& AGeomMng,
			 const std::string& 	   AFile)
{
	// we build a mesh from the file
	VTKReader<IGMesh> r(AGeomMng.getMeshView());
	r.read(AFile);
	AGeomMng.updateFromMesh();

}
/*----------------------------------------------------------------------------*/
void VTKFacetedGeomReadAndWrite::
exportVTK( geom::FacetedGeomManager& AGeomMng, const std::string& AFile)
{
	VTKWriter<IGMesh> w(AGeomMng.getMeshView());
	w.write(AFile,N|F);
}

/*----------------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------------*/
