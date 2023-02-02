/*
 * IGMeshWriterImpl.cpp
 *
 *  Created on: 21 mai 2014
 *      Author: ledouxf
 */




#include <GMDS/IOVTK/VTKWriter.h>
#include <GMDS/IOVTK/VTKReader.h>
#include <GMDS/IG/IGMesh.h>

namespace gmds{

template class VTKWriter<IGMesh>;
template class VTKReader<IGMesh>;
}
