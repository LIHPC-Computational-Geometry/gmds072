/*----------------------------------------------------------------------------*/
/* Copyright: CEA
 * contributors: F. Ledoux and N. Le Goff (2015)
 *
 * franck.ledoux@cea.fr
 * nicolas.le-goff@cea.fr
 *
 * The GMDS library is a computer program whose purpose is to provide a set of
 * functionnalities to represent and handle any type of meshes (2D, 3D,
 * triangles, tetrahedra, quad, hexa, polygons, polyhedra, etc.) and write
 * meshing algorithms. So it gathers many mathematical objects like points,
 * segment, quaternions, etc. and basic algorithms useful to build more evolved
 * ones.
 *
 * This software is governed by the CeCILL-C license under French law and
 * abiding by the rules of distribution of free software.  You can  use, 
 * modify and/ or redistribute the software under the terms of the CeCILL-C
 * license as circulated by CEA, CNRS and INRIA at the following URL
 * "http://www.cecill.info". 
 *
 * As a counterpart to the access to the source code and  rights to copy,
 * modify and redistribute granted by the license, users are provided only
 * with a limited warranty  and the software's author,  the holder of the
 * economic rights,  and the successive licensors  have only  limited
 * liability. 
 *
 * In this respect, the user's attention is drawn to the risks associated
 * with loading,  using,  modifying and/or developing or reproducing the
 * software by the user in light of its specific status of free software,
 * that may mean  that it is complicated to manipulate,  and  that  also
 * therefore means  that it is reserved for developers  and  experienced
 * professionals having in-depth computer knowledge. Users are therefore
 * encouraged to load and test the software's suitability as regards their
 * requirements in conditions enabling the security of their systems and/or 
 * data to be ensured and,  more generally, to use and operate it in the 
 * same conditions as regards security. 
 *
 * The fact that you are presently reading this means that you have had
 * knowledge of the CeCILL-C license and that you accept its terms.
 */
/*----------------------------------------------------------------------------*/
/** \file    IReader.h
 *  \author  F. LEDOUX
 *  \date    03/17/2009
 */
/*----------------------------------------------------------------------------*/
#ifndef GMDS_IREADER_H_
#define GMDS_IREADER_H_
/*----------------------------------------------------------------------------*/
// GMDS header files
#include <GMDS/Utils/CommonTypes.h>
//#include <GMDS/IG/IGMesh.h>
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
#include <GMDS/IO/IReader_def.h>
/*----------------------------------------------------------------------------*/
template<typename TMesh>
IReader<TMesh>::IReader(TMesh& AMesh)
:mesh_(AMesh)
{}
/*----------------------------------------------------------------------------*/
template<typename TMesh>
IReader<TMesh>::~IReader()
{}
/*----------------------------------------------------------------------------*/
template<typename TMesh>
void IReader<TMesh>::updateMeshIDContainers()
{
	mesh_.updateIDContainers();
}
/*----------------------------------------------------------------------------*/
template<typename TMesh>
void IReader<TMesh>::specifyMaxNodeID(const TInt& AID)
{
	mesh_.clearAndResizeNodeIDContainer(AID);
}
/*----------------------------------------------------------------------------*/
template<typename TMesh>
void IReader<TMesh>::specifyMaxEdgeID(const TInt& AID)
{
	mesh_.clearAndResizeEdgeIDContainer(AID);
}
/*----------------------------------------------------------------------------*/
template<typename TMesh>
void IReader<TMesh>::specifyMaxFaceID(const TInt& AID)
{
	mesh_.clearAndResizeFaceIDContainer(AID);
}
/*----------------------------------------------------------------------------*/
template<typename TMesh>
void IReader<TMesh>::specifyMaxRegionID(const TInt& AID)
{
	mesh_.clearAndResizeRegionIDContainer(AID);
}
/*----------------------------------------------------------------------------*/
template<typename TMesh>
void IReader<TMesh>::newNode(const TCoord& AX, const TCoord& AY,
							 const TCoord& AZ, const TCellID& AGID)
{
	mesh_.newNodeWithID(AX,AY,AZ,AGID);
}
/*----------------------------------------------------------------------------*/
template<typename TMesh>
void IReader<TMesh>::newNode(const TCoord& AX, const TCoord& AY,
							 const TCellID& AGID)
{
	mesh_.newNodeWithID(AX,AY,AGID);
}
/*----------------------------------------------------------------------------*/
template<typename TMesh>
void IReader<TMesh>::newEdge(const TCellID& AN1, const TCellID& AN2,
							 const TCellID& AGID)
{
	mesh_.newEdgeWithID(AN1,AN2,AGID);
}
/*----------------------------------------------------------------------------*/
template<typename TMesh>
void IReader<TMesh>::newTriangle(const TCellID& AN1, const TCellID& AN2,
								 const TCellID& AN3, const TCellID& AGID)
{
	mesh_.newTriangleWithID(AN1,AN2,AN3,AGID);
}
/*----------------------------------------------------------------------------*/
template<typename TMesh>
void IReader<TMesh>::newQuad(const TCellID& AN1, const TCellID& AN2,
							 const TCellID& AN3, const TCellID& AN4,
							 const TCellID& AGID)
{
	mesh_.newQuadWithID(AN1,AN2,AN3,AN4,AGID);
}
/*----------------------------------------------------------------------------*/
template<typename TMesh>
void IReader<TMesh>::newPolygon(std::vector<TCellID>& AIDs, const TCellID& AGID)
{
	mesh_.newPolygonWithID(AIDs,AGID);
}
/*----------------------------------------------------------------------------*/
template<typename TMesh>
void IReader<TMesh>::newFace(std::vector<TCellID>& AIDs, const TCellID& AGID)
{
	mesh_.newFaceWithID(AIDs,AGID);
}
/*----------------------------------------------------------------------------*/
template<typename TMesh>
void IReader<TMesh>::newTet(const TCellID& AN1, const TCellID& AN2,
							const TCellID& AN3, const TCellID& AN4,
							const TCellID& AGID)
{
	mesh_.newTetWithID(AN1,AN2,AN3,AN4,AGID);
}
/*----------------------------------------------------------------------------*/
template<typename TMesh>
void IReader<TMesh>::newPyramid(const TCellID& AN1, const TCellID& AN2,
								const TCellID& AN3, const TCellID& AN4,
								const TCellID& AN5, const TCellID& AGID)
{
	mesh_.newPyramidWithID(AN1,AN2,AN3,AN4,AN5,AGID);
}

/*----------------------------------------------------------------------------*/
template<typename TMesh>
void IReader<TMesh>::newHex(const TCellID& AN1, const TCellID& AN2,
							const TCellID& AN3, const TCellID& AN4,
							const TCellID& AN5, const TCellID& AN6,
							const TCellID& AN7, const TCellID& AN8,
							const TCellID& AGID)
{
	mesh_.newHexWithID(AN1,AN2,AN3,AN4,AN5,AN6,AN7,AN8,AGID);
}
/*----------------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------------*/
#endif /* GMDS_IREADER_H_ */
/*----------------------------------------------------------------------------*/
