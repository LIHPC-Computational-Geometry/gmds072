/*----------------------------------------------------------------------------*/
/* Copyright: CEA
 * contributors: F. Ledoux and N. Le Goff (2015)
 *
 * franck.ledoux@cea.fr
 * nicolas.le-goff@cea.fr
 *
 * This software is a computer program whose purpose is to provide a set of
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
 * data to be ensured and, more generally, to use and operate it in the 
 * same conditions as regards security. 
 *
 * The fact that you are presently reading this means that you have had
 * knowledge of the CeCILL-C license and that you accept its terms.
 */
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/** \file    MeshSplitter3Refinement.cpp
 *  \author  legoff
 *  \date    23/01/2015
 */
/*----------------------------------------------------------------------------*/
// External File Headers
/*----------------------------------------------------------------------------*/
#include <GMDS/MeshSplitter3Refinement.h>
#include <GMDS/IO/VTKWriter.h>
/*----------------------------------------------------------------------------*/
namespace gmds {
/*----------------------------------------------------------------------------*/
MeshSplitter3Refinement::MeshSplitter3Refinement(
		gmds::IGMesh& AMesh,
		gmds::geom::GeomManager& AManager,
		gmds::geom::GeomMeshIntersectionService& AService)
: MeshSplitter(AMesh,AManager,AService)
{
	tableMarkedNodes<8> key;
	tableMarkedNodes<8> table;

	// all
	key = tableMarkedNodes<8> (0,0,0,0,0,0,0,0); table = tableMarkedNodes<8> (0,0,0,0,0,0,0,0); lookupNodes_[key] = table;
	key = tableMarkedNodes<8> (1,1,1,1,1,1,1,1); table = tableMarkedNodes<8> (1,1,1,1,1,1,1,1); lookupNodes_[key] = table;

	// corners
	key = tableMarkedNodes<8> (1,0,0,0,0,0,0,0); table = tableMarkedNodes<8> (1,0,0,0,0,0,0,0); lookupNodes_[key] = table;
	key = tableMarkedNodes<8> (0,1,0,0,0,0,0,0); table = tableMarkedNodes<8> (0,1,0,0,0,0,0,0); lookupNodes_[key] = table;
	key = tableMarkedNodes<8> (0,0,1,0,0,0,0,0); table = tableMarkedNodes<8> (0,0,1,0,0,0,0,0); lookupNodes_[key] = table;
	key = tableMarkedNodes<8> (0,0,0,1,0,0,0,0); table = tableMarkedNodes<8> (0,0,0,1,0,0,0,0); lookupNodes_[key] = table;
	key = tableMarkedNodes<8> (0,0,0,0,1,0,0,0); table = tableMarkedNodes<8> (0,0,0,0,1,0,0,0); lookupNodes_[key] = table;
	key = tableMarkedNodes<8> (0,0,0,0,0,1,0,0); table = tableMarkedNodes<8> (0,0,0,0,0,1,0,0); lookupNodes_[key] = table;
	key = tableMarkedNodes<8> (0,0,0,0,0,0,1,0); table = tableMarkedNodes<8> (0,0,0,0,0,0,1,0); lookupNodes_[key] = table;
	key = tableMarkedNodes<8> (0,0,0,0,0,0,0,1); table = tableMarkedNodes<8> (0,0,0,0,0,0,0,1); lookupNodes_[key] = table;

	// edges
	key = tableMarkedNodes<8> (1,1,0,0,0,0,0,0); table = tableMarkedNodes<8> (1,1,0,0,0,0,0,0); lookupNodes_[key] = table; // top
	key = tableMarkedNodes<8> (0,1,1,0,0,0,0,0); table = tableMarkedNodes<8> (0,1,1,0,0,0,0,0); lookupNodes_[key] = table;
	key = tableMarkedNodes<8> (0,0,1,1,0,0,0,0); table = tableMarkedNodes<8> (0,0,1,1,0,0,0,0); lookupNodes_[key] = table;
	key = tableMarkedNodes<8> (1,0,0,1,0,0,0,0); table = tableMarkedNodes<8> (1,0,0,1,0,0,0,0); lookupNodes_[key] = table;
	key = tableMarkedNodes<8> (1,0,0,0,1,0,0,0); table = tableMarkedNodes<8> (1,0,0,0,1,0,0,0); lookupNodes_[key] = table; // middle
	key = tableMarkedNodes<8> (0,1,0,0,0,1,0,0); table = tableMarkedNodes<8> (0,1,0,0,0,1,0,0); lookupNodes_[key] = table;
	key = tableMarkedNodes<8> (0,0,1,0,0,0,1,0); table = tableMarkedNodes<8> (0,0,1,0,0,0,1,0); lookupNodes_[key] = table;
	key = tableMarkedNodes<8> (0,0,0,1,0,0,0,1); table = tableMarkedNodes<8> (0,0,0,1,0,0,0,1); lookupNodes_[key] = table;
	key = tableMarkedNodes<8> (0,0,0,0,1,1,0,0); table = tableMarkedNodes<8> (0,0,0,0,1,1,0,0); lookupNodes_[key] = table; // bottom
	key = tableMarkedNodes<8> (0,0,0,0,0,1,1,0); table = tableMarkedNodes<8> (0,0,0,0,0,1,1,0); lookupNodes_[key] = table;
	key = tableMarkedNodes<8> (0,0,0,0,0,0,1,1); table = tableMarkedNodes<8> (0,0,0,0,0,0,1,1); lookupNodes_[key] = table;
	key = tableMarkedNodes<8> (0,0,0,0,1,0,0,1); table = tableMarkedNodes<8> (0,0,0,0,1,0,0,1); lookupNodes_[key] = table;

	// diagonals
	key = tableMarkedNodes<8> (1,0,1,0,0,0,0,0); table = tableMarkedNodes<8> (1,1,1,1,0,0,0,0); lookupNodes_[key] = table; // top
	key = tableMarkedNodes<8> (0,1,0,1,0,0,0,0); table = tableMarkedNodes<8> (1,1,1,1,0,0,0,0); lookupNodes_[key] = table;
	key = tableMarkedNodes<8> (0,0,0,0,1,0,1,0); table = tableMarkedNodes<8> (0,0,0,0,1,1,1,1); lookupNodes_[key] = table; // bottom
	key = tableMarkedNodes<8> (0,0,0,0,0,1,0,1); table = tableMarkedNodes<8> (0,0,0,0,1,1,1,1); lookupNodes_[key] = table;
	key = tableMarkedNodes<8> (1,0,0,0,0,1,0,0); table = tableMarkedNodes<8> (1,1,0,0,1,1,0,0); lookupNodes_[key] = table; // left
	key = tableMarkedNodes<8> (0,1,0,0,1,0,0,0); table = tableMarkedNodes<8> (1,1,0,0,1,1,0,0); lookupNodes_[key] = table;
	key = tableMarkedNodes<8> (0,0,1,0,0,0,0,1); table = tableMarkedNodes<8> (0,0,1,1,0,0,1,1); lookupNodes_[key] = table; // right
	key = tableMarkedNodes<8> (0,0,0,1,0,0,1,0); table = tableMarkedNodes<8> (0,0,1,1,0,0,1,1); lookupNodes_[key] = table;
	key = tableMarkedNodes<8> (0,0,0,1,1,0,0,0); table = tableMarkedNodes<8> (1,0,0,1,1,0,0,1); lookupNodes_[key] = table; // front
	key = tableMarkedNodes<8> (1,0,0,0,0,0,0,1); table = tableMarkedNodes<8> (1,0,0,1,1,0,0,1); lookupNodes_[key] = table;
	key = tableMarkedNodes<8> (0,0,1,0,0,1,0,0); table = tableMarkedNodes<8> (0,1,1,0,0,1,1,0); lookupNodes_[key] = table; // back
	key = tableMarkedNodes<8> (0,1,0,0,0,0,1,0); table = tableMarkedNodes<8> (0,1,1,0,0,1,1,0); lookupNodes_[key] = table;

	// faces
	key = tableMarkedNodes<8> (1,1,1,1,0,0,0,0); table = tableMarkedNodes<8> (1,1,1,1,0,0,0,0); lookupNodes_[key] = table; // top
	key = tableMarkedNodes<8> (0,0,0,0,1,1,1,1); table = tableMarkedNodes<8> (0,0,0,0,1,1,1,1); lookupNodes_[key] = table; // bottom
	key = tableMarkedNodes<8> (1,1,0,0,1,1,0,0); table = tableMarkedNodes<8> (1,1,0,0,1,1,0,0); lookupNodes_[key] = table; // left
	key = tableMarkedNodes<8> (0,0,1,1,0,0,1,1); table = tableMarkedNodes<8> (0,0,1,1,0,0,1,1); lookupNodes_[key] = table; // right
	key = tableMarkedNodes<8> (1,0,0,1,1,0,0,1); table = tableMarkedNodes<8> (1,0,0,1,1,0,0,1); lookupNodes_[key] = table; // front
	key = tableMarkedNodes<8> (0,1,1,0,0,1,1,0); table = tableMarkedNodes<8> (0,1,1,0,0,1,1,0); lookupNodes_[key] = table; // back

	// faces 3 nodes
	key = tableMarkedNodes<8> (0,1,1,1,0,0,0,0); table = tableMarkedNodes<8> (1,1,1,1,0,0,0,0); lookupNodes_[key] = table; // top
	key = tableMarkedNodes<8> (1,0,1,1,0,0,0,0); table = tableMarkedNodes<8> (1,1,1,1,0,0,0,0); lookupNodes_[key] = table;
	key = tableMarkedNodes<8> (1,1,0,1,0,0,0,0); table = tableMarkedNodes<8> (1,1,1,1,0,0,0,0); lookupNodes_[key] = table;
	key = tableMarkedNodes<8> (1,1,1,0,0,0,0,0); table = tableMarkedNodes<8> (1,1,1,1,0,0,0,0); lookupNodes_[key] = table;
	key = tableMarkedNodes<8> (0,0,0,0,0,1,1,1); table = tableMarkedNodes<8> (0,0,0,0,1,1,1,1); lookupNodes_[key] = table; // bottom
	key = tableMarkedNodes<8> (0,0,0,0,1,0,1,1); table = tableMarkedNodes<8> (0,0,0,0,1,1,1,1); lookupNodes_[key] = table;
	key = tableMarkedNodes<8> (0,0,0,0,1,1,0,1); table = tableMarkedNodes<8> (0,0,0,0,1,1,1,1); lookupNodes_[key] = table;
	key = tableMarkedNodes<8> (0,0,0,0,1,1,1,0); table = tableMarkedNodes<8> (0,0,0,0,1,1,1,1); lookupNodes_[key] = table;
	key = tableMarkedNodes<8> (0,1,0,0,1,1,0,0); table = tableMarkedNodes<8> (1,1,0,0,1,1,0,0); lookupNodes_[key] = table; // left
	key = tableMarkedNodes<8> (1,0,0,0,1,1,0,0); table = tableMarkedNodes<8> (1,1,0,0,1,1,0,0); lookupNodes_[key] = table;
	key = tableMarkedNodes<8> (1,1,0,0,0,1,0,0); table = tableMarkedNodes<8> (1,1,0,0,1,1,0,0); lookupNodes_[key] = table;
	key = tableMarkedNodes<8> (1,1,0,0,1,0,0,0); table = tableMarkedNodes<8> (1,1,0,0,1,1,0,0); lookupNodes_[key] = table;
	key = tableMarkedNodes<8> (0,0,0,1,0,0,1,1); table = tableMarkedNodes<8> (0,0,1,1,0,0,1,1); lookupNodes_[key] = table; // right
	key = tableMarkedNodes<8> (0,0,1,0,0,0,1,1); table = tableMarkedNodes<8> (0,0,1,1,0,0,1,1); lookupNodes_[key] = table;
	key = tableMarkedNodes<8> (0,0,1,1,0,0,0,1); table = tableMarkedNodes<8> (0,0,1,1,0,0,1,1); lookupNodes_[key] = table;
	key = tableMarkedNodes<8> (0,0,1,1,0,0,1,0); table = tableMarkedNodes<8> (0,0,1,1,0,0,1,1); lookupNodes_[key] = table;
	key = tableMarkedNodes<8> (0,0,0,1,1,0,0,1); table = tableMarkedNodes<8> (1,0,0,1,1,0,0,1); lookupNodes_[key] = table; // front
	key = tableMarkedNodes<8> (1,0,0,0,1,0,0,1); table = tableMarkedNodes<8> (1,0,0,1,1,0,0,1); lookupNodes_[key] = table;
	key = tableMarkedNodes<8> (1,0,0,1,0,0,0,1); table = tableMarkedNodes<8> (1,0,0,1,1,0,0,1); lookupNodes_[key] = table;
	key = tableMarkedNodes<8> (1,0,0,1,1,0,0,0); table = tableMarkedNodes<8> (1,0,0,1,1,0,0,1); lookupNodes_[key] = table;
	key = tableMarkedNodes<8> (0,0,1,0,0,1,1,0); table = tableMarkedNodes<8> (0,1,1,0,0,1,1,0); lookupNodes_[key] = table; // back
	key = tableMarkedNodes<8> (0,1,0,0,0,1,1,0); table = tableMarkedNodes<8> (0,1,1,0,0,1,1,0); lookupNodes_[key] = table;
	key = tableMarkedNodes<8> (0,1,1,0,0,0,1,0); table = tableMarkedNodes<8> (0,1,1,0,0,1,1,0); lookupNodes_[key] = table;
	key = tableMarkedNodes<8> (0,1,1,0,0,1,0,0); table = tableMarkedNodes<8> (0,1,1,0,0,1,1,0); lookupNodes_[key] = table;



	// faces
	tableMarkedNodes<4> keyFace;

	// all
	keyFace = tableMarkedNodes<4> (0,0,0,0); nbNodesOnFaceToBuild_[keyFace] = 1;
											 nbFacesOnFaceToBuild_[keyFace] = 0;
	keyFace = tableMarkedNodes<4> (1,1,1,1); nbNodesOnFaceToBuild_[keyFace] = 12;
	 	 	 	 	 	 	 	 	 	 	 nbFacesOnFaceToBuild_[keyFace] = 9;

	// corners
	keyFace = tableMarkedNodes<4> (1,0,0,0); nbNodesOnFaceToBuild_[keyFace] = 3;
											 nbFacesOnFaceToBuild_[keyFace] = 3;
	keyFace = tableMarkedNodes<4> (0,1,0,0); nbNodesOnFaceToBuild_[keyFace] = 3;
	 	 	 	 	 	 	 	 	 	 	 nbFacesOnFaceToBuild_[keyFace] = 3;
	keyFace = tableMarkedNodes<4> (0,0,1,0); nbNodesOnFaceToBuild_[keyFace] = 3;
	 	 	 	 	 	 	 	 	 	 	 nbFacesOnFaceToBuild_[keyFace] = 3;
	keyFace = tableMarkedNodes<4> (0,0,0,1); nbNodesOnFaceToBuild_[keyFace] = 3;
	 	 	 	 	 	 	 	 	 	 	 nbFacesOnFaceToBuild_[keyFace] = 3;

	// edges
	keyFace = tableMarkedNodes<4> (1,1,0,0); nbNodesOnFaceToBuild_[keyFace] = 8;
	 	 	 	 	 	 	 	 	 	 	 nbFacesOnFaceToBuild_[keyFace] = 7;
	keyFace = tableMarkedNodes<4> (0,1,1,0); nbNodesOnFaceToBuild_[keyFace] = 8;
	 	 	 	 	 	 	 	 	 	 	 nbFacesOnFaceToBuild_[keyFace] = 7;
	keyFace = tableMarkedNodes<4> (0,0,1,1); nbNodesOnFaceToBuild_[keyFace] = 8;
	 	 	 	 	 	 	 	 	 	 	 nbFacesOnFaceToBuild_[keyFace] = 7;
	keyFace = tableMarkedNodes<4> (1,0,0,1); nbNodesOnFaceToBuild_[keyFace] = 8;
	 	 	 	 	 	 	 	 	 	 	 nbFacesOnFaceToBuild_[keyFace] = 7;



	// all
	std::vector<TCellID> nodesIDs;
	keyFace = tableMarkedNodes<4> (0,0,0,0); nodesOnFaceToBuild_[keyFace] = nodesIDs;

	keyFace = tableMarkedNodes<4> (1,1,1,1); nodesIDs.resize(nbNodesOnFaceToBuild_[keyFace]);
	nodesIDs[0] =  1; nodesIDs[1] =  2; nodesIDs[2]  =  4; nodesIDs[3]  =  5;
	nodesIDs[4] =  6; nodesIDs[5] =  7; nodesIDs[6]  =  8; nodesIDs[7]  =  9;
	nodesIDs[8] = 10; nodesIDs[9] = 11; nodesIDs[10] = 13; nodesIDs[11] = 14;
	nodesOnFaceToBuild_[keyFace] = nodesIDs;

	// corners
	keyFace = tableMarkedNodes<4> (1,0,0,0); nodesIDs.resize(nbNodesOnFaceToBuild_[keyFace]);
	nodesIDs[0] =  1; nodesIDs[1] =  4; nodesIDs[2]  =  5;
	nodesOnFaceToBuild_[keyFace] = nodesIDs;
	keyFace = tableMarkedNodes<4> (0,1,0,0); nodesIDs.resize(nbNodesOnFaceToBuild_[keyFace]);
	nodesIDs[0] =  2; nodesIDs[1] =  6; nodesIDs[2]  =  7;
	nodesOnFaceToBuild_[keyFace] = nodesIDs;
	keyFace = tableMarkedNodes<4> (0,0,1,0); nodesIDs.resize(nbNodesOnFaceToBuild_[keyFace]);
	nodesIDs[0] = 10; nodesIDs[1] = 11; nodesIDs[2]  = 14;
	nodesOnFaceToBuild_[keyFace] = nodesIDs;
	keyFace = tableMarkedNodes<4> (0,0,0,1); nodesIDs.resize(nbNodesOnFaceToBuild_[keyFace]);
	nodesIDs[0] =  8; nodesIDs[1] =  9; nodesIDs[2]  = 13;
	nodesOnFaceToBuild_[keyFace] = nodesIDs;

	// edges
	keyFace = tableMarkedNodes<4> (1,1,0,0); nodesIDs.resize(nbNodesOnFaceToBuild_[keyFace]);
	nodesIDs[0] =  1; nodesIDs[1] =  2; nodesIDs[2] =  4; nodesIDs[3] =  5;
	nodesIDs[4] =  6; nodesIDs[5] =  7; nodesIDs[6] =  9; nodesIDs[7] = 10;
	nodesOnFaceToBuild_[keyFace] = nodesIDs;
	keyFace = tableMarkedNodes<4> (0,1,1,0); nodesIDs.resize(nbNodesOnFaceToBuild_[keyFace]);
	nodesIDs[0] =  2; nodesIDs[1] =  5; nodesIDs[2] =  6; nodesIDs[3] =  7;
	nodesIDs[4] =  9; nodesIDs[5] = 10; nodesIDs[6] = 11; nodesIDs[7] = 14;
	nodesOnFaceToBuild_[keyFace] = nodesIDs;
	keyFace = tableMarkedNodes<4> (0,0,1,1); nodesIDs.resize(nbNodesOnFaceToBuild_[keyFace]);
	nodesIDs[0] =  5; nodesIDs[1] =  6; nodesIDs[2] =  8; nodesIDs[3] =  9;
	nodesIDs[4] = 10; nodesIDs[5] = 11; nodesIDs[6] = 13; nodesIDs[7] = 14;
	nodesOnFaceToBuild_[keyFace] = nodesIDs;
	keyFace = tableMarkedNodes<4> (1,0,0,1); nodesIDs.resize(nbNodesOnFaceToBuild_[keyFace]);
	nodesIDs[0] =  1; nodesIDs[1] =  4; nodesIDs[2] =  5; nodesIDs[3] =  6;
	nodesIDs[4] =  8; nodesIDs[5] =  9; nodesIDs[6] = 10; nodesIDs[7] = 13;
	nodesOnFaceToBuild_[keyFace] = nodesIDs;

	std::vector<FaceIDs> facesIDsArray;
	struct FaceIDs faceIDs;

	// all
	keyFace = tableMarkedNodes<4> (0,0,0,0);
	facesIDsArray.clear();
	facesOnFaceToBuild_[keyFace] = facesIDsArray;

	keyFace = tableMarkedNodes<4> (1,1,1,1);
	facesIDsArray.clear();
	faceIDs.ids[0] =  0; faceIDs.ids[1] =  1; faceIDs.ids[2] =  5; faceIDs.ids[3] =  4;
	facesIDsArray.push_back(faceIDs);
	faceIDs.ids[0] =  1; faceIDs.ids[1] =  2; faceIDs.ids[2] =  6; faceIDs.ids[3] =  5;
	facesIDsArray.push_back(faceIDs);
	faceIDs.ids[0] =  2; faceIDs.ids[1] =  3; faceIDs.ids[2] =  7; faceIDs.ids[3] =  6;
	facesIDsArray.push_back(faceIDs);
	faceIDs.ids[0] =  4; faceIDs.ids[1] =  5; faceIDs.ids[2] =  9; faceIDs.ids[3] =  8;
	facesIDsArray.push_back(faceIDs);
	faceIDs.ids[0] =  5; faceIDs.ids[1] =  6; faceIDs.ids[2] = 10; faceIDs.ids[3] =  9;
	facesIDsArray.push_back(faceIDs);
	faceIDs.ids[0] =  6; faceIDs.ids[1] =  7; faceIDs.ids[2] = 11; faceIDs.ids[3] = 10;
	facesIDsArray.push_back(faceIDs);
	faceIDs.ids[0] =  8; faceIDs.ids[1] =  9; faceIDs.ids[2] = 13; faceIDs.ids[3] = 12;
	facesIDsArray.push_back(faceIDs);
	faceIDs.ids[0] =  9; faceIDs.ids[1] = 10; faceIDs.ids[2] = 14; faceIDs.ids[3] = 13;
	facesIDsArray.push_back(faceIDs);
	faceIDs.ids[0] = 10; faceIDs.ids[1] = 11; faceIDs.ids[2] = 15; faceIDs.ids[3] = 14;
	facesIDsArray.push_back(faceIDs);
	facesOnFaceToBuild_[keyFace] = facesIDsArray;

	// corners
	keyFace = tableMarkedNodes<4> (1,0,0,0);
	facesIDsArray.clear();
	faceIDs.ids[0] =  0; faceIDs.ids[1] =  1; faceIDs.ids[2] =  5; faceIDs.ids[3] =  4;
	facesIDsArray.push_back(faceIDs);
	faceIDs.ids[0] =  1; faceIDs.ids[1] =  3; faceIDs.ids[2] = 15; faceIDs.ids[3] =  5;
	facesIDsArray.push_back(faceIDs);
	faceIDs.ids[0] =  4; faceIDs.ids[1] =  5; faceIDs.ids[2] = 15; faceIDs.ids[3] = 12;
	facesIDsArray.push_back(faceIDs);
	facesOnFaceToBuild_[keyFace] = facesIDsArray;
	keyFace = tableMarkedNodes<4> (0,1,0,0);
	facesIDsArray.clear();
	faceIDs.ids[0] =  0; faceIDs.ids[1] =  2; faceIDs.ids[2] =  6; faceIDs.ids[3] = 12;
	facesIDsArray.push_back(faceIDs);
	faceIDs.ids[0] =  2; faceIDs.ids[1] =  3; faceIDs.ids[2] =  7; faceIDs.ids[3] =  6;
	facesIDsArray.push_back(faceIDs);
	faceIDs.ids[0] =  6; faceIDs.ids[1] =  7; faceIDs.ids[2] = 15; faceIDs.ids[3] = 12;
	facesIDsArray.push_back(faceIDs);
	facesOnFaceToBuild_[keyFace] = facesIDsArray;
	keyFace = tableMarkedNodes<4> (0,0,1,0);
	facesIDsArray.clear();
	faceIDs.ids[0] =  0; faceIDs.ids[1] =  3; faceIDs.ids[2] = 11; faceIDs.ids[3] = 10;
	facesIDsArray.push_back(faceIDs);
	faceIDs.ids[0] =  0; faceIDs.ids[1] = 10; faceIDs.ids[2] = 14; faceIDs.ids[3] = 12;
	facesIDsArray.push_back(faceIDs);
	faceIDs.ids[0] = 10; faceIDs.ids[1] = 11; faceIDs.ids[2] = 15; faceIDs.ids[3] = 14;
	facesIDsArray.push_back(faceIDs);
	facesOnFaceToBuild_[keyFace] = facesIDsArray;
	keyFace = tableMarkedNodes<4> (0,0,0,1);
	facesIDsArray.clear();
	faceIDs.ids[0] =  0; faceIDs.ids[1] =  3; faceIDs.ids[2] =  9; faceIDs.ids[3] =  8;
	facesIDsArray.push_back(faceIDs);
	faceIDs.ids[0] =  9; faceIDs.ids[1] =  3; faceIDs.ids[2] = 15; faceIDs.ids[3] = 13;
	facesIDsArray.push_back(faceIDs);
	faceIDs.ids[0] =  8; faceIDs.ids[1] =  9; faceIDs.ids[2] = 13; faceIDs.ids[3] = 12;
	facesIDsArray.push_back(faceIDs);
	facesOnFaceToBuild_[keyFace] = facesIDsArray;

	// edges
	keyFace = tableMarkedNodes<4> (1,1,0,0);
	facesIDsArray.clear();
	faceIDs.ids[0] =  0; faceIDs.ids[1] =  1; faceIDs.ids[2] =  5; faceIDs.ids[3] =  4;
	facesIDsArray.push_back(faceIDs);
	faceIDs.ids[0] =  1; faceIDs.ids[1] =  2; faceIDs.ids[2] =  6; faceIDs.ids[3] =  5;
	facesIDsArray.push_back(faceIDs);
	faceIDs.ids[0] =  2; faceIDs.ids[1] =  3; faceIDs.ids[2] =  7; faceIDs.ids[3] =  6;
	facesIDsArray.push_back(faceIDs);
	faceIDs.ids[0] =  4; faceIDs.ids[1] =  5; faceIDs.ids[2] =  9; faceIDs.ids[3] = 12;
	facesIDsArray.push_back(faceIDs);
	faceIDs.ids[0] =  5; faceIDs.ids[1] =  6; faceIDs.ids[2] = 10; faceIDs.ids[3] =  9;
	facesIDsArray.push_back(faceIDs);
	faceIDs.ids[0] =  6; faceIDs.ids[1] =  7; faceIDs.ids[2] = 15; faceIDs.ids[3] = 10;
	facesIDsArray.push_back(faceIDs);
	faceIDs.ids[0] =  9; faceIDs.ids[1] = 10; faceIDs.ids[2] = 15; faceIDs.ids[3] = 12;
	facesIDsArray.push_back(faceIDs);
	facesOnFaceToBuild_[keyFace] = facesIDsArray;
	keyFace = tableMarkedNodes<4> (0,1,1,0);
	facesIDsArray.clear();
	faceIDs.ids[0] =  0; faceIDs.ids[1] =  5; faceIDs.ids[2] =  9; faceIDs.ids[3] = 12;
	facesIDsArray.push_back(faceIDs);
	faceIDs.ids[0] =  0; faceIDs.ids[1] =  2; faceIDs.ids[2] =  6; faceIDs.ids[3] =  5;
	facesIDsArray.push_back(faceIDs);
	faceIDs.ids[0] =  2; faceIDs.ids[1] =  3; faceIDs.ids[2] =  7; faceIDs.ids[3] =  6;
	facesIDsArray.push_back(faceIDs);
	faceIDs.ids[0] =  5; faceIDs.ids[1] =  6; faceIDs.ids[2] = 10; faceIDs.ids[3] =  9;
	facesIDsArray.push_back(faceIDs);
	faceIDs.ids[0] =  6; faceIDs.ids[1] =  7; faceIDs.ids[2] = 11; faceIDs.ids[3] = 10;
	facesIDsArray.push_back(faceIDs);
	faceIDs.ids[0] =  9; faceIDs.ids[1] = 10; faceIDs.ids[2] = 14; faceIDs.ids[3] = 12;
	facesIDsArray.push_back(faceIDs);
	faceIDs.ids[0] = 10; faceIDs.ids[1] = 11; faceIDs.ids[2] = 15; faceIDs.ids[3] = 14;
	facesIDsArray.push_back(faceIDs);
	facesOnFaceToBuild_[keyFace] = facesIDsArray;
	keyFace = tableMarkedNodes<4> (0,0,1,1);
	facesIDsArray.clear();
	faceIDs.ids[0] =  0; faceIDs.ids[1] =  3; faceIDs.ids[2] =  6; faceIDs.ids[3] =  5;
	facesIDsArray.push_back(faceIDs);
	faceIDs.ids[0] =  0; faceIDs.ids[1] =  5; faceIDs.ids[2] =  9; faceIDs.ids[3] =  8;
	facesIDsArray.push_back(faceIDs);
	faceIDs.ids[0] =  5; faceIDs.ids[1] =  6; faceIDs.ids[2] = 10; faceIDs.ids[3] =  9;
	facesIDsArray.push_back(faceIDs);
	faceIDs.ids[0] =  6; faceIDs.ids[1] =  3; faceIDs.ids[2] = 11; faceIDs.ids[3] = 10;
	facesIDsArray.push_back(faceIDs);
	faceIDs.ids[0] =  8; faceIDs.ids[1] =  9; faceIDs.ids[2] = 13; faceIDs.ids[3] = 12;
	facesIDsArray.push_back(faceIDs);
	faceIDs.ids[0] =  9; faceIDs.ids[1] = 10; faceIDs.ids[2] = 14; faceIDs.ids[3] = 13;
	facesIDsArray.push_back(faceIDs);
	faceIDs.ids[0] = 10; faceIDs.ids[1] = 11; faceIDs.ids[2] = 15; faceIDs.ids[3] = 14;
	facesIDsArray.push_back(faceIDs);
	facesOnFaceToBuild_[keyFace] = facesIDsArray;
	keyFace = tableMarkedNodes<4> (1,0,0,1);
	facesIDsArray.clear();
	faceIDs.ids[0] =  0; faceIDs.ids[1] =  1; faceIDs.ids[2] =  5; faceIDs.ids[3] =  4;
	facesIDsArray.push_back(faceIDs);
	faceIDs.ids[0] =  1; faceIDs.ids[1] =  3; faceIDs.ids[2] =  6; faceIDs.ids[3] =  5;
	facesIDsArray.push_back(faceIDs);
	faceIDs.ids[0] =  4; faceIDs.ids[1] =  5; faceIDs.ids[2] =  9; faceIDs.ids[3] =  8;
	facesIDsArray.push_back(faceIDs);
	faceIDs.ids[0] =  5; faceIDs.ids[1] =  6; faceIDs.ids[2] = 10; faceIDs.ids[3] =  9;
	facesIDsArray.push_back(faceIDs);
	faceIDs.ids[0] =  6; faceIDs.ids[1] =  3; faceIDs.ids[2] = 15; faceIDs.ids[3] = 10;
	facesIDsArray.push_back(faceIDs);
	faceIDs.ids[0] =  8; faceIDs.ids[1] =  9; faceIDs.ids[2] = 13; faceIDs.ids[3] = 12;
	facesIDsArray.push_back(faceIDs);
	faceIDs.ids[0] =  9; faceIDs.ids[1] = 10; faceIDs.ids[2] = 15; faceIDs.ids[3] = 13;
	facesIDsArray.push_back(faceIDs);
	facesOnFaceToBuild_[keyFace] = facesIDsArray;


	tableMarkedNodes<8> keyRegion;
	struct RegionIDs regionPermutIDs;

	// all
	keyRegion = tableMarkedNodes<8> (0,0,0,0,0,0,0,0);
	regionPermutIDs.ids[0] = 0; regionPermutIDs.ids[1] = 1; regionPermutIDs.ids[2] = 2; regionPermutIDs.ids[3] = 3;
	regionPermutIDs.ids[4] = 4; regionPermutIDs.ids[5] = 5; regionPermutIDs.ids[6] = 6; regionPermutIDs.ids[7] = 7;
	nodePermutToUnitHex_[keyRegion] = regionPermutIDs;
	keyRegion = tableMarkedNodes<8> (1,1,1,1,1,1,1,1);
	regionPermutIDs.ids[0] = 0; regionPermutIDs.ids[1] = 1; regionPermutIDs.ids[2] = 2; regionPermutIDs.ids[3] = 3;
	regionPermutIDs.ids[4] = 4; regionPermutIDs.ids[5] = 5; regionPermutIDs.ids[6] = 6; regionPermutIDs.ids[7] = 7;
	nodePermutToUnitHex_[keyRegion] = regionPermutIDs;

	// corners
	keyRegion = tableMarkedNodes<8> (1,0,0,0,0,0,0,0);
	regionPermutIDs.ids[0] = 0; regionPermutIDs.ids[1] = 1; regionPermutIDs.ids[2] = 2; regionPermutIDs.ids[3] = 3;
	regionPermutIDs.ids[4] = 4; regionPermutIDs.ids[5] = 5; regionPermutIDs.ids[6] = 6; regionPermutIDs.ids[7] = 7;
	nodePermutToUnitHex_[keyRegion] = regionPermutIDs;
	keyRegion = tableMarkedNodes<8> (0,1,0,0,0,0,0,0);
	regionPermutIDs.ids[0] = 3; regionPermutIDs.ids[1] = 0; regionPermutIDs.ids[2] = 1; regionPermutIDs.ids[3] = 2;
	regionPermutIDs.ids[4] = 7; regionPermutIDs.ids[5] = 4; regionPermutIDs.ids[6] = 5; regionPermutIDs.ids[7] = 6;
	nodePermutToUnitHex_[keyRegion] = regionPermutIDs;
	keyRegion = tableMarkedNodes<8> (0,0,1,0,0,0,0,0);
	regionPermutIDs.ids[0] = 2; regionPermutIDs.ids[1] = 3; regionPermutIDs.ids[2] = 0; regionPermutIDs.ids[3] = 1;
	regionPermutIDs.ids[4] = 6; regionPermutIDs.ids[5] = 7; regionPermutIDs.ids[6] = 4; regionPermutIDs.ids[7] = 5;
	nodePermutToUnitHex_[keyRegion] = regionPermutIDs;
	keyRegion = tableMarkedNodes<8> (0,0,0,1,0,0,0,0);
	regionPermutIDs.ids[0] = 1; regionPermutIDs.ids[1] = 2; regionPermutIDs.ids[2] = 3; regionPermutIDs.ids[3] = 0;
	regionPermutIDs.ids[4] = 5; regionPermutIDs.ids[5] = 6; regionPermutIDs.ids[6] = 7; regionPermutIDs.ids[7] = 4;
	nodePermutToUnitHex_[keyRegion] = regionPermutIDs;
	keyRegion = tableMarkedNodes<8> (0,0,0,0,1,0,0,0);
	regionPermutIDs.ids[0] = 1; regionPermutIDs.ids[1] = 5; regionPermutIDs.ids[2] = 6; regionPermutIDs.ids[3] = 2;
	regionPermutIDs.ids[4] = 0; regionPermutIDs.ids[5] = 4; regionPermutIDs.ids[6] = 7; regionPermutIDs.ids[7] = 3;
	nodePermutToUnitHex_[keyRegion] = regionPermutIDs;
	keyRegion = tableMarkedNodes<8> (0,0,0,0,0,1,0,0);
	regionPermutIDs.ids[0] = 5; regionPermutIDs.ids[1] = 4; regionPermutIDs.ids[2] = 7; regionPermutIDs.ids[3] = 6;
	regionPermutIDs.ids[4] = 1; regionPermutIDs.ids[5] = 0; regionPermutIDs.ids[6] = 3; regionPermutIDs.ids[7] = 2;
	nodePermutToUnitHex_[keyRegion] = regionPermutIDs;
	keyRegion = tableMarkedNodes<8> (0,0,0,0,0,0,1,0);
	regionPermutIDs.ids[0] = 6; regionPermutIDs.ids[1] = 2; regionPermutIDs.ids[2] = 1; regionPermutIDs.ids[3] = 5;
	regionPermutIDs.ids[4] = 7; regionPermutIDs.ids[5] = 3; regionPermutIDs.ids[6] = 0; regionPermutIDs.ids[7] = 4;
	nodePermutToUnitHex_[keyRegion] = regionPermutIDs;
	keyRegion = tableMarkedNodes<8> (0,0,0,0,0,0,0,1);
	regionPermutIDs.ids[0] = 7; regionPermutIDs.ids[1] = 6; regionPermutIDs.ids[2] = 5; regionPermutIDs.ids[3] = 4;
	regionPermutIDs.ids[4] = 3; regionPermutIDs.ids[5] = 2; regionPermutIDs.ids[6] = 1; regionPermutIDs.ids[7] = 0;
	nodePermutToUnitHex_[keyRegion] = regionPermutIDs;

	// edges
	keyRegion = tableMarkedNodes<8> (1,1,0,0,0,0,0,0); // top
	regionPermutIDs.ids[0] = 0; regionPermutIDs.ids[1] = 1; regionPermutIDs.ids[2] = 2; regionPermutIDs.ids[3] = 3;
	regionPermutIDs.ids[4] = 4; regionPermutIDs.ids[5] = 5; regionPermutIDs.ids[6] = 6; regionPermutIDs.ids[7] = 7;
	nodePermutToUnitHex_[keyRegion] = regionPermutIDs;
	keyRegion = tableMarkedNodes<8> (0,1,1,0,0,0,0,0);
	regionPermutIDs.ids[0] = 3; regionPermutIDs.ids[1] = 0; regionPermutIDs.ids[2] = 1; regionPermutIDs.ids[3] = 2;
	regionPermutIDs.ids[4] = 7; regionPermutIDs.ids[5] = 4; regionPermutIDs.ids[6] = 5; regionPermutIDs.ids[7] = 6;
	nodePermutToUnitHex_[keyRegion] = regionPermutIDs;
	keyRegion = tableMarkedNodes<8> (0,0,1,1,0,0,0,0);
	regionPermutIDs.ids[0] = 2; regionPermutIDs.ids[1] = 3; regionPermutIDs.ids[2] = 0; regionPermutIDs.ids[3] = 1;
	regionPermutIDs.ids[4] = 6; regionPermutIDs.ids[5] = 7; regionPermutIDs.ids[6] = 4; regionPermutIDs.ids[7] = 5;
	nodePermutToUnitHex_[keyRegion] = regionPermutIDs;
	keyRegion = tableMarkedNodes<8> (1,0,0,1,0,0,0,0);
	regionPermutIDs.ids[0] = 1; regionPermutIDs.ids[1] = 2; regionPermutIDs.ids[2] = 3; regionPermutIDs.ids[3] = 0;
	regionPermutIDs.ids[4] = 5; regionPermutIDs.ids[5] = 6; regionPermutIDs.ids[6] = 7; regionPermutIDs.ids[7] = 4;
	nodePermutToUnitHex_[keyRegion] = regionPermutIDs;
	keyRegion = tableMarkedNodes<8> (1,0,0,0,1,0,0,0); // middle
	regionPermutIDs.ids[0] = 1; regionPermutIDs.ids[1] = 5; regionPermutIDs.ids[2] = 6; regionPermutIDs.ids[3] = 2;
	regionPermutIDs.ids[4] = 0; regionPermutIDs.ids[5] = 4; regionPermutIDs.ids[6] = 7; regionPermutIDs.ids[7] = 3;
	nodePermutToUnitHex_[keyRegion] = regionPermutIDs;
	keyRegion = tableMarkedNodes<8> (0,1,0,0,0,1,0,0);
	regionPermutIDs.ids[0] = 4; regionPermutIDs.ids[1] = 0; regionPermutIDs.ids[2] = 3; regionPermutIDs.ids[3] = 7;
	regionPermutIDs.ids[4] = 5; regionPermutIDs.ids[5] = 1; regionPermutIDs.ids[6] = 2; regionPermutIDs.ids[7] = 6;
	nodePermutToUnitHex_[keyRegion] = regionPermutIDs;
	keyRegion = tableMarkedNodes<8> (0,0,1,0,0,0,1,0);
	regionPermutIDs.ids[0] = 7; regionPermutIDs.ids[1] = 4; regionPermutIDs.ids[2] = 0; regionPermutIDs.ids[3] = 3;
	regionPermutIDs.ids[4] = 6; regionPermutIDs.ids[5] = 5; regionPermutIDs.ids[6] = 1; regionPermutIDs.ids[7] = 2;
	nodePermutToUnitHex_[keyRegion] = regionPermutIDs;
	keyRegion = tableMarkedNodes<8> (0,0,0,1,0,0,0,1);
	regionPermutIDs.ids[0] = 5; regionPermutIDs.ids[1] = 6; regionPermutIDs.ids[2] = 2; regionPermutIDs.ids[3] = 1;
	regionPermutIDs.ids[4] = 4; regionPermutIDs.ids[5] = 7; regionPermutIDs.ids[6] = 3; regionPermutIDs.ids[7] = 0;
	nodePermutToUnitHex_[keyRegion] = regionPermutIDs;
	keyRegion = tableMarkedNodes<8> (0,0,0,0,1,1,0,0); // bottom
	regionPermutIDs.ids[0] = 3; regionPermutIDs.ids[1] = 2; regionPermutIDs.ids[2] = 6; regionPermutIDs.ids[3] = 7;
	regionPermutIDs.ids[4] = 0; regionPermutIDs.ids[5] = 1; regionPermutIDs.ids[6] = 5; regionPermutIDs.ids[7] = 4;
	nodePermutToUnitHex_[keyRegion] = regionPermutIDs;
	keyRegion = tableMarkedNodes<8> (0,0,0,0,0,1,1,0);
	regionPermutIDs.ids[0] = 6; regionPermutIDs.ids[1] = 5; regionPermutIDs.ids[2] = 4; regionPermutIDs.ids[3] = 7;
	regionPermutIDs.ids[4] = 2; regionPermutIDs.ids[5] = 1; regionPermutIDs.ids[6] = 0; regionPermutIDs.ids[7] = 3;
	nodePermutToUnitHex_[keyRegion] = regionPermutIDs;
	keyRegion = tableMarkedNodes<8> (0,0,0,0,0,0,1,1);
	regionPermutIDs.ids[0] = 7; regionPermutIDs.ids[1] = 6; regionPermutIDs.ids[2] = 5; regionPermutIDs.ids[3] = 4;
	regionPermutIDs.ids[4] = 3; regionPermutIDs.ids[5] = 2; regionPermutIDs.ids[6] = 1; regionPermutIDs.ids[7] = 0;
	nodePermutToUnitHex_[keyRegion] = regionPermutIDs;
	keyRegion = tableMarkedNodes<8> (0,0,0,0,1,0,0,1);
	regionPermutIDs.ids[0] = 4; regionPermutIDs.ids[1] = 7; regionPermutIDs.ids[2] = 6; regionPermutIDs.ids[3] = 5;
	regionPermutIDs.ids[4] = 0; regionPermutIDs.ids[5] = 3; regionPermutIDs.ids[6] = 2; regionPermutIDs.ids[7] = 1;
	nodePermutToUnitHex_[keyRegion] = regionPermutIDs;


	// faces
	keyRegion = tableMarkedNodes<8> (1,1,1,1,0,0,0,0); // top
	regionPermutIDs.ids[0] = 0; regionPermutIDs.ids[1] = 1; regionPermutIDs.ids[2] = 2; regionPermutIDs.ids[3] = 3;
	regionPermutIDs.ids[4] = 4; regionPermutIDs.ids[5] = 5; regionPermutIDs.ids[6] = 6; regionPermutIDs.ids[7] = 7;
	nodePermutToUnitHex_[keyRegion] = regionPermutIDs;
	keyRegion = tableMarkedNodes<8> (0,0,0,0,1,1,1,1); // bottom
	regionPermutIDs.ids[0] = 4; regionPermutIDs.ids[1] = 7; regionPermutIDs.ids[2] = 6; regionPermutIDs.ids[3] = 5;
	regionPermutIDs.ids[4] = 0; regionPermutIDs.ids[5] = 3; regionPermutIDs.ids[6] = 2; regionPermutIDs.ids[7] = 1;
	nodePermutToUnitHex_[keyRegion] = regionPermutIDs;
	keyRegion = tableMarkedNodes<8> (1,1,0,0,1,1,0,0); // left
	regionPermutIDs.ids[0] = 3; regionPermutIDs.ids[1] = 2; regionPermutIDs.ids[2] = 6; regionPermutIDs.ids[3] = 7;
	regionPermutIDs.ids[4] = 0; regionPermutIDs.ids[5] = 1; regionPermutIDs.ids[6] = 5; regionPermutIDs.ids[7] = 4;
	nodePermutToUnitHex_[keyRegion] = regionPermutIDs;
	keyRegion = tableMarkedNodes<8> (0,0,1,1,0,0,1,1); // right
	regionPermutIDs.ids[0] = 4; regionPermutIDs.ids[1] = 5; regionPermutIDs.ids[2] = 1; regionPermutIDs.ids[3] = 0;
	regionPermutIDs.ids[4] = 7; regionPermutIDs.ids[5] = 6; regionPermutIDs.ids[6] = 2; regionPermutIDs.ids[7] = 3;
	nodePermutToUnitHex_[keyRegion] = regionPermutIDs;
	keyRegion = tableMarkedNodes<8> (1,0,0,1,1,0,0,1); // front
	regionPermutIDs.ids[0] = 1; regionPermutIDs.ids[1] = 5; regionPermutIDs.ids[2] = 6; regionPermutIDs.ids[3] = 2;
	regionPermutIDs.ids[4] = 0; regionPermutIDs.ids[5] = 4; regionPermutIDs.ids[6] = 7; regionPermutIDs.ids[7] = 3;
	nodePermutToUnitHex_[keyRegion] = regionPermutIDs;
	keyRegion = tableMarkedNodes<8> (0,1,1,0,0,1,1,0); // back
	regionPermutIDs.ids[0] = 4; regionPermutIDs.ids[1] = 0; regionPermutIDs.ids[2] = 3; regionPermutIDs.ids[3] = 7;
	regionPermutIDs.ids[4] = 5; regionPermutIDs.ids[5] = 1; regionPermutIDs.ids[6] = 2; regionPermutIDs.ids[7] = 6;
	nodePermutToUnitHex_[keyRegion] = regionPermutIDs;



	keyRegion = tableMarkedNodes<8> (0,0,0,0,0,0,0,0); nbNodesOnRegionToBuild_[keyRegion] =  0; // nothing
	keyRegion = tableMarkedNodes<8> (1,0,0,0,0,0,0,0); nbNodesOnRegionToBuild_[keyRegion] =  7; // corner
	keyRegion = tableMarkedNodes<8> (1,1,0,0,0,0,0,0); nbNodesOnRegionToBuild_[keyRegion] = 20; // edge
	keyRegion = tableMarkedNodes<8> (1,1,1,1,0,0,0,0); nbNodesOnRegionToBuild_[keyRegion] = 40; // face
	keyRegion = tableMarkedNodes<8> (1,1,1,1,1,1,1,1); nbNodesOnRegionToBuild_[keyRegion] = 56; // region
//	keyRegion = tableMarkedNodes<8> (0,0,0,0,0,0,0,0); nbNodesOnRegionToBuild_[keyRegion] =  8; // nothing
//	keyRegion = tableMarkedNodes<8> (1,0,0,0,0,0,0,0); nbNodesOnRegionToBuild_[keyRegion] = 15; // corner
//	keyRegion = tableMarkedNodes<8> (1,1,0,0,0,0,0,0); nbNodesOnRegionToBuild_[keyRegion] = 28; // edge
//	keyRegion = tableMarkedNodes<8> (1,1,1,1,0,0,0,0); nbNodesOnRegionToBuild_[keyRegion] = 48; // face
//	keyRegion = tableMarkedNodes<8> (1,1,1,1,1,1,1,1); nbNodesOnRegionToBuild_[keyRegion] = 64; // region

	keyRegion = tableMarkedNodes<8> (0,0,0,0,0,0,0,0); nbRegionsOnRegionToBuild_[keyRegion] =  1; // nothing
	keyRegion = tableMarkedNodes<8> (1,0,0,0,0,0,0,0); nbRegionsOnRegionToBuild_[keyRegion] =  4; // corner
	keyRegion = tableMarkedNodes<8> (1,1,0,0,0,0,0,0); nbRegionsOnRegionToBuild_[keyRegion] = 11; // edge
	keyRegion = tableMarkedNodes<8> (1,1,1,1,0,0,0,0); nbRegionsOnRegionToBuild_[keyRegion] = 22; // face
	keyRegion = tableMarkedNodes<8> (1,1,1,1,1,1,1,1); nbRegionsOnRegionToBuild_[keyRegion] = 27; // region

	std::vector<NodeHexIJK> nodesIJKs;
	keyRegion = tableMarkedNodes<8> (0,0,0,0,0,0,0,0); // nothing
	nodesIJKs.clear(); nodesIJKs.resize(nbNodesOnRegionToBuild_[keyRegion]);
	nodesOnRegionToBuild_[keyRegion] = nodesIJKs;
	keyRegion = tableMarkedNodes<8> (1,0,0,0,0,0,0,0); // corner
	nodesIJKs.clear(); nodesIJKs.resize(nbNodesOnRegionToBuild_[keyRegion]);
	nodesIJKs[ 0].ids[0]=0; nodesIJKs[ 0].ids[1]=0; nodesIJKs[ 0].ids[2]=4;
	nodesIJKs[ 1].ids[0]=2; nodesIJKs[ 1].ids[1]=0; nodesIJKs[ 1].ids[2]=4;
	nodesIJKs[ 2].ids[0]=0; nodesIJKs[ 2].ids[1]=2; nodesIJKs[ 2].ids[2]=4;
	nodesIJKs[ 3].ids[0]=2; nodesIJKs[ 3].ids[1]=2; nodesIJKs[ 3].ids[2]=4;

	nodesIJKs[ 4].ids[0]=2; nodesIJKs[ 4].ids[1]=0; nodesIJKs[ 4].ids[2]=6;
	nodesIJKs[ 5].ids[0]=0; nodesIJKs[ 5].ids[1]=2; nodesIJKs[ 5].ids[2]=6;
	nodesIJKs[ 6].ids[0]=2; nodesIJKs[ 6].ids[1]=2; nodesIJKs[ 6].ids[2]=6;

	nodesOnRegionToBuild_[keyRegion] = nodesIJKs;

	keyRegion = tableMarkedNodes<8> (1,1,0,0,0,0,0,0); // edge
	nodesIJKs.clear(); nodesIJKs.resize(nbNodesOnRegionToBuild_[keyRegion]);
	nodesIJKs[ 0].ids[0]=0; nodesIJKs[ 0].ids[1]=2; nodesIJKs[ 0].ids[2]=2;
	nodesIJKs[ 1].ids[0]=4; nodesIJKs[ 1].ids[1]=2; nodesIJKs[ 1].ids[2]=2;
	nodesIJKs[ 2].ids[0]=0; nodesIJKs[ 2].ids[1]=4; nodesIJKs[ 2].ids[2]=2;
	nodesIJKs[ 3].ids[0]=4; nodesIJKs[ 3].ids[1]=4; nodesIJKs[ 3].ids[2]=2;

	nodesIJKs[ 4].ids[0]=0; nodesIJKs[ 4].ids[1]=0; nodesIJKs[ 4].ids[2]=4;
	nodesIJKs[ 5].ids[0]=2; nodesIJKs[ 5].ids[1]=0; nodesIJKs[ 5].ids[2]=4;
	nodesIJKs[ 6].ids[0]=0; nodesIJKs[ 6].ids[1]=2; nodesIJKs[ 6].ids[2]=4;
	nodesIJKs[ 7].ids[0]=2; nodesIJKs[ 7].ids[1]=2; nodesIJKs[ 7].ids[2]=4;
	nodesIJKs[ 8].ids[0]=0; nodesIJKs[ 8].ids[1]=4; nodesIJKs[ 8].ids[2]=4;
	nodesIJKs[ 9].ids[0]=2; nodesIJKs[ 9].ids[1]=4; nodesIJKs[ 9].ids[2]=4;
	nodesIJKs[10].ids[0]=0; nodesIJKs[10].ids[1]=6; nodesIJKs[10].ids[2]=4;
	nodesIJKs[11].ids[0]=2; nodesIJKs[11].ids[1]=6; nodesIJKs[11].ids[2]=4;

	nodesIJKs[12].ids[0]=2; nodesIJKs[12].ids[1]=0; nodesIJKs[12].ids[2]=6;
	nodesIJKs[13].ids[0]=0; nodesIJKs[13].ids[1]=2; nodesIJKs[13].ids[2]=6;
	nodesIJKs[14].ids[0]=2; nodesIJKs[14].ids[1]=2; nodesIJKs[14].ids[2]=6;
	nodesIJKs[15].ids[0]=4; nodesIJKs[15].ids[1]=2; nodesIJKs[15].ids[2]=6;
	nodesIJKs[16].ids[0]=0; nodesIJKs[16].ids[1]=4; nodesIJKs[16].ids[2]=6;
	nodesIJKs[17].ids[0]=2; nodesIJKs[17].ids[1]=4; nodesIJKs[17].ids[2]=6;
	nodesIJKs[18].ids[0]=4; nodesIJKs[18].ids[1]=4; nodesIJKs[18].ids[2]=6;
	nodesIJKs[19].ids[0]=2; nodesIJKs[19].ids[1]=6; nodesIJKs[19].ids[2]=6;

	nodesOnRegionToBuild_[keyRegion] = nodesIJKs;

	keyRegion = tableMarkedNodes<8> (1,1,1,1,0,0,0,0); // face
	nodesIJKs.clear(); nodesIJKs.resize(nbNodesOnRegionToBuild_[keyRegion]);
	nodesIJKs[ 0].ids[0]=2; nodesIJKs[ 0].ids[1]=0; nodesIJKs[ 0].ids[2]=2;
	nodesIJKs[ 1].ids[0]=4; nodesIJKs[ 1].ids[1]=0; nodesIJKs[ 1].ids[2]=2;
	nodesIJKs[ 2].ids[0]=0; nodesIJKs[ 2].ids[1]=2; nodesIJKs[ 2].ids[2]=2;
	nodesIJKs[ 3].ids[0]=6; nodesIJKs[ 3].ids[1]=2; nodesIJKs[ 3].ids[2]=2;
	nodesIJKs[ 4].ids[0]=0; nodesIJKs[ 4].ids[1]=4; nodesIJKs[ 4].ids[2]=2;
	nodesIJKs[ 5].ids[0]=6; nodesIJKs[ 5].ids[1]=4; nodesIJKs[ 5].ids[2]=2;
	nodesIJKs[ 6].ids[0]=2; nodesIJKs[ 6].ids[1]=6; nodesIJKs[ 6].ids[2]=2;
	nodesIJKs[ 7].ids[0]=4; nodesIJKs[ 7].ids[1]=6; nodesIJKs[ 7].ids[2]=2;

	nodesIJKs[ 8].ids[0]=2; nodesIJKs[ 8].ids[1]=2; nodesIJKs[ 8].ids[2]=3;
	nodesIJKs[ 9].ids[0]=4; nodesIJKs[ 9].ids[1]=2; nodesIJKs[ 9].ids[2]=3;
	nodesIJKs[10].ids[0]=2; nodesIJKs[10].ids[1]=4; nodesIJKs[10].ids[2]=3;
	nodesIJKs[11].ids[0]=4; nodesIJKs[11].ids[1]=4; nodesIJKs[11].ids[2]=3;

	nodesIJKs[12].ids[0]=0; nodesIJKs[12].ids[1]=0; nodesIJKs[12].ids[2]=4;
	nodesIJKs[13].ids[0]=2; nodesIJKs[13].ids[1]=0; nodesIJKs[13].ids[2]=4;
	nodesIJKs[14].ids[0]=4; nodesIJKs[14].ids[1]=0; nodesIJKs[14].ids[2]=4;
	nodesIJKs[15].ids[0]=6; nodesIJKs[15].ids[1]=0; nodesIJKs[15].ids[2]=4;
	nodesIJKs[16].ids[0]=0; nodesIJKs[16].ids[1]=2; nodesIJKs[16].ids[2]=4;
	nodesIJKs[17].ids[0]=2; nodesIJKs[17].ids[1]=2; nodesIJKs[17].ids[2]=4;
	nodesIJKs[18].ids[0]=4; nodesIJKs[18].ids[1]=2; nodesIJKs[18].ids[2]=4;
	nodesIJKs[19].ids[0]=6; nodesIJKs[19].ids[1]=2; nodesIJKs[19].ids[2]=4;
	nodesIJKs[20].ids[0]=0; nodesIJKs[20].ids[1]=4; nodesIJKs[20].ids[2]=4;
	nodesIJKs[21].ids[0]=2; nodesIJKs[21].ids[1]=4; nodesIJKs[21].ids[2]=4;
	nodesIJKs[22].ids[0]=4; nodesIJKs[22].ids[1]=4; nodesIJKs[22].ids[2]=4;
	nodesIJKs[23].ids[0]=6; nodesIJKs[23].ids[1]=4; nodesIJKs[23].ids[2]=4;
	nodesIJKs[24].ids[0]=0; nodesIJKs[24].ids[1]=6; nodesIJKs[24].ids[2]=4;
	nodesIJKs[25].ids[0]=2; nodesIJKs[25].ids[1]=6; nodesIJKs[25].ids[2]=4;
	nodesIJKs[26].ids[0]=4; nodesIJKs[26].ids[1]=6; nodesIJKs[26].ids[2]=4;
	nodesIJKs[27].ids[0]=6; nodesIJKs[27].ids[1]=6; nodesIJKs[27].ids[2]=4;

	nodesIJKs[28].ids[0]=2; nodesIJKs[28].ids[1]=0; nodesIJKs[28].ids[2]=6;
	nodesIJKs[29].ids[0]=4; nodesIJKs[29].ids[1]=0; nodesIJKs[29].ids[2]=6;
	nodesIJKs[30].ids[0]=0; nodesIJKs[30].ids[1]=2; nodesIJKs[30].ids[2]=6;
	nodesIJKs[31].ids[0]=2; nodesIJKs[31].ids[1]=2; nodesIJKs[31].ids[2]=6;
	nodesIJKs[32].ids[0]=4; nodesIJKs[32].ids[1]=2; nodesIJKs[32].ids[2]=6;
	nodesIJKs[33].ids[0]=6; nodesIJKs[33].ids[1]=2; nodesIJKs[33].ids[2]=6;
	nodesIJKs[34].ids[0]=0; nodesIJKs[34].ids[1]=4; nodesIJKs[34].ids[2]=6;
	nodesIJKs[35].ids[0]=2; nodesIJKs[35].ids[1]=4; nodesIJKs[35].ids[2]=6;
	nodesIJKs[36].ids[0]=4; nodesIJKs[36].ids[1]=4; nodesIJKs[36].ids[2]=6;
	nodesIJKs[37].ids[0]=6; nodesIJKs[37].ids[1]=4; nodesIJKs[37].ids[2]=6;
	nodesIJKs[38].ids[0]=2; nodesIJKs[38].ids[1]=6; nodesIJKs[38].ids[2]=6;
	nodesIJKs[39].ids[0]=4; nodesIJKs[39].ids[1]=6; nodesIJKs[39].ids[2]=6;

	nodesOnRegionToBuild_[keyRegion] = nodesIJKs;

	keyRegion = tableMarkedNodes<8> (1,1,1,1,1,1,1,1); // region
	nodesIJKs.clear(); nodesIJKs.resize(nbNodesOnRegionToBuild_[keyRegion]);
	nodesIJKs[ 0].ids[0]=2; nodesIJKs[ 0].ids[1]=0; nodesIJKs[ 0].ids[2]=0;
	nodesIJKs[ 1].ids[0]=4; nodesIJKs[ 1].ids[1]=0; nodesIJKs[ 1].ids[2]=0;
	nodesIJKs[ 2].ids[0]=0; nodesIJKs[ 2].ids[1]=2; nodesIJKs[ 2].ids[2]=0;
	nodesIJKs[ 3].ids[0]=2; nodesIJKs[ 3].ids[1]=2; nodesIJKs[ 3].ids[2]=0;
	nodesIJKs[ 4].ids[0]=4; nodesIJKs[ 4].ids[1]=2; nodesIJKs[ 4].ids[2]=0;
	nodesIJKs[ 5].ids[0]=6; nodesIJKs[ 5].ids[1]=2; nodesIJKs[ 5].ids[2]=0;
	nodesIJKs[ 6].ids[0]=0; nodesIJKs[ 6].ids[1]=4; nodesIJKs[ 6].ids[2]=0;
	nodesIJKs[ 7].ids[0]=2; nodesIJKs[ 7].ids[1]=4; nodesIJKs[ 7].ids[2]=0;
	nodesIJKs[ 8].ids[0]=4; nodesIJKs[ 8].ids[1]=4; nodesIJKs[ 8].ids[2]=0;
	nodesIJKs[ 9].ids[0]=6; nodesIJKs[ 9].ids[1]=4; nodesIJKs[ 9].ids[2]=0;
	nodesIJKs[10].ids[0]=2; nodesIJKs[10].ids[1]=6; nodesIJKs[10].ids[2]=0;
	nodesIJKs[11].ids[0]=4; nodesIJKs[11].ids[1]=6; nodesIJKs[11].ids[2]=0;

	nodesIJKs[12].ids[0]=0; nodesIJKs[12].ids[1]=0; nodesIJKs[12].ids[2]=2;
	nodesIJKs[13].ids[0]=2; nodesIJKs[13].ids[1]=0; nodesIJKs[13].ids[2]=2;
	nodesIJKs[14].ids[0]=4; nodesIJKs[14].ids[1]=0; nodesIJKs[14].ids[2]=2;
	nodesIJKs[15].ids[0]=6; nodesIJKs[15].ids[1]=0; nodesIJKs[15].ids[2]=2;
	nodesIJKs[16].ids[0]=0; nodesIJKs[16].ids[1]=2; nodesIJKs[16].ids[2]=2;
	nodesIJKs[17].ids[0]=2; nodesIJKs[17].ids[1]=2; nodesIJKs[17].ids[2]=2;
	nodesIJKs[18].ids[0]=4; nodesIJKs[18].ids[1]=2; nodesIJKs[18].ids[2]=2;
	nodesIJKs[19].ids[0]=6; nodesIJKs[19].ids[1]=2; nodesIJKs[19].ids[2]=2;
	nodesIJKs[20].ids[0]=0; nodesIJKs[20].ids[1]=4; nodesIJKs[20].ids[2]=2;
	nodesIJKs[21].ids[0]=2; nodesIJKs[21].ids[1]=4; nodesIJKs[21].ids[2]=2;
	nodesIJKs[22].ids[0]=4; nodesIJKs[22].ids[1]=4; nodesIJKs[22].ids[2]=2;
	nodesIJKs[23].ids[0]=6; nodesIJKs[23].ids[1]=4; nodesIJKs[23].ids[2]=2;
	nodesIJKs[24].ids[0]=0; nodesIJKs[24].ids[1]=6; nodesIJKs[24].ids[2]=2;
	nodesIJKs[25].ids[0]=2; nodesIJKs[25].ids[1]=6; nodesIJKs[25].ids[2]=2;
	nodesIJKs[26].ids[0]=4; nodesIJKs[26].ids[1]=6; nodesIJKs[26].ids[2]=2;
	nodesIJKs[27].ids[0]=6; nodesIJKs[27].ids[1]=6; nodesIJKs[27].ids[2]=2;

	nodesIJKs[28].ids[0]=0; nodesIJKs[28].ids[1]=0; nodesIJKs[28].ids[2]=4;
	nodesIJKs[29].ids[0]=2; nodesIJKs[29].ids[1]=0; nodesIJKs[29].ids[2]=4;
	nodesIJKs[30].ids[0]=4; nodesIJKs[30].ids[1]=0; nodesIJKs[30].ids[2]=4;
	nodesIJKs[31].ids[0]=6; nodesIJKs[31].ids[1]=0; nodesIJKs[31].ids[2]=4;
	nodesIJKs[32].ids[0]=0; nodesIJKs[32].ids[1]=2; nodesIJKs[32].ids[2]=4;
	nodesIJKs[33].ids[0]=2; nodesIJKs[33].ids[1]=2; nodesIJKs[33].ids[2]=4;
	nodesIJKs[34].ids[0]=4; nodesIJKs[34].ids[1]=2; nodesIJKs[34].ids[2]=4;
	nodesIJKs[35].ids[0]=6; nodesIJKs[35].ids[1]=2; nodesIJKs[35].ids[2]=4;
	nodesIJKs[36].ids[0]=0; nodesIJKs[36].ids[1]=4; nodesIJKs[36].ids[2]=4;
	nodesIJKs[37].ids[0]=2; nodesIJKs[37].ids[1]=4; nodesIJKs[37].ids[2]=4;
	nodesIJKs[38].ids[0]=4; nodesIJKs[38].ids[1]=4; nodesIJKs[38].ids[2]=4;
	nodesIJKs[39].ids[0]=6; nodesIJKs[39].ids[1]=4; nodesIJKs[39].ids[2]=4;
	nodesIJKs[40].ids[0]=0; nodesIJKs[40].ids[1]=6; nodesIJKs[40].ids[2]=4;
	nodesIJKs[41].ids[0]=2; nodesIJKs[41].ids[1]=6; nodesIJKs[41].ids[2]=4;
	nodesIJKs[42].ids[0]=4; nodesIJKs[42].ids[1]=6; nodesIJKs[42].ids[2]=4;
	nodesIJKs[43].ids[0]=6; nodesIJKs[43].ids[1]=6; nodesIJKs[43].ids[2]=4;

	nodesIJKs[44].ids[0]=2; nodesIJKs[44].ids[1]=0; nodesIJKs[44].ids[2]=6;
	nodesIJKs[45].ids[0]=4; nodesIJKs[45].ids[1]=0; nodesIJKs[45].ids[2]=6;
	nodesIJKs[46].ids[0]=0; nodesIJKs[46].ids[1]=2; nodesIJKs[46].ids[2]=6;
	nodesIJKs[47].ids[0]=2; nodesIJKs[47].ids[1]=2; nodesIJKs[47].ids[2]=6;
	nodesIJKs[48].ids[0]=4; nodesIJKs[48].ids[1]=2; nodesIJKs[48].ids[2]=6;
	nodesIJKs[49].ids[0]=6; nodesIJKs[49].ids[1]=2; nodesIJKs[49].ids[2]=6;
	nodesIJKs[50].ids[0]=0; nodesIJKs[50].ids[1]=4; nodesIJKs[50].ids[2]=6;
	nodesIJKs[51].ids[0]=2; nodesIJKs[51].ids[1]=4; nodesIJKs[51].ids[2]=6;
	nodesIJKs[52].ids[0]=4; nodesIJKs[52].ids[1]=4; nodesIJKs[52].ids[2]=6;
	nodesIJKs[53].ids[0]=6; nodesIJKs[53].ids[1]=4; nodesIJKs[53].ids[2]=6;
	nodesIJKs[54].ids[0]=2; nodesIJKs[54].ids[1]=6; nodesIJKs[54].ids[2]=6;
	nodesIJKs[55].ids[0]=4; nodesIJKs[55].ids[1]=6; nodesIJKs[55].ids[2]=6;

	nodesOnRegionToBuild_[keyRegion] = nodesIJKs;

//	nodesIJKs[ 0].ids[0]=0; nodesIJKs[ 0].ids[1]=0; nodesIJKs[ 0].ids[2]=0;
//	nodesIJKs[ 1].ids[0]=0; nodesIJKs[ 1].ids[1]=0; nodesIJKs[ 1].ids[2]=0;
//	nodesIJKs[ 2].ids[0]=0; nodesIJKs[ 2].ids[1]=0; nodesIJKs[ 2].ids[2]=0;
//	nodesIJKs[ 3].ids[0]=0; nodesIJKs[ 3].ids[1]=0; nodesIJKs[ 3].ids[2]=0;
//	nodesIJKs[ 4].ids[0]=0; nodesIJKs[ 4].ids[1]=0; nodesIJKs[ 4].ids[2]=0;
//	nodesIJKs[ 5].ids[0]=0; nodesIJKs[ 5].ids[1]=0; nodesIJKs[ 5].ids[2]=0;
//	nodesIJKs[ 6].ids[0]=0; nodesIJKs[ 6].ids[1]=0; nodesIJKs[ 6].ids[2]=0;
//	nodesIJKs[ 7].ids[0]=0; nodesIJKs[ 7].ids[1]=0; nodesIJKs[ 7].ids[2]=0;
//	nodesIJKs[ 8].ids[0]=0; nodesIJKs[ 8].ids[1]=0; nodesIJKs[ 8].ids[2]=0;
//	nodesIJKs[ 9].ids[0]=0; nodesIJKs[ 9].ids[1]=0; nodesIJKs[ 9].ids[2]=0;


	std::vector<HexIJK> hexIJKs;
	struct HexIJK hexIJK;

	keyRegion = tableMarkedNodes<8> (0,0,0,0,0,0,0,0); // nothing
	hexIJKs.clear(); //hexIJKs.resize(nbRegionsOnRegionToBuild_[keyRegion]);
	regionsOnRegionToBuild_[keyRegion] = hexIJKs;

	keyRegion = tableMarkedNodes<8> (1,0,0,0,0,0,0,0); // corner
	hexIJKs.clear(); //hexIJKs.resize(nbRegionsOnRegionToBuild_[keyRegion]);

	hexIJK.ijk[0].ids[0] = 0; hexIJK.ijk[0].ids[1] = 0; hexIJK.ijk[0].ids[2] = 4;
	hexIJK.ijk[1].ids[0] = 0; hexIJK.ijk[1].ids[1] = 2; hexIJK.ijk[1].ids[2] = 4;
	hexIJK.ijk[2].ids[0] = 2; hexIJK.ijk[2].ids[1] = 2; hexIJK.ijk[2].ids[2] = 4;
	hexIJK.ijk[3].ids[0] = 2; hexIJK.ijk[3].ids[1] = 0; hexIJK.ijk[3].ids[2] = 4;
	hexIJK.ijk[4].ids[0] = 0; hexIJK.ijk[4].ids[1] = 0; hexIJK.ijk[4].ids[2] = 0;
	hexIJK.ijk[5].ids[0] = 0; hexIJK.ijk[5].ids[1] = 6; hexIJK.ijk[5].ids[2] = 0;
	hexIJK.ijk[6].ids[0] = 6; hexIJK.ijk[6].ids[1] = 6; hexIJK.ijk[6].ids[2] = 0;
	hexIJK.ijk[7].ids[0] = 6; hexIJK.ijk[7].ids[1] = 0; hexIJK.ijk[7].ids[2] = 0; hexIJKs.push_back(hexIJK);
	hexIJK.ijk[0].ids[0] = 0; hexIJK.ijk[0].ids[1] = 0; hexIJK.ijk[0].ids[2] = 6;
	hexIJK.ijk[1].ids[0] = 0; hexIJK.ijk[1].ids[1] = 2; hexIJK.ijk[1].ids[2] = 6;
	hexIJK.ijk[2].ids[0] = 2; hexIJK.ijk[2].ids[1] = 2; hexIJK.ijk[2].ids[2] = 6;
	hexIJK.ijk[3].ids[0] = 2; hexIJK.ijk[3].ids[1] = 0; hexIJK.ijk[3].ids[2] = 6;
	hexIJK.ijk[4].ids[0] = 0; hexIJK.ijk[4].ids[1] = 0; hexIJK.ijk[4].ids[2] = 4;
	hexIJK.ijk[5].ids[0] = 0; hexIJK.ijk[5].ids[1] = 2; hexIJK.ijk[5].ids[2] = 4;
	hexIJK.ijk[6].ids[0] = 2; hexIJK.ijk[6].ids[1] = 2; hexIJK.ijk[6].ids[2] = 4;
	hexIJK.ijk[7].ids[0] = 2; hexIJK.ijk[7].ids[1] = 0; hexIJK.ijk[7].ids[2] = 4; hexIJKs.push_back(hexIJK);
	hexIJK.ijk[0].ids[0] = 0; hexIJK.ijk[0].ids[1] = 2; hexIJK.ijk[0].ids[2] = 6;
	hexIJK.ijk[1].ids[0] = 0; hexIJK.ijk[1].ids[1] = 6; hexIJK.ijk[1].ids[2] = 6;
	hexIJK.ijk[2].ids[0] = 6; hexIJK.ijk[2].ids[1] = 6; hexIJK.ijk[2].ids[2] = 6;
	hexIJK.ijk[3].ids[0] = 2; hexIJK.ijk[3].ids[1] = 2; hexIJK.ijk[3].ids[2] = 6;
	hexIJK.ijk[4].ids[0] = 0; hexIJK.ijk[4].ids[1] = 2; hexIJK.ijk[4].ids[2] = 4;
	hexIJK.ijk[5].ids[0] = 0; hexIJK.ijk[5].ids[1] = 6; hexIJK.ijk[5].ids[2] = 0;
	hexIJK.ijk[6].ids[0] = 6; hexIJK.ijk[6].ids[1] = 6; hexIJK.ijk[6].ids[2] = 0;
	hexIJK.ijk[7].ids[0] = 2; hexIJK.ijk[7].ids[1] = 2; hexIJK.ijk[7].ids[2] = 4; hexIJKs.push_back(hexIJK);
	hexIJK.ijk[0].ids[0] = 2; hexIJK.ijk[0].ids[1] = 0; hexIJK.ijk[0].ids[2] = 6;
	hexIJK.ijk[1].ids[0] = 2; hexIJK.ijk[1].ids[1] = 2; hexIJK.ijk[1].ids[2] = 6;
	hexIJK.ijk[2].ids[0] = 6; hexIJK.ijk[2].ids[1] = 6; hexIJK.ijk[2].ids[2] = 6;
	hexIJK.ijk[3].ids[0] = 6; hexIJK.ijk[3].ids[1] = 0; hexIJK.ijk[3].ids[2] = 6;
	hexIJK.ijk[4].ids[0] = 2; hexIJK.ijk[4].ids[1] = 0; hexIJK.ijk[4].ids[2] = 4;
	hexIJK.ijk[5].ids[0] = 2; hexIJK.ijk[5].ids[1] = 2; hexIJK.ijk[5].ids[2] = 4;
	hexIJK.ijk[6].ids[0] = 6; hexIJK.ijk[6].ids[1] = 6; hexIJK.ijk[6].ids[2] = 0;
	hexIJK.ijk[7].ids[0] = 6; hexIJK.ijk[7].ids[1] = 0; hexIJK.ijk[7].ids[2] = 0; hexIJKs.push_back(hexIJK);

	regionsOnRegionToBuild_[keyRegion] = hexIJKs;

	keyRegion = tableMarkedNodes<8> (1,1,0,0,0,0,0,0); // edge
	hexIJKs.clear(); //hexIJKs.resize(nbRegionsOnRegionToBuild_[keyRegion]);

	hexIJK.ijk[0].ids[0] = 0; hexIJK.ijk[0].ids[1] = 0; hexIJK.ijk[0].ids[2] = 4;
	hexIJK.ijk[1].ids[0] = 0; hexIJK.ijk[1].ids[1] = 2; hexIJK.ijk[1].ids[2] = 4;
	hexIJK.ijk[2].ids[0] = 2; hexIJK.ijk[2].ids[1] = 2; hexIJK.ijk[2].ids[2] = 4;
	hexIJK.ijk[3].ids[0] = 2; hexIJK.ijk[3].ids[1] = 0; hexIJK.ijk[3].ids[2] = 4;
	hexIJK.ijk[4].ids[0] = 0; hexIJK.ijk[4].ids[1] = 0; hexIJK.ijk[4].ids[2] = 0;
	hexIJK.ijk[5].ids[0] = 0; hexIJK.ijk[5].ids[1] = 2; hexIJK.ijk[5].ids[2] = 2;
	hexIJK.ijk[6].ids[0] = 4; hexIJK.ijk[6].ids[1] = 2; hexIJK.ijk[6].ids[2] = 2;
	hexIJK.ijk[7].ids[0] = 6; hexIJK.ijk[7].ids[1] = 0; hexIJK.ijk[7].ids[2] = 0; hexIJKs.push_back(hexIJK); //  1
	hexIJK.ijk[0].ids[0] = 0; hexIJK.ijk[0].ids[1] = 2; hexIJK.ijk[0].ids[2] = 2;
	hexIJK.ijk[1].ids[0] = 0; hexIJK.ijk[1].ids[1] = 4; hexIJK.ijk[1].ids[2] = 2;
	hexIJK.ijk[2].ids[0] = 4; hexIJK.ijk[2].ids[1] = 4; hexIJK.ijk[2].ids[2] = 2;
	hexIJK.ijk[3].ids[0] = 4; hexIJK.ijk[3].ids[1] = 2; hexIJK.ijk[3].ids[2] = 2;
	hexIJK.ijk[4].ids[0] = 0; hexIJK.ijk[4].ids[1] = 0; hexIJK.ijk[4].ids[2] = 0;
	hexIJK.ijk[5].ids[0] = 0; hexIJK.ijk[5].ids[1] = 6; hexIJK.ijk[5].ids[2] = 0;
	hexIJK.ijk[6].ids[0] = 6; hexIJK.ijk[6].ids[1] = 6; hexIJK.ijk[6].ids[2] = 0;
	hexIJK.ijk[7].ids[0] = 6; hexIJK.ijk[7].ids[1] = 0; hexIJK.ijk[7].ids[2] = 0; hexIJKs.push_back(hexIJK); //  2
	hexIJK.ijk[0].ids[0] = 0; hexIJK.ijk[0].ids[1] = 2; hexIJK.ijk[0].ids[2] = 4;
	hexIJK.ijk[1].ids[0] = 0; hexIJK.ijk[1].ids[1] = 4; hexIJK.ijk[1].ids[2] = 4;
	hexIJK.ijk[2].ids[0] = 2; hexIJK.ijk[2].ids[1] = 4; hexIJK.ijk[2].ids[2] = 4;
	hexIJK.ijk[3].ids[0] = 2; hexIJK.ijk[3].ids[1] = 2; hexIJK.ijk[3].ids[2] = 4;
	hexIJK.ijk[4].ids[0] = 0; hexIJK.ijk[4].ids[1] = 2; hexIJK.ijk[4].ids[2] = 2;
	hexIJK.ijk[5].ids[0] = 0; hexIJK.ijk[5].ids[1] = 4; hexIJK.ijk[5].ids[2] = 2;
	hexIJK.ijk[6].ids[0] = 4; hexIJK.ijk[6].ids[1] = 4; hexIJK.ijk[6].ids[2] = 2;
	hexIJK.ijk[7].ids[0] = 4; hexIJK.ijk[7].ids[1] = 2; hexIJK.ijk[7].ids[2] = 2; hexIJKs.push_back(hexIJK); //  3
	hexIJK.ijk[0].ids[0] = 0; hexIJK.ijk[0].ids[1] = 4; hexIJK.ijk[0].ids[2] = 4;
	hexIJK.ijk[1].ids[0] = 0; hexIJK.ijk[1].ids[1] = 6; hexIJK.ijk[1].ids[2] = 4;
	hexIJK.ijk[2].ids[0] = 2; hexIJK.ijk[2].ids[1] = 6; hexIJK.ijk[2].ids[2] = 4;
	hexIJK.ijk[3].ids[0] = 2; hexIJK.ijk[3].ids[1] = 4; hexIJK.ijk[3].ids[2] = 4;
	hexIJK.ijk[4].ids[0] = 0; hexIJK.ijk[4].ids[1] = 4; hexIJK.ijk[4].ids[2] = 2;
	hexIJK.ijk[5].ids[0] = 0; hexIJK.ijk[5].ids[1] = 6; hexIJK.ijk[5].ids[2] = 0;
	hexIJK.ijk[6].ids[0] = 6; hexIJK.ijk[6].ids[1] = 6; hexIJK.ijk[6].ids[2] = 0;
	hexIJK.ijk[7].ids[0] = 4; hexIJK.ijk[7].ids[1] = 4; hexIJK.ijk[7].ids[2] = 2; hexIJKs.push_back(hexIJK); //  4
	hexIJK.ijk[0].ids[0] = 0; hexIJK.ijk[0].ids[1] = 0; hexIJK.ijk[0].ids[2] = 6;
	hexIJK.ijk[1].ids[0] = 0; hexIJK.ijk[1].ids[1] = 2; hexIJK.ijk[1].ids[2] = 6;
	hexIJK.ijk[2].ids[0] = 2; hexIJK.ijk[2].ids[1] = 2; hexIJK.ijk[2].ids[2] = 6;
	hexIJK.ijk[3].ids[0] = 2; hexIJK.ijk[3].ids[1] = 0; hexIJK.ijk[3].ids[2] = 6;
	hexIJK.ijk[4].ids[0] = 0; hexIJK.ijk[4].ids[1] = 0; hexIJK.ijk[4].ids[2] = 4;
	hexIJK.ijk[5].ids[0] = 0; hexIJK.ijk[5].ids[1] = 2; hexIJK.ijk[5].ids[2] = 4;
	hexIJK.ijk[6].ids[0] = 2; hexIJK.ijk[6].ids[1] = 2; hexIJK.ijk[6].ids[2] = 4;
	hexIJK.ijk[7].ids[0] = 2; hexIJK.ijk[7].ids[1] = 0; hexIJK.ijk[7].ids[2] = 4; hexIJKs.push_back(hexIJK); //  5
	hexIJK.ijk[0].ids[0] = 2; hexIJK.ijk[0].ids[1] = 0; hexIJK.ijk[0].ids[2] = 6;
	hexIJK.ijk[1].ids[0] = 2; hexIJK.ijk[1].ids[1] = 2; hexIJK.ijk[1].ids[2] = 6;
	hexIJK.ijk[2].ids[0] = 4; hexIJK.ijk[2].ids[1] = 2; hexIJK.ijk[2].ids[2] = 6;
	hexIJK.ijk[3].ids[0] = 6; hexIJK.ijk[3].ids[1] = 0; hexIJK.ijk[3].ids[2] = 6;
	hexIJK.ijk[4].ids[0] = 2; hexIJK.ijk[4].ids[1] = 0; hexIJK.ijk[4].ids[2] = 4;
	hexIJK.ijk[5].ids[0] = 2; hexIJK.ijk[5].ids[1] = 2; hexIJK.ijk[5].ids[2] = 4;
	hexIJK.ijk[6].ids[0] = 4; hexIJK.ijk[6].ids[1] = 2; hexIJK.ijk[6].ids[2] = 2;
	hexIJK.ijk[7].ids[0] = 6; hexIJK.ijk[7].ids[1] = 0; hexIJK.ijk[7].ids[2] = 0; hexIJKs.push_back(hexIJK); //  6
	hexIJK.ijk[0].ids[0] = 0; hexIJK.ijk[0].ids[1] = 2; hexIJK.ijk[0].ids[2] = 6;
	hexIJK.ijk[1].ids[0] = 0; hexIJK.ijk[1].ids[1] = 4; hexIJK.ijk[1].ids[2] = 6;
	hexIJK.ijk[2].ids[0] = 2; hexIJK.ijk[2].ids[1] = 4; hexIJK.ijk[2].ids[2] = 6;
	hexIJK.ijk[3].ids[0] = 2; hexIJK.ijk[3].ids[1] = 2; hexIJK.ijk[3].ids[2] = 6;
	hexIJK.ijk[4].ids[0] = 0; hexIJK.ijk[4].ids[1] = 2; hexIJK.ijk[4].ids[2] = 4;
	hexIJK.ijk[5].ids[0] = 0; hexIJK.ijk[5].ids[1] = 4; hexIJK.ijk[5].ids[2] = 4;
	hexIJK.ijk[6].ids[0] = 2; hexIJK.ijk[6].ids[1] = 4; hexIJK.ijk[6].ids[2] = 4;
	hexIJK.ijk[7].ids[0] = 2; hexIJK.ijk[7].ids[1] = 2; hexIJK.ijk[7].ids[2] = 4; hexIJKs.push_back(hexIJK); //  7
	hexIJK.ijk[0].ids[0] = 2; hexIJK.ijk[0].ids[1] = 2; hexIJK.ijk[0].ids[2] = 6;
	hexIJK.ijk[1].ids[0] = 2; hexIJK.ijk[1].ids[1] = 4; hexIJK.ijk[1].ids[2] = 6;
	hexIJK.ijk[2].ids[0] = 4; hexIJK.ijk[2].ids[1] = 4; hexIJK.ijk[2].ids[2] = 6;
	hexIJK.ijk[3].ids[0] = 4; hexIJK.ijk[3].ids[1] = 2; hexIJK.ijk[3].ids[2] = 6;
	hexIJK.ijk[4].ids[0] = 2; hexIJK.ijk[4].ids[1] = 2; hexIJK.ijk[4].ids[2] = 4;
	hexIJK.ijk[5].ids[0] = 2; hexIJK.ijk[5].ids[1] = 4; hexIJK.ijk[5].ids[2] = 4;
	hexIJK.ijk[6].ids[0] = 4; hexIJK.ijk[6].ids[1] = 4; hexIJK.ijk[6].ids[2] = 2;
	hexIJK.ijk[7].ids[0] = 4; hexIJK.ijk[7].ids[1] = 2; hexIJK.ijk[7].ids[2] = 2; hexIJKs.push_back(hexIJK); //  8
	hexIJK.ijk[0].ids[0] = 4; hexIJK.ijk[0].ids[1] = 2; hexIJK.ijk[0].ids[2] = 6;
	hexIJK.ijk[1].ids[0] = 4; hexIJK.ijk[1].ids[1] = 4; hexIJK.ijk[1].ids[2] = 6;
	hexIJK.ijk[2].ids[0] = 6; hexIJK.ijk[2].ids[1] = 6; hexIJK.ijk[2].ids[2] = 6;
	hexIJK.ijk[3].ids[0] = 6; hexIJK.ijk[3].ids[1] = 0; hexIJK.ijk[3].ids[2] = 6;
	hexIJK.ijk[4].ids[0] = 4; hexIJK.ijk[4].ids[1] = 2; hexIJK.ijk[4].ids[2] = 2;
	hexIJK.ijk[5].ids[0] = 4; hexIJK.ijk[5].ids[1] = 4; hexIJK.ijk[5].ids[2] = 2;
	hexIJK.ijk[6].ids[0] = 6; hexIJK.ijk[6].ids[1] = 6; hexIJK.ijk[6].ids[2] = 0;
	hexIJK.ijk[7].ids[0] = 6; hexIJK.ijk[7].ids[1] = 0; hexIJK.ijk[7].ids[2] = 0; hexIJKs.push_back(hexIJK); //  9
	hexIJK.ijk[0].ids[0] = 0; hexIJK.ijk[0].ids[1] = 4; hexIJK.ijk[0].ids[2] = 6;
	hexIJK.ijk[1].ids[0] = 0; hexIJK.ijk[1].ids[1] = 6; hexIJK.ijk[1].ids[2] = 6;
	hexIJK.ijk[2].ids[0] = 2; hexIJK.ijk[2].ids[1] = 6; hexIJK.ijk[2].ids[2] = 6;
	hexIJK.ijk[3].ids[0] = 2; hexIJK.ijk[3].ids[1] = 4; hexIJK.ijk[3].ids[2] = 6;
	hexIJK.ijk[4].ids[0] = 0; hexIJK.ijk[4].ids[1] = 4; hexIJK.ijk[4].ids[2] = 4;
	hexIJK.ijk[5].ids[0] = 0; hexIJK.ijk[5].ids[1] = 6; hexIJK.ijk[5].ids[2] = 4;
	hexIJK.ijk[6].ids[0] = 2; hexIJK.ijk[6].ids[1] = 6; hexIJK.ijk[6].ids[2] = 4;
	hexIJK.ijk[7].ids[0] = 2; hexIJK.ijk[7].ids[1] = 4; hexIJK.ijk[7].ids[2] = 4; hexIJKs.push_back(hexIJK); // 10
	hexIJK.ijk[0].ids[0] = 2; hexIJK.ijk[0].ids[1] = 4; hexIJK.ijk[0].ids[2] = 6;
	hexIJK.ijk[1].ids[0] = 2; hexIJK.ijk[1].ids[1] = 6; hexIJK.ijk[1].ids[2] = 6;
	hexIJK.ijk[2].ids[0] = 6; hexIJK.ijk[2].ids[1] = 6; hexIJK.ijk[2].ids[2] = 6;
	hexIJK.ijk[3].ids[0] = 4; hexIJK.ijk[3].ids[1] = 4; hexIJK.ijk[3].ids[2] = 6;
	hexIJK.ijk[4].ids[0] = 2; hexIJK.ijk[4].ids[1] = 4; hexIJK.ijk[4].ids[2] = 4;
	hexIJK.ijk[5].ids[0] = 2; hexIJK.ijk[5].ids[1] = 6; hexIJK.ijk[5].ids[2] = 4;
	hexIJK.ijk[6].ids[0] = 6; hexIJK.ijk[6].ids[1] = 6; hexIJK.ijk[6].ids[2] = 0;
	hexIJK.ijk[7].ids[0] = 4; hexIJK.ijk[7].ids[1] = 4; hexIJK.ijk[7].ids[2] = 2; hexIJKs.push_back(hexIJK); // 11

	regionsOnRegionToBuild_[keyRegion] = hexIJKs;

	keyRegion = tableMarkedNodes<8> (1,1,1,1,0,0,0,0); // face
	hexIJKs.clear(); //hexIJKs.resize(nbRegionsOnRegionToBuild_[keyRegion]);

	hexIJK.ijk[0].ids[0] = 2; hexIJK.ijk[0].ids[1] = 0; hexIJK.ijk[0].ids[2] = 2;
	hexIJK.ijk[1].ids[0] = 2; hexIJK.ijk[1].ids[1] = 2; hexIJK.ijk[1].ids[2] = 3;
	hexIJK.ijk[2].ids[0] = 4; hexIJK.ijk[2].ids[1] = 2; hexIJK.ijk[2].ids[2] = 3;
	hexIJK.ijk[3].ids[0] = 4; hexIJK.ijk[3].ids[1] = 0; hexIJK.ijk[3].ids[2] = 2;
	hexIJK.ijk[4].ids[0] = 0; hexIJK.ijk[4].ids[1] = 0; hexIJK.ijk[4].ids[2] = 0;
	hexIJK.ijk[5].ids[0] = 0; hexIJK.ijk[5].ids[1] = 2; hexIJK.ijk[5].ids[2] = 2;
	hexIJK.ijk[6].ids[0] = 6; hexIJK.ijk[6].ids[1] = 2; hexIJK.ijk[6].ids[2] = 2;
	hexIJK.ijk[7].ids[0] = 6; hexIJK.ijk[7].ids[1] = 0; hexIJK.ijk[7].ids[2] = 0; hexIJKs.push_back(hexIJK); //  1
	hexIJK.ijk[0].ids[0] = 0; hexIJK.ijk[0].ids[1] = 2; hexIJK.ijk[0].ids[2] = 2;
	hexIJK.ijk[1].ids[0] = 0; hexIJK.ijk[1].ids[1] = 4; hexIJK.ijk[1].ids[2] = 2;
	hexIJK.ijk[2].ids[0] = 6; hexIJK.ijk[2].ids[1] = 4; hexIJK.ijk[2].ids[2] = 2;
	hexIJK.ijk[3].ids[0] = 6; hexIJK.ijk[3].ids[1] = 2; hexIJK.ijk[3].ids[2] = 2;
	hexIJK.ijk[4].ids[0] = 0; hexIJK.ijk[4].ids[1] = 0; hexIJK.ijk[4].ids[2] = 0;
	hexIJK.ijk[5].ids[0] = 0; hexIJK.ijk[5].ids[1] = 6; hexIJK.ijk[5].ids[2] = 0;
	hexIJK.ijk[6].ids[0] = 6; hexIJK.ijk[6].ids[1] = 6; hexIJK.ijk[6].ids[2] = 0;
	hexIJK.ijk[7].ids[0] = 6; hexIJK.ijk[7].ids[1] = 0; hexIJK.ijk[7].ids[2] = 0; hexIJKs.push_back(hexIJK); //  2
	hexIJK.ijk[0].ids[0] = 2; hexIJK.ijk[0].ids[1] = 4; hexIJK.ijk[0].ids[2] = 3;
	hexIJK.ijk[1].ids[0] = 2; hexIJK.ijk[1].ids[1] = 6; hexIJK.ijk[1].ids[2] = 2;
	hexIJK.ijk[2].ids[0] = 4; hexIJK.ijk[2].ids[1] = 6; hexIJK.ijk[2].ids[2] = 2;
	hexIJK.ijk[3].ids[0] = 4; hexIJK.ijk[3].ids[1] = 4; hexIJK.ijk[3].ids[2] = 3;
	hexIJK.ijk[4].ids[0] = 0; hexIJK.ijk[4].ids[1] = 4; hexIJK.ijk[4].ids[2] = 2;
	hexIJK.ijk[5].ids[0] = 0; hexIJK.ijk[5].ids[1] = 6; hexIJK.ijk[5].ids[2] = 0;
	hexIJK.ijk[6].ids[0] = 6; hexIJK.ijk[6].ids[1] = 6; hexIJK.ijk[6].ids[2] = 0;
	hexIJK.ijk[7].ids[0] = 6; hexIJK.ijk[7].ids[1] = 4; hexIJK.ijk[7].ids[2] = 2; hexIJKs.push_back(hexIJK); //  3
	hexIJK.ijk[0].ids[0] = 0; hexIJK.ijk[0].ids[1] = 0; hexIJK.ijk[0].ids[2] = 4;
	hexIJK.ijk[1].ids[0] = 0; hexIJK.ijk[1].ids[1] = 2; hexIJK.ijk[1].ids[2] = 4;
	hexIJK.ijk[2].ids[0] = 2; hexIJK.ijk[2].ids[1] = 2; hexIJK.ijk[2].ids[2] = 4;
	hexIJK.ijk[3].ids[0] = 2; hexIJK.ijk[3].ids[1] = 0; hexIJK.ijk[3].ids[2] = 4;
	hexIJK.ijk[4].ids[0] = 0; hexIJK.ijk[4].ids[1] = 0; hexIJK.ijk[4].ids[2] = 0;
	hexIJK.ijk[5].ids[0] = 0; hexIJK.ijk[5].ids[1] = 2; hexIJK.ijk[5].ids[2] = 2;
	hexIJK.ijk[6].ids[0] = 2; hexIJK.ijk[6].ids[1] = 2; hexIJK.ijk[6].ids[2] = 3;
	hexIJK.ijk[7].ids[0] = 2; hexIJK.ijk[7].ids[1] = 0; hexIJK.ijk[7].ids[2] = 2; hexIJKs.push_back(hexIJK); //  4
	hexIJK.ijk[0].ids[0] = 2; hexIJK.ijk[0].ids[1] = 0; hexIJK.ijk[0].ids[2] = 4;
	hexIJK.ijk[1].ids[0] = 2; hexIJK.ijk[1].ids[1] = 2; hexIJK.ijk[1].ids[2] = 4;
	hexIJK.ijk[2].ids[0] = 4; hexIJK.ijk[2].ids[1] = 2; hexIJK.ijk[2].ids[2] = 4;
	hexIJK.ijk[3].ids[0] = 4; hexIJK.ijk[3].ids[1] = 0; hexIJK.ijk[3].ids[2] = 4;
	hexIJK.ijk[4].ids[0] = 2; hexIJK.ijk[4].ids[1] = 0; hexIJK.ijk[4].ids[2] = 2;
	hexIJK.ijk[5].ids[0] = 2; hexIJK.ijk[5].ids[1] = 2; hexIJK.ijk[5].ids[2] = 3;
	hexIJK.ijk[6].ids[0] = 4; hexIJK.ijk[6].ids[1] = 2; hexIJK.ijk[6].ids[2] = 3;
	hexIJK.ijk[7].ids[0] = 4; hexIJK.ijk[7].ids[1] = 0; hexIJK.ijk[7].ids[2] = 2; hexIJKs.push_back(hexIJK); //  5
	hexIJK.ijk[0].ids[0] = 4; hexIJK.ijk[0].ids[1] = 0; hexIJK.ijk[0].ids[2] = 4;
	hexIJK.ijk[1].ids[0] = 4; hexIJK.ijk[1].ids[1] = 2; hexIJK.ijk[1].ids[2] = 4;
	hexIJK.ijk[2].ids[0] = 6; hexIJK.ijk[2].ids[1] = 2; hexIJK.ijk[2].ids[2] = 4;
	hexIJK.ijk[3].ids[0] = 6; hexIJK.ijk[3].ids[1] = 0; hexIJK.ijk[3].ids[2] = 4;
	hexIJK.ijk[4].ids[0] = 4; hexIJK.ijk[4].ids[1] = 0; hexIJK.ijk[4].ids[2] = 2;
	hexIJK.ijk[5].ids[0] = 4; hexIJK.ijk[5].ids[1] = 2; hexIJK.ijk[5].ids[2] = 3;
	hexIJK.ijk[6].ids[0] = 6; hexIJK.ijk[6].ids[1] = 2; hexIJK.ijk[6].ids[2] = 2;
	hexIJK.ijk[7].ids[0] = 6; hexIJK.ijk[7].ids[1] = 0; hexIJK.ijk[7].ids[2] = 0; hexIJKs.push_back(hexIJK); //  6
	hexIJK.ijk[0].ids[0] = 0; hexIJK.ijk[0].ids[1] = 2; hexIJK.ijk[0].ids[2] = 4;
	hexIJK.ijk[1].ids[0] = 0; hexIJK.ijk[1].ids[1] = 4; hexIJK.ijk[1].ids[2] = 4;
	hexIJK.ijk[2].ids[0] = 2; hexIJK.ijk[2].ids[1] = 4; hexIJK.ijk[2].ids[2] = 4;
	hexIJK.ijk[3].ids[0] = 2; hexIJK.ijk[3].ids[1] = 2; hexIJK.ijk[3].ids[2] = 4;
	hexIJK.ijk[4].ids[0] = 0; hexIJK.ijk[4].ids[1] = 2; hexIJK.ijk[4].ids[2] = 2;
	hexIJK.ijk[5].ids[0] = 0; hexIJK.ijk[5].ids[1] = 4; hexIJK.ijk[5].ids[2] = 2;
	hexIJK.ijk[6].ids[0] = 2; hexIJK.ijk[6].ids[1] = 4; hexIJK.ijk[6].ids[2] = 3;
	hexIJK.ijk[7].ids[0] = 2; hexIJK.ijk[7].ids[1] = 2; hexIJK.ijk[7].ids[2] = 3; hexIJKs.push_back(hexIJK); //  7
	hexIJK.ijk[0].ids[0] = 2; hexIJK.ijk[0].ids[1] = 2; hexIJK.ijk[0].ids[2] = 3;
	hexIJK.ijk[1].ids[0] = 2; hexIJK.ijk[1].ids[1] = 4; hexIJK.ijk[1].ids[2] = 3;
	hexIJK.ijk[2].ids[0] = 4; hexIJK.ijk[2].ids[1] = 4; hexIJK.ijk[2].ids[2] = 3;
	hexIJK.ijk[3].ids[0] = 4; hexIJK.ijk[3].ids[1] = 2; hexIJK.ijk[3].ids[2] = 3;
	hexIJK.ijk[4].ids[0] = 0; hexIJK.ijk[4].ids[1] = 2; hexIJK.ijk[4].ids[2] = 2;
	hexIJK.ijk[5].ids[0] = 0; hexIJK.ijk[5].ids[1] = 4; hexIJK.ijk[5].ids[2] = 2;
	hexIJK.ijk[6].ids[0] = 6; hexIJK.ijk[6].ids[1] = 4; hexIJK.ijk[6].ids[2] = 2;
	hexIJK.ijk[7].ids[0] = 6; hexIJK.ijk[7].ids[1] = 2; hexIJK.ijk[7].ids[2] = 2; hexIJKs.push_back(hexIJK); //  8
	hexIJK.ijk[0].ids[0] = 2; hexIJK.ijk[0].ids[1] = 2; hexIJK.ijk[0].ids[2] = 4;
	hexIJK.ijk[1].ids[0] = 2; hexIJK.ijk[1].ids[1] = 4; hexIJK.ijk[1].ids[2] = 4;
	hexIJK.ijk[2].ids[0] = 4; hexIJK.ijk[2].ids[1] = 4; hexIJK.ijk[2].ids[2] = 4;
	hexIJK.ijk[3].ids[0] = 4; hexIJK.ijk[3].ids[1] = 2; hexIJK.ijk[3].ids[2] = 4;
	hexIJK.ijk[4].ids[0] = 2; hexIJK.ijk[4].ids[1] = 2; hexIJK.ijk[4].ids[2] = 3;
	hexIJK.ijk[5].ids[0] = 2; hexIJK.ijk[5].ids[1] = 4; hexIJK.ijk[5].ids[2] = 3;
	hexIJK.ijk[6].ids[0] = 4; hexIJK.ijk[6].ids[1] = 4; hexIJK.ijk[6].ids[2] = 3;
	hexIJK.ijk[7].ids[0] = 4; hexIJK.ijk[7].ids[1] = 2; hexIJK.ijk[7].ids[2] = 3; hexIJKs.push_back(hexIJK); //  9
	hexIJK.ijk[0].ids[0] = 4; hexIJK.ijk[0].ids[1] = 2; hexIJK.ijk[0].ids[2] = 4;
	hexIJK.ijk[1].ids[0] = 4; hexIJK.ijk[1].ids[1] = 4; hexIJK.ijk[1].ids[2] = 4;
	hexIJK.ijk[2].ids[0] = 6; hexIJK.ijk[2].ids[1] = 4; hexIJK.ijk[2].ids[2] = 4;
	hexIJK.ijk[3].ids[0] = 6; hexIJK.ijk[3].ids[1] = 2; hexIJK.ijk[3].ids[2] = 4;
	hexIJK.ijk[4].ids[0] = 4; hexIJK.ijk[4].ids[1] = 2; hexIJK.ijk[4].ids[2] = 3;
	hexIJK.ijk[5].ids[0] = 4; hexIJK.ijk[5].ids[1] = 4; hexIJK.ijk[5].ids[2] = 3;
	hexIJK.ijk[6].ids[0] = 6; hexIJK.ijk[6].ids[1] = 4; hexIJK.ijk[6].ids[2] = 2;
	hexIJK.ijk[7].ids[0] = 6; hexIJK.ijk[7].ids[1] = 2; hexIJK.ijk[7].ids[2] = 2; hexIJKs.push_back(hexIJK); // 10
	hexIJK.ijk[0].ids[0] = 0; hexIJK.ijk[0].ids[1] = 4; hexIJK.ijk[0].ids[2] = 4;
	hexIJK.ijk[1].ids[0] = 0; hexIJK.ijk[1].ids[1] = 6; hexIJK.ijk[1].ids[2] = 4;
	hexIJK.ijk[2].ids[0] = 2; hexIJK.ijk[2].ids[1] = 6; hexIJK.ijk[2].ids[2] = 4;
	hexIJK.ijk[3].ids[0] = 2; hexIJK.ijk[3].ids[1] = 4; hexIJK.ijk[3].ids[2] = 4;
	hexIJK.ijk[4].ids[0] = 0; hexIJK.ijk[4].ids[1] = 4; hexIJK.ijk[4].ids[2] = 2;
	hexIJK.ijk[5].ids[0] = 0; hexIJK.ijk[5].ids[1] = 6; hexIJK.ijk[5].ids[2] = 0;
	hexIJK.ijk[6].ids[0] = 2; hexIJK.ijk[6].ids[1] = 6; hexIJK.ijk[6].ids[2] = 2;
	hexIJK.ijk[7].ids[0] = 2; hexIJK.ijk[7].ids[1] = 4; hexIJK.ijk[7].ids[2] = 3; hexIJKs.push_back(hexIJK); // 11
	hexIJK.ijk[0].ids[0] = 2; hexIJK.ijk[0].ids[1] = 4; hexIJK.ijk[0].ids[2] = 4;
	hexIJK.ijk[1].ids[0] = 2; hexIJK.ijk[1].ids[1] = 6; hexIJK.ijk[1].ids[2] = 4;
	hexIJK.ijk[2].ids[0] = 4; hexIJK.ijk[2].ids[1] = 6; hexIJK.ijk[2].ids[2] = 4;
	hexIJK.ijk[3].ids[0] = 4; hexIJK.ijk[3].ids[1] = 4; hexIJK.ijk[3].ids[2] = 4;
	hexIJK.ijk[4].ids[0] = 2; hexIJK.ijk[4].ids[1] = 4; hexIJK.ijk[4].ids[2] = 3;
	hexIJK.ijk[5].ids[0] = 2; hexIJK.ijk[5].ids[1] = 6; hexIJK.ijk[5].ids[2] = 2;
	hexIJK.ijk[6].ids[0] = 4; hexIJK.ijk[6].ids[1] = 6; hexIJK.ijk[6].ids[2] = 2;
	hexIJK.ijk[7].ids[0] = 4; hexIJK.ijk[7].ids[1] = 4; hexIJK.ijk[7].ids[2] = 3; hexIJKs.push_back(hexIJK); // 12
	hexIJK.ijk[0].ids[0] = 4; hexIJK.ijk[0].ids[1] = 4; hexIJK.ijk[0].ids[2] = 4;
	hexIJK.ijk[1].ids[0] = 4; hexIJK.ijk[1].ids[1] = 6; hexIJK.ijk[1].ids[2] = 4;
	hexIJK.ijk[2].ids[0] = 6; hexIJK.ijk[2].ids[1] = 6; hexIJK.ijk[2].ids[2] = 4;
	hexIJK.ijk[3].ids[0] = 6; hexIJK.ijk[3].ids[1] = 4; hexIJK.ijk[3].ids[2] = 4;
	hexIJK.ijk[4].ids[0] = 4; hexIJK.ijk[4].ids[1] = 4; hexIJK.ijk[4].ids[2] = 3;
	hexIJK.ijk[5].ids[0] = 4; hexIJK.ijk[5].ids[1] = 6; hexIJK.ijk[5].ids[2] = 2;
	hexIJK.ijk[6].ids[0] = 6; hexIJK.ijk[6].ids[1] = 6; hexIJK.ijk[6].ids[2] = 0;
	hexIJK.ijk[7].ids[0] = 6; hexIJK.ijk[7].ids[1] = 4; hexIJK.ijk[7].ids[2] = 2; hexIJKs.push_back(hexIJK); // 13
	hexIJK.ijk[0].ids[0] = 0; hexIJK.ijk[0].ids[1] = 0; hexIJK.ijk[0].ids[2] = 6;
	hexIJK.ijk[1].ids[0] = 0; hexIJK.ijk[1].ids[1] = 2; hexIJK.ijk[1].ids[2] = 6;
	hexIJK.ijk[2].ids[0] = 2; hexIJK.ijk[2].ids[1] = 2; hexIJK.ijk[2].ids[2] = 6;
	hexIJK.ijk[3].ids[0] = 2; hexIJK.ijk[3].ids[1] = 0; hexIJK.ijk[3].ids[2] = 6;
	hexIJK.ijk[4].ids[0] = 0; hexIJK.ijk[4].ids[1] = 0; hexIJK.ijk[4].ids[2] = 4;
	hexIJK.ijk[5].ids[0] = 0; hexIJK.ijk[5].ids[1] = 2; hexIJK.ijk[5].ids[2] = 4;
	hexIJK.ijk[6].ids[0] = 2; hexIJK.ijk[6].ids[1] = 2; hexIJK.ijk[6].ids[2] = 4;
	hexIJK.ijk[7].ids[0] = 2; hexIJK.ijk[7].ids[1] = 0; hexIJK.ijk[7].ids[2] = 4; hexIJKs.push_back(hexIJK); // 14
	hexIJK.ijk[0].ids[0] = 2; hexIJK.ijk[0].ids[1] = 0; hexIJK.ijk[0].ids[2] = 6;
	hexIJK.ijk[1].ids[0] = 2; hexIJK.ijk[1].ids[1] = 2; hexIJK.ijk[1].ids[2] = 6;
	hexIJK.ijk[2].ids[0] = 4; hexIJK.ijk[2].ids[1] = 2; hexIJK.ijk[2].ids[2] = 6;
	hexIJK.ijk[3].ids[0] = 4; hexIJK.ijk[3].ids[1] = 0; hexIJK.ijk[3].ids[2] = 6;
	hexIJK.ijk[4].ids[0] = 2; hexIJK.ijk[4].ids[1] = 0; hexIJK.ijk[4].ids[2] = 4;
	hexIJK.ijk[5].ids[0] = 2; hexIJK.ijk[5].ids[1] = 2; hexIJK.ijk[5].ids[2] = 4;
	hexIJK.ijk[6].ids[0] = 4; hexIJK.ijk[6].ids[1] = 2; hexIJK.ijk[6].ids[2] = 4;
	hexIJK.ijk[7].ids[0] = 4; hexIJK.ijk[7].ids[1] = 0; hexIJK.ijk[7].ids[2] = 4; hexIJKs.push_back(hexIJK); // 15
	hexIJK.ijk[0].ids[0] = 4; hexIJK.ijk[0].ids[1] = 0; hexIJK.ijk[0].ids[2] = 6;
	hexIJK.ijk[1].ids[0] = 4; hexIJK.ijk[1].ids[1] = 2; hexIJK.ijk[1].ids[2] = 6;
	hexIJK.ijk[2].ids[0] = 6; hexIJK.ijk[2].ids[1] = 2; hexIJK.ijk[2].ids[2] = 6;
	hexIJK.ijk[3].ids[0] = 6; hexIJK.ijk[3].ids[1] = 0; hexIJK.ijk[3].ids[2] = 6;
	hexIJK.ijk[4].ids[0] = 4; hexIJK.ijk[4].ids[1] = 0; hexIJK.ijk[4].ids[2] = 4;
	hexIJK.ijk[5].ids[0] = 4; hexIJK.ijk[5].ids[1] = 2; hexIJK.ijk[5].ids[2] = 4;
	hexIJK.ijk[6].ids[0] = 6; hexIJK.ijk[6].ids[1] = 2; hexIJK.ijk[6].ids[2] = 4;
	hexIJK.ijk[7].ids[0] = 6; hexIJK.ijk[7].ids[1] = 0; hexIJK.ijk[7].ids[2] = 4; hexIJKs.push_back(hexIJK); // 16
	hexIJK.ijk[0].ids[0] = 0; hexIJK.ijk[0].ids[1] = 2; hexIJK.ijk[0].ids[2] = 6;
	hexIJK.ijk[1].ids[0] = 0; hexIJK.ijk[1].ids[1] = 4; hexIJK.ijk[1].ids[2] = 6;
	hexIJK.ijk[2].ids[0] = 2; hexIJK.ijk[2].ids[1] = 4; hexIJK.ijk[2].ids[2] = 6;
	hexIJK.ijk[3].ids[0] = 2; hexIJK.ijk[3].ids[1] = 2; hexIJK.ijk[3].ids[2] = 6;
	hexIJK.ijk[4].ids[0] = 0; hexIJK.ijk[4].ids[1] = 2; hexIJK.ijk[4].ids[2] = 4;
	hexIJK.ijk[5].ids[0] = 0; hexIJK.ijk[5].ids[1] = 4; hexIJK.ijk[5].ids[2] = 4;
	hexIJK.ijk[6].ids[0] = 2; hexIJK.ijk[6].ids[1] = 4; hexIJK.ijk[6].ids[2] = 4;
	hexIJK.ijk[7].ids[0] = 2; hexIJK.ijk[7].ids[1] = 2; hexIJK.ijk[7].ids[2] = 4; hexIJKs.push_back(hexIJK); // 17
	hexIJK.ijk[0].ids[0] = 2; hexIJK.ijk[0].ids[1] = 2; hexIJK.ijk[0].ids[2] = 6;
	hexIJK.ijk[1].ids[0] = 2; hexIJK.ijk[1].ids[1] = 4; hexIJK.ijk[1].ids[2] = 6;
	hexIJK.ijk[2].ids[0] = 4; hexIJK.ijk[2].ids[1] = 4; hexIJK.ijk[2].ids[2] = 6;
	hexIJK.ijk[3].ids[0] = 4; hexIJK.ijk[3].ids[1] = 2; hexIJK.ijk[3].ids[2] = 6;
	hexIJK.ijk[4].ids[0] = 2; hexIJK.ijk[4].ids[1] = 2; hexIJK.ijk[4].ids[2] = 4;
	hexIJK.ijk[5].ids[0] = 2; hexIJK.ijk[5].ids[1] = 4; hexIJK.ijk[5].ids[2] = 4;
	hexIJK.ijk[6].ids[0] = 4; hexIJK.ijk[6].ids[1] = 4; hexIJK.ijk[6].ids[2] = 4;
	hexIJK.ijk[7].ids[0] = 4; hexIJK.ijk[7].ids[1] = 2; hexIJK.ijk[7].ids[2] = 4; hexIJKs.push_back(hexIJK); // 18
	hexIJK.ijk[0].ids[0] = 4; hexIJK.ijk[0].ids[1] = 2; hexIJK.ijk[0].ids[2] = 6;
	hexIJK.ijk[1].ids[0] = 4; hexIJK.ijk[1].ids[1] = 4; hexIJK.ijk[1].ids[2] = 6;
	hexIJK.ijk[2].ids[0] = 6; hexIJK.ijk[2].ids[1] = 4; hexIJK.ijk[2].ids[2] = 6;
	hexIJK.ijk[3].ids[0] = 6; hexIJK.ijk[3].ids[1] = 2; hexIJK.ijk[3].ids[2] = 6;
	hexIJK.ijk[4].ids[0] = 4; hexIJK.ijk[4].ids[1] = 2; hexIJK.ijk[4].ids[2] = 4;
	hexIJK.ijk[5].ids[0] = 4; hexIJK.ijk[5].ids[1] = 4; hexIJK.ijk[5].ids[2] = 4;
	hexIJK.ijk[6].ids[0] = 6; hexIJK.ijk[6].ids[1] = 4; hexIJK.ijk[6].ids[2] = 4;
	hexIJK.ijk[7].ids[0] = 6; hexIJK.ijk[7].ids[1] = 2; hexIJK.ijk[7].ids[2] = 4; hexIJKs.push_back(hexIJK); // 19
	hexIJK.ijk[0].ids[0] = 0; hexIJK.ijk[0].ids[1] = 4; hexIJK.ijk[0].ids[2] = 6;
	hexIJK.ijk[1].ids[0] = 0; hexIJK.ijk[1].ids[1] = 6; hexIJK.ijk[1].ids[2] = 6;
	hexIJK.ijk[2].ids[0] = 2; hexIJK.ijk[2].ids[1] = 6; hexIJK.ijk[2].ids[2] = 6;
	hexIJK.ijk[3].ids[0] = 2; hexIJK.ijk[3].ids[1] = 4; hexIJK.ijk[3].ids[2] = 6;
	hexIJK.ijk[4].ids[0] = 0; hexIJK.ijk[4].ids[1] = 4; hexIJK.ijk[4].ids[2] = 4;
	hexIJK.ijk[5].ids[0] = 0; hexIJK.ijk[5].ids[1] = 6; hexIJK.ijk[5].ids[2] = 4;
	hexIJK.ijk[6].ids[0] = 2; hexIJK.ijk[6].ids[1] = 6; hexIJK.ijk[6].ids[2] = 4;
	hexIJK.ijk[7].ids[0] = 2; hexIJK.ijk[7].ids[1] = 4; hexIJK.ijk[7].ids[2] = 4; hexIJKs.push_back(hexIJK); // 20
	hexIJK.ijk[0].ids[0] = 2; hexIJK.ijk[0].ids[1] = 4; hexIJK.ijk[0].ids[2] = 6;
	hexIJK.ijk[1].ids[0] = 2; hexIJK.ijk[1].ids[1] = 6; hexIJK.ijk[1].ids[2] = 6;
	hexIJK.ijk[2].ids[0] = 4; hexIJK.ijk[2].ids[1] = 6; hexIJK.ijk[2].ids[2] = 6;
	hexIJK.ijk[3].ids[0] = 4; hexIJK.ijk[3].ids[1] = 4; hexIJK.ijk[3].ids[2] = 6;
	hexIJK.ijk[4].ids[0] = 2; hexIJK.ijk[4].ids[1] = 4; hexIJK.ijk[4].ids[2] = 4;
	hexIJK.ijk[5].ids[0] = 2; hexIJK.ijk[5].ids[1] = 6; hexIJK.ijk[5].ids[2] = 4;
	hexIJK.ijk[6].ids[0] = 4; hexIJK.ijk[6].ids[1] = 6; hexIJK.ijk[6].ids[2] = 4;
	hexIJK.ijk[7].ids[0] = 4; hexIJK.ijk[7].ids[1] = 4; hexIJK.ijk[7].ids[2] = 4; hexIJKs.push_back(hexIJK); // 21
	hexIJK.ijk[0].ids[0] = 4; hexIJK.ijk[0].ids[1] = 4; hexIJK.ijk[0].ids[2] = 6;
	hexIJK.ijk[1].ids[0] = 4; hexIJK.ijk[1].ids[1] = 6; hexIJK.ijk[1].ids[2] = 6;
	hexIJK.ijk[2].ids[0] = 6; hexIJK.ijk[2].ids[1] = 6; hexIJK.ijk[2].ids[2] = 6;
	hexIJK.ijk[3].ids[0] = 6; hexIJK.ijk[3].ids[1] = 4; hexIJK.ijk[3].ids[2] = 6;
	hexIJK.ijk[4].ids[0] = 4; hexIJK.ijk[4].ids[1] = 4; hexIJK.ijk[4].ids[2] = 4;
	hexIJK.ijk[5].ids[0] = 4; hexIJK.ijk[5].ids[1] = 6; hexIJK.ijk[5].ids[2] = 4;
	hexIJK.ijk[6].ids[0] = 6; hexIJK.ijk[6].ids[1] = 6; hexIJK.ijk[6].ids[2] = 4;
	hexIJK.ijk[7].ids[0] = 6; hexIJK.ijk[7].ids[1] = 4; hexIJK.ijk[7].ids[2] = 4; hexIJKs.push_back(hexIJK); // 22

	regionsOnRegionToBuild_[keyRegion] = hexIJKs;

	keyRegion = tableMarkedNodes<8> (1,1,1,1,1,1,1,1); // all
	hexIJKs.clear(); //hexIJKs.resize(nbRegionsOnRegionToBuild_[keyRegion]);

	hexIJK.ijk[0].ids[0] = 0; hexIJK.ijk[0].ids[1] = 0; hexIJK.ijk[0].ids[2] = 2;
	hexIJK.ijk[1].ids[0] = 0; hexIJK.ijk[1].ids[1] = 2; hexIJK.ijk[1].ids[2] = 2;
	hexIJK.ijk[2].ids[0] = 2; hexIJK.ijk[2].ids[1] = 2; hexIJK.ijk[2].ids[2] = 2;
	hexIJK.ijk[3].ids[0] = 2; hexIJK.ijk[3].ids[1] = 0; hexIJK.ijk[3].ids[2] = 2;
	hexIJK.ijk[4].ids[0] = 0; hexIJK.ijk[4].ids[1] = 0; hexIJK.ijk[4].ids[2] = 0;
	hexIJK.ijk[5].ids[0] = 0; hexIJK.ijk[5].ids[1] = 2; hexIJK.ijk[5].ids[2] = 0;
	hexIJK.ijk[6].ids[0] = 2; hexIJK.ijk[6].ids[1] = 2; hexIJK.ijk[6].ids[2] = 0;
	hexIJK.ijk[7].ids[0] = 2; hexIJK.ijk[7].ids[1] = 0; hexIJK.ijk[7].ids[2] = 0; hexIJKs.push_back(hexIJK); //  1
	hexIJK.ijk[0].ids[0] = 2; hexIJK.ijk[0].ids[1] = 0; hexIJK.ijk[0].ids[2] = 2;
	hexIJK.ijk[1].ids[0] = 2; hexIJK.ijk[1].ids[1] = 2; hexIJK.ijk[1].ids[2] = 2;
	hexIJK.ijk[2].ids[0] = 4; hexIJK.ijk[2].ids[1] = 2; hexIJK.ijk[2].ids[2] = 2;
	hexIJK.ijk[3].ids[0] = 4; hexIJK.ijk[3].ids[1] = 0; hexIJK.ijk[3].ids[2] = 2;
	hexIJK.ijk[4].ids[0] = 2; hexIJK.ijk[4].ids[1] = 0; hexIJK.ijk[4].ids[2] = 0;
	hexIJK.ijk[5].ids[0] = 2; hexIJK.ijk[5].ids[1] = 2; hexIJK.ijk[5].ids[2] = 0;
	hexIJK.ijk[6].ids[0] = 4; hexIJK.ijk[6].ids[1] = 2; hexIJK.ijk[6].ids[2] = 0;
	hexIJK.ijk[7].ids[0] = 4; hexIJK.ijk[7].ids[1] = 0; hexIJK.ijk[7].ids[2] = 0; hexIJKs.push_back(hexIJK); //  2
	hexIJK.ijk[0].ids[0] = 4; hexIJK.ijk[0].ids[1] = 0; hexIJK.ijk[0].ids[2] = 2;
	hexIJK.ijk[1].ids[0] = 4; hexIJK.ijk[1].ids[1] = 2; hexIJK.ijk[1].ids[2] = 2;
	hexIJK.ijk[2].ids[0] = 6; hexIJK.ijk[2].ids[1] = 2; hexIJK.ijk[2].ids[2] = 2;
	hexIJK.ijk[3].ids[0] = 6; hexIJK.ijk[3].ids[1] = 0; hexIJK.ijk[3].ids[2] = 2;
	hexIJK.ijk[4].ids[0] = 4; hexIJK.ijk[4].ids[1] = 0; hexIJK.ijk[4].ids[2] = 0;
	hexIJK.ijk[5].ids[0] = 4; hexIJK.ijk[5].ids[1] = 2; hexIJK.ijk[5].ids[2] = 0;
	hexIJK.ijk[6].ids[0] = 6; hexIJK.ijk[6].ids[1] = 2; hexIJK.ijk[6].ids[2] = 0;
	hexIJK.ijk[7].ids[0] = 6; hexIJK.ijk[7].ids[1] = 0; hexIJK.ijk[7].ids[2] = 0; hexIJKs.push_back(hexIJK); //  3
	hexIJK.ijk[0].ids[0] = 0; hexIJK.ijk[0].ids[1] = 2; hexIJK.ijk[0].ids[2] = 2;
	hexIJK.ijk[1].ids[0] = 0; hexIJK.ijk[1].ids[1] = 4; hexIJK.ijk[1].ids[2] = 2;
	hexIJK.ijk[2].ids[0] = 2; hexIJK.ijk[2].ids[1] = 4; hexIJK.ijk[2].ids[2] = 2;
	hexIJK.ijk[3].ids[0] = 2; hexIJK.ijk[3].ids[1] = 2; hexIJK.ijk[3].ids[2] = 2;
	hexIJK.ijk[4].ids[0] = 0; hexIJK.ijk[4].ids[1] = 2; hexIJK.ijk[4].ids[2] = 0;
	hexIJK.ijk[5].ids[0] = 0; hexIJK.ijk[5].ids[1] = 4; hexIJK.ijk[5].ids[2] = 0;
	hexIJK.ijk[6].ids[0] = 2; hexIJK.ijk[6].ids[1] = 4; hexIJK.ijk[6].ids[2] = 0;
	hexIJK.ijk[7].ids[0] = 2; hexIJK.ijk[7].ids[1] = 2; hexIJK.ijk[7].ids[2] = 0; hexIJKs.push_back(hexIJK); //  4
	hexIJK.ijk[0].ids[0] = 2; hexIJK.ijk[0].ids[1] = 2; hexIJK.ijk[0].ids[2] = 2;
	hexIJK.ijk[1].ids[0] = 2; hexIJK.ijk[1].ids[1] = 4; hexIJK.ijk[1].ids[2] = 2;
	hexIJK.ijk[2].ids[0] = 4; hexIJK.ijk[2].ids[1] = 4; hexIJK.ijk[2].ids[2] = 2;
	hexIJK.ijk[3].ids[0] = 4; hexIJK.ijk[3].ids[1] = 2; hexIJK.ijk[3].ids[2] = 2;
	hexIJK.ijk[4].ids[0] = 2; hexIJK.ijk[4].ids[1] = 2; hexIJK.ijk[4].ids[2] = 0;
	hexIJK.ijk[5].ids[0] = 2; hexIJK.ijk[5].ids[1] = 4; hexIJK.ijk[5].ids[2] = 0;
	hexIJK.ijk[6].ids[0] = 4; hexIJK.ijk[6].ids[1] = 4; hexIJK.ijk[6].ids[2] = 0;
	hexIJK.ijk[7].ids[0] = 4; hexIJK.ijk[7].ids[1] = 2; hexIJK.ijk[7].ids[2] = 0; hexIJKs.push_back(hexIJK); //  5
	hexIJK.ijk[0].ids[0] = 4; hexIJK.ijk[0].ids[1] = 2; hexIJK.ijk[0].ids[2] = 2;
	hexIJK.ijk[1].ids[0] = 4; hexIJK.ijk[1].ids[1] = 4; hexIJK.ijk[1].ids[2] = 2;
	hexIJK.ijk[2].ids[0] = 6; hexIJK.ijk[2].ids[1] = 4; hexIJK.ijk[2].ids[2] = 2;
	hexIJK.ijk[3].ids[0] = 6; hexIJK.ijk[3].ids[1] = 2; hexIJK.ijk[3].ids[2] = 2;
	hexIJK.ijk[4].ids[0] = 4; hexIJK.ijk[4].ids[1] = 2; hexIJK.ijk[4].ids[2] = 0;
	hexIJK.ijk[5].ids[0] = 4; hexIJK.ijk[5].ids[1] = 4; hexIJK.ijk[5].ids[2] = 0;
	hexIJK.ijk[6].ids[0] = 6; hexIJK.ijk[6].ids[1] = 4; hexIJK.ijk[6].ids[2] = 0;
	hexIJK.ijk[7].ids[0] = 6; hexIJK.ijk[7].ids[1] = 2; hexIJK.ijk[7].ids[2] = 0; hexIJKs.push_back(hexIJK); //  6
	hexIJK.ijk[0].ids[0] = 0; hexIJK.ijk[0].ids[1] = 4; hexIJK.ijk[0].ids[2] = 2;
	hexIJK.ijk[1].ids[0] = 0; hexIJK.ijk[1].ids[1] = 6; hexIJK.ijk[1].ids[2] = 2;
	hexIJK.ijk[2].ids[0] = 2; hexIJK.ijk[2].ids[1] = 6; hexIJK.ijk[2].ids[2] = 2;
	hexIJK.ijk[3].ids[0] = 2; hexIJK.ijk[3].ids[1] = 4; hexIJK.ijk[3].ids[2] = 2;
	hexIJK.ijk[4].ids[0] = 0; hexIJK.ijk[4].ids[1] = 4; hexIJK.ijk[4].ids[2] = 0;
	hexIJK.ijk[5].ids[0] = 0; hexIJK.ijk[5].ids[1] = 6; hexIJK.ijk[5].ids[2] = 0;
	hexIJK.ijk[6].ids[0] = 2; hexIJK.ijk[6].ids[1] = 6; hexIJK.ijk[6].ids[2] = 0;
	hexIJK.ijk[7].ids[0] = 2; hexIJK.ijk[7].ids[1] = 4; hexIJK.ijk[7].ids[2] = 0; hexIJKs.push_back(hexIJK); //  7
	hexIJK.ijk[0].ids[0] = 2; hexIJK.ijk[0].ids[1] = 4; hexIJK.ijk[0].ids[2] = 2;
	hexIJK.ijk[1].ids[0] = 2; hexIJK.ijk[1].ids[1] = 6; hexIJK.ijk[1].ids[2] = 2;
	hexIJK.ijk[2].ids[0] = 4; hexIJK.ijk[2].ids[1] = 6; hexIJK.ijk[2].ids[2] = 2;
	hexIJK.ijk[3].ids[0] = 4; hexIJK.ijk[3].ids[1] = 4; hexIJK.ijk[3].ids[2] = 2;
	hexIJK.ijk[4].ids[0] = 2; hexIJK.ijk[4].ids[1] = 4; hexIJK.ijk[4].ids[2] = 0;
	hexIJK.ijk[5].ids[0] = 2; hexIJK.ijk[5].ids[1] = 6; hexIJK.ijk[5].ids[2] = 0;
	hexIJK.ijk[6].ids[0] = 4; hexIJK.ijk[6].ids[1] = 6; hexIJK.ijk[6].ids[2] = 0;
	hexIJK.ijk[7].ids[0] = 4; hexIJK.ijk[7].ids[1] = 4; hexIJK.ijk[7].ids[2] = 0; hexIJKs.push_back(hexIJK); //  8
	hexIJK.ijk[0].ids[0] = 4; hexIJK.ijk[0].ids[1] = 4; hexIJK.ijk[0].ids[2] = 2;
	hexIJK.ijk[1].ids[0] = 4; hexIJK.ijk[1].ids[1] = 6; hexIJK.ijk[1].ids[2] = 2;
	hexIJK.ijk[2].ids[0] = 6; hexIJK.ijk[2].ids[1] = 6; hexIJK.ijk[2].ids[2] = 2;
	hexIJK.ijk[3].ids[0] = 6; hexIJK.ijk[3].ids[1] = 4; hexIJK.ijk[3].ids[2] = 2;
	hexIJK.ijk[4].ids[0] = 4; hexIJK.ijk[4].ids[1] = 4; hexIJK.ijk[4].ids[2] = 0;
	hexIJK.ijk[5].ids[0] = 4; hexIJK.ijk[5].ids[1] = 6; hexIJK.ijk[5].ids[2] = 0;
	hexIJK.ijk[6].ids[0] = 6; hexIJK.ijk[6].ids[1] = 6; hexIJK.ijk[6].ids[2] = 0;
	hexIJK.ijk[7].ids[0] = 6; hexIJK.ijk[7].ids[1] = 4; hexIJK.ijk[7].ids[2] = 0; hexIJKs.push_back(hexIJK); //  9

	hexIJK.ijk[0].ids[0] = 0; hexIJK.ijk[0].ids[1] = 0; hexIJK.ijk[0].ids[2] = 4;
	hexIJK.ijk[1].ids[0] = 0; hexIJK.ijk[1].ids[1] = 2; hexIJK.ijk[1].ids[2] = 4;
	hexIJK.ijk[2].ids[0] = 2; hexIJK.ijk[2].ids[1] = 2; hexIJK.ijk[2].ids[2] = 4;
	hexIJK.ijk[3].ids[0] = 2; hexIJK.ijk[3].ids[1] = 0; hexIJK.ijk[3].ids[2] = 4;
	hexIJK.ijk[4].ids[0] = 0; hexIJK.ijk[4].ids[1] = 0; hexIJK.ijk[4].ids[2] = 2;
	hexIJK.ijk[5].ids[0] = 0; hexIJK.ijk[5].ids[1] = 2; hexIJK.ijk[5].ids[2] = 2;
	hexIJK.ijk[6].ids[0] = 2; hexIJK.ijk[6].ids[1] = 2; hexIJK.ijk[6].ids[2] = 2;
	hexIJK.ijk[7].ids[0] = 2; hexIJK.ijk[7].ids[1] = 0; hexIJK.ijk[7].ids[2] = 2; hexIJKs.push_back(hexIJK); // 10
	hexIJK.ijk[0].ids[0] = 2; hexIJK.ijk[0].ids[1] = 0; hexIJK.ijk[0].ids[2] = 4;
	hexIJK.ijk[1].ids[0] = 2; hexIJK.ijk[1].ids[1] = 2; hexIJK.ijk[1].ids[2] = 4;
	hexIJK.ijk[2].ids[0] = 4; hexIJK.ijk[2].ids[1] = 2; hexIJK.ijk[2].ids[2] = 4;
	hexIJK.ijk[3].ids[0] = 4; hexIJK.ijk[3].ids[1] = 0; hexIJK.ijk[3].ids[2] = 4;
	hexIJK.ijk[4].ids[0] = 2; hexIJK.ijk[4].ids[1] = 0; hexIJK.ijk[4].ids[2] = 2;
	hexIJK.ijk[5].ids[0] = 2; hexIJK.ijk[5].ids[1] = 2; hexIJK.ijk[5].ids[2] = 2;
	hexIJK.ijk[6].ids[0] = 4; hexIJK.ijk[6].ids[1] = 2; hexIJK.ijk[6].ids[2] = 2;
	hexIJK.ijk[7].ids[0] = 4; hexIJK.ijk[7].ids[1] = 0; hexIJK.ijk[7].ids[2] = 2; hexIJKs.push_back(hexIJK); // 11
	hexIJK.ijk[0].ids[0] = 4; hexIJK.ijk[0].ids[1] = 0; hexIJK.ijk[0].ids[2] = 4;
	hexIJK.ijk[1].ids[0] = 4; hexIJK.ijk[1].ids[1] = 2; hexIJK.ijk[1].ids[2] = 4;
	hexIJK.ijk[2].ids[0] = 6; hexIJK.ijk[2].ids[1] = 2; hexIJK.ijk[2].ids[2] = 4;
	hexIJK.ijk[3].ids[0] = 6; hexIJK.ijk[3].ids[1] = 0; hexIJK.ijk[3].ids[2] = 4;
	hexIJK.ijk[4].ids[0] = 4; hexIJK.ijk[4].ids[1] = 0; hexIJK.ijk[4].ids[2] = 2;
	hexIJK.ijk[5].ids[0] = 4; hexIJK.ijk[5].ids[1] = 2; hexIJK.ijk[5].ids[2] = 2;
	hexIJK.ijk[6].ids[0] = 6; hexIJK.ijk[6].ids[1] = 2; hexIJK.ijk[6].ids[2] = 2;
	hexIJK.ijk[7].ids[0] = 6; hexIJK.ijk[7].ids[1] = 0; hexIJK.ijk[7].ids[2] = 2; hexIJKs.push_back(hexIJK); // 12
	hexIJK.ijk[0].ids[0] = 0; hexIJK.ijk[0].ids[1] = 2; hexIJK.ijk[0].ids[2] = 4;
	hexIJK.ijk[1].ids[0] = 0; hexIJK.ijk[1].ids[1] = 4; hexIJK.ijk[1].ids[2] = 4;
	hexIJK.ijk[2].ids[0] = 2; hexIJK.ijk[2].ids[1] = 4; hexIJK.ijk[2].ids[2] = 4;
	hexIJK.ijk[3].ids[0] = 2; hexIJK.ijk[3].ids[1] = 2; hexIJK.ijk[3].ids[2] = 4;
	hexIJK.ijk[4].ids[0] = 0; hexIJK.ijk[4].ids[1] = 2; hexIJK.ijk[4].ids[2] = 2;
	hexIJK.ijk[5].ids[0] = 0; hexIJK.ijk[5].ids[1] = 4; hexIJK.ijk[5].ids[2] = 2;
	hexIJK.ijk[6].ids[0] = 2; hexIJK.ijk[6].ids[1] = 4; hexIJK.ijk[6].ids[2] = 2;
	hexIJK.ijk[7].ids[0] = 2; hexIJK.ijk[7].ids[1] = 2; hexIJK.ijk[7].ids[2] = 2; hexIJKs.push_back(hexIJK); // 13
	hexIJK.ijk[0].ids[0] = 2; hexIJK.ijk[0].ids[1] = 2; hexIJK.ijk[0].ids[2] = 4;
	hexIJK.ijk[1].ids[0] = 2; hexIJK.ijk[1].ids[1] = 4; hexIJK.ijk[1].ids[2] = 4;
	hexIJK.ijk[2].ids[0] = 4; hexIJK.ijk[2].ids[1] = 4; hexIJK.ijk[2].ids[2] = 4;
	hexIJK.ijk[3].ids[0] = 4; hexIJK.ijk[3].ids[1] = 2; hexIJK.ijk[3].ids[2] = 4;
	hexIJK.ijk[4].ids[0] = 2; hexIJK.ijk[4].ids[1] = 2; hexIJK.ijk[4].ids[2] = 2;
	hexIJK.ijk[5].ids[0] = 2; hexIJK.ijk[5].ids[1] = 4; hexIJK.ijk[5].ids[2] = 2;
	hexIJK.ijk[6].ids[0] = 4; hexIJK.ijk[6].ids[1] = 4; hexIJK.ijk[6].ids[2] = 2;
	hexIJK.ijk[7].ids[0] = 4; hexIJK.ijk[7].ids[1] = 2; hexIJK.ijk[7].ids[2] = 2; hexIJKs.push_back(hexIJK); // 14
	hexIJK.ijk[0].ids[0] = 4; hexIJK.ijk[0].ids[1] = 2; hexIJK.ijk[0].ids[2] = 4;
	hexIJK.ijk[1].ids[0] = 4; hexIJK.ijk[1].ids[1] = 4; hexIJK.ijk[1].ids[2] = 4;
	hexIJK.ijk[2].ids[0] = 6; hexIJK.ijk[2].ids[1] = 4; hexIJK.ijk[2].ids[2] = 4;
	hexIJK.ijk[3].ids[0] = 6; hexIJK.ijk[3].ids[1] = 2; hexIJK.ijk[3].ids[2] = 4;
	hexIJK.ijk[4].ids[0] = 4; hexIJK.ijk[4].ids[1] = 2; hexIJK.ijk[4].ids[2] = 2;
	hexIJK.ijk[5].ids[0] = 4; hexIJK.ijk[5].ids[1] = 4; hexIJK.ijk[5].ids[2] = 2;
	hexIJK.ijk[6].ids[0] = 6; hexIJK.ijk[6].ids[1] = 4; hexIJK.ijk[6].ids[2] = 2;
	hexIJK.ijk[7].ids[0] = 6; hexIJK.ijk[7].ids[1] = 2; hexIJK.ijk[7].ids[2] = 2; hexIJKs.push_back(hexIJK); // 15
	hexIJK.ijk[0].ids[0] = 0; hexIJK.ijk[0].ids[1] = 4; hexIJK.ijk[0].ids[2] = 4;
	hexIJK.ijk[1].ids[0] = 0; hexIJK.ijk[1].ids[1] = 6; hexIJK.ijk[1].ids[2] = 4;
	hexIJK.ijk[2].ids[0] = 2; hexIJK.ijk[2].ids[1] = 6; hexIJK.ijk[2].ids[2] = 4;
	hexIJK.ijk[3].ids[0] = 2; hexIJK.ijk[3].ids[1] = 4; hexIJK.ijk[3].ids[2] = 4;
	hexIJK.ijk[4].ids[0] = 0; hexIJK.ijk[4].ids[1] = 4; hexIJK.ijk[4].ids[2] = 2;
	hexIJK.ijk[5].ids[0] = 0; hexIJK.ijk[5].ids[1] = 6; hexIJK.ijk[5].ids[2] = 2;
	hexIJK.ijk[6].ids[0] = 2; hexIJK.ijk[6].ids[1] = 6; hexIJK.ijk[6].ids[2] = 2;
	hexIJK.ijk[7].ids[0] = 2; hexIJK.ijk[7].ids[1] = 4; hexIJK.ijk[7].ids[2] = 2; hexIJKs.push_back(hexIJK); // 16
	hexIJK.ijk[0].ids[0] = 2; hexIJK.ijk[0].ids[1] = 4; hexIJK.ijk[0].ids[2] = 4;
	hexIJK.ijk[1].ids[0] = 2; hexIJK.ijk[1].ids[1] = 6; hexIJK.ijk[1].ids[2] = 4;
	hexIJK.ijk[2].ids[0] = 4; hexIJK.ijk[2].ids[1] = 6; hexIJK.ijk[2].ids[2] = 4;
	hexIJK.ijk[3].ids[0] = 4; hexIJK.ijk[3].ids[1] = 4; hexIJK.ijk[3].ids[2] = 4;
	hexIJK.ijk[4].ids[0] = 2; hexIJK.ijk[4].ids[1] = 4; hexIJK.ijk[4].ids[2] = 2;
	hexIJK.ijk[5].ids[0] = 2; hexIJK.ijk[5].ids[1] = 6; hexIJK.ijk[5].ids[2] = 2;
	hexIJK.ijk[6].ids[0] = 4; hexIJK.ijk[6].ids[1] = 6; hexIJK.ijk[6].ids[2] = 2;
	hexIJK.ijk[7].ids[0] = 4; hexIJK.ijk[7].ids[1] = 4; hexIJK.ijk[7].ids[2] = 2; hexIJKs.push_back(hexIJK); // 17
	hexIJK.ijk[0].ids[0] = 4; hexIJK.ijk[0].ids[1] = 4; hexIJK.ijk[0].ids[2] = 4;
	hexIJK.ijk[1].ids[0] = 4; hexIJK.ijk[1].ids[1] = 6; hexIJK.ijk[1].ids[2] = 4;
	hexIJK.ijk[2].ids[0] = 6; hexIJK.ijk[2].ids[1] = 6; hexIJK.ijk[2].ids[2] = 4;
	hexIJK.ijk[3].ids[0] = 6; hexIJK.ijk[3].ids[1] = 4; hexIJK.ijk[3].ids[2] = 4;
	hexIJK.ijk[4].ids[0] = 4; hexIJK.ijk[4].ids[1] = 4; hexIJK.ijk[4].ids[2] = 2;
	hexIJK.ijk[5].ids[0] = 4; hexIJK.ijk[5].ids[1] = 6; hexIJK.ijk[5].ids[2] = 2;
	hexIJK.ijk[6].ids[0] = 6; hexIJK.ijk[6].ids[1] = 6; hexIJK.ijk[6].ids[2] = 2;
	hexIJK.ijk[7].ids[0] = 6; hexIJK.ijk[7].ids[1] = 4; hexIJK.ijk[7].ids[2] = 2; hexIJKs.push_back(hexIJK); // 18

	hexIJK.ijk[0].ids[0] = 0; hexIJK.ijk[0].ids[1] = 0; hexIJK.ijk[0].ids[2] = 6;
	hexIJK.ijk[1].ids[0] = 0; hexIJK.ijk[1].ids[1] = 2; hexIJK.ijk[1].ids[2] = 6;
	hexIJK.ijk[2].ids[0] = 2; hexIJK.ijk[2].ids[1] = 2; hexIJK.ijk[2].ids[2] = 6;
	hexIJK.ijk[3].ids[0] = 2; hexIJK.ijk[3].ids[1] = 0; hexIJK.ijk[3].ids[2] = 6;
	hexIJK.ijk[4].ids[0] = 0; hexIJK.ijk[4].ids[1] = 0; hexIJK.ijk[4].ids[2] = 4;
	hexIJK.ijk[5].ids[0] = 0; hexIJK.ijk[5].ids[1] = 2; hexIJK.ijk[5].ids[2] = 4;
	hexIJK.ijk[6].ids[0] = 2; hexIJK.ijk[6].ids[1] = 2; hexIJK.ijk[6].ids[2] = 4;
	hexIJK.ijk[7].ids[0] = 2; hexIJK.ijk[7].ids[1] = 0; hexIJK.ijk[7].ids[2] = 4; hexIJKs.push_back(hexIJK); // 19
	hexIJK.ijk[0].ids[0] = 2; hexIJK.ijk[0].ids[1] = 0; hexIJK.ijk[0].ids[2] = 6;
	hexIJK.ijk[1].ids[0] = 2; hexIJK.ijk[1].ids[1] = 2; hexIJK.ijk[1].ids[2] = 6;
	hexIJK.ijk[2].ids[0] = 4; hexIJK.ijk[2].ids[1] = 2; hexIJK.ijk[2].ids[2] = 6;
	hexIJK.ijk[3].ids[0] = 4; hexIJK.ijk[3].ids[1] = 0; hexIJK.ijk[3].ids[2] = 6;
	hexIJK.ijk[4].ids[0] = 2; hexIJK.ijk[4].ids[1] = 0; hexIJK.ijk[4].ids[2] = 4;
	hexIJK.ijk[5].ids[0] = 2; hexIJK.ijk[5].ids[1] = 2; hexIJK.ijk[5].ids[2] = 4;
	hexIJK.ijk[6].ids[0] = 4; hexIJK.ijk[6].ids[1] = 2; hexIJK.ijk[6].ids[2] = 4;
	hexIJK.ijk[7].ids[0] = 4; hexIJK.ijk[7].ids[1] = 0; hexIJK.ijk[7].ids[2] = 4; hexIJKs.push_back(hexIJK); // 20
	hexIJK.ijk[0].ids[0] = 4; hexIJK.ijk[0].ids[1] = 0; hexIJK.ijk[0].ids[2] = 6;
	hexIJK.ijk[1].ids[0] = 4; hexIJK.ijk[1].ids[1] = 2; hexIJK.ijk[1].ids[2] = 6;
	hexIJK.ijk[2].ids[0] = 6; hexIJK.ijk[2].ids[1] = 2; hexIJK.ijk[2].ids[2] = 6;
	hexIJK.ijk[3].ids[0] = 6; hexIJK.ijk[3].ids[1] = 0; hexIJK.ijk[3].ids[2] = 6;
	hexIJK.ijk[4].ids[0] = 4; hexIJK.ijk[4].ids[1] = 0; hexIJK.ijk[4].ids[2] = 4;
	hexIJK.ijk[5].ids[0] = 4; hexIJK.ijk[5].ids[1] = 2; hexIJK.ijk[5].ids[2] = 4;
	hexIJK.ijk[6].ids[0] = 6; hexIJK.ijk[6].ids[1] = 2; hexIJK.ijk[6].ids[2] = 4;
	hexIJK.ijk[7].ids[0] = 6; hexIJK.ijk[7].ids[1] = 0; hexIJK.ijk[7].ids[2] = 4; hexIJKs.push_back(hexIJK); // 21
	hexIJK.ijk[0].ids[0] = 0; hexIJK.ijk[0].ids[1] = 2; hexIJK.ijk[0].ids[2] = 6;
	hexIJK.ijk[1].ids[0] = 0; hexIJK.ijk[1].ids[1] = 4; hexIJK.ijk[1].ids[2] = 6;
	hexIJK.ijk[2].ids[0] = 2; hexIJK.ijk[2].ids[1] = 4; hexIJK.ijk[2].ids[2] = 6;
	hexIJK.ijk[3].ids[0] = 2; hexIJK.ijk[3].ids[1] = 2; hexIJK.ijk[3].ids[2] = 6;
	hexIJK.ijk[4].ids[0] = 0; hexIJK.ijk[4].ids[1] = 2; hexIJK.ijk[4].ids[2] = 4;
	hexIJK.ijk[5].ids[0] = 0; hexIJK.ijk[5].ids[1] = 4; hexIJK.ijk[5].ids[2] = 4;
	hexIJK.ijk[6].ids[0] = 2; hexIJK.ijk[6].ids[1] = 4; hexIJK.ijk[6].ids[2] = 4;
	hexIJK.ijk[7].ids[0] = 2; hexIJK.ijk[7].ids[1] = 2; hexIJK.ijk[7].ids[2] = 4; hexIJKs.push_back(hexIJK); // 22
	hexIJK.ijk[0].ids[0] = 2; hexIJK.ijk[0].ids[1] = 2; hexIJK.ijk[0].ids[2] = 6;
	hexIJK.ijk[1].ids[0] = 2; hexIJK.ijk[1].ids[1] = 4; hexIJK.ijk[1].ids[2] = 6;
	hexIJK.ijk[2].ids[0] = 4; hexIJK.ijk[2].ids[1] = 4; hexIJK.ijk[2].ids[2] = 6;
	hexIJK.ijk[3].ids[0] = 4; hexIJK.ijk[3].ids[1] = 2; hexIJK.ijk[3].ids[2] = 6;
	hexIJK.ijk[4].ids[0] = 2; hexIJK.ijk[4].ids[1] = 2; hexIJK.ijk[4].ids[2] = 4;
	hexIJK.ijk[5].ids[0] = 2; hexIJK.ijk[5].ids[1] = 4; hexIJK.ijk[5].ids[2] = 4;
	hexIJK.ijk[6].ids[0] = 4; hexIJK.ijk[6].ids[1] = 4; hexIJK.ijk[6].ids[2] = 4;
	hexIJK.ijk[7].ids[0] = 4; hexIJK.ijk[7].ids[1] = 2; hexIJK.ijk[7].ids[2] = 4; hexIJKs.push_back(hexIJK); // 23
	hexIJK.ijk[0].ids[0] = 4; hexIJK.ijk[0].ids[1] = 2; hexIJK.ijk[0].ids[2] = 6;
	hexIJK.ijk[1].ids[0] = 4; hexIJK.ijk[1].ids[1] = 4; hexIJK.ijk[1].ids[2] = 6;
	hexIJK.ijk[2].ids[0] = 6; hexIJK.ijk[2].ids[1] = 4; hexIJK.ijk[2].ids[2] = 6;
	hexIJK.ijk[3].ids[0] = 6; hexIJK.ijk[3].ids[1] = 2; hexIJK.ijk[3].ids[2] = 6;
	hexIJK.ijk[4].ids[0] = 4; hexIJK.ijk[4].ids[1] = 2; hexIJK.ijk[4].ids[2] = 4;
	hexIJK.ijk[5].ids[0] = 4; hexIJK.ijk[5].ids[1] = 4; hexIJK.ijk[5].ids[2] = 4;
	hexIJK.ijk[6].ids[0] = 6; hexIJK.ijk[6].ids[1] = 4; hexIJK.ijk[6].ids[2] = 4;
	hexIJK.ijk[7].ids[0] = 6; hexIJK.ijk[7].ids[1] = 2; hexIJK.ijk[7].ids[2] = 4; hexIJKs.push_back(hexIJK); // 24
	hexIJK.ijk[0].ids[0] = 0; hexIJK.ijk[0].ids[1] = 4; hexIJK.ijk[0].ids[2] = 6;
	hexIJK.ijk[1].ids[0] = 0; hexIJK.ijk[1].ids[1] = 6; hexIJK.ijk[1].ids[2] = 6;
	hexIJK.ijk[2].ids[0] = 2; hexIJK.ijk[2].ids[1] = 6; hexIJK.ijk[2].ids[2] = 6;
	hexIJK.ijk[3].ids[0] = 2; hexIJK.ijk[3].ids[1] = 4; hexIJK.ijk[3].ids[2] = 6;
	hexIJK.ijk[4].ids[0] = 0; hexIJK.ijk[4].ids[1] = 4; hexIJK.ijk[4].ids[2] = 4;
	hexIJK.ijk[5].ids[0] = 0; hexIJK.ijk[5].ids[1] = 6; hexIJK.ijk[5].ids[2] = 4;
	hexIJK.ijk[6].ids[0] = 2; hexIJK.ijk[6].ids[1] = 6; hexIJK.ijk[6].ids[2] = 4;
	hexIJK.ijk[7].ids[0] = 2; hexIJK.ijk[7].ids[1] = 4; hexIJK.ijk[7].ids[2] = 4; hexIJKs.push_back(hexIJK); // 25
	hexIJK.ijk[0].ids[0] = 2; hexIJK.ijk[0].ids[1] = 4; hexIJK.ijk[0].ids[2] = 6;
	hexIJK.ijk[1].ids[0] = 2; hexIJK.ijk[1].ids[1] = 6; hexIJK.ijk[1].ids[2] = 6;
	hexIJK.ijk[2].ids[0] = 4; hexIJK.ijk[2].ids[1] = 6; hexIJK.ijk[2].ids[2] = 6;
	hexIJK.ijk[3].ids[0] = 4; hexIJK.ijk[3].ids[1] = 4; hexIJK.ijk[3].ids[2] = 6;
	hexIJK.ijk[4].ids[0] = 2; hexIJK.ijk[4].ids[1] = 4; hexIJK.ijk[4].ids[2] = 4;
	hexIJK.ijk[5].ids[0] = 2; hexIJK.ijk[5].ids[1] = 6; hexIJK.ijk[5].ids[2] = 4;
	hexIJK.ijk[6].ids[0] = 4; hexIJK.ijk[6].ids[1] = 6; hexIJK.ijk[6].ids[2] = 4;
	hexIJK.ijk[7].ids[0] = 4; hexIJK.ijk[7].ids[1] = 4; hexIJK.ijk[7].ids[2] = 4; hexIJKs.push_back(hexIJK); // 26
	hexIJK.ijk[0].ids[0] = 4; hexIJK.ijk[0].ids[1] = 4; hexIJK.ijk[0].ids[2] = 6;
	hexIJK.ijk[1].ids[0] = 4; hexIJK.ijk[1].ids[1] = 6; hexIJK.ijk[1].ids[2] = 6;
	hexIJK.ijk[2].ids[0] = 6; hexIJK.ijk[2].ids[1] = 6; hexIJK.ijk[2].ids[2] = 6;
	hexIJK.ijk[3].ids[0] = 6; hexIJK.ijk[3].ids[1] = 4; hexIJK.ijk[3].ids[2] = 6;
	hexIJK.ijk[4].ids[0] = 4; hexIJK.ijk[4].ids[1] = 4; hexIJK.ijk[4].ids[2] = 4;
	hexIJK.ijk[5].ids[0] = 4; hexIJK.ijk[5].ids[1] = 6; hexIJK.ijk[5].ids[2] = 4;
	hexIJK.ijk[6].ids[0] = 6; hexIJK.ijk[6].ids[1] = 6; hexIJK.ijk[6].ids[2] = 4;
	hexIJK.ijk[7].ids[0] = 6; hexIJK.ijk[7].ids[1] = 4; hexIJK.ijk[7].ids[2] = 4; hexIJKs.push_back(hexIJK); // 27

	regionsOnRegionToBuild_[keyRegion] = hexIJKs;


//	hexIJK.ijk[0].ids[0] = 0; hexIJK.ijk[0].ids[1] = 0; hexIJK.ijk[0].ids[2] = 0;
//	hexIJK.ijk[1].ids[0] = 0; hexIJK.ijk[1].ids[1] = 0; hexIJK.ijk[1].ids[2] = 0;
//	hexIJK.ijk[2].ids[0] = 0; hexIJK.ijk[2].ids[1] = 0; hexIJK.ijk[2].ids[2] = 0;
//	hexIJK.ijk[3].ids[0] = 0; hexIJK.ijk[3].ids[1] = 0; hexIJK.ijk[3].ids[2] = 0;
//	hexIJK.ijk[4].ids[0] = 0; hexIJK.ijk[4].ids[1] = 0; hexIJK.ijk[4].ids[2] = 0;
//	hexIJK.ijk[5].ids[0] = 0; hexIJK.ijk[5].ids[1] = 0; hexIJK.ijk[5].ids[2] = 0;
//	hexIJK.ijk[6].ids[0] = 0; hexIJK.ijk[6].ids[1] = 0; hexIJK.ijk[6].ids[2] = 0;
//	hexIJK.ijk[7].ids[0] = 0; hexIJK.ijk[7].ids[1] = 0; hexIJK.ijk[7].ids[2] = 0; hexIJKs.push_back(hexIJK);


//	std::vector<FaceIDs> tableFace;
//	struct FaceIDs faceIDs;
//	faceIDs.ids[0] = 0; faceIDs.ids[1] = 1; faceIDs.ids[2] = 2; faceIDs.ids[3] = 3;
//	tableFace.push_back(faceIDs);
//
//	// all
//	keyFace = tableMarkedNodes<4> (0,0,0,0);
//	tableFace.clear();
//	faceIDs.ids[0] = 0; faceIDs.ids[1] = 3; faceIDs.ids[2] = 12; faceIDs.ids[3] = 15; tableFace.push_back(faceIDs);
//	facesToBuild_[keyFace] = tableFace;
//
//	keyFace = tableMarkedNodes<4> (1,1,1,1);
//	tableFace.clear();
//	faceIDs.ids[0] =  0; faceIDs.ids[1] =  1; faceIDs.ids[2] =  5; faceIDs.ids[3] =  4; tableFace.push_back(faceIDs);
//	faceIDs.ids[0] =  1; faceIDs.ids[1] =  2; faceIDs.ids[2] =  6; faceIDs.ids[3] =  5; tableFace.push_back(faceIDs);
//	faceIDs.ids[0] =  2; faceIDs.ids[1] =  3; faceIDs.ids[2] =  7; faceIDs.ids[3] =  6; tableFace.push_back(faceIDs);
//	faceIDs.ids[0] =  4; faceIDs.ids[1] =  5; faceIDs.ids[2] =  9; faceIDs.ids[3] =  8; tableFace.push_back(faceIDs);
//	faceIDs.ids[0] =  5; faceIDs.ids[1] =  6; faceIDs.ids[2] = 10; faceIDs.ids[3] =  9; tableFace.push_back(faceIDs);
//	faceIDs.ids[0] =  6; faceIDs.ids[1] =  7; faceIDs.ids[2] = 11; faceIDs.ids[3] = 10; tableFace.push_back(faceIDs);
//	faceIDs.ids[0] =  8; faceIDs.ids[1] =  9; faceIDs.ids[2] = 13; faceIDs.ids[3] = 12; tableFace.push_back(faceIDs);
//	faceIDs.ids[0] =  9; faceIDs.ids[1] = 10; faceIDs.ids[2] = 14; faceIDs.ids[3] = 13; tableFace.push_back(faceIDs);
//	faceIDs.ids[0] = 10; faceIDs.ids[1] = 11; faceIDs.ids[2] = 15; faceIDs.ids[3] = 14; tableFace.push_back(faceIDs);
//	facesToBuild_[keyFace] = tableFace;
//
//	// corners
//	keyFace = tableMarkedNodes<4> (1,0,0,0);
//	tableFace.clear();
//	faceIDs.ids[0] =  0; faceIDs.ids[1] =  1; faceIDs.ids[2] =  5; faceIDs.ids[3] =  4; tableFace.push_back(faceIDs);
//	faceIDs.ids[0] =  5; faceIDs.ids[1] =  3; faceIDs.ids[2] = 15; faceIDs.ids[3] = 13; tableFace.push_back(faceIDs);
//	faceIDs.ids[0] =  4; faceIDs.ids[1] =  5; faceIDs.ids[2] = 13; faceIDs.ids[3] = 12; tableFace.push_back(faceIDs);
//	facesToBuild_[keyFace] = tableFace;

	buildSurfacesNeighbors();
}
/*----------------------------------------------------------------------------*/
MeshSplitter3Refinement::~MeshSplitter3Refinement()
{}
/*----------------------------------------------------------------------------*/
void
MeshSplitter3Refinement::markRegionsToSplit(
		const int& AMarkRegionToSplit,
		const int& AMarkFaceToSplit,
		const int& AMarkEdgeToSplit,
		const int& AMarkNodeToSplit)
{
	this->initialization();

	// criteria for splitting
	markNeighborRegionsInOutModelAABBTree(AMarkRegionToSplit);

	markRegionsIntersectedTwiceByModelAABBTree(AMarkRegionToSplit);

	// mark nodes owned by marked regions
	{
		IGMesh::region_iterator it  = this->mesh_.regions_begin();
		for(;!it.isDone(); it.next()) {
			Region current_region = it.value();

			if(this->mesh_.isMarked(current_region,AMarkRegionToSplit)) {
				std::vector<Node> nodes = current_region.get<Node>();
				for(unsigned int iNode=0; iNode<nodes.size(); iNode++) {
					this->mesh_.mark(nodes[iNode],AMarkNodeToSplit);
				}
			}
		}
	}

	// count number of marked regions and nodes
	{
		unsigned int nbMarkedRegions = 0;

		IGMesh::region_iterator itr  = this->mesh_.regions_begin();
		for(;!itr.isDone(); itr.next()) {
			Region current_region = itr.value();

			if(this->mesh_.isMarked(current_region,AMarkRegionToSplit)) {
				nbMarkedRegions++;
			}

		}

		unsigned int nbMarkedNodes = 0;

		IGMesh::node_iterator itn  = this->mesh_.nodes_begin();
		for(;!itn.isDone(); itn.next()) {
			Node current_node = itn.value();

			if(this->mesh_.isMarked(current_node,AMarkNodeToSplit)) {
				nbMarkedNodes++;
			}

		}

		std::cout<<"markRegionsToSplit : number of marked regions "<< nbMarkedRegions << std::endl;
		std::cout<<"markRegionsToSplit : number of marked nodes "<< nbMarkedNodes << std::endl;
	}

}
/*----------------------------------------------------------------------------*/
void
MeshSplitter3Refinement::markEntitiesToSplitForAddEdges(
		gmds::Node ANode,
		gmds::Face AFace,
                const int& AMarkRegionToSplit,
                const int& AMarkFaceToSplit,
                const int& AMarkEdgeToSplit,
                const int& AMarkNodeToSplit)
{
	// select the node opposite to ANode in AFace
	std::vector<gmds::Node> nodes = AFace.get<gmds::Node>();

	gmds::Node node2Refine;

	for(unsigned int iNode=0; iNode<nodes.size(); iNode++) {
		
		if(ANode == nodes[iNode]) {
			node2Refine = nodes[(iNode+2)%nodes.size()];
			break;
		}	
	}
	
	this->mesh_.mark(node2Refine,AMarkNodeToSplit);

	//std::vector<gmds::Edge> edges = node2Refine.get<gmds::Edge>();
	//for(unsigned int iEdge=0; iEdge<edges.size(); iEdge++) {
	//	this->mesh_.mark(edges[iEdge],AMarkEdgeToSplit);
	//}
	//std::vector<gmds::Face> faces = node2Refine.get<gmds::Face>();
        //for(unsigned int iFace=0; iFace<faces.size(); iFace++) {
        //        this->mesh_.mark(faces[iFace],AMarkFaceToSplit);
        //}
        
	IGMesh::region_iterator itr  = this->mesh_.regions_begin();
        for(;!itr.isDone(); itr.next()) {
                Region current_region = itr.value();

		std::vector<gmds::Node> nodes = current_region.get<Node>();
		for(unsigned int iNode=0; iNode<nodes.size(); iNode++) {
			if(this->mesh_.isMarked(nodes[iNode],AMarkNodeToSplit)) {
				this->mesh_.mark(current_region,AMarkRegionToSplit);
			}
		}
	}

	//std::vector<gmds::Region> regions = node2Refine.get<gmds::Region>();
        //for(unsigned int iRegion=0; iRegion<regions.size(); iRegion++) {
        //        this->mesh_.mark(regions[iRegion],AMarkRegionToSplit);
        //}
}
/*----------------------------------------------------------------------------*/
void
MeshSplitter3Refinement::avoidBadConfigurations(
		const int& AMarkRegionToSplit,
		const int& AMarkFaceToSplit,
		const int& AMarkEdgeToSplit,
		const int& AMarkNodeToSplit)
{
	unsigned int nbModifiedBadCell = 0;

	// the bad configuration resolution propagates across the mesh.
	// It will result in obtaining a "topological "convex area.
	do {
		nbModifiedBadCell = 0;

		IGMesh::region_iterator it  = this->mesh_.regions_begin();
		for(;!it.isDone(); it.next()) {
			Region current_region = it.value();

			std::vector<Node> nodes = current_region.get<Node>();

			std::vector<bool> tableBool(nodes.size());

			bool isSplit = false;
			for(unsigned int iNode=0; iNode<nodes.size(); iNode++) {
				tableBool[iNode] = this->mesh_.isMarked(nodes[iNode],AMarkNodeToSplit);

				if(this->mesh_.isMarked(nodes[iNode],AMarkNodeToSplit)) {
					isSplit = true;
				}
			}

			// no node of this region has to be split, hence we do nothing
			if(!isSplit) {
				continue;
			}

			// at this point we have a region that has to be split, and this algorithm
			// currently only works on hexahedra.
			if(current_region.getType() != gmds::GMDS_HEX) {
				throw GMDSException("MeshSplitter3Refinement::avoidBadConfigurations "
						"can not work on non hex cells");
			}

			// this region has at least one node that must be split, hence the region is marked
			this->mesh_.mark(current_region,AMarkRegionToSplit);

			tableMarkedNodes<8> tableMark(tableBool[0],tableBool[1],tableBool[2],tableBool[3],
										  tableBool[4],tableBool[5],tableBool[6],tableBool[7]);

			tableMarkedNodes<8> tableMarkTarget;
			if(lookupNodes_.find(tableMark) != lookupNodes_.end()) {
				tableMarkTarget = lookupNodes_[tableMark];
			} else {
				tableMarkTarget = lookupNodes_[tableMarkedNodes<8> (1,1,1,1,1,1,1,1)];
			}

			//
			if(tableMarkTarget != tableMark) {

				for(unsigned int iNode=0; iNode<nodes.size(); iNode++) {
					if(tableMarkTarget.isMarked(iNode)) {
						this->mesh_.mark(nodes[iNode],AMarkNodeToSplit);
					}
				}
				nbModifiedBadCell++;
			}

		} // for(;!it->isDone(); it->next()) {
	}
	while(nbModifiedBadCell != 0);

	// count number of marked regions and nodes
	// mark faces and edges that need to be split
	{
		unsigned int nbMarkedRegions = 0;

		IGMesh::region_iterator itr  = this->mesh_.regions_begin();
		for(;!itr.isDone(); itr.next()) {
			Region current_region = itr.value();

			if(this->mesh_.isMarked(current_region,AMarkRegionToSplit)) {
				nbMarkedRegions++;
			}
		}

		unsigned int nbMarkedFaces = 0;

		IGMesh::face_iterator itf  = this->mesh_.faces_begin();
		for(;!itf.isDone(); itf.next()) {
			Face current_face = itf.value();

			std::vector<Node> nodes = current_face.get<Node>();
			for(unsigned int iNode=0; iNode<nodes.size(); iNode++) {
				if(this->mesh_.isMarked(nodes[iNode],AMarkNodeToSplit)) {
					this->mesh_.mark(current_face,AMarkFaceToSplit);
					nbMarkedFaces++;
					break;
				}
			}
		}

		unsigned int nbMarkedEdges = 0;

		IGMesh::edge_iterator ite  = this->mesh_.edges_begin();
		for(;!ite.isDone(); ite.next()) {
			Edge current_edge = ite.value();

			std::vector<Node> nodes = current_edge.get<Node>();
			for(unsigned int iNode=0; iNode<nodes.size(); iNode++) {
				if(this->mesh_.isMarked(nodes[iNode],AMarkNodeToSplit)) {
					this->mesh_.mark(current_edge,AMarkEdgeToSplit);
					nbMarkedEdges++;
					break;
				}
			}
		}

		unsigned int nbMarkedNodes = 0;

		IGMesh::node_iterator itn  = this->mesh_.nodes_begin();
		for(;!itn.isDone(); itn.next()) {
			Node current_node = itn.value();

			// criteria for splitting?
			if(this->mesh_.isMarked(current_node,AMarkNodeToSplit)) {
				nbMarkedNodes++;
			}
		}

		std::cout<<"markRegionsToSplit : number of marked regions "<< nbMarkedRegions << std::endl;
		std::cout<<"markRegionsToSplit : number of marked faces "<< nbMarkedFaces << std::endl;
		std::cout<<"markRegionsToSplit : number of marked edges "<< nbMarkedEdges << std::endl;
		std::cout<<"markRegionsToSplit : number of marked nodes "<< nbMarkedNodes << std::endl;
	}
}
/*----------------------------------------------------------------------------*/
void
MeshSplitter3Refinement::split(
		const int& AMarkRegionToSplit,
		const int& AMarkFaceToSplit,
		const int& AMarkEdgeToSplit,
		const int& AMarkNodeToSplit)
{
//	// mark entities to treat, meaning that entities created during the
//	// spliting phase will not be split
//	int markOrigineEntities = this->mesh_.getNewMark();
//	{
//		typename Mesh<TMask>::edges_iterator it  = this->mesh_.edges_begin();
//		for(;!it->isDone(); it->next()) {
//			Edge* current_edge = it->currentItem();
//
//			this->mesh_.mark(current_edge,markOrigineEntities);
//		}
//	}
//	{
//		typename Mesh<TMask>::faces_iterator it  = this->mesh_.faces_begin();
//		for(;!it->isDone(); it->next()) {
//			Face* current_face = it->currentItem();
//
//			this->mesh_.mark(current_face,markOrigineEntities);
//		}
//	}

	// split the edges
	{
		std::cout<<"begin split the edges"<<std::endl;
		unsigned int nbSplitEdges = 0;

		IGMesh::edge_iterator it  = this->mesh_.edges_begin();
		for(;!it.isDone(); it.next()) {
			Edge current_edge = it.value();

			if(this->mesh_.isMarked(current_edge,AMarkEdgeToSplit)) {
				if(splitEdge(AMarkNodeToSplit,AMarkEdgeToSplit,current_edge)) {
					nbSplitEdges++;
				} else {
					throw GMDSException("MeshSplitter3Refinement::split : edge should be split but is not... problem!");
				}
			}
		}
		std::cout<<"number of split edges "<< nbSplitEdges <<std::endl;
		std::cout<<"end split the edges"<<std::endl;
	}

//	exportVTK("afterEdge.unf");

	// split the faces
	{
		std::cout<<"begin split the faces"<<std::endl;
		unsigned int nbSplitFaces = 0;

		IGMesh::face_iterator it  = this->mesh_.faces_begin();
		for(;!it.isDone(); it.next()) {
			Face current_face = it.value();

			if(this->mesh_.isMarked(current_face,AMarkFaceToSplit)) {
				if(splitFace(AMarkNodeToSplit,AMarkFaceToSplit,current_face)) {
					nbSplitFaces++;
				} else {
					throw GMDSException("MeshSplitter3Refinement::split : face should be split but is not... problem!");
				}
			}
		}
		std::cout<<"number of split faces "<< nbSplitFaces <<std::endl;
		std::cout<<"end split the faces"<<std::endl;
	}

//	exportVTK("afterFace.unf");

	// split the regions
	{
		std::vector<Region> regions2Split;

		std::cout<<"begin split the regions"<<std::endl;

		IGMesh::region_iterator it  = this->mesh_.regions_begin();
		for(;!it.isDone(); it.next()) {
			Region current_region = it.value();

			if(this->mesh_.isMarked(current_region,AMarkRegionToSplit)) {
				regions2Split.push_back(current_region);
			}
		}

		for(unsigned int iRegion=0; iRegion<regions2Split.size(); iRegion++) {
			Region current_region = regions2Split[iRegion];
			std::cout<<"split region id "<<current_region.getID()<<std::endl;
			splitRegion(AMarkNodeToSplit,current_region);
		}

		std::cout<<"end split the regions"<<std::endl;
	}

//	this->mesh_.unmarkAll(markOrigineEntities);
//	this->mesh_.freeMark(markOrigineEntities);

//	exportVTK("afterRegion.unf");

}
/*----------------------------------------------------------------------------*/
bool
MeshSplitter3Refinement::splitEdge(const int& AMarkNodeToSplit, const int& AMarkEdgeToSplit, Edge& AEdge)
{
	std::vector<Node> edgeNodes = AEdge.get<Node>();

	if(!this->mesh_.isMarked(edgeNodes[0],AMarkNodeToSplit) && !this->mesh_.isMarked(edgeNodes[1],AMarkNodeToSplit)) {
		// nothing to do : this edge is not split
		return false;
	}

	this->mesh_.mark(AEdge,AMarkEdgeToSplit);

	// get classification in order to propagate it to the newly
	// created edges if it exists
	Variable<geom::GeomEntity*>* edgeClassification;

	if(this->mesh_.doesGeometricClassificationExist(1)) {
		edgeClassification = this->mesh_.getGeometricClassification(1);
	}

	std::vector<Node> edgeSplitNodes;
	edgeSplitNodes.resize(2,Node());

	if(this->mesh_.isMarked(edgeNodes[0],AMarkNodeToSplit) && this->mesh_.isMarked(edgeNodes[1],AMarkNodeToSplit)) {

		gmds::math::Point point1 = edgeNodes[0].getPoint();
		gmds::math::Point point2 = edgeNodes[1].getPoint();

		gmds::math::Point point_tmp1 = ((point1 * 2.) + point2) * (1./3.);
		gmds::math::Point point_tmp2 = (point1 + (point2 * 2.)) * (1./3.);

		Node newNode1 = this->mesh_.newNode(point_tmp1);
		this->splitNodesCreated_.insert(newNode1);
		Node newNode2 = this->mesh_.newNode(point_tmp2);
		this->splitNodesCreated_.insert(newNode2);

		edgeSplitNodes[0] = newNode1;
		edgeSplitNodes[1] = newNode2;

		Edge edge0 = this->mesh_.newEdge(edgeNodes[0],newNode1);
		this->splitEdgesCreated_.insert(edge0);
		Edge edge1 = this->mesh_.newEdge(newNode1,newNode2);
		this->splitEdgesCreated_.insert(edge1);
		Edge edge2 = this->mesh_.newEdge(newNode2,edgeNodes[1]);
		this->splitEdgesCreated_.insert(edge2);

		if(this->mesh_.doesGeometricClassificationExist(1)) {
			(*edgeClassification)[edge0.getID()] = (*edgeClassification)[AEdge.getID()];
			(*edgeClassification)[edge1.getID()] = (*edgeClassification)[AEdge.getID()];
			(*edgeClassification)[edge2.getID()] = (*edgeClassification)[AEdge.getID()];
		}
	}

	if(this->mesh_.isMarked(edgeNodes[0],AMarkNodeToSplit) && !this->mesh_.isMarked(edgeNodes[1],AMarkNodeToSplit)) {

		gmds::math::Point point1 = edgeNodes[0].getPoint();
		gmds::math::Point point2 = edgeNodes[1].getPoint();

		gmds::math::Point point_tmp1 = ((point1 * 2.) + point2) * (1./3.);

		Node newNode1 = this->mesh_.newNode(point_tmp1);
		this->splitNodesCreated_.insert(newNode1);

		edgeSplitNodes[0] = newNode1;

		Edge edge0 = this->mesh_.newEdge(edgeNodes[0],newNode1);
		this->splitEdgesCreated_.insert(edge0);
		Edge edge1 = this->mesh_.newEdge(newNode1,edgeNodes[1]);
		this->splitEdgesCreated_.insert(edge1);

		if(this->mesh_.doesGeometricClassificationExist(1)) {
			(*edgeClassification)[edge0.getID()] = (*edgeClassification)[AEdge.getID()];
			(*edgeClassification)[edge1.getID()] = (*edgeClassification)[AEdge.getID()];
		}
	}

	if(!this->mesh_.isMarked(edgeNodes[0],AMarkNodeToSplit) && this->mesh_.isMarked(edgeNodes[1],AMarkNodeToSplit)) {

		gmds::math::Point point1 = edgeNodes[0].getPoint();
		gmds::math::Point point2 = edgeNodes[1].getPoint();

		gmds::math::Point point_tmp2 = (point1 + (point2 * 2.)) * (1./3.);

		Node newNode2 = this->mesh_.newNode(point_tmp2);
		this->splitNodesCreated_.insert(newNode2);

		edgeSplitNodes[1] = newNode2;

		Edge edge0 = this->mesh_.newEdge(edgeNodes[0],newNode2);
		this->splitEdgesCreated_.insert(edge0);
		Edge edge1 = this->mesh_.newEdge(newNode2,edgeNodes[1]);
		this->splitEdgesCreated_.insert(edge1);

		if(this->mesh_.doesGeometricClassificationExist(1)) {
			(*edgeClassification)[edge0.getID()] = (*edgeClassification)[AEdge.getID()];
			(*edgeClassification)[edge1.getID()] = (*edgeClassification)[AEdge.getID()];
		}

//		gmds::Marks32 marks = this->mesh_.getMarks(AEdge);
//		this->mesh_.mark<Edge>(edge0,marks);
//		this->mesh_.mark<Edge>(edge1,marks);
//		this->mesh_.unmark(edge0,AMarkEdgeToSplit);
//		this->mesh_.unmark(edge1,AMarkEdgeToSplit);
	}

	edgeSplitNodes_[AEdge] = edgeSplitNodes;

	return true;
}
/*----------------------------------------------------------------------------*/
bool
MeshSplitter3Refinement::splitFace(const int& AMarkNodeToSplit, const int& AMarkFaceToSplit, Face& AFace)
{
	std::vector<Node> faceNodes = AFace.get<Node>();

	bool tableBool[4];

	bool isSplit = false;
	for(unsigned int iNode=0; iNode<faceNodes.size(); iNode++) {
		tableBool[iNode] = this->mesh_.isMarked(faceNodes[iNode],AMarkNodeToSplit);

		if(this->mesh_.isMarked(faceNodes[iNode],AMarkNodeToSplit)) {
			isSplit = true;
		}
	}

	// no node of this face has to be split, hence we do nothing
	if(!isSplit) {
		return false;
	}

	this->mesh_.mark(AFace,AMarkFaceToSplit);

	// get classification in order to propagate it to the newly
	// created faces if it exists
	Variable<geom::GeomEntity*>* faceClassification;

	if(this->mesh_.doesGeometricClassificationExist(2)) {
		faceClassification = this->mesh_.getGeometricClassification(2);
	}

	// this face has at least one node that must be split
	tableMarkedNodes<4> tableMark(tableBool[0],tableBool[1],tableBool[2],tableBool[3]);

	std::vector<Node> faceSplitNodes;
	faceSplitNodes.resize(12,Node());
	gmds::math::Point cornersPoints[4];

	for(unsigned int iNode=0; iNode<faceNodes.size(); iNode++) {
		cornersPoints[iNode] = faceNodes[iNode].getPoint();
	}

	// place each edge in the face
	std::vector<Edge> faceEdges = AFace.get<Edge>();
	std::vector<Edge> faceEdgeOrdered;
	bool egdeIsOrdered[4];
	for(unsigned int iNode=0; iNode<faceNodes.size(); iNode++) {
		Node current_node = faceNodes[iNode];
		Node next_node = faceNodes[(iNode+1)%faceNodes.size()];

		for(unsigned int iEdge=0; iEdge<faceEdges.size(); iEdge++) {
			Edge current_edge = faceEdges[iEdge];

			std::vector<Node> edgeNodes = current_edge.get<Node>();
			if((edgeNodes[0] == current_node) && (edgeNodes[1] == next_node)) {
				faceEdgeOrdered.push_back(current_edge);
				egdeIsOrdered[iNode] = true;
				break;
			}
			if((edgeNodes[1] == current_node) && (edgeNodes[0] == next_node)) {
				faceEdgeOrdered.push_back(current_edge);
				egdeIsOrdered[iNode] = false;
				break;
			}

		}
	}
	if(faceEdgeOrdered.size() != faceEdges.size())
		throw GMDSException("MeshSplitter3Refinement::splitFace could not find every edge in face");

	for(unsigned int iNode=0; iNode<nbNodesOnFaceToBuild_[tableMark]; iNode++) {
		TCellID nodeID = (nodesOnFaceToBuild_[tableMark])[iNode];

		Node newNode;
		unsigned int index = 999;

		switch(nodeID) {
		case  0: {
			throw GMDSException("MeshSplitter3Refinement::splitFace bad id of node");
		}
		case  1: {
			index = 0;
			if(egdeIsOrdered[0])
				newNode = getNodeOfEdge(faceEdgeOrdered[0],1);
			else
				newNode = getNodeOfEdge(faceEdgeOrdered[0],2);
			break;
		}
		case  2: {
			index = 1;
			if(egdeIsOrdered[0])
				newNode = getNodeOfEdge(faceEdgeOrdered[0],2);
			else
				newNode = getNodeOfEdge(faceEdgeOrdered[0],1);
			break;
		}
		case  3: {
			throw GMDSException("MeshSplitter3Refinement::splitFace bad id of node");
		}
		case  4: {
			index = 2;
			if(egdeIsOrdered[3])
				newNode = getNodeOfEdge(faceEdgeOrdered[3],2);
			else
				newNode = getNodeOfEdge(faceEdgeOrdered[3],1);
			break;
		}
		case  5: {
			index = 3;
			gmds::math::Point point = (2.*((2.*cornersPoints[0] +    cornersPoints[1]) * (1./3.)) +
											   ((2.*cornersPoints[3] +    cornersPoints[2]) * (1./3.))) * (1./3.);
			newNode = this->mesh_.newNode(point);
			this->splitNodesCreated_.insert(newNode);
			break;
		}
		case  6: {
			index = 4;
			gmds::math::Point point = (2.*((   cornersPoints[0] + 2.*cornersPoints[1]) * (1./3.)) +
											   ((   cornersPoints[3] + 2.*cornersPoints[2]) * (1./3.))) * (1./3.);
			newNode = this->mesh_.newNode(point);
			this->splitNodesCreated_.insert(newNode);
			break;
		}
		case  7: {
			index = 5;
			if(egdeIsOrdered[1])
				newNode = getNodeOfEdge(faceEdgeOrdered[1],1);
			else
				newNode = getNodeOfEdge(faceEdgeOrdered[1],2);
			break;
		}
		case  8: {
			index = 6;
			if(egdeIsOrdered[3])
				newNode = getNodeOfEdge(faceEdgeOrdered[3],1);
			else
				newNode = getNodeOfEdge(faceEdgeOrdered[3],2);
			break;
		}
		case  9: {
			index = 7;
			gmds::math::Point point = (   ((2.*cornersPoints[0] +    cornersPoints[1]) * (1./3.)) +
											2.*((2.*cornersPoints[3] +    cornersPoints[2]) * (1./3.))) * (1./3.);
			newNode = this->mesh_.newNode(point);
			this->splitNodesCreated_.insert(newNode);
			break;
		}
		case 10: {
			index = 8;
			gmds::math::Point point = (   ((   cornersPoints[0] + 2.*cornersPoints[1]) * (1./3.)) +
											2.*((   cornersPoints[3] + 2.*cornersPoints[2]) * (1./3.))) * (1./3.);
			newNode = this->mesh_.newNode(point);
			this->splitNodesCreated_.insert(newNode);
			break;
		}
		case 11: {
			index = 9;
			if(egdeIsOrdered[1])
				newNode = getNodeOfEdge(faceEdgeOrdered[1],2);
			else
				newNode = getNodeOfEdge(faceEdgeOrdered[1],1);
			break;
		}
		case 12: {
			throw GMDSException("MeshSplitter3Refinement::splitFace bad id of node");
		}
		case 13: {
			index = 10;
			if(egdeIsOrdered[2])
				newNode = getNodeOfEdge(faceEdgeOrdered[2],2);
			else
				newNode = getNodeOfEdge(faceEdgeOrdered[2],1);
			break;
		}
		case 14: {
			index = 11;
			if(egdeIsOrdered[2])
				newNode = getNodeOfEdge(faceEdgeOrdered[2],1);
			else
				newNode = getNodeOfEdge(faceEdgeOrdered[2],2);
			break;
		}
		case 15: {
			throw GMDSException("MeshSplitter3Refinement::splitFace bad id of node");
		}
		default:
			throw GMDSException("MeshSplitter3Refinement::splitFace bad id of node");
		}

		faceSplitNodes[index] = newNode;

	} // for(unsigned int iNode=0; iNode<nbNodesOnFaceToBuild_[tableMark]; iNode++) {

	faceSplitNodes_[AFace] = faceSplitNodes;

	gmds::Marks32 marks = this->mesh_.getMarks<Face>(AFace);

	for(unsigned int iFace=0; iFace<nbFacesOnFaceToBuild_[tableMark]; iFace++) {
		Node node0 = getNodeOfFace(AFace,(facesOnFaceToBuild_[tableMark])[iFace].ids[0]);
		Node node1 = getNodeOfFace(AFace,(facesOnFaceToBuild_[tableMark])[iFace].ids[1]);
		Node node2 = getNodeOfFace(AFace,(facesOnFaceToBuild_[tableMark])[iFace].ids[2]);
		Node node3 = getNodeOfFace(AFace,(facesOnFaceToBuild_[tableMark])[iFace].ids[3]);

		Face newFace_tmp = this->mesh_.newQuad(node0,node1,node2,node3);
		this->splitFacesCreated_.insert(newFace_tmp);

		if(this->mesh_.doesGeometricClassificationExist(2)) {
			(*faceClassification)[newFace_tmp.getID()] = (*faceClassification)[AFace.getID()];
		}

//		this->mesh_.mark(newFace_tmp,marks);
//		this->mesh_.unmark(newFace_tmp,AMarkFaceToSplit);
	}

	return true;
}
/*----------------------------------------------------------------------------*/
void
MeshSplitter3Refinement::splitRegion(const int& AMarkNodeToSplit, Region& ARegion)
{
	std::vector<Node> regionNodes = ARegion.get<Node>();

	// gather the nodes that will be needed
	// Do we use dummy nodes? That would fill a 4*4*4 node array.
	Node nodes[343];
	for(unsigned int iNode=0; iNode<343; iNode++) {
		nodes[iNode] = Node();
	}

	bool tableBool[8];

	// determine the template this hex should be split with
	for(unsigned int iNode=0; iNode<regionNodes.size(); iNode++) {
		tableBool[iNode] = this->mesh_.isMarked(regionNodes[iNode],AMarkNodeToSplit);
	}

	// get classification in order to propagate it to the newly
	// created faces if it exists
	Variable<geom::GeomEntity*>* regionClassification;

	if(this->mesh_.doesGeometricClassificationExist(3)) {
		regionClassification = this->mesh_.getGeometricClassification(3);
	}

	// this region has at least one node that must be split
	tableMarkedNodes<8> tableMark(tableBool[0],tableBool[1],tableBool[2],tableBool[3],
								  tableBool[4],tableBool[5],tableBool[6],tableBool[7]);

	// permut the nodes of the region in order to match one of the five templates
	if(nodePermutToUnitHex_.find(tableMark) == nodePermutToUnitHex_.end()) {
		throw GMDSException("tableMark not found in nodePermutToUnitHex_");
	}

	RegionIDs permutNodes = nodePermutToUnitHex_[tableMark];

//	if (permutNodes == nodePermutToUnitHex_.end()) {
//		throw GMDSException("tableMark nor found in nodePermutToUnitHex_");
//	}


	std::vector<Node> regionNodesOrdered(regionNodes.size());
	for(unsigned int iNode=0; iNode<regionNodes.size(); iNode++) {
		regionNodesOrdered[permutNodes.ids[iNode]] = regionNodes[iNode];
	}
	ARegion.set<Node>(regionNodesOrdered);
	regionNodes = ARegion.get<Node>();

	// build a new tableMark with new nodes ordering
	bool tableOrderedBool[8];

	for(unsigned int iNode=0; iNode<regionNodes.size(); iNode++) {
		tableOrderedBool[permutNodes.ids[iNode]] = tableBool[iNode];
	}

	tableMarkedNodes<8> tableOrdered(tableOrderedBool[0],tableOrderedBool[1],tableOrderedBool[2],tableOrderedBool[3],
			tableOrderedBool[4],tableOrderedBool[5],tableOrderedBool[6],tableOrderedBool[7]);

	// order faces inside regions
	std::vector<Face> faces = ARegion.get<Face>();
	std::vector<Face> faceOrdered; // order is top bottom left right front back
	unsigned int faceOffset[faces.size()];
	unsigned int faceIsOrdered[faces.size()];

	for(unsigned int iFace=0; iFace<faces.size(); iFace++) {
		for(unsigned int iFace2=0; iFace2<faces.size(); iFace2++) {
			std::vector<Node> faceNodes = faces[iFace2].get<Node>();
			std::vector<Node> faceRegionNodes;
			faceRegionNodes.push_back(regionNodes[MESH_SPLITTER_3_REFINEMENT_HEX_ORDERED_FACES[iFace][0]]);
			faceRegionNodes.push_back(regionNodes[MESH_SPLITTER_3_REFINEMENT_HEX_ORDERED_FACES[iFace][1]]);
			faceRegionNodes.push_back(regionNodes[MESH_SPLITTER_3_REFINEMENT_HEX_ORDERED_FACES[iFace][2]]);
			faceRegionNodes.push_back(regionNodes[MESH_SPLITTER_3_REFINEMENT_HEX_ORDERED_FACES[iFace][3]]);

			unsigned int nbCommonNodes = 0;

			for(unsigned int iNode1=0; iNode1<faceNodes.size(); iNode1++) {
				for(unsigned int iNode2=0; iNode2<faceRegionNodes.size(); iNode2++) {

					if(faceNodes[iNode1] == faceRegionNodes[iNode2]) {
						nbCommonNodes++;
					}
				}
			}
			if(nbCommonNodes == 1) {
				throw GMDSException("MeshSplitter3Refinement::splitRegion cannot match only 1 nodes");
			}
			if(nbCommonNodes == 3) {
				throw GMDSException("MeshSplitter3Refinement::splitRegion cannot match only 3 nodes");
			}
			if(nbCommonNodes == 4) {
				faceOrdered.push_back(faces[iFace2]);

				// we have matched this face so now we determine whether it is ordered or not and the offset
				for(unsigned int iNode1=0; iNode1<faceNodes.size(); iNode1++) {
					if(faceNodes[iNode1] == faceRegionNodes[0]) {
						faceOffset[iFace] = iNode1;
						if(faceNodes[(iNode1+1)%faceNodes.size()] == faceRegionNodes[1]) {
							faceIsOrdered[iFace] = true;
						} else {
							faceIsOrdered[iFace] = false;
						}
						break;
					}
				}
				break;
			}
		}
	}

	// we will fill the nodes pool for tapping in when building the sub-hexahedra.
	for(unsigned int iNode=0; iNode<nbNodesOnRegionToBuild_[tableOrdered]; iNode++) {
		NodeHexIJK nodeHexIJK = (nodesOnRegionToBuild_[tableOrdered])[iNode];

		Node newNode;

		if(nodeHexIJK.ids[0] == 0) {
			// get node previously created on left face
//			blabbla
			unsigned int iFace = 2;
			Face current_face = faceOrdered[iFace];

			std::vector<Node> nodesOfFace = current_face.get<Node>();
//			std::cout<<"face "<<nodesOfFace[0]->getID()<<" "<<nodesOfFace[1]->getID()<<" "<<nodesOfFace[2]->getID()<<" "<<nodesOfFace[3]->getID()<<std::endl;

			TCellID indexI = nodeHexIJK.ids[2] / 2;
			TCellID indexJ = nodeHexIJK.ids[1] / 2;
			TCellID faceLocalId = getNodeIdOfFace(faceOffset[iFace],faceIsOrdered[iFace],indexI,indexJ);
			newNode = getNodeOfFace(current_face,faceLocalId);

		} else if(nodeHexIJK.ids[1] == 0) {
			// get node previously created on front face
			unsigned int iFace = 4;
			Face current_face = faceOrdered[iFace];
			TCellID indexI = nodeHexIJK.ids[0] / 2;
			TCellID indexJ = nodeHexIJK.ids[2] / 2;
			TCellID faceLocalId = getNodeIdOfFace(faceOffset[iFace],faceIsOrdered[iFace],indexI,indexJ);
			newNode = getNodeOfFace(current_face,faceLocalId);

		} else if(nodeHexIJK.ids[2] == 0) {
			// get node previously created on bottom face
			unsigned int iFace = 1;
			Face current_face = faceOrdered[iFace];
			TCellID indexI = nodeHexIJK.ids[0] / 2;
			TCellID indexJ = nodeHexIJK.ids[1] / 2;
			TCellID faceLocalId = getNodeIdOfFace(faceOffset[iFace],faceIsOrdered[iFace],indexI,indexJ);
			newNode = getNodeOfFace(current_face,faceLocalId);

		} else if(nodeHexIJK.ids[0] == 6) {
			// get node previously created on right face
			unsigned int iFace = 3;
			Face current_face = faceOrdered[iFace];
			TCellID indexI = nodeHexIJK.ids[2] / 2;
			TCellID indexJ = nodeHexIJK.ids[1] / 2;
			TCellID faceLocalId = getNodeIdOfFace(faceOffset[iFace],faceIsOrdered[iFace],indexI,indexJ);
			newNode = getNodeOfFace(current_face,faceLocalId);

		} else if(nodeHexIJK.ids[1] == 6) {
			// get node previously created on back face
			unsigned int iFace = 5;
			Face current_face = faceOrdered[iFace];
			TCellID indexI = nodeHexIJK.ids[0] / 2;
			TCellID indexJ = nodeHexIJK.ids[2] / 2;
			TCellID faceLocalId = getNodeIdOfFace(faceOffset[iFace],faceIsOrdered[iFace],indexI,indexJ);
			newNode = getNodeOfFace(current_face,faceLocalId);

		} else if(nodeHexIJK.ids[2] == 6) {
			// get node previously created on top face
			unsigned int iFace = 0;
			Face current_face = faceOrdered[iFace];
			TCellID indexI = nodeHexIJK.ids[0] / 2;
			TCellID indexJ = nodeHexIJK.ids[1] / 2;
			TCellID faceLocalId = getNodeIdOfFace(faceOffset[iFace],faceIsOrdered[iFace],indexI,indexJ);
			newNode = getNodeOfFace(current_face,faceLocalId);

		} else {
			// it is an inside node, hence it has to be created
			TCoord coords[3];

			// we will interpolate the coordinates t=for the new nodes
			// method is :
			// using for example all of 4 nodes of front face and 4 nodes
			// of back face for x;
			// all of 4 nodes of right face and 4 nodes
			// of left face for y;
			// all of 4 nodes of bottom face and 4 nodes
			// of top face for z;

//			std::cout<<"coords"<<std::endl;
//			std::cout<<regionNodes[4]->getID()<<" "<<regionNodes[7]->getID()<<" "<<regionNodes[3]->getID()<<" "<<regionNodes[0]->getID()<<std::endl;
//			std::cout<<regionNodes[5]->getID()<<" "<<regionNodes[6]->getID()<<" "<<regionNodes[2]->getID()<<" "<<regionNodes[1]->getID()<<std::endl;
//			std::cout<<regionNodes[4]->getY()<<" "<<regionNodes[7]->getY()<<" "<<regionNodes[3]->getY()<<" "<<regionNodes[0]->getY()<<std::endl;
//			std::cout<<regionNodes[5]->getY()<<" "<<regionNodes[6]->getY()<<" "<<regionNodes[2]->getY()<<" "<<regionNodes[1]->getY()<<std::endl;
//			std::cout<<(6 - nodeHexIJK.ids[1])<<" "<<nodeHexIJK.ids[1]<<std::endl;
//			std::cout<<(regionNodes[4]->getY() + regionNodes[7]->getY() + regionNodes[3]->getY() + regionNodes[0]->getY())<<std::endl;
//			std::cout<<(regionNodes[5]->getY() + regionNodes[6]->getY() + regionNodes[2]->getY() + regionNodes[1]->getY())<<std::endl;

//			coords[0] = ((regionNodes[4]->getX() + regionNodes[0]->getX() + regionNodes[1]->getX() + regionNodes[5]->getX()) * (6 - nodeHexIJK.ids[0]) +
//					(regionNodes[7]->getX() + regionNodes[3]->getX() + regionNodes[2]->getX() + regionNodes[6]->getX()) * nodeHexIJK.ids[0]) / (4.*6.);
//			coords[1] = ((regionNodes[4]->getY() + regionNodes[7]->getY() + regionNodes[3]->getY() + regionNodes[0]->getY()) * (6 - nodeHexIJK.ids[1]) +
//					(regionNodes[5]->getY() + regionNodes[6]->getY() + regionNodes[2]->getY() + regionNodes[1]->getY()) * nodeHexIJK.ids[1]) / (4.*6.);
//			coords[2] = ((regionNodes[4]->getZ() + regionNodes[7]->getZ() + regionNodes[6]->getZ() + regionNodes[5]->getZ()) * (6 - nodeHexIJK.ids[2]) +
//					(regionNodes[0]->getZ() + regionNodes[3]->getZ() + regionNodes[2]->getZ() + regionNodes[1]->getZ()) * nodeHexIJK.ids[2]) / (4.*6.);


			gmds::math::Point oneminusx = (6 - nodeHexIJK.ids[0])/6. * (6 - nodeHexIJK.ids[1])/6. * (6 - nodeHexIJK.ids[2])/6. * regionNodes[4].getPoint()
					+ (nodeHexIJK.ids[0])/6. * (6 - nodeHexIJK.ids[1])/6. * (6 - nodeHexIJK.ids[2])/6. * regionNodes[7].getPoint()
					+ (6 - nodeHexIJK.ids[0])/6. * (nodeHexIJK.ids[1])/6. * (6 - nodeHexIJK.ids[2])/6. * regionNodes[5].getPoint()
					+ (6 - nodeHexIJK.ids[0])/6. * (6 - nodeHexIJK.ids[1])/6. * (nodeHexIJK.ids[2])/6. * regionNodes[0].getPoint()
					+ (nodeHexIJK.ids[0])/6. * (6 - nodeHexIJK.ids[1])/6. * (nodeHexIJK.ids[2])/6. * regionNodes[3].getPoint()
					+ (6 - nodeHexIJK.ids[0])/6. * (nodeHexIJK.ids[1])/6. * (nodeHexIJK.ids[2])/6. * regionNodes[1].getPoint()
					+ (nodeHexIJK.ids[0])/6. * (nodeHexIJK.ids[1])/6. * (6 - nodeHexIJK.ids[2])/6. * regionNodes[6].getPoint()
					+ (nodeHexIJK.ids[0])/6. * (nodeHexIJK.ids[1])/6. * (nodeHexIJK.ids[2])/6. * regionNodes[2].getPoint();

//			blablabla
//			GEPETO::Point<3,TBase> newPoint(0.,0.,0.);
//			GEPETO::Point<3,TBase> newPoint07(0.,0.,0.);
//			GEPETO::Point<3,TBase> newPoint16(0.,0.,0.);
//			GEPETO::Point<3,TBase> newPoint25(0.,0.,0.);
//			GEPETO::Point<3,TBase> newPoint34(0.,0.,0.);




//			newNode = this->mesh_.newNode(coords[0],coords[1],coords[2]);
			newNode = this->mesh_.newNode(oneminusx);
			this->splitNodesCreated_.insert(newNode);
		}

		unsigned int index = nodeHexIJK.ids[2]*7*7 + nodeHexIJK.ids[1]*7 + nodeHexIJK.ids[0];
		nodes[index] = newNode;

	} // for(unsigned int iNode=0; iNode<nbNodesOnRegionToBuild_[tableMark]; iNode++) {

	// we have to add the vertices of the region
	unsigned int index = 6*7*7 + 0*7 + 0;
	nodes[index] = regionNodes[0];
	index = 6*7*7 + 6*7 + 0;
	nodes[index] = regionNodes[1];
	index = 6*7*7 + 6*7 + 6;
	nodes[index] = regionNodes[2];
	index = 6*7*7 + 0*7 + 6;
	nodes[index] = regionNodes[3];
	index = 0*7*7 + 0*7 + 0;
	nodes[index] = regionNodes[4];
	index = 0*7*7 + 6*7 + 0;
	nodes[index] = regionNodes[5];
	index = 0*7*7 + 6*7 + 6;
	nodes[index] = regionNodes[6];
	index = 0*7*7 + 0*7 + 6;
	nodes[index] = regionNodes[7];


//	exportVTK("afterNodeInRegion.unf");

	// hex building
	{

		if(nbRegionsOnRegionToBuild_.find(tableOrdered) == nbRegionsOnRegionToBuild_.end()) {
			throw GMDSException("MeshSplitter3Refinement::splitRegion tableOrdered not found in nbRegionsOnRegionToBuild_");
		}

		unsigned int nbRegion2Create = nbRegionsOnRegionToBuild_[tableOrdered];

		if(regionsOnRegionToBuild_.find(tableOrdered) == regionsOnRegionToBuild_.end()) {
			throw GMDSException("MeshSplitter3Refinement::splitRegion tableOrdered not found in regionsOnRegionToBuild_");
		}

		std::vector<HexIJK> hexIJKs = regionsOnRegionToBuild_[tableOrdered];

		gmds::Marks32 marks = this->mesh_.getMarks(ARegion);

		for(unsigned int iRegion=0; iRegion<nbRegion2Create; iRegion++) {

			HexIJK hexijk = hexIJKs[iRegion];

			unsigned int index[8];
			Node hexNodes[8];
			for(unsigned int iNode=0; iNode<8; iNode++) {
				index[iNode] = hexijk.ijk[iNode].ids[2]*7*7 + hexijk.ijk[iNode].ids[1]*7 + hexijk.ijk[iNode].ids[0];
				hexNodes[iNode] = nodes[index[iNode]];
			}

			Region newRegion_tmp = this->mesh_.newHex(hexNodes[0],hexNodes[1],hexNodes[2],hexNodes[3],
					hexNodes[4],hexNodes[5],hexNodes[6],hexNodes[7]);

			if(this->mesh_.doesGeometricClassificationExist(3)) {
				(*regionClassification)[newRegion_tmp.getID()] = (*regionClassification)[ARegion.getID()];
			}

//			this->mesh_.mark(newRegion_tmp,marks);
//			this->mesh_.unmark(newRegion_tmp,AMarkRegionToSplit);
		}
	}



	// delete the region
	this->splitRegions_.insert(ARegion);
	this->mesh_.deleteRegion(ARegion);
}
/*----------------------------------------------------------------------------*/
Node
MeshSplitter3Refinement::getNodeOfEdge(Edge& AEdge, const TCellID& AId)
{
	std::vector<Node> nodes = AEdge.get<Node>();

	switch(AId) {
	case 0:
		return nodes[0];
	case 1:
		return edgeSplitNodes_[AEdge][0];
	case 2:
		return edgeSplitNodes_[AEdge][1];
	case 3:
		return nodes[1];
	default:
		throw GMDSException("MeshSplitter3Refinement::getNodeOfEdge AId is greater than 3");
	}

}
/*----------------------------------------------------------------------------*/
Node
MeshSplitter3Refinement::getNodeOfFace(Face& AFace, const TCellID& AId)
{
	std::vector<Node> nodes = AFace.get<Node>();

	switch(AId) {
	case 0:
		return nodes[0];
	case  1:
		return faceSplitNodes_[AFace][0];
	case  2:
		return faceSplitNodes_[AFace][1];
	case  3:
		return nodes[1];
	case  4:
		return faceSplitNodes_[AFace][2];
	case  5:
		return faceSplitNodes_[AFace][3];
	case  6:
		return faceSplitNodes_[AFace][4];
	case  7:
		return faceSplitNodes_[AFace][5];
	case  8:
		return faceSplitNodes_[AFace][6];
	case  9:
		return faceSplitNodes_[AFace][7];
	case 10:
		return faceSplitNodes_[AFace][8];
	case 11:
		return faceSplitNodes_[AFace][9];
	case 12:
		return nodes[3];
	case 13:
		return faceSplitNodes_[AFace][10];
	case 14:
		return faceSplitNodes_[AFace][11];
	case 15:
		return nodes[2];
	default:
		throw GMDSException("MeshSplitter3Refinement::getNodeOfFace AId is greater than 16");
	}

}
/*----------------------------------------------------------------------------*/
TCellID
MeshSplitter3Refinement::getNodeIdOfFace(const TCellID& AOffset, const bool& AIsDirect, const TCellID& AIndexI, const TCellID& AIndexJ)
{
	TCellID localId = 0;

	if(AIsDirect) {
		if(AOffset == 0) {
			localId = AIndexJ*4 + AIndexI;
		}
		if(AOffset == 1) {
			localId = AIndexI*4 + (3-AIndexJ);
		}
		if(AOffset == 2) {
			localId = (3-AIndexJ)*4 + (3-AIndexI);
		}
		if(AOffset == 3) {
			localId = (3-AIndexI)*4 + AIndexJ;
		}
	} else {
		if(AOffset == 0) {
			localId = AIndexI*4 + AIndexJ ;
		}
		if(AOffset == 1) {
			localId = AIndexJ*4 + (3-AIndexI);
		}
		if(AOffset == 2) {
			localId = (3-AIndexI)*4 + (3-AIndexJ);
		}
		if(AOffset == 3) {
			localId = (3-AIndexJ)*4 + AIndexI;
		}
	}

	return localId;
}
/*----------------------------------------------------------------------------*/
void
MeshSplitter3Refinement::exportVTK(const std::string& AFile)
{
	gmds::VTKWriter<gmds::IGMesh> w(this->mesh_);
	w.write(AFile,gmds::N|gmds::E|gmds::F|gmds::R);
}
/*----------------------------------------------------------------------------*/
bool
MeshSplitter3Refinement::regionIsIntersectedTwiceByModel(Region& ARegion)
{
	// building the cell hexahedron
	std::vector<Node> nodes = ARegion.get<Node>();

	if(nodes.size() != 8)
		throw GMDSException("Not yet implemented!");

	gmds::math::Hexahedron hexa(
			nodes[4].getPoint(),nodes[7].getPoint(),
			nodes[6].getPoint(),nodes[5].getPoint(),
			nodes[0].getPoint(),nodes[3].getPoint(),
			nodes[2].getPoint(),nodes[1].getPoint());

	// check whether the hex is inside, outside or intersects the model
	std::vector<gmds::geom::GeomSurface*> surfaces;
	this->manager_.getSurfaces(surfaces);

//	if(this->service_.intersectsTwiceTheModel(surfaces,this->surfacesNeighbors_,hexa))
	if(this->service_.intersectsTwiceTheModel(surfaces,this->surfacesNeighbors_,ARegion))
	{
		return true;
	}

	return false;
}
/*----------------------------------------------------------------------------*/
bool
MeshSplitter3Refinement::regionIsIntersectedByModel(Region& ARegion)
{
	// building the cell hexahedron
	std::vector<Node> nodes = ARegion.get<Node>();

	std::vector<gmds::geom::GeomSurface*> surfaces;
	this->manager_.getSurfaces(surfaces);

	if(this->service_.intersects(surfaces,ARegion))
	{
		return true;
	}

	return false;
}
/*----------------------------------------------------------------------------*/
void
MeshSplitter3Refinement::markNeighborRegionsInOutModelAABBTree(
		const int& AMarkRegionToSplit)
{
	IGMesh::region_iterator it  = this->mesh_.regions_begin();

	std::vector<Region> intersectedRegions;
	std::map<gmds::TCellID, std::vector<gmds::math::Triangle> > in_out_triangles;

	std::vector<gmds::geom::GeomSurface* > surfaces;
	this->manager_.getSurfaces(surfaces);

	this->service_.initialization(surfaces);

	// first get the intersecting cells IDs bbox / triangles bbox pairs.
	std::map<gmds::TCellID, std::vector<gmds::math::Triangle*> > inOutTriangles;

	GNode* boxTree = this->aabbRegionsTree_;
	this->service_.intersects(boxTree,inOutTriangles);

	std::map<TCellID,bool> nodesAreInside;
	std::map<TCellID,bool> regionsAreInside;
	std::map<TCellID,bool> regionsAreIntersected;

	// used for progress display
	int iPercent = 0;
	unsigned int iCell = 0;

	std::map<gmds::TCellID, std::vector<gmds::math::Triangle*> >::iterator itCellID = inOutTriangles.begin();


	// marking intersecting cells
	for(; itCellID != inOutTriangles.end(); itCellID++) {

                gmds::Region current_region = this->mesh_.get<gmds::Region>(itCellID->first);
		intersectedRegions.push_back(current_region);
		regionsAreIntersected[current_region.getID()] = true;

		for(int iTri=0; iTri<itCellID->second.size(); iTri++ ) {
		  in_out_triangles[itCellID->first].push_back(*(itCellID->second[iTri]));
		}

	} // for(; itCellID != inOutSurfTriangles.end(); itCellID++)

	for(unsigned int iRegion=0; iRegion<intersectedRegions.size(); iRegion++) {
		Region current_region = intersectedRegions[iRegion];
		double ratio;
		bool isInsideRegion = this->service_.isMainlyInsideMC(
				in_out_triangles[current_region.getID()],current_region,ratio);

		regionsAreInside[current_region.getID()] = isInsideRegion;

		std::vector<Node> nodes = current_region.get<Node>();
		for(unsigned int iNode=0; iNode<nodes.size(); iNode++) {
			gmds::math::Point point = nodes[iNode].getPoint();
			bool isInsidePoint = this->service_.isInside(
					in_out_triangles[current_region.getID()],point);

			nodesAreInside[nodes[iNode].getID()] = isInsidePoint;
		}

	} // for(;it!=ite;it++){

	// for each couple of neighboring intersected regions
	// we will identify potential problems and mark them as refinement candidates
	for(unsigned int iRegion=0; iRegion<intersectedRegions.size(); iRegion++) {
		Region current_region = intersectedRegions[iRegion];

		std::vector<Node> nodes = current_region.get<Node>();

		std::vector<Face> faces = current_region.get<Face>();

		for(int iFace=0; iFace<faces.size(); iFace++) {

			Face current_face = faces[iFace];

			std::vector<Region> regions = current_face.get<Region>();
			if(regions.size() == 1) {
				continue;
			}

			Region neighborRegion;
			if(regions[0] == current_region) {
				neighborRegion = regions[1];
			} else {
				neighborRegion = regions[0];
			}

			// both regions are intersected
			if(regionsAreIntersected.find(neighborRegion.getID()) != regionsAreIntersected.end()) {

				// common face must be either wholly in or out
				std::vector<Node> nodes = current_face.get<Node>();

				bool isWhollyInOut = true;
				bool isIn = nodesAreInside[nodes[0].getID()];
				for(unsigned int iNode=1; iNode<nodes.size(); iNode++) {
					if(nodesAreInside[nodes[iNode].getID()] != isIn) {
						isWhollyInOut = false;
						break;
					}
				}

				if(!isWhollyInOut) {
					continue;
				} else {

					// face IN and both regions OUT
					if(isIn && (!regionsAreInside[current_region.getID()] && !regionsAreInside[neighborRegion.getID()])) {
						this->mesh_.mark(current_region,AMarkRegionToSplit);
						this->mesh_.mark(neighborRegion,AMarkRegionToSplit);
					}

					// face OUT and both regions IN
					if(!isIn && (regionsAreInside[current_region.getID()] && regionsAreInside[neighborRegion.getID()])) {
						this->mesh_.mark(current_region,AMarkRegionToSplit);
						this->mesh_.mark(neighborRegion,AMarkRegionToSplit);
					}

					// any other idea?

				}
			}

		} // for(int iFace=0; iFace<faces.size(); iFace++)
	}

}
/*----------------------------------------------------------------------------*/
void
MeshSplitter3Refinement::markNeighborRegionsInOutModel(
		const int& AMarkRegionToSplit)
{
	IGMesh::region_iterator it  = this->mesh_.regions_begin();

	std::vector<Region> intersectedRegions;
	std::map<gmds::TCellID, std::vector<gmds::math::Triangle> > in_out_triangles;

	std::vector<gmds::geom::GeomSurface* > surfaces;
	this->manager_.getSurfaces(surfaces);

	this->service_.initialization(surfaces);

	std::map<TCellID,bool> nodesAreInside;
	std::map<TCellID,bool> regionsAreInside;
	std::map<TCellID,bool> regionsAreIntersected;

	// used for progress display
	int iPercent = 0;
	unsigned int iCell = 0;

	// marking intersecting cells
	for(;!it.isDone();it.next()) {

		gmds::Region current_region = it.value();

		// display progression
		if(iCell == ((this->mesh_.getNbRegions()/10)*iPercent)) {
			std::cout<<"markNeighborRegionsInOutModel "<<iCell<<" of "<<this->mesh_.getNbRegions()<<" "<<iPercent*10<<"%"<<std::endl;
			iPercent++;
		}
		iCell++;

		if(this->mesh_.isMarked(current_region,AMarkRegionToSplit)) {
			continue;
		}

	    // check whether the hex is inside, outside or intersects the model
		std::vector<gmds::math::Triangle> intersectedTriangles;
		bool regionIsIntersected = this->service_.intersects(surfaces, current_region, intersectedTriangles);
		if(regionIsIntersected)
		{
			intersectedRegions.push_back(current_region);
			in_out_triangles[current_region.getID()] = intersectedTriangles;

			regionsAreIntersected[current_region.getID()] = true;
		} else {
			continue;
		}

		double ratio;
		bool isInsideRegion = this->service_.isMainlyInsideMC(
				in_out_triangles[current_region.getID()],current_region,ratio);

		regionsAreInside[current_region.getID()] = isInsideRegion;

		std::vector<Node> nodes = current_region.get<Node>();
		for(unsigned int iNode=0; iNode<nodes.size(); iNode++) {
			gmds::math::Point point = nodes[iNode].getPoint();
			bool isInsidePoint = this->service_.isInside(
					in_out_triangles[current_region.getID()],point);

			nodesAreInside[nodes[iNode].getID()] = isInsidePoint;
		}

	} // for(;it!=ite;it++){

	// for each couple of neighboring intersected regions
	// we will identify potential problems and mark them as refinement candidates
	for(unsigned int iRegion=0; iRegion<intersectedRegions.size(); iRegion++) {
		Region current_region = intersectedRegions[iRegion];

		std::vector<Node> nodes = current_region.get<Node>();

		std::vector<Face> faces = current_region.get<Face>();

		for(int iFace=0; iFace<faces.size(); iFace++) {

			Face current_face = faces[iFace];

			std::vector<Region> regions = current_face.get<Region>();
			if(regions.size() == 1) {
				continue;
			}

			Region neighborRegion;
			if(regions[0] == current_region) {
				neighborRegion = regions[1];
			} else {
				neighborRegion = regions[0];
			}

			// both regions are intersected
			if(regionsAreIntersected.find(neighborRegion.getID()) != regionsAreIntersected.end()) {

				// common face must be either wholly in or out
				std::vector<Node> nodes = current_face.get<Node>();

				bool isWhollyInOut = true;
				bool isIn = nodesAreInside[nodes[0].getID()];
				for(unsigned int iNode=1; iNode<nodes.size(); iNode++) {
					if(nodesAreInside[nodes[iNode].getID()] != isIn) {
						isWhollyInOut = false;
						break;
					}
				}

				if(!isWhollyInOut) {
					continue;
				} else {

					// face IN and both regions OUT
					if(isIn && (!regionsAreInside[current_region.getID()] && !regionsAreInside[neighborRegion.getID()])) {
						this->mesh_.mark(current_region,AMarkRegionToSplit);
						this->mesh_.mark(neighborRegion,AMarkRegionToSplit);
					}

					// face OUT and both regions IN
					if(!isIn && (regionsAreInside[current_region.getID()] && regionsAreInside[neighborRegion.getID()])) {
						this->mesh_.mark(current_region,AMarkRegionToSplit);
						this->mesh_.mark(neighborRegion,AMarkRegionToSplit);
					}

					// any other idea?

				}
			}

		} // for(int iFace=0; iFace<faces.size(); iFace++)
	}

}
/*----------------------------------------------------------------------------*/
void
MeshSplitter3Refinement::markRegionsIntersectedTwiceByModelAABBTree(
		const int& AMarkRegionToSplit)
{
	IGMesh::region_iterator it  = this->mesh_.regions_begin();

	std::vector<gmds::geom::GeomSurface* > surfaces;
	this->manager_.getSurfaces(surfaces);

	this->service_.initialization(surfaces);

	// first get the intersecting cells IDs bbox / triangles bbox pairs.
	std::map<gmds::TCellID, std::map<gmds::geom::GeomSurface*,std::vector<gmds::math::Triangle*> > > inOutSurfTriangles;

	GNode* boxTree = this->aabbRegionsTree_;
	this->service_.intersects(boxTree,inOutSurfTriangles);

	unsigned int nbIsOnIntersection = 0;

	std::map<gmds::TCellID, std::map<gmds::geom::GeomSurface*,std::vector<gmds::math::Triangle*> > >::iterator itCellID = inOutSurfTriangles.begin();

	for(; itCellID != inOutSurfTriangles.end(); itCellID++) {

		std::vector<gmds::geom::GeomSurface*> surfacesFound;

		gmds::Region current_region = this->mesh_.get<gmds::Region>(itCellID->first);
		bool isRegionIntersected = false;

		std::map<gmds::geom::GeomSurface*,std::vector<gmds::math::Triangle*> >::iterator itSurf = itCellID->second.begin();
		for(; itSurf!=itCellID->second.end(); itSurf++) {

			for(unsigned int iTriangle=0; iTriangle<itCellID->second[itSurf->first].size(); iTriangle++) {

				if(this->service_.intersects(*(itCellID->second[itSurf->first][iTriangle]),current_region)) {
					surfacesFound.push_back(itSurf->first);
					break;
				}
			}
		}

		for(unsigned int iSurface_found=0; iSurface_found<surfacesFound.size(); iSurface_found++) {
			for(unsigned int iSurface_found2=iSurface_found+1; iSurface_found2<surfacesFound.size(); iSurface_found2++) {
				if((this->surfacesNeighbors_[surfacesFound[iSurface_found]]).find(surfacesFound[iSurface_found2]) == (this->surfacesNeighbors_[surfacesFound[iSurface_found]]).end()) {
					this->mesh_.mark(current_region,AMarkRegionToSplit);
					break;
				}
			}
		}
	} // for(; itCellID != inOutSurfTriangles.end(); itCellID++)



//	// used for progress display
//	int iPercent = 0;
//	unsigned int iCell = 0;
//
//	// marking intersecting cells
//	for(;!it.isDone();it.next()) {
//
//		gmds::Region current_region = it.value();
//
//		// display progression
//		if(iCell == ((this->mesh_.getNbRegions()/10)*iPercent)) {
//			std::cout<<"markRegionsIsIntersectedTwiceByModel "<<iCell<<" of "<<this->mesh_.getNbRegions()<<" "<<iPercent*10<<"%"<<std::endl;
//			iPercent++;
//		}
//		iCell++;
//
//		if(this->mesh_.isMarked(current_region,AMarkRegionToSplit)) {
//			continue;
//		}
//
//		if(regionIsIntersectedTwiceByModel(current_region)) {
//			this->mesh_.mark(current_region,AMarkRegionToSplit);
//		}
//
//	} // for(;!it.isDone();it.next())



}
/*----------------------------------------------------------------------------*/
void
MeshSplitter3Refinement::markRegionsIntersectedTwiceByModel(
		const int& AMarkRegionToSplit)
{
	// used for progress display
	int iPercent = 0;
	unsigned int iCell = 0;

	// marking intersecting cells
	IGMesh::region_iterator it  = this->mesh_.regions_begin();

	for(;!it.isDone();it.next()) {

		gmds::Region current_region = it.value();

		// display progression
		if(iCell == ((this->mesh_.getNbRegions()/10)*iPercent)) {
			std::cout<<"markRegionsIsIntersectedTwiceByModel "<<iCell<<" of "<<this->mesh_.getNbRegions()<<" "<<iPercent*10<<"%"<<std::endl;
			iPercent++;
		}
		iCell++;

		if(this->mesh_.isMarked(current_region,AMarkRegionToSplit)) {
			continue;
		}

		if(regionIsIntersectedTwiceByModel(current_region)) {
			this->mesh_.mark(current_region,AMarkRegionToSplit);
		}

	} // for(;!it.isDone();it.next())

}
/*----------------------------------------------------------------------------*/
void
MeshSplitter3Refinement::markRegionsIntersectedByModel(
		const int& AMarkRegionToSplit)
{
	IGMesh::region_iterator it  = this->mesh_.regions_begin();

	std::vector<gmds::geom::GeomSurface* > surfaces;
	this->manager_.getSurfaces(surfaces);

	this->service_.initialization(surfaces);

	// used for progress display
	int iPercent = 0;
	unsigned int iCell = 0;

	// marking intersecting cells
	for(;!it.isDone();it.next()) {

		gmds::Region current_region = it.value();

		// display progression
		if(iCell == ((this->mesh_.getNbRegions()/10)*iPercent)) {
			std::cout<<"markRegionsIsIntersectedTwiceByModel "<<iCell<<" of "<<this->mesh_.getNbRegions()<<" "<<iPercent*10<<"%"<<std::endl;
			iPercent++;
		}
		iCell++;

		if(this->mesh_.isMarked(current_region,AMarkRegionToSplit)) {
			continue;
		}

		if(regionIsIntersectedByModel(current_region)) {
			this->mesh_.mark(current_region,AMarkRegionToSplit);
		}

	} // for(;!it.isDone();it.next())

}
/*----------------------------------------------------------------------------*/
void
MeshSplitter3Refinement::initialization()
{
    GSList* list = NULL;
	gmds::IGMesh::region_iterator itr  = this->mesh_.regions_begin();

	for(;!itr.isDone();itr.next()) {
		double minXYZ[3];
		double maxXYZ[3];

		gmds::Region current_region = itr.value();
		current_region.computeBoundingBox(minXYZ,maxXYZ);

		gpointer pointer = GINT_TO_POINTER(current_region.getID());
		GtsBBox* bbox = gts_bbox_new(
								gts_bbox_class (),
								pointer,
								minXYZ[0],minXYZ[1],minXYZ[2],
								maxXYZ[0],maxXYZ[1],maxXYZ[2]);

		list = g_slist_prepend(list,bbox);
	}
	this->aabbRegionsTree_ = gts_bb_tree_new(list);
}
/*----------------------------------------------------------------------------*/
void
MeshSplitter3Refinement::buildSurfacesNeighbors()
{
	std::vector<gmds::geom::GeomSurface*> surfaces;
	this->manager_.getSurfaces(surfaces);

	for(unsigned int iSurface=0; iSurface<surfaces.size(); iSurface++) {
		std::vector<gmds::geom::GeomCurve*> curves;
		surfaces[iSurface]->get(curves);

		for(unsigned int iCurve=0; iCurve<curves.size(); iCurve++) {
			std::vector<gmds::geom::GeomSurface*> surfaces_tmp;
			curves[iCurve]->get(surfaces_tmp);

			for(unsigned int iSurface_tmp=0; iSurface_tmp<surfaces_tmp.size(); iSurface_tmp++) {
				(surfacesNeighbors_[surfaces[iSurface]]).insert(surfaces_tmp[iSurface_tmp]);
			}
		}

	}

}
/*----------------------------------------------------------------------------*/
std::set<Region>
MeshSplitter3Refinement::getSplitRegions() const
{
	return splitRegions_;
}
/*----------------------------------------------------------------------------*/
std::set<Face>
MeshSplitter3Refinement::getSplitFaces() const
{
	return splitFaces_;
}
/*----------------------------------------------------------------------------*/
std::set<Edge>
MeshSplitter3Refinement::getSplitEdges() const
{
	return splitEdges_;
}
/*----------------------------------------------------------------------------*/
std::set<Region>
MeshSplitter3Refinement::getCreatedRegions() const
{
	return splitRegionsCreated_;
}
/*----------------------------------------------------------------------------*/
std::set<Face>
MeshSplitter3Refinement::getCreatedFaces() const
{
	return splitFacesCreated_;
}
/*----------------------------------------------------------------------------*/
std::set<Edge>
MeshSplitter3Refinement::getCreatedEdges() const
{
	return splitEdgesCreated_;
}
/*----------------------------------------------------------------------------*/
std::set<Node>
MeshSplitter3Refinement::getCreatedNodes() const
{
	return splitNodesCreated_;
}
/*----------------------------------------------------------------------------*/
} // end namespace gmds
/*----------------------------------------------------------------------------*/
