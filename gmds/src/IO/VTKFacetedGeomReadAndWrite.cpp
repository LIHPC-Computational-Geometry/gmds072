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
/*
 * VTKFacetedGeomReader.cpp
 *
 *  Created on: 29 août 2014
 *      Author: ledouxf
 */
/*----------------------------------------------------------------------------*/
#include "GMDS/IO/VTKFacetedGeomReadAndWrite.h"
#include "GMDS/IO/VTKWriter.h"
#include "GMDS/IO/VTKReader.h"
#include "GMDS/Algo/BoundaryOperator.h"
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
void VTKFacetedGeomReadAndWrite::build(
                         geom::FacetedGeomManager& AGeomMng,
                         const std::string&        AFile,
			 const bool ASingleSurface)
{
	std::cout<<"VTKFacetedGeomReadAndWrite::build begin"<<std::endl;

        // we build a mesh from the file
	gmds::MeshModel mod = gmds::DIM3|gmds::N|gmds::F|gmds::R|gmds::F2N|gmds::R2N
					|gmds::F2R|gmds::R2F;
	gmds::IGMesh mesh(mod);

        VTKReader<IGMesh> r(mesh);
        r.read(AFile);

	// we only keep the skin of the mesh
	std::cout<<"VTKFacetedGeomReadAndWrite::build begin IGMeshDoctor"<<std::endl;
	gmds::IGMeshDoctor doc(&mesh);
	doc.buildFacesAndR2F();
	//doc.buildEdgesAndX2E();
	doc.updateUpwardConnectivity();	
	std::cout<<"VTKFacetedGeomReadAndWrite::build end IGMeshDoctor"<<std::endl;

	//BoundaryOperator boundaryOp(&mesh);
	//if (!boundaryOp.isValid()) {
	//	std::cout << "Invalid model for boundary operations" << std::endl;
	//	throw GMDSException("VTKFacetedGeomReadAndWrite::build : Invalid model for boundary operations");
	//}
	//int markFaceSurf;
	//int markEdgeSurf;
	//int markNodeSurf;
	//boundaryOp.markCellsOnSurfaces(markFaceSurf,markEdgeSurf,markNodeSurf);

	int markFaceSurf = mesh.getNewMark<gmds::Face>();
	int markNodeSurf = mesh.getNewMark<gmds::Node>();
    
    gmds::IGMesh::face_iterator itf  = mesh.faces_begin();
    for(;!itf.isDone();itf.next()) {
        gmds::Face current_face = itf.value();
        
        std::vector<gmds::Region> regions = current_face.get<gmds::Region>();
        if(regions.size() == 1) {
            std::vector<gmds::Node> nodes = current_face.get<gmds::Node>();
            
            for(unsigned int iNode=0; iNode<nodes.size(); iNode++) {
                mesh.mark(nodes[iNode],markNodeSurf);
            }
            
            bool isFaceOutward = regions[0].isFaceOrientedOutward(nodes);
            if(isFaceOutward) {
                if(current_face.getType() == GMDS_TRIANGLE) {
                    mesh.mark(current_face,markFaceSurf);
                } else {
                    if(current_face.getType() == GMDS_QUAD) {
                        gmds::Face tri1 = mesh.newTriangle(nodes[0],nodes[1],nodes[2]);
                        gmds::Face tri2 = mesh.newTriangle(nodes[0],nodes[2],nodes[3]);
                        mesh.mark(tri1,markFaceSurf);
                        mesh.mark(tri2,markFaceSurf);
                    } else {
                        throw GMDSException("VTKFacetedGeomReadAndWrite::build : not implemented for non tri/quad.");
                    }
                }
            } else {
                if(current_face.getType() == GMDS_TRIANGLE) {
                    std::vector<gmds::Node> nodes_tmp;
                    nodes_tmp.push_back(nodes[2]);
                    nodes_tmp.push_back(nodes[1]);
                    nodes_tmp.push_back(nodes[0]);
                    current_face.set<gmds::Node>(nodes_tmp);
                    
                    mesh.mark(current_face,markFaceSurf);
                } else {
                    if(current_face.getType() == GMDS_QUAD) {
                        gmds::Face tri1 = mesh.newTriangle(nodes[0],nodes[3],nodes[2]);
                        gmds::Face tri2 = mesh.newTriangle(nodes[0],nodes[2],nodes[1]);
                        mesh.mark(tri1,markFaceSurf);
                        mesh.mark(tri2,markFaceSurf);
                    } else {
                        throw GMDSException("VTKFacetedGeomReadAndWrite::build : not implemented for non tri/quad.");
                    }
                }
                
            }
        }
    }


	//gmds::MeshModel mod2 = gmds::DIM3|gmds::N|gmds::E|gmds::F|gmds::F2N|gmds::E2N|gmds::F2E;
	//gmds::IGMesh mesh2(mod2);
	gmds::IGMesh& mesh2 = AGeomMng.getMeshView();
	
	// copy relevant (boundary) nodes and faces (triangles) to mesh2
	{
		std::map<TCellID,TCellID> oldNodde2NewNode;

		gmds::IGMesh::node_iterator itn  = mesh.nodes_begin();
        	for(;!itn.isDone();itn.next()) {
                	gmds::Node current_node = itn.value();
		
			if(mesh.isMarked(current_node,markNodeSurf)) {
				gmds::Node newNode = mesh2.newNode(current_node.getPoint());
				oldNodde2NewNode[current_node.getID()] = newNode.getID();
			}
		}

		gmds::IGMesh::face_iterator itf  = mesh.faces_begin();
                for(;!itf.isDone();itf.next()) {
                        gmds::Face current_face = itf.value();

			if(mesh.isMarked(current_face,markFaceSurf)) {
				std::vector<TCellID> nodesIDs = current_face.getIDs<gmds::Node>();

				mesh2.newTriangle(oldNodde2NewNode[nodesIDs[0]],oldNodde2NewNode[nodesIDs[1]],oldNodde2NewNode[nodesIDs[2]]);
			}
		}

		std::cout<<"mesh2.getNbNodes() "<<mesh2.getNbNodes()<<std::endl;
		std::cout<<"mesh2.getNbTriangles() "<<mesh2.getNbTriangles()<<std::endl;
	}

	mesh.unmarkAll<gmds::Face>(markFaceSurf);
	mesh.unmarkAll<gmds::Node>(markNodeSurf);
	mesh.freeMark<gmds::Face>(markFaceSurf);
	mesh.freeMark<gmds::Node>(markNodeSurf);

        AGeomMng.buildFromMesh(ASingleSurface);
}
/*----------------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------------*/
