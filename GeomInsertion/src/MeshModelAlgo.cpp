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
/** \file    MeshModelAlgo.cpp
 *  \author  legoff
 *  \date    19/05/2015
 */
/*----------------------------------------------------------------------------*/
#include <GMDS/MeshModelAlgo.h>
/*----------------------------------------------------------------------------*/
// stl File Headers
/*----------------------------------------------------------------------------*/
#include <sstream>
#include <set>
#include <map>
/*----------------------------------------------------------------------------*/
// External File Headers
/*----------------------------------------------------------------------------*/
#include <GMDS/IG/IGMesh.h>
#include <GMDS/IO/VTKWriter.h>
#include <GMDS/CAD/GeomVolume.h>
#include <GMDS/CAD/GeomSurface.h>
#include <GMDS/CAD/GeomCurve.h>
#include <GMDS/CAD/GeomPoint.h>
#include <GMDS/IG/IGMeshQualityEvaluation.h>
/*----------------------------------------------------------------------------*/
namespace gmds {
/*----------------------------------------------------------------------------*/
MeshModelAlgo::MeshModelAlgo(
		IGMesh& AMesh,
		gmds::geom::GeomManager& AManager)
:mesh_(AMesh), manager_(AManager)
{}
/*----------------------------------------------------------------------------*/
MeshModelAlgo::~MeshModelAlgo()
{}
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
void
MeshModelAlgo::removeUnassociatedFaces()
{
	// first list all the faces to remove
	Variable<geom::GeomEntity* >* surfaceClassification= this->mesh_.getGeometricClassification(2);	

	std::set<gmds::Face> faces2Remove;
	gmds::IGMesh::face_iterator itf  = this->mesh_.faces_begin();
        for(;!itf.isDone();itf.next()) {
        	gmds::Face current_face = itf.value();
		
		if((*surfaceClassification)[current_face.getID()] == NULL) {
			faces2Remove.insert(current_face);
		}
	}


	// remove adjacencies that contain any of those faces
	if(this->mesh_.getModel().has(R2F)) {
		gmds::IGMesh::region_iterator itr  = this->mesh_.regions_begin();
	        for(;!itr.isDone();itr.next()) {
        	        gmds::Region current_region = itr.value();
			
			std::vector<gmds::Face> faces = current_region.get<gmds::Face>();
			for(int iFace=0; iFace<faces.size(); iFace++) {
				if(faces2Remove.find(faces[iFace]) != faces2Remove.end()) {
					current_region.remove<gmds::Face>(faces[iFace]);
				}
			}
		}
	}	
	
	if(this->mesh_.getModel().has(F2F)) {
                gmds::IGMesh::face_iterator itf  = this->mesh_.faces_begin();
                for(;!itf.isDone();itf.next()) {
                        gmds::Face current_face = itf.value();

                        std::vector<gmds::Face> faces = current_face.get<gmds::Face>();
                        for(int iFace=0; iFace<faces.size(); iFace++) {
                                if(faces2Remove.find(faces[iFace]) != faces2Remove.end()) {
                                        current_face.remove<gmds::Face>(faces[iFace]);
                                }
                        }
                }
        }

	if(this->mesh_.getModel().has(E2F)) {
                gmds::IGMesh::edge_iterator ite  = this->mesh_.edges_begin();
                for(;!ite.isDone();ite.next()) {
                        gmds::Edge current_edge = ite.value();

                        std::vector<gmds::Face> faces = current_edge.get<gmds::Face>();
                        for(int iFace=0; iFace<faces.size(); iFace++) {
                                if(faces2Remove.find(faces[iFace]) != faces2Remove.end()) {
                                        current_edge.remove<gmds::Face>(faces[iFace]);
                                }
                        }
                }
        }	

	if(this->mesh_.getModel().has(N2F)) {
                gmds::IGMesh::node_iterator itn  = this->mesh_.nodes_begin();
                for(;!itn.isDone();itn.next()) {
                        gmds::Node current_node = itn.value();

                        std::vector<gmds::Face> faces = current_node.get<gmds::Face>();
                        for(int iFace=0; iFace<faces.size(); iFace++) {
                                if(faces2Remove.find(faces[iFace]) != faces2Remove.end()) {
                                        current_node.remove<gmds::Face>(faces[iFace]);
                                }
                        }
                }
        }

	// finally delete the faces
	std::set<gmds::Face>::iterator it = faces2Remove.begin();
	for(; it!=faces2Remove.end(); it++) {
		this->mesh_.deleteFace(*it);
	}
}
/*----------------------------------------------------------------------------*/
void
MeshModelAlgo::displayMeshQuality()
{
	IGMeshQualityEvaluation qualEval;

	// scaled jacobian
	Variable<double>* scaledJacobianRegionVar;
        if(this->mesh_.doesVariableExist(GMDS_REGION,"scaledJacobianRegionVar")) {
                scaledJacobianRegionVar = this->mesh_.getVariable<double>(GMDS_REGION,"scaledJacobianRegionVar");
        } else {
                scaledJacobianRegionVar = this->mesh_.newVariable<double>(GMDS_REGION,"scaledJacobianRegionVar");
        }
	
	double minScaledJacobian =  HUGE_VALF;
	double maxScaledJacobian = -HUGE_VALF;
	double meanScaledJacobian = 0;

        gmds::IGMesh::region_iterator itr  = this->mesh_.regions_begin();
        for(;!itr.isDone();itr.next()) {
        	gmds::Region current_region = itr.value();

		double scaledJacobian = current_region.computeScaledJacobian();
		//double scaledJacobian = qualEval.scaledJacobian(current_region);
		if(scaledJacobian < minScaledJacobian) minScaledJacobian = scaledJacobian;
		if(scaledJacobian > maxScaledJacobian) maxScaledJacobian = scaledJacobian;
		meanScaledJacobian += scaledJacobian;

		(*scaledJacobianRegionVar)[current_region.getID()] = scaledJacobian;
	}
	meanScaledJacobian = meanScaledJacobian / this->mesh_.getNbRegions();

	std::cout<<"MeshQuality scaled jacobian min "<<minScaledJacobian<<" max "<<maxScaledJacobian<<" mean "<<meanScaledJacobian<<std::endl;
}
/*----------------------------------------------------------------------------*/
void
MeshModelAlgo::associateNodes()
{
	// mesh model requirements N|E|F|R|E2N|F2N|R2N
	if(! (mesh_.getModel().has(N) && mesh_.getModel().has(E) && mesh_.getModel().has(F) && mesh_.getModel().has(R)
	&& mesh_.getModel().has(E2N) && mesh_.getModel().has(F2N) && mesh_.getModel().has(R2N))) {
		throw GMDSException("MeshModelAlgo::associateNodes mesh model does meet the requirements N|E|F|R|E2N|F2N|R2N.");
	}

	std::vector<gmds::geom::GeomPoint*> points;
	std::vector<gmds::geom::GeomCurve*> curves;
	std::vector<gmds::geom::GeomSurface*> surfaces;
	std::vector<gmds::geom::GeomVolume*> volumes;

	//AVol->get(points);
	//AVol->get(curves);
	//AVol->get(surfaces);

	Variable<geom::GeomEntity* >* nodeClassification   = this->mesh_.getGeometricClassification(0);
	Variable<geom::GeomEntity* >* edgeClassification   = this->mesh_.getGeometricClassification(1);
	Variable<geom::GeomEntity* >* faceClassification   = this->mesh_.getGeometricClassification(2);
	Variable<geom::GeomEntity* >* regionClassification = this->mesh_.getGeometricClassification(3);

	gmds::IGMesh::region_iterator itr  = this->mesh_.regions_begin();
	for(;!itr.isDone(); itr.next()) {
		gmds::Region current_region = itr.value();

		if((*regionClassification)[current_region.getID()]  != NULL) {

			std::vector<gmds::Node> nodes = current_region.get<Node>();

			for(int iNode=0; iNode<nodes.size(); iNode++) {
				if((*nodeClassification)[nodes[iNode].getID()] == NULL || (*nodeClassification)[nodes[iNode].getID()]->getDim()  != 0) {
					(*nodeClassification)[nodes[iNode].getID()] = (*regionClassification)[current_region.getID()];
				}
			}
		}
	}

	gmds::IGMesh::face_iterator itf  = this->mesh_.faces_begin();
        for(;!itf.isDone(); itf.next()) {
                gmds::Face current_face = itf.value();

                if((*faceClassification)[current_face.getID()]  != NULL) {

                        std::vector<gmds::Node> nodes = current_face.get<Node>();

                        for(int iNode=0; iNode<nodes.size(); iNode++) {
                                if((*nodeClassification)[nodes[iNode].getID()] == NULL || (*nodeClassification)[nodes[iNode].getID()]->getDim()  != 0) {
                                        (*nodeClassification)[nodes[iNode].getID()] = (*faceClassification)[current_face.getID()];
                                }
                        }
                }
        }

	gmds::IGMesh::edge_iterator ite  = this->mesh_.edges_begin();
        for(;!ite.isDone(); ite.next()) {
                gmds::Edge current_edge = ite.value();

                if((*edgeClassification)[current_edge.getID()]  != NULL) {

                        std::vector<gmds::Node> nodes = current_edge.get<Node>();

                        for(int iNode=0; iNode<nodes.size(); iNode++) {
                                if((*nodeClassification)[nodes[iNode].getID()] == NULL || (*nodeClassification)[nodes[iNode].getID()]->getDim()  != 0) {
                                        (*nodeClassification)[nodes[iNode].getID()] = (*edgeClassification)[current_edge.getID()];
                                }
                        }
                }
        }
}
/*----------------------------------------------------------------------------*/
void
MeshModelAlgo::markBoundaryNodes(int AMarkBoundaryNodes)
{
	gmds::IGMesh::face_iterator itf  = this->mesh_.faces_begin();
        for(;!itf.isDone(); itf.next()) {
                gmds::Face current_face = itf.value();

		if((current_face.get<Region>()).size() == 1) {
			std::vector<Node> nodes = current_face.get<Node>();
			for(int iNode=0; iNode<nodes.size(); iNode++) {
				this->mesh_.mark(nodes[iNode],AMarkBoundaryNodes);
			}
		} else {
			if((current_face.get<Region>()).size() == 0) {
				throw GMDSException("MeshModelAlgo::markBoundaryNodes missing F2R relation.");
			}
		}
		
	}
}
/*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/
