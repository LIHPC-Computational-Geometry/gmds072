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
/** \file    LaplacianSmoothing.cpp
 *  \author  legoff
 *  \date    02/12/2014
 */
/*----------------------------------------------------------------------------*/
// stl File Headers
/*----------------------------------------------------------------------------*/
#include <set>
/*----------------------------------------------------------------------------*/
// GMDS File Headers
/*----------------------------------------------------------------------------*/
#include <GMDS/LaplacianSmoothing.h>
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
LaplacianSmoothing::LaplacianSmoothing(IGMesh& AMesh,
		gmds::geom::GeomManager& AManager,
		gmds::geom::GeomMeshIntersectionService& AService)
:mesh_(AMesh), manager_(AManager), service_(AService)
{

}
/*----------------------------------------------------------------------------*/
LaplacianSmoothing::LaplacianSmoothing(const LaplacianSmoothing& ALaplacianSmoothing)
:mesh_(ALaplacianSmoothing.mesh_), manager_(ALaplacianSmoothing.manager_), service_(ALaplacianSmoothing.service_)
{

}
/*----------------------------------------------------------------------------*/
LaplacianSmoothing::~LaplacianSmoothing()
{

}
/*----------------------------------------------------------------------------*/
LaplacianSmoothing&
LaplacianSmoothing::operator=(const LaplacianSmoothing& ALaplacianSmoothing)
{
	return *this;
}
/*----------------------------------------------------------------------------*/
void
LaplacianSmoothing::initNodesAdjacencies()
{
	node2nodes_.clear();

	// we will build every fake edges from the regions in order to build node adjacency
	std::set<gmds::FakeEdge> allFakeEdges;

	gmds::IGMesh::region_iterator itr  = this->mesh_.regions_begin();

	for(;!itr.isDone();itr.next()) {
		Region current_region = itr.value();

		std::vector<gmds::FakeEdge> fakeEdges = current_region.getFakeEdges();
		for(int iEdge=0; iEdge<fakeEdges.size(); iEdge++) {
			allFakeEdges.insert(fakeEdges[iEdge]);
		}
	}

	std::set<gmds::FakeEdge>::iterator ite = allFakeEdges.begin();
	for(;ite != allFakeEdges.end(); ite++) {

		gmds::Node node0 = this->mesh_.get<gmds::Node>(ite->getID().getID1());
		gmds::Node node1 = this->mesh_.get<gmds::Node>(ite->getID().getID2());

		if(node2nodes_.find(node0) == node2nodes_.end()) {
			node2nodes_[node0] = std::vector<gmds::Node>();
		}
		node2nodes_[node0].push_back(node1);

		if(node2nodes_.find(node1) == node2nodes_.end()) {
			node2nodes_[node1] = std::vector<gmds::Node>();
		}
		node2nodes_[node1].push_back(node0);
	}
}
/*----------------------------------------------------------------------------*/
void
LaplacianSmoothing::exec(int AMarkFixedNodes)
{
	Variable<geom::GeomEntity* >* nodesClassification = this->mesh_.getGeometricClassification(0);
	Variable<geom::GeomEntity* >* edgesClassification = this->mesh_.getGeometricClassification(1);
	Variable<geom::GeomEntity* >* facesClassification = this->mesh_.getGeometricClassification(2);

	std::vector<gmds::geom::GeomPoint* > points;
	this->manager_.getPoints(points);
	std::vector<gmds::geom::GeomCurve* > curves;
	this->manager_.getCurves(curves);
	std::vector<gmds::geom::GeomSurface* > surfaces;
	this->manager_.getSurfaces(surfaces);

	gmds::IGMesh::node_iterator itn  = this->mesh_.nodes_begin();

	for(;!itn.isDone();itn.next()) {
		Node current_node = itn.value();

		// check whether the node can be moved or not
		if(AMarkFixedNodes != -1 && this->mesh_.isMarked(current_node,AMarkFixedNodes)) {
			continue;
		}

		if((*nodesClassification)[current_node.getID()] == NULL || 
		((*nodesClassification)[current_node.getID()] != NULL && (*nodesClassification)[current_node.getID()]->getDim() != 0)) {
			std::vector<gmds::Node> nodes = node2nodes_[current_node];

			gmds::math::Point point(0.,0.,0.);

			double sumFactor = 0.;
			for(int iNode=0; iNode<nodes.size(); iNode++) {

				double factor = 1.;
                                if((*nodesClassification)[nodes[iNode].getID()] != NULL && (*nodesClassification)[nodes[iNode].getID()]->getDim() == 0) {
                                        factor = 0.1;
                                } else if((*nodesClassification)[nodes[iNode].getID()] != NULL && (*nodesClassification)[nodes[iNode].getID()]->getDim() == 1) {
                                        factor = 0.1;
} else if((*nodesClassification)[nodes[iNode].getID()] != NULL && (*nodesClassification)[nodes[iNode].getID()]->getDim() == 2) {
                                        factor = 0.1;
}
                                        sumFactor += factor;

				//point = point + nodes[iNode].getPoint();
				point = point + factor * nodes[iNode].getPoint();
			}
			//point = point * (1. / nodes.size());
			point = point * (1. / sumFactor);

			if((*nodesClassification)[current_node.getID()] != NULL) {

				// we project if the node is associated to a curve or surface
				if((*nodesClassification)[current_node.getID()]->getDim() == 1) {
					gmds::geom::GeomCurve* curve = dynamic_cast<gmds::geom::GeomCurve*>((*nodesClassification)[current_node.getID()]);
					curve->project(point);
				}
				if((*nodesClassification)[current_node.getID()]->getDim() == 2) {
					gmds::geom::GeomSurface* surface = dynamic_cast<gmds::geom::GeomSurface*>((*nodesClassification)[current_node.getID()]);
					//surface->project(point);
                                        this->service_.project(surface,point);
				}
			}
			current_node.setPoint(point);

		} else {
			// it is on a GeomPoint

		}
	} // for(;!itn.isDone();itn.next()) {


}
/*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/
