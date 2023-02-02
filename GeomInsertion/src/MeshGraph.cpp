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
/** \file    MeshGraph.cpp
 *  \author  legoff
 *  \date    06/11/2014
 */
/*----------------------------------------------------------------------------*/
#include <set>
/*----------------------------------------------------------------------------*/
#include <GMDS/MeshGraph.h>
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
MeshGraph::MeshGraph(IGMesh& AMesh, const int& AMark)
:mesh_(AMesh)
{
	mark_ = AMark;
}
/*----------------------------------------------------------------------------*/
MeshGraph::~MeshGraph()
{}
/*----------------------------------------------------------------------------*/
bool MeshGraph::findShortestPath
(Node& ANodeBegin, Node& ANodeEnd, std::vector<Edge>& APath)
{
	std::cout<<"findShortestPath"<<std::endl;

	APath.clear();

	// first check if a path exists
//	int markNodeIsCovered = this->mesh_.getNewMark();
//
//	this->mesh_.unmarkAll(markNodeIsCovered);
//	this->mesh_.freeMark(markNodeIsCovered);

	Variable<double>* nodeDistance = this->mesh_.newVariable<double>(GMDS_NODE,"nodeDistance");
	Variable<int>* nodeDone = this->mesh_.newVariable<int>(GMDS_NODE,"nodeDone");
	Variable<int>* nodeVisited = this->mesh_.newVariable<int>(GMDS_NODE,"nodeVisited");
	Variable<Node>* nodeAncestor = this->mesh_.newVariable<Node>(GMDS_NODE,"nodeAncestor");

	IGMesh::node_iterator it  = this->mesh_.nodes_begin();

	for(;!it.isDone();it.next()) {
		(*nodeDistance)[it.value().getID()] = HUGE_VALF;
		(*nodeDone)[it.value().getID()] = 0;
		(*nodeVisited)[it.value().getID()] = 0;
		(*nodeAncestor)[it.value().getID()] = Node();
	}

	Node current_node = ANodeBegin;
	(*nodeDistance)[ANodeBegin.getID()] = 0.;
	(*nodeDone)[ANodeBegin.getID()] = 1;
	(*nodeVisited)[ANodeBegin.getID()] = 1;

	std::set<Node> nodes_pool;

	while(true) {

		std::vector<Edge> edges = current_node.get<Edge>();
		std::set<Node> nodes_neighbors;

		for(unsigned int iEdge=0; iEdge<edges.size(); iEdge++) {
			// we only allow marked edges
			if(this->mesh_.isMarked(edges[iEdge],mark_)) {
				std::vector<Node> nodes = edges[iEdge].get<Node>();
				if(nodes[0] != current_node) {
					nodes_neighbors.insert(nodes[0]);
				} else {
					nodes_neighbors.insert(nodes[1]);
				}
			}
		}

		// each edge carries a value of 1. for the path cost; so we do not concern ourselves
		// with edges.
		std::set<Node>::iterator itn  = nodes_neighbors.begin();
		std::set<Node>::iterator itne = nodes_neighbors.end();
		for(; itn!=itne; itn++) {
			if((*nodeDone)[(*itn).getID()]) {

			} else {
				double node_distance_old = (*nodeDistance)[(*itn).getID()];
//				double node_distance_new = (*nodeDistance)[current_node->getID()] + 1.;
				double dist = current_node.getPoint().distance((*itn).getPoint());
				double node_distance_new = (*nodeDistance)[current_node.getID()] + dist;

				if(node_distance_new < node_distance_old) {
					(*nodeDistance)[(*itn).getID()] = node_distance_new;
					(*nodeAncestor)[(*itn).getID()] = current_node;
					(*nodeVisited)[(*itn).getID()] = 1;

					nodes_pool.insert(*itn);
				}
			}

		}

		(*nodeDone)[current_node.getID()] = 1;

		// select next node
		// we take the nearest from origin
		double minDist = HUGE_VALF;
		bool found = false;

		std::set<Node>::iterator itnn  = nodes_pool.begin();
		std::set<Node>::iterator itnne = nodes_pool.end();

		for(; itnn!=itnne; itnn++) {
			if(!(*nodeDone)[(*itnn).getID()]) {
				if((*nodeDistance)[(*itnn).getID()] < minDist) {
					current_node = *itnn;
					minDist = (*nodeDistance)[(*itnn).getID()];
					found = true;
				}
				//break;
			}
		}

//		if(itnn == itnne)
//			throw GMDSException("path not found.");
		if(!found)
			throw GMDSException("path not found.");

		if(current_node == ANodeEnd)
			break;
		else
			continue;
	}

	// here current_node is equal to ANodeEnd

	do {
		Node next_node = (*nodeAncestor)[current_node.getID()];

		std::vector<Edge> edges = current_node.get<Edge>();
		for(unsigned int iEdge=0; iEdge<edges.size(); iEdge++) {
			std::vector<Node> nodes = edges[iEdge].get<Node>();
			if((nodes[0] == next_node) || (nodes[1] == next_node)) {
				APath.push_back(edges[iEdge]);
			}

			current_node = next_node;
		}
	}
	while (current_node != ANodeBegin);

	this->mesh_.deleteVariable(GMDS_NODE,"nodeDistance");
	this->mesh_.deleteVariable(GMDS_NODE,"nodeDone");
	this->mesh_.deleteVariable(GMDS_NODE,"nodeVisited");
	this->mesh_.deleteVariable(GMDS_NODE,"nodeAncestor");

	return true;;
}
/*----------------------------------------------------------------------------*/
bool
MeshGraph::findShortestPathHaussdorf (
		double& ADistance,
		Node& ANodeBegin,
		Node& ANodeEnd,
		std::vector<Edge>& APath,
		gmds::geom::GeomCurve* ACurve)
{
	std::cout<<"findShortestPathHaussdorf"<<std::endl;

	APath.clear();

	Variable<double>* nodeDistance = this->mesh_.newVariable<double>(GMDS_NODE,"nodeDistance");
	Variable<int>* nodeDone = this->mesh_.newVariable<int>(GMDS_NODE,"nodeDone");
	Variable<int>* nodeVisited = this->mesh_.newVariable<int>(GMDS_NODE,"nodeVisited");
	Variable<Node>* nodeAncestor = this->mesh_.newVariable<Node>(GMDS_NODE,"nodeAncestor");

	IGMesh::node_iterator it  = this->mesh_.nodes_begin();

	for(;!it.isDone();it.next()) {
		(*nodeDistance)[it.value().getID()] = HUGE_VALF;
		(*nodeDone)[it.value().getID()] = 0;
		(*nodeVisited)[it.value().getID()] = 0;
		(*nodeAncestor)[it.value().getID()] = Node();
	}

	Node current_node = ANodeBegin;
	(*nodeDistance)[ANodeBegin.getID()] = 0.;
	(*nodeDone)[ANodeBegin.getID()] = 1;
	(*nodeVisited)[ANodeBegin.getID()] = 1;

	std::set<Node> nodes_pool;

	while(true) {

		std::vector<Edge> edges = current_node.get<Edge>();
		std::set<Node> nodes_neighbors;

		for(unsigned int iEdge=0; iEdge<edges.size(); iEdge++) {
			// we only allow marked edges
			if(this->mesh_.isMarked(edges[iEdge],mark_)) {
				std::vector<Node> nodes = edges[iEdge].get<Node>();
				if(nodes[0] != current_node) {
					nodes_neighbors.insert(nodes[0]);
				} else {
					nodes_neighbors.insert(nodes[1]);
				}
			}
		}

		// each edge carries a value of 1. for the path cost; so we do not concern ourselves
		// with edges.
		std::set<Node>::iterator itn  = nodes_neighbors.begin();
		std::set<Node>::iterator itne = nodes_neighbors.end();
		for(; itn!=itne; itn++) {
			if((*nodeDone)[(*itn).getID()]) {

			} else {
				double node_distance_old = (*nodeDistance)[(*itn).getID()];

				gmds::math::Segment segment(current_node.getPoint(),(*itn).getPoint());
				double dist = ACurve->computeDistanceHaussdorf(segment);

				double node_distance_new = (*nodeDistance)[current_node.getID()] + dist;
				//double node_distance_new = (*nodeDistance)[current_node.getID()];
				//if(node_distance_new<dist) {
				//	node_distance_new = dist;
				//}

				if(node_distance_new < node_distance_old) {
					(*nodeDistance)[(*itn).getID()] = node_distance_new;
					(*nodeAncestor)[(*itn).getID()] = current_node;
					(*nodeVisited)[(*itn).getID()] = 1;

					nodes_pool.insert(*itn);
				}
			}

		}

		(*nodeDone)[current_node.getID()] = 1;

		// select next node
		// we take the nearest from origin
		double minDist = HUGE_VALF;
		bool found = false;

		std::set<Node>::iterator itnn  = nodes_pool.begin();
		std::set<Node>::iterator itnne = nodes_pool.end();

		for(; itnn!=itnne; itnn++) {
			if(!(*nodeDone)[(*itnn).getID()]) {
				if((*nodeDistance)[(*itnn).getID()] < minDist) {
					current_node = *itnn;
					minDist = (*nodeDistance)[(*itnn).getID()];
					found = true;
				}
				//break;
			}
		}

		if(!found) {
//			throw GMDSException("path not found.");
			this->mesh_.deleteVariable(GMDS_NODE,"nodeDistance");
			this->mesh_.deleteVariable(GMDS_NODE,"nodeDone");
			this->mesh_.deleteVariable(GMDS_NODE,"nodeVisited");
			this->mesh_.deleteVariable(GMDS_NODE,"nodeAncestor");
			return false;
		}

		if(current_node == ANodeEnd) {
			break;
		} else {
			continue;
		}
	}

	// here current_node is equal to ANodeEnd

	ADistance = (*nodeDistance)[ANodeEnd.getID()];

	do {
		Node next_node = (*nodeAncestor)[current_node.getID()];

		std::vector<Edge> edges = current_node.get<Edge>();
		for(unsigned int iEdge=0; iEdge<edges.size(); iEdge++) {
			std::vector<Node> nodes = edges[iEdge].get<Node>();
			if((nodes[0] == next_node) || (nodes[1] == next_node)) {
				APath.push_back(edges[iEdge]);
			}

			current_node = next_node;
		}
	}
	while (current_node != ANodeBegin);

	this->mesh_.deleteVariable(GMDS_NODE,"nodeDistance");
	this->mesh_.deleteVariable(GMDS_NODE,"nodeDone");
	this->mesh_.deleteVariable(GMDS_NODE,"nodeVisited");
	this->mesh_.deleteVariable(GMDS_NODE,"nodeAncestor");

	return true;
}
/*----------------------------------------------------------------------------*/
bool
MeshGraph::findShortestPathHaussdorfSubCurve (
		double& ADistance,
		const int& ANbPoints,
		const int& AIPoint,
		Node& ANodeBegin,
		Node& ANodeEnd,
		std::vector<Edge>& APath,
		gmds::geom::GeomCurve* ACurve)
{
	std::cout<<"findShortestPathHaussdorfSubCurve"<<std::endl;

	APath.clear();

	Variable<double>* nodeDistance = this->mesh_.newVariable<double>(GMDS_NODE,"nodeDistance");
	Variable<int>* nodeDone = this->mesh_.newVariable<int>(GMDS_NODE,"nodeDone");
	Variable<int>* nodeVisited = this->mesh_.newVariable<int>(GMDS_NODE,"nodeVisited");
	Variable<Node>* nodeAncestor = this->mesh_.newVariable<Node>(GMDS_NODE,"nodeAncestor");

	IGMesh::node_iterator it  = this->mesh_.nodes_begin();

	for(;!it.isDone();it.next()) {
		(*nodeDistance)[it.value().getID()] = HUGE_VALF;
		(*nodeDone)[it.value().getID()] = 0;
		(*nodeVisited)[it.value().getID()] = 0;
		(*nodeAncestor)[it.value().getID()] = Node();
	}

	Node current_node = ANodeBegin;
	(*nodeDistance)[ANodeBegin.getID()] = 0.;
	(*nodeDone)[ANodeBegin.getID()] = 1;
	(*nodeVisited)[ANodeBegin.getID()] = 1;

	std::set<Node> nodes_pool;

	while(true) {

		std::vector<Edge> edges = current_node.get<Edge>();
		std::set<Node> nodes_neighbors;

		for(unsigned int iEdge=0; iEdge<edges.size(); iEdge++) {
			// we only allow marked edges
			if(this->mesh_.isMarked(edges[iEdge],mark_)) {
				std::vector<Node> nodes = edges[iEdge].get<Node>();
				if(nodes[0] != current_node) {
					nodes_neighbors.insert(nodes[0]);
				} else {
					nodes_neighbors.insert(nodes[1]);
				}
			}
		}

		// each edge carries a value of 1. for the path cost; so we do not concern ourselves
		// with edges.
		std::set<Node>::iterator itn  = nodes_neighbors.begin();
		std::set<Node>::iterator itne = nodes_neighbors.end();
		for(; itn!=itne; itn++) {
			if((*nodeDone)[(*itn).getID()]) {

			} else {
				double node_distance_old = (*nodeDistance)[(*itn).getID()];

				gmds::math::Segment segment(current_node.getPoint(),(*itn).getPoint());
				double dist = ACurve->computeDistanceHaussdorf(segment);
				//double dist = ACurve->computeDistanceHaussdorfSubCurve(ANbPoints,AIPoint,segment);

				//double node_distance_new = (*nodeDistance)[current_node.getID()] + dist;
				double node_distance_new = (*nodeDistance)[current_node.getID()];
				if(node_distance_new<dist)
					node_distance_new = dist;

				if(node_distance_new < node_distance_old) {
					(*nodeDistance)[(*itn).getID()] = node_distance_new;
					(*nodeAncestor)[(*itn).getID()] = current_node;
					(*nodeVisited)[(*itn).getID()] = 1;

					nodes_pool.insert(*itn);
				}
			}

		}

		(*nodeDone)[current_node.getID()] = 1;

		// select next node
		// we take the nearest from origin
		double minDist = HUGE_VALF;
		bool found = false;

		std::set<Node>::iterator itnn  = nodes_pool.begin();
		std::set<Node>::iterator itnne = nodes_pool.end();

		for(; itnn!=itnne; itnn++) {
			if(!(*nodeDone)[(*itnn).getID()]) {
				if((*nodeDistance)[(*itnn).getID()] < minDist) {
					current_node = *itnn;
					minDist = (*nodeDistance)[(*itnn).getID()];
					found = true;
				}
				//break;
			}
		}

		if(!found) {
//			throw GMDSException("path not found.");
			this->mesh_.deleteVariable(GMDS_NODE,"nodeDistance");
			this->mesh_.deleteVariable(GMDS_NODE,"nodeDone");
			this->mesh_.deleteVariable(GMDS_NODE,"nodeVisited");
			this->mesh_.deleteVariable(GMDS_NODE,"nodeAncestor");
			return false;
		}

		if(current_node == ANodeEnd) {
			break;
		} else {
			continue;
		}
	}

	// here current_node is equal to ANodeEnd
	ADistance = (*nodeDistance)[ANodeEnd.getID()];

	do {
		Node next_node = (*nodeAncestor)[current_node.getID()];

		std::vector<Edge> edges = current_node.get<Edge>();
		for(unsigned int iEdge=0; iEdge<edges.size(); iEdge++) {
			std::vector<Node> nodes = edges[iEdge].get<Node>();
			if((nodes[0] == next_node) || (nodes[1] == next_node)) {
				APath.push_back(edges[iEdge]);
			}

			current_node = next_node;
		}
	}
	while (current_node != ANodeBegin);

	this->mesh_.deleteVariable(GMDS_NODE,"nodeDistance");
	this->mesh_.deleteVariable(GMDS_NODE,"nodeDone");
	this->mesh_.deleteVariable(GMDS_NODE,"nodeVisited");
	this->mesh_.deleteVariable(GMDS_NODE,"nodeAncestor");

	return true;
}
/*----------------------------------------------------------------------------*/
//template<int TMask, typename TBase>
//bool MeshGraph<TMask,TBase>::findNearlyShortestPathLoop
//(Node* ANode, std::vector<Edge*>& APath, const TBase& AThreshold)
//{
//	std::cout<<"findShortestPath"<<std::endl;
//
//	APath.clear();
//
//	// first check if a path exists
////	int markNodeIsCovered = this->mesh_.getNewMark();
////
////	this->mesh_.unmarkAll(markNodeIsCovered);
////	this->mesh_.freeMark(markNodeIsCovered);
//
//	Variable<double>* nodeDistance = this->mesh_.template newVariable<double>(GMDS_NODE,"nodeDistance");
//	Variable<int>* nodeDone = this->mesh_.template newVariable<int>(GMDS_NODE,"nodeDone");
//	Variable<int>* nodeVisited = this->mesh_.template newVariable<int>(GMDS_NODE,"nodeVisited");
//	Variable<Node*>* nodeAncestor = this->mesh_.template newVariable<Node*>(GMDS_NODE,"nodeAncestor");
//
//	typename Mesh<TMask>::nodes_iterator it  = this->mesh_.nodes_begin();
//	typename Mesh<TMask>::nodes_iterator ite = this->mesh_.nodes_end();
//
//	for(;it!=ite;it++) {
//		(*nodeDistance)[(*it)->getID()] = HUGE_VALF;
//		(*nodeDone)[(*it)->getID()] = 0;
//		(*nodeVisited)[(*it)->getID()] = 0;
//		(*nodeAncestor)[(*it)->getID()] = NULL;
//	}
//
//	Node* current_node = ANode;
//	(*nodeDistance)[ANode->getID()] = 0.;
//	(*nodeDone)[ANode->getID()] = 1;
//	(*nodeVisited)[ANode->getID()] = 1;
//
//	std::set<Node*> nodes_pool;
//
//	while(true) {
//
//		std::cout<<current_node->getID()<<" "<<(*nodeDistance)[current_node->getID()]<<std::endl;
//
//		std::vector<Edge*> edges = current_node->getEdges();
//		std::set<Node*> nodes_neighbors;
//
//		for(unsigned int iEdge=0; iEdge<edges.size(); iEdge++) {
//			// we only allow marked edges
//			if(this->mesh_.isMarked(edges[iEdge],mark_)) {
//				std::vector<Node*> nodes = edges[iEdge]->getNodes();
//				if(nodes[0] != current_node) {
//					nodes_neighbors.insert(nodes[0]);
//				} else {
//					nodes_neighbors.insert(nodes[1]);
//				}
//			}
//		}
//
//		// each edge carries a value of 1. for the path cost; so we do not concern ourselves
//		// with edges.
//		std::set<Node*>::iterator itn  = nodes_neighbors.begin();
//		std::set<Node*>::iterator itne = nodes_neighbors.end();
//		for(; itn!=itne; itn++) {
//			if((*nodeDone)[(*itn)->getID()]) {
//
//			} else {
//				double node_distance_old = (*nodeDistance)[(*itn)->getID()];
//				double node_distance_new = (*nodeDistance)[current_node->getID()] + 1.;
//
//				if(node_distance_new < node_distance_old) {
//					(*nodeDistance)[(*itn)->getID()] = node_distance_new;
//					(*nodeAncestor)[(*itn)->getID()] = current_node;
//					(*nodeVisited)[(*itn)->getID()] = 1;
//
//					nodes_pool.insert(*itn);
//				}
//			}
//
//		}
//
//		(*nodeDone)[current_node->getID()] = 1;
//
//		// select next node
//		// we take the nearest from origin
//		double minDist = HUGE_VALF;
//		bool found = false;
//
//		std::set<Node*>::iterator itnn  = nodes_pool.begin();
//		std::set<Node*>::iterator itnne = nodes_pool.end();
//
//		for(; itnn!=itnne; itnn++) {
//			if(!(*nodeDone)[(*itnn)->getID()]) {
//				if((*nodeDistance)[(*itnn)->getID()] < minDist) {
//					current_node = *itnn;
//					minDist = (*nodeDistance)[(*itnn)->getID()];
//					found = true;
//				}
//				//break;
//			}
//		}
//
////		if(itnn == itnne)
////			throw GMDSException("path not found.");
//		if(!found)
//			throw GMDSException("path not found.");
//
//		if(current_node == ANode)
//			break;
//		else
//			continue;
//	}
//
//	// here current_node is equal to ANodeEnd
//
//	do {
//		Node* next_node = (*nodeAncestor)[current_node->getID()];
//
//		std::vector<Edge*> edges = current_node->getEdges();
//		for(unsigned int iEdge=0; iEdge<edges.size(); iEdge++) {
//			std::vector<Node*> nodes = edges[iEdge]->getNodes();
//			if((nodes[0] == next_node) || (nodes[1] == next_node)) {
//				APath.push_back(edges[iEdge]);
//			}
//
//			current_node = next_node;
//		}
//	}
//	while (current_node != ANode);
//
//	this->mesh_.deleteVariable(GMDS_NODE,"nodeDistance");
//	this->mesh_.deleteVariable(GMDS_NODE,"nodeDone");
//	this->mesh_.deleteVariable(GMDS_NODE,"nodeVisited");
//	this->mesh_.deleteVariable(GMDS_NODE,"nodeAncestor");
//
//	return true;;
//}
///*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/
