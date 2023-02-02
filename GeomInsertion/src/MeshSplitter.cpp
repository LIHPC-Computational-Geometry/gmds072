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
/** \file    MeshSplitter.cpp
 *  \author  legoff
 *  \date    23/01/2015
 */
/*----------------------------------------------------------------------------*/
#include <GMDS/IG/IG.h>
#include <GMDS/MeshSplitter.h>

#include <GMDS/MeshModelAlgo.h>
/*----------------------------------------------------------------------------*/
namespace gmds {
/*----------------------------------------------------------------------------*/
MeshSplitter::MeshSplitter(
		gmds::IGMesh& AMesh,
		gmds::geom::GeomManager& AManager,
		gmds::geom::GeomMeshIntersectionService& AService)
:mesh_(AMesh), manager_(AManager), service_(AService)
{}
/*----------------------------------------------------------------------------*/
MeshSplitter::~MeshSplitter()
{}
/*----------------------------------------------------------------------------*/
void
MeshSplitter::splitMesh()
{
	std::cout<<"begin meshDoctor"<<std::endl;
	IGMeshDoctor doc(&(this->mesh_));
	doc.buildFacesAndR2F();
	doc.buildEdgesAndX2E();
	doc.updateUpwardConnectivity();
	std::cout<<"end meshDoctor"<<std::endl;

	int markRegionToSplit = this->mesh_.getNewMark<Region>();
	int markFaceToSplit = this->mesh_.getNewMark<Face>();
	int markEdgeToSplit = this->mesh_.getNewMark<Edge>();
	int markNodeToSplit = this->mesh_.getNewMark<Node>();

	std::cout<<"begin markRegionsToSplit"<<std::endl;
	markRegionsToSplit(markRegionToSplit,markFaceToSplit,markEdgeToSplit,markNodeToSplit);
	std::cout<<"end markRegionsToSplit"<<std::endl;

	std::cout<<"begin avoidBadConfigurations"<<std::endl;
	avoidBadConfigurations(markRegionToSplit,markFaceToSplit,markEdgeToSplit,markNodeToSplit);
	std::cout<<"end avoidBadConfigurations"<<std::endl;

	std::cout<<"begin split"<<std::endl;
	split(markRegionToSplit,markFaceToSplit,markEdgeToSplit,markNodeToSplit);
	std::cout<<"end split"<<std::endl;

	std::cout<<"begin clean mesh"<<std::endl;
	cleanMesh(markRegionToSplit,markFaceToSplit,markEdgeToSplit,markNodeToSplit);
	std::cout<<"end clean mesh"<<std::endl;

	this->mesh_.unmarkAll<Region>(markRegionToSplit);
	this->mesh_.unmarkAll<Face>(markFaceToSplit);
	this->mesh_.unmarkAll<Edge>(markEdgeToSplit);
	this->mesh_.unmarkAll<Node>(markNodeToSplit);
	this->mesh_.freeMark<Region>(markRegionToSplit);
	this->mesh_.freeMark<Face>(markFaceToSplit);
	this->mesh_.freeMark<Edge>(markEdgeToSplit);
	this->mesh_.freeMark<Node>(markNodeToSplit);
}
/*----------------------------------------------------------------------------*/
void
MeshSplitter::addEdges2Node(gmds::Node ANode, gmds::Face AFace)
{
	std::cout<<"addEdges2Node node "<<ANode.getID()<<" face "<<AFace.getID()<<std::endl; 

	std::cout<<"begin meshDoctor"<<std::endl;
        IGMeshDoctor doc(&(this->mesh_));
        doc.buildFacesAndR2F();
        doc.buildEdgesAndX2E();
        doc.updateUpwardConnectivity();
        std::cout<<"end meshDoctor"<<std::endl;

        int markRegionToSplit = this->mesh_.getNewMark<Region>();
        int markFaceToSplit = this->mesh_.getNewMark<Face>();
        int markEdgeToSplit = this->mesh_.getNewMark<Edge>();
        int markNodeToSplit = this->mesh_.getNewMark<Node>();

        std::cout<<"begin markRegionsToSplit"<<std::endl;
        //markRegionsToSplit(markRegionToSplit,markFaceToSplit,markEdgeToSplit,markNodeToSplit);
        markEntitiesToSplitForAddEdges(ANode,AFace,
			markRegionToSplit,markFaceToSplit,markEdgeToSplit,markNodeToSplit
	);
        std::cout<<"end markRegionsToSplit"<<std::endl;
	
	std::cout<<"begin avoidBadConfigurations"<<std::endl;
        avoidBadConfigurations(markRegionToSplit,markFaceToSplit,markEdgeToSplit,markNodeToSplit);
        std::cout<<"end avoidBadConfigurations"<<std::endl;

        std::cout<<"begin split"<<std::endl;
        split(markRegionToSplit,markFaceToSplit,markEdgeToSplit,markNodeToSplit);
        std::cout<<"end split"<<std::endl;

        std::cout<<"begin clean mesh"<<std::endl;
        cleanMesh(markRegionToSplit,markFaceToSplit,markEdgeToSplit,markNodeToSplit);
        std::cout<<"end clean mesh"<<std::endl;

        this->mesh_.unmarkAll<Region>(markRegionToSplit);
        this->mesh_.unmarkAll<Face>(markFaceToSplit);
        this->mesh_.unmarkAll<Edge>(markEdgeToSplit);
        this->mesh_.unmarkAll<Node>(markNodeToSplit);
        this->mesh_.freeMark<Region>(markRegionToSplit);
        this->mesh_.freeMark<Face>(markFaceToSplit);
        this->mesh_.freeMark<Edge>(markEdgeToSplit);
        this->mesh_.freeMark<Node>(markNodeToSplit);
}
/*----------------------------------------------------------------------------*/
void
MeshSplitter::refineQuads2EdgesOnCurve()
{
	std::cout<<"begin refineQuads2EdgesOnCurve" <<std::endl;

	Variable<geom::GeomEntity* >* nodeClassification   = this->mesh_.getGeometricClassification(0);
        Variable<geom::GeomEntity* >* edgeClassification   = this->mesh_.getGeometricClassification(1);
        Variable<geom::GeomEntity* >* faceClassification   = this->mesh_.getGeometricClassification(2);
        Variable<geom::GeomEntity* >* regionClassification = this->mesh_.getGeometricClassification(3);

	// first detect bad quads/nodes
	bool wasRefined = false;

	do {	
	wasRefined = false;
	gmds::IGMesh::face_iterator itf  = this->mesh_.faces_begin();
        for(;!itf.isDone() && !wasRefined; itf.next()) {

		gmds::Face current_face = itf.value();

		std::vector<TCellID> nodesIDs = current_face.getIDs<Node>();
		for(int iNode=0; iNode<nodesIDs.size(); iNode++) {
			if(((*nodeClassification)[nodesIDs[iNode]] != NULL) && ((*nodeClassification)[nodesIDs[iNode]]->getDim() == 1)) {
				if(((*nodeClassification)[nodesIDs[(iNode+1)%nodesIDs.size()]] == (*nodeClassification)[nodesIDs[iNode]]) &&  ((*nodeClassification)[nodesIDs[(iNode+2)%nodesIDs.size()]] == (*nodeClassification)[nodesIDs[iNode]])) {
					
					// check whether the face/regions are of bad quality
					std::vector<Region> regions = current_face.get<Region>();
					double minQuality = HUGE_VALF;
					for(int iRegion=0; iRegion<regions.size(); iRegion++) {
						double quality = regions[iRegion].computeScaledJacobian();
						if(quality<minQuality) {
							minQuality = quality;
						}
					}	
					if(minQuality<=0.) {
						Node node2Refine = this->mesh_.get<Node>(nodesIDs[(iNode+1)%nodesIDs.size()]);
						addEdges2Node(node2Refine,current_face);
						wasRefined = true;

						IGMeshDoctor doc(&(this->mesh_));
						doc.buildFacesAndR2F();
						doc.buildEdgesAndX2E();
    						doc.updateUpwardConnectivity();

						MeshModelAlgo meshModelAlgo(this->mesh_,this->manager_);
						meshModelAlgo.associateNodes();
						break;
					}
				}
			}
		}
	}
	std::cout<<"wasRefined "<<wasRefined<<std::endl;
	} while(wasRefined);
	std::cout<<"end refineQuads2EdgesOnCurve"<<std::endl;
}
/*----------------------------------------------------------------------------*/
void
MeshSplitter::refineQuads2EdgesOnCurveImproved()
{
	std::cout<<"begin refineQuads2EdgesOnCurveImproved" <<std::endl;

	std::cout<<"begin meshDoctor"<<std::endl;
        IGMeshDoctor doc(&(this->mesh_));
        doc.buildFacesAndR2F();
        doc.buildEdgesAndX2E();
        doc.updateUpwardConnectivity();
        std::cout<<"end meshDoctor"<<std::endl;

	int markRegionToSplit = this->mesh_.getNewMark<Region>();
	int markFaceToSplit = this->mesh_.getNewMark<Face>();
	int markEdgeToSplit = this->mesh_.getNewMark<Edge>();
	int markNodeToSplit = this->mesh_.getNewMark<Node>();

	Variable<geom::GeomEntity* >* nodeClassification   = this->mesh_.getGeometricClassification(0);
        Variable<geom::GeomEntity* >* edgeClassification   = this->mesh_.getGeometricClassification(1);
        Variable<geom::GeomEntity* >* faceClassification   = this->mesh_.getGeometricClassification(2);
        Variable<geom::GeomEntity* >* regionClassification = this->mesh_.getGeometricClassification(3);

	//bool wasRefined = false;
	bool nodesMarked2Refine = false;

	do {
	  //wasRefined = false;
	nodesMarked2Refine = false;
	//bool isConflictFree = true;

	gmds::IGMesh::face_iterator itf  = this->mesh_.faces_begin();
        for(; !itf.isDone(); itf.next()) {

		gmds::Face current_face = itf.value();

		// looking for a node to add edges to
		bool nodeFound = false;
		Node node2addEdges;
		Node node2Refine;
		
		std::vector<Edge> orderedEdges;
		current_face.getOrderedEdges(orderedEdges);
		
		for(int iEdge=0; iEdge<orderedEdges.size(); iEdge++) {
			if(((*edgeClassification)[orderedEdges[iEdge].getID()] != NULL) && ((*edgeClassification)[orderedEdges[iEdge].getID()]->getDim() == 1)) {
				if((*edgeClassification)[orderedEdges[(iEdge+1)%orderedEdges.size()].getID()] == (*edgeClassification)[orderedEdges[iEdge].getID()]) {
					// check the angle between the two edges; if near Pi then split
					std::vector<Node> nodes1 = orderedEdges[iEdge].get<Node>();
					std::vector<Node> nodes2 = orderedEdges[(iEdge+1)%orderedEdges.size()].get<Node>();
					
					gmds::Node node_middle;
					gmds::Node node_left;
					gmds::Node node_right;

					if(nodes1[0] == nodes2[0]) {
						node_middle = nodes1[0];
						node_left = nodes1[1];
						node_right = nodes2[1];
					} else if(nodes1[0] == nodes2[1]) {
                                                node_middle = nodes1[0];
                                                node_left = nodes1[1];
                                                node_right = nodes2[0];
                                        } else if(nodes1[1] == nodes2[0]) {
                                                node_middle = nodes1[1];
                                                node_left = nodes1[0];
                                                node_right = nodes2[1];
                                        } else {
                                                node_middle = nodes1[1];
                                                node_left = nodes1[0];
                                                node_right = nodes2[0];
					}
/*
					gmds::math::Vector v1(node_middle.getPoint(),node_left.getPoint());
					gmds::math::Vector v2(node_middle.getPoint(),node_right.getPoint());

					double angle = v1.angle(v2);
                                        if(fabs(angle) > gmds::math::Constants::PI*(7./8.)) {
                                                node2addEdges = node_middle;

						std::vector<Node> nodes = current_face.get<Node>();
						for(int iNode=0; iNode<nodes.size(); iNode++) {
							if(nodes[iNode] == node_middle) {
								node2Refine = nodes[(iNode+2)%nodes.size()];
								break;
							}
						}
                                                nodeFound = true;
                                                break;
                                        }			
*/
					// in this case we always split, whatever the angle
					node2addEdges = node_middle;
					std::vector<Node> nodes = current_face.get<Node>();
                                        for(int iNode=0; iNode<nodes.size(); iNode++) {
                                                if(nodes[iNode] == node_middle) {
                                                        node2Refine = nodes[(iNode+2)%nodes.size()];
                                                        break;
						}
					}
					nodeFound = true;
					break;
				}
			}
		} // for(int iEdge=0; iEdge<orderedEdges.size(); iEdge++ {

		if(nodeFound) {
		  // check that marking the node will not force us to mark other nodes to refine
		  // in order to meet the requirements of the template 
		  // cause that could force us to not add edges where we want them
		  /*	
		  // this would work if we had the N2R connectivity   
		  std::vector<Region> regions = node2Refine.get<Region>();
		  for(int iRegion=0; iRegion<regions.size(); iRegion++) {
		    if(this->mesh_.isMarked(regions[iRegion],markRegionToSplit)) {
		      isConflictFree = false;
		    }
		    }
		  */
		  bool isConflictFree = true;

		  gmds::IGMesh::region_iterator itr  = this->mesh_.regions_begin();
		  for(; !itr.isDone(); itr.next()) {
		    gmds::Region current_region = itr.value();

		    if(!isConflictFree) {
		      break;
		    }

		    if(this->mesh_.isMarked(current_region,markRegionToSplit)) {

		      std::vector<Node> nodes = current_region.get<Node>();

		      for(int iNode=0; iNode<nodes.size(); iNode++) {
			if(nodes[iNode] == node2Refine){
			  isConflictFree = false;
			}
		      }
		    }
		  }
		  if(isConflictFree) {
		    markEntitiesToSplitForAddEdges(node2addEdges,current_face,
						   markRegionToSplit,markFaceToSplit,markEdgeToSplit,markNodeToSplit
						   );
		    nodesMarked2Refine = true;
		  } else {
		    // in this case we stop looking for nodes to refine and launch the refinement
		    // in order to avoid conflict
		    //break;
		  }
		}
	} // for(; !itf.isDone(); itf.next())
	
	if(nodesMarked2Refine) {

	  std::cout<<"call to split"<<std::endl;
	  avoidBadConfigurations(markRegionToSplit,markFaceToSplit,markEdgeToSplit,markNodeToSplit);
	  split(markRegionToSplit,markFaceToSplit,markEdgeToSplit,markNodeToSplit);
	  cleanMesh(markRegionToSplit,markFaceToSplit,markEdgeToSplit,markNodeToSplit);
	  this->mesh_.unmarkAll<Region>(markRegionToSplit);
	  this->mesh_.unmarkAll<Face>(markFaceToSplit);
	  this->mesh_.unmarkAll<Edge>(markEdgeToSplit);
	  this->mesh_.unmarkAll<Node>(markNodeToSplit);
	  IGMeshDoctor doc(&(this->mesh_));
	  doc.buildFacesAndR2F();
	  doc.buildEdgesAndX2E();
	  doc.updateUpwardConnectivity();
	   
	  MeshModelAlgo meshModelAlgo(this->mesh_,this->manager_);
	  meshModelAlgo.associateNodes();
	}
			
	} while(nodesMarked2Refine);
	
	this->mesh_.unmarkAll<Region>(markRegionToSplit);
        this->mesh_.unmarkAll<Face>(markFaceToSplit);
        this->mesh_.unmarkAll<Edge>(markEdgeToSplit);
        this->mesh_.unmarkAll<Node>(markNodeToSplit);
        this->mesh_.freeMark<Region>(markRegionToSplit);
        this->mesh_.freeMark<Face>(markFaceToSplit);
        this->mesh_.freeMark<Edge>(markEdgeToSplit);
        this->mesh_.freeMark<Node>(markNodeToSplit);

	std::cout<<"end refineQuads2EdgesOnCurveImproved"<<std::endl;
}
/*----------------------------------------------------------------------------*/
void
MeshSplitter::cleanMesh(const int& AMarkRegionToSplit,
		const int& AMarkFaceToSplit,
		const int& AMarkEdgeToSplit,
		const int& AMarkNodeToSplit)
{
	{
		IGMesh::node_iterator itn  = this->mesh_.nodes_begin();
		for(;!itn.isDone(); itn.next()) {
			Node current_node = itn.value();
//
//			std::vector<Edge> edges = current_node.get<Edge>();
//			for(unsigned int iEdge=0; iEdge<edges.size(); iEdge++) {
//				if(this->mesh_.isMarked(edges[iEdge],AMarkEdgeToSplit)) {
//					current_node.remove<Edge>(edges[iEdge]);
//				}
//			}
//
//			std::vector<Face> faces = current_node.get<Face>();
//			for(unsigned int iFace=0; iFace<faces.size(); iFace++) {
//				if(this->mesh_.isMarked(faces[iFace],AMarkFaceToSplit)) {
//					current_node.remove<Face>(faces[iFace]);
//				}
//			}
			//current_node.removeAll<Edge>();
			//current_node.removeAll<Face>();
			//current_node.removeAll<Region>();
		}
	}

	{
		IGMesh::edge_iterator ite  = this->mesh_.edges_begin();
		for(;!ite.isDone(); ite.next()) {
			Edge current_edge = ite.value();
//
//			std::vector<Face> faces = current_edge.get<Face>();
//			for(unsigned int iFace=0; iFace<faces.size(); iFace++) {
//				if(this->mesh_.isMarked(faces[iFace],AMarkFaceToSplit)) {
//					current_edge.remove<Face>(faces[iFace]);
//				}
//			}
//
			//current_edge.removeAll<Face>();
			//current_edge.removeAll<Region>();
		}
	}

	{
		IGMesh::face_iterator itf  = this->mesh_.faces_begin();
		for(;!itf.isDone(); itf.next()) {
			Face current_face = itf.value();
//
//			std::vector<Edge> edges = current_face.get<Edge>();
//			for(unsigned int iEdge=0; iEdge<edges.size(); iEdge++) {
//				if(this->mesh_.isMarked(edges[iEdge],AMarkEdgeToSplit)) {
//					current_face.remove<Edge>(edges[iEdge]);
//				}
//			}
//
			current_face.removeAll<Region>();
		}
	}

//	{
//		IGMesh::edge_iterator ite  = this->mesh_.edges_begin();
//		for(;!ite.isDone(); ite.next()) {
//			Edge current_edge = ite.value();
//
//			if(this->mesh_.isMarked(current_edge,AMarkEdgeToSplit)) {
//				this->mesh_.deleteEdge(current_edge);
//			}
//		}
//	}

	{
		std::vector<TCellID> edges2DeleteIDs;

		IGMesh::edge_iterator ite  = this->mesh_.edges_begin();
		for(;!ite.isDone(); ite.next()) {
			Edge current_edge = ite.value();

			if(this->mesh_.isMarked(current_edge,AMarkEdgeToSplit)) {
//				this->mesh_.deleteEdge(current_edge);
				edges2DeleteIDs.push_back(current_edge.getID());
			}
		}

		for(int iEdge=0; iEdge<edges2DeleteIDs.size(); iEdge++) {
			this->mesh_.deleteEdge(edges2DeleteIDs[iEdge]);
		}
	}

	{
		std::vector<TCellID> faces2DeleteIDs;

		IGMesh::face_iterator itf  = this->mesh_.faces_begin();
		for(;!itf.isDone(); itf.next()) {
			Face current_face = itf.value();

			if(this->mesh_.isMarked(current_face,AMarkFaceToSplit)) {
//				this->mesh_.deleteFace(current_face);
				faces2DeleteIDs.push_back(current_face.getID());
			}
		}

		for(int iFace=0; iFace<faces2DeleteIDs.size(); iFace++) {
			this->mesh_.deleteFace(faces2DeleteIDs[iFace]);
		}
	}

//	// regions that had connectivities to removed edges/faces were removed, so
//	// we don't need to reset these connectivities.
//	{
//		IGMesh::region_iterator itr  = this->mesh_.regions_begin();
//		for(;!itr.isDone(); itr.next()) {
//			Region current_region = itr.value();
//
//			current_region.removeAll<Face>();
////			current_region->removeEdges();
//		}
//	}
//	{
//		typename Mesh<TMask>::nodes_iterator itn  = this->mesh_.nodes_begin();
//		for(;!itn->isDone(); itn->next()) {
//			Node* current_node = itn->currentItem();
//		}
//	}

}
/*----------------------------------------------------------------------------*/
} // end namespace gmds
/*----------------------------------------------------------------------------*/
