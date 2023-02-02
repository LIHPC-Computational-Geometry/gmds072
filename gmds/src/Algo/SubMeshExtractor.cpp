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
/** \file    SubMeshExtractor.cpp
 *  \author  F. LEDOUX
 *  \date    08/31/2010
 */
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
#include <GMDS/Algo/SubMeshExtractor.h>
/*----------------------------------------------------------------------------*/
//
//SubMeshExtractor::Extractor(IGMesh& AMesh)
//:mesh_(AMesh)
//{}
///*----------------------------------------------------------------------------*/
//
//SubMeshExtractor::~Extractor()
//{}
///*----------------------------------------------------------------------------*/
//
//void SubMeshExtractor::
//getSubMeshRepartitionRemove(const int AMark, const TInt PartID,
//		Variable<TCellID>* AVarNodeLID, Variable<TCellID>* AVarEdgeLID,
//		Variable<TCellID>* AVarFaceLID, Variable<TCellID>* AVarRegionLID,
//		const Variable<TInt>* AVarNodeOwner,
//		const Variable<TInt>* AVarEdgeOwner,
//		const Variable<TInt>* AVarFaceOwner,
//		const Variable<TInt>* AVarRegionOwner,
//		const int AMarkShared)
//{
//	// TODO here we need to fill the masters table for those which migrated to submeshes.
//	std::map<id,std::pair<TInt,TCellID> > nodesMasterMigration;
//	std::map<id,std::pair<TInt,TCellID> > edgesMasterMigration;
//	std::map<id,std::pair<TInt,TCellID> > facesMasterMigration;
//	std::map<id,std::pair<TInt,TCellID> > regionsMasterMigration;
//
//	if(MeshDescriptor<M>::hasRegions){
//		typename IGMesh::regions_iterator itRegion  = mesh_.regions_begin();
//		typename IGMesh::regions_iterator itRegione = mesh_.regions_end();
//
//		for(;itRegion!=itRegione;itRegion++){
//			if(mesh_.isMaster(*itRegion) && (*AVarRegionOwner)[(*itRegion)->getID()] != PartID){
//				std::pair<TInt,TCellID> p((*AVarRegionOwner)[(*itRegion)->getID()],(*AVarRegionLID)[(*itRegion)->getID()]);
//				regionsMasterMigration[(*itRegion)->getID()] = p;
//			}
//		}
//	}
//
//	if(MeshDescriptor<M>::hasFaces){
//		typename IGMesh::faces_iterator itFace  = mesh_.faces_begin();
//		typename IGMesh::faces_iterator itFacee = mesh_.faces_end();
//
//		for(;itFace!=itFacee;itFace++){
//			if(mesh_.isMaster(*itFace) && (*AVarFaceOwner)[(*itFace)->getID()] != PartID){
//				std::pair<TInt,TCellID> p((*AVarFaceOwner)[(*itFace)->getID()],(*AVarFaceLID)[(*itFace)->getID()]);
//				facesMasterMigration[(*itFace)->getID()] = p;
//			}
//		}
//	}
//
//	if(MeshDescriptor<M>::hasEdges){
//		typename IGMesh::edges_iterator itEdge  = mesh_.edges_begin();
//		typename IGMesh::edges_iterator itEdgee = mesh_.edges_end();
//
//		for(;itEdge!=itEdgee;itEdge++){
//			if(mesh_.isMaster(*itEdge) && (*AVarEdgeOwner)[(*itEdge)->getID()] != PartID){
//				std::pair<TInt,TCellID> p((*AVarEdgeOwner)[(*itEdge)->getID()],(*AVarEdgeLID)[(*itEdge)->getID()]);
//				edgesMasterMigration[(*itEdge)->getID()] = p;
//			}
//		}
//	}
//
//	if(MeshDescriptor<M>::hasNodes){
//		typename IGMesh::nodes_iterator itNode  = mesh_.nodes_begin();
//		typename IGMesh::nodes_iterator itNodee = mesh_.nodes_end();
//
//		for(;itNode!=itNodee;itNode++){
//			if(mesh_.isMaster(*itNode) && (*AVarNodeOwner)[(*itNode)->getID()] != PartID){
//				std::pair<TInt,TCellID> p((*AVarNodeOwner)[(*itNode)->getID()],(*AVarNodeLID)[(*itNode)->getID()]);
//				nodesMasterMigration[(*itNode)->getID()] = p;
//			}
//		}
//	}
//
//	// first pass, in order to update the connectivities between entities, removing those
//	// that will be removed.
//	if(MeshDescriptor<M>::hasRegions){
//		typename IGMesh::regions_iterator itRegion  = mesh_.regions_begin();
//		typename IGMesh::regions_iterator itRegione = mesh_.regions_end();
//
//		for(;itRegion!=itRegione;itRegion++){
//			if(mesh_.isMarked(*itRegion,AMark)){
//
//				/* connection R->R*/
//				if(MeshDescriptor<M>::RKnowR){
//					std::vector<Region*> regions = (*itRegion)->getRegions();
//					std::vector<Region*>::iterator itR  = regions.begin();
//					std::vector<Region*>::iterator itRe = regions.end();
//					for(;itR!=itRe;itR++){
//						if(!mesh_.isMarked(*itR,AMark)){
//							(*itRegion)->removeRegion(*itR);
//						}
//					}
//				}
//
//				/* connection R->F*/
//				if(MeshDescriptor<M>::RKnowF){
//					std::vector<Face*> faces = (*itRegion)->getFaces();
//					std::vector<Face*>::iterator itF  = faces.begin();
//					std::vector<Face*>::iterator itFe = faces.end();
//					for(;itF!=itFe;itF++){
//						if(!mesh_.isMarked(*itF,AMark)){
//							(*itRegion)->removeFace(*itF);
//						}
//					}
//				}
//
//				/* connection R->E*/
//				if(MeshDescriptor<M>::RKnowE){
//					std::vector<Edge*> edges = (*itRegion)->getEdges();
//					std::vector<Edge*>::iterator itE  = edges.begin();
//					std::vector<Edge*>::iterator itEe = edges.end();
//					for(;itE!=itEe;itE++){
//						if(!mesh_.isMarked(*itE,AMark)){
//							(*itRegion)->removeEdge(*itE);
//						}
//					}
//				}
//
//				/* connection R->N*/
//				if(MeshDescriptor<M>::RKnowN){
//					std::vector<Node*> nodes = (*itRegion)->getNodes();
//					std::vector<Node*>::iterator itN  = nodes.begin();
//					std::vector<Node*>::iterator itNe = nodes.end();
//					for(;itN!=itNe;itN++){
//						if(!mesh_.isMarked(*itN,AMark)){
//							(*itRegion)->removeNode(*itN);
//						}
//					}
//				}
//
//			}
//		}
//	} // if(MeshDescriptor<M>::hasRegions)
//
//	if(MeshDescriptor<M>::hasFaces){
//		typename IGMesh::faces_iterator itFace  = mesh_.faces_begin();
//		typename IGMesh::faces_iterator itFacee = mesh_.faces_end();
//
//		for(;itFace!=itFacee;itFace++){
//			if(mesh_.isMarked(*itFace,AMark)){
//
//				/* connection F->R*/
//				if(MeshDescriptor<M>::FKnowR){
//					std::vector<Region*> regions = (*itFace)->getRegions();
//					std::vector<Region*>::iterator itR  = regions.begin();
//					std::vector<Region*>::iterator itRe = regions.end();
//					for(;itR!=itRe;itR++){
//						if(!mesh_.isMarked(*itR,AMark)){
//							(*itFace)->removeRegion(*itR);
//						}
//					}
//				}
//
//				/* connection F->F*/
//				if(MeshDescriptor<M>::FKnowF){
//					std::vector<Face*> faces = (*itFace)->getFaces();
//					std::vector<Face*>::iterator itF  = faces.begin();
//					std::vector<Face*>::iterator itFe = faces.end();
//					for(;itF!=itFe;itF++){
//						if(!mesh_.isMarked(*itF,AMark)){
//							(*itFace)->removeFace(*itF);
//						}
//					}
//				}
//
//				/* connection F->E*/
//				if(MeshDescriptor<M>::FKnowE){
//					std::vector<Edge*> edges = (*itFace)->getEdges();
//					std::vector<Edge*>::iterator itE  = edges.begin();
//					std::vector<Edge*>::iterator itEe = edges.end();
//					for(;itE!=itEe;itE++){
//						if(!mesh_.isMarked(*itE,AMark)){
//							(*itFace)->removeEdge(*itE);
//						}
//					}
//				}
//
//				/* connection F->N*/
//				if(MeshDescriptor<M>::FKnowN){
//					std::vector<Node*> nodes = (*itFace)->getNodes();
//					std::vector<Node*>::iterator itN  = nodes.begin();
//					std::vector<Node*>::iterator itNe = nodes.end();
//					for(;itN!=itNe;itN++){
//						if(!mesh_.isMarked(*itN,AMark)){
//							(*itFace)->removeNode(*itN);
//						}
//					}
//				}
//
//			}
//		}
//	} // if(MeshDescriptor<M>::hasFaces)
//
//	if(MeshDescriptor<M>::hasEdges){
//		typename IGMesh::edges_iterator itEdge  = mesh_.edges_begin();
//		typename IGMesh::edges_iterator itEdgee = mesh_.edges_end();
//
//		for(;itEdge!=itEdgee;itEdge++){
//			if(mesh_.isMarked(*itEdge,AMark)){
//
//				/* connection E->R*/
//				if(MeshDescriptor<M>::EKnowR){
//					std::vector<Region*> regions = (*itEdge)->getRegions();
//					std::vector<Region*>::iterator itR  = regions.begin();
//					std::vector<Region*>::iterator itRe = regions.end();
//					for(;itR!=itRe;itR++){
//						if(!mesh_.isMarked(*itR,AMark)){
//							(*itEdge)->removeRegion(*itR);
//						}
//					}
//				}
//
//				/* connection E->F*/
//				if(MeshDescriptor<M>::EKnowF){
//					std::vector<Face*> faces = (*itEdge)->getFaces();
//					std::vector<Face*>::iterator itF  = faces.begin();
//					std::vector<Face*>::iterator itFe = faces.end();
//					for(;itF!=itFe;itF++){
//						if(!mesh_.isMarked(*itF,AMark)){
//							(*itEdge)->removeFace(*itF);
//						}
//					}
//				}
//
//				/* connection E->E*/
//				if(MeshDescriptor<M>::EKnowE){
//					std::vector<Edge*> edges = (*itEdge)->getEdges();
//					std::vector<Edge*>::iterator itE  = edges.begin();
//					std::vector<Edge*>::iterator itEe = edges.end();
//					for(;itE!=itEe;itE++){
//						if(!mesh_.isMarked(*itE,AMark)){
//							(*itEdge)->removeEdge(*itE);
//						}
//					}
//				}
//
//				/* connection E->N*/
//				if(MeshDescriptor<M>::EKnowN){
//					std::vector<Node*> nodes = (*itEdge)->getNodes();
//					std::vector<Node*>::iterator itN  = nodes.begin();
//					std::vector<Node*>::iterator itNe = nodes.end();
//					for(;itN!=itNe;itN++){
//						if(!mesh_.isMarked(*itN,AMark)){
//							(*itEdge)->removeNode(*itN);
//						}
//					}
//				}
//
//			}
//		}
//	} // if(MeshDescriptor<M>::hasEdges)
//
//	if(MeshDescriptor<M>::hasNodes){
//		typename IGMesh::nodes_iterator itNode  = mesh_.nodes_begin();
//		typename IGMesh::nodes_iterator itNodee = mesh_.nodes_end();
//
//		for(;itNode!=itNodee;itNode++){
//			if(mesh_.isMarked(*itNode,AMark)){
//
//				/* connection N->R*/
//				if(MeshDescriptor<M>::NKnowR){
//					std::vector<Region*> regions = (*itNode)->getRegions();
//					std::vector<Region*>::iterator itR  = regions.begin();
//					std::vector<Region*>::iterator itRe = regions.end();
//					for(;itR!=itRe;itR++){
//						if(!mesh_.isMarked(*itR,AMark)){
//							(*itNode)->removeRegion(*itR);
//						}
//					}
//				}
//
//				/* connection N->F*/
//				if(MeshDescriptor<M>::NKnowF){
//					std::vector<Face*> faces = (*itNode)->getFaces();
//					std::vector<Face*>::iterator itF  = faces.begin();
//					std::vector<Face*>::iterator itFe = faces.end();
//					for(;itF!=itFe;itF++){
//						if(!mesh_.isMarked(*itF,AMark)){
//							(*itNode)->removeFace(*itF);
//						}
//					}
//				}
//
//				/* connection N->E*/
//				if(MeshDescriptor<M>::NKnowE){
//					std::vector<Edge*> edges = (*itNode)->getEdges();
//					std::vector<Edge*>::iterator itE  = edges.begin();
//					std::vector<Edge*>::iterator itEe = edges.end();
//					for(;itE!=itEe;itE++){
//						if(!mesh_.isMarked(*itE,AMark)){
//							(*itNode)->removeEdge(*itE);
//						}
//					}
//				}
//
//				/* connection N->N*/
//				if(MeshDescriptor<M>::NKnowN){
//					std::vector<Node*> nodes = (*itNode)->getNodes();
//					std::vector<Node*>::iterator itN  = nodes.begin();
//					std::vector<Node*>::iterator itNe = nodes.end();
//					for(;itN!=itNe;itN++){
//						if(!mesh_.isMarked(*itN,AMark)){
//							(*itNode)->removeNode(*itN);
//						}
//					}
//				}
//
//			}
//		}
//	} // if(MeshDescriptor<M>::hasNodes)
//
//	// actually removing entities.
//	if (MeshDescriptor<M>::hasRegions){
//		typename IGMesh::regions_iterator itRegion  = mesh_.regions_begin();
//		typename IGMesh::regions_iterator itRegione = mesh_.regions_end();
//
//		for(;itRegion!=itRegione;itRegion++){
//			if(!mesh_.isMarked(*itRegion,AMark)){
//				mesh_.deleteRegion(*itRegion);
//			}
//			else{
//				if(mesh_.isMarked(*itRegion,AMarkShared)){
//					if(!mesh_.isShared(*itRegion)){
//						if((*AVarRegionOwner)[(*itRegion)->getID()]==PartID){
//							mesh_.setMaster(*itRegion);
//						}
//						else{
//							mesh_.setSlave(*itRegion,(*AVarRegionLID)[(*itRegion)->getID()],(*AVarRegionOwner)[(*itRegion)->getID()]);
//						}
//					}
//				}
//			}
//		}
//	}
//
//	if (MeshDescriptor<M>::hasFaces){
//		typename IGMesh::faces_iterator itFace  = mesh_.faces_begin();
//		typename IGMesh::faces_iterator itFacee = mesh_.faces_end();
//
//		for(;itFace!=itFacee;itFace++){
//			if(!mesh_.isMarked(*itFace,AMark)){
//				mesh_.deleteFace(*itFace);
//			}
//			else{
//				if(mesh_.isMarked(*itFace,AMarkShared)){
//					if(!mesh_.isShared(*itFace)){
//						if((*AVarFaceOwner)[(*itFace)->getID()]==PartID){
//							mesh_.setMaster(*itFace);
//						}
//						else{
//							mesh_.setSlave(*itFace,(*AVarFaceLID)[(*itFace)->getID()],(*AVarFaceOwner)[(*itFace)->getID()]);
//						}
//					}
//				}
//			}
//		}
//	}
//
//	if (MeshDescriptor<M>::hasEdges){
//		typename IGMesh::edges_iterator itEdge  = mesh_.edges_begin();
//		typename IGMesh::edges_iterator itEdgee = mesh_.edges_end();
//
//		for(;itEdge!=itEdgee;itEdge++){
//			if(!mesh_.isMarked(*itEdge,AMark)){
//				mesh_.deleteEdge(*itEdge);
//			}
//			else{
//				if(mesh_.isMarked(*itEdge,AMarkShared)){
//					if(!mesh_.isShared(*itEdge)){
//						if((*AVarEdgeOwner)[(*itEdge)->getID()]==PartID){
//							mesh_.setMaster(*itEdge);
//						}
//						else{
//							mesh_.setSlave(*itEdge,(*AVarEdgeLID)[(*itEdge)->getID()],(*AVarEdgeOwner)[(*itEdge)->getID()]);
//						}
//					}
//				}
//			}
//		}
//	}
//
//	if (MeshDescriptor<M>::hasNodes){
//		typename IGMesh::nodes_iterator itNode  = mesh_.nodes_begin();
//		typename IGMesh::nodes_iterator itNodee = mesh_.nodes_end();
//
//		for(;itNode!=itNodee;itNode++){
//			if(!mesh_.isMarked(*itNode,AMark)){
//				mesh_.deleteNode(*itNode);
//			}
//			else{
//				if(mesh_.isMarked(*itNode,AMarkShared)){
//					if(!mesh_.isShared(*itNode)){
//						if((*AVarNodeOwner)[(*itNode)->getID()]==PartID){
//							mesh_.setMaster(*itNode);
//						}
//						else{
//							mesh_.setSlave(*itNode,(*AVarNodeLID)[(*itNode)->getID()],(*AVarNodeOwner)[(*itNode)->getID()]);
//						}
//					}
//				}
//			}
//		}
//	}
//
//	return;
//}
///*----------------------------------------------------------------------------*/
//
//void SubMeshExtractor::
//getSubMeshRepartition(IGMesh& ADestMesh, const int AMark, const TInt PartID,
//		Variable<TCellID>* AVarNodeLID, Variable<TCellID>* AVarEdgeLID, Variable<TCellID>* AVarFaceLID, Variable<TCellID>* AVarRegionLID,
//		const Variable<TInt>* AVarNodeOwner, const Variable<TInt>* AVarEdgeOwner, const Variable<TInt>* AVarFaceOwner, const Variable<TInt>* AVarRegionOwner,
//		const int AMarkShared)
//{
//	// clear the destination mesh.
//	ADestMesh.clear();
//
//	/* ADestMesh is built starting from nodes to regions. During this process,
//	 * downward connections can be built too. A second process from regions to
//	 * nodes is done to build upward connections. */
//	std::map<id,TCellID> mapID2Nodes;
//	std::map<id,TCellID> mapID2Edges;
//	std::map<id,TCellID> mapID2Faces;
//	std::map<id,TCellID> mapID2Regions;
//
//	if(MeshDescriptor<M>::hasNodes){
//		typename IGMesh::nodes_iterator itNode  = mesh_.nodes_begin();
//		typename IGMesh::nodes_iterator itNodee = mesh_.nodes_end();
//
//		for(;itNode!=itNodee;itNode++){
//			if(mesh_.isMarked(*itNode,AMark)){
//
//				Node* nodeDest = ADestMesh.newNode((*itNode)->getX(),(*itNode)->getY(),(*itNode)->getZ());
//
//				mapID2Nodes[(*itNode)->getID()] = nodeDest->getID();
//
//				if ((*AVarNodeOwner)[(*itNode)->getID()]==PartID){
//					(*AVarNodeLID)[(*itNode)->getID()] = nodeDest->getID();
//				}
//
//				if(mesh_.isMaster(*itNode)){
//					if((*AVarNodeOwner)[(*itNode)->getID()]==PartID){
//						ADestMesh.setMaster(nodeDest);
//					}
//					else{
//						ADestMesh.setSlave(nodeDest,(*AVarNodeLID)[(*itNode)->getID()],(*AVarNodeOwner)[(*itNode)->getID()]);
//					}
//				}
//				else{
//					if(mesh_.isSlave(*itNode)){
//						TInt masterPart;
//						id masterLID;
//						mesh_.getMasterData(*itNode,masterPart,masterLID);
//						ADestMesh.setSlave(nodeDest,masterLID,masterPart);
//					}
//					else{
//						if(mesh_.isMarked(*itNode,AMarkShared)){
//							if ((*AVarNodeOwner)[(*itNode)->getID()]==PartID){
//								ADestMesh.setMaster(nodeDest);
//							}
//							else{
//								ADestMesh.setSlave(nodeDest,(*AVarNodeLID)[(*itNode)->getID()],(*AVarNodeOwner)[(*itNode)->getID()]);
//							}
//						}
//					}
//				}
//			}
//		}
//
//	} // if(MeshDescriptor<M>::hasNodes)
//
//	if (MeshDescriptor<M>::hasEdges){
//		if(MeshDescriptor<M>::EKnowN){
//			typename IGMesh::edges_iterator itEdge  = mesh_.edges_begin();
//			typename IGMesh::edges_iterator itEdgee = mesh_.edges_end();
//
//			for(;itEdge!=itEdgee;itEdge++){
//				if(mesh_.isMarked(*itEdge,AMark)){
//					std::vector<TCellID> nodesID = (*itEdge)->getNodeIDs();
//
//					Edge* edgeDest = ADestMesh.newEdge(mapID2Nodes[nodesID[0]],mapID2Nodes[nodesID[1]]);
//
//					mapID2Edges[(*itEdge)->getID()] = edgeDest->getID();
//
//					if ((*AVarEdgeOwner)[(*itEdge)->getID()]==PartID){
//						(*AVarEdgeLID)[(*itEdge)->getID()] = edgeDest->getID();
//					}
//
//					if(mesh_.isMaster(*itEdge)){
//						if((*AVarEdgeOwner)[(*itEdge)->getID()]==PartID){
//							ADestMesh.setMaster(edgeDest);
//						}
//						else{
//							ADestMesh.setSlave(edgeDest,(*AVarEdgeLID)[(*itEdge)->getID()],(*AVarEdgeOwner)[(*itEdge)->getID()]);
//						}
//					}
//					else{
//						if(mesh_.isSlave(*itEdge)){
//							TInt masterPart;
//							id masterLID;
//							mesh_.getMasterData(*itEdge,masterPart,masterLID);
//							ADestMesh.setSlave(edgeDest,masterLID,masterPart);
//						}
//						else{
//							if(mesh_.isMarked(*itEdge,AMarkShared)){
//								if ((*AVarEdgeOwner)[(*itEdge)->getID()]==PartID){
//									ADestMesh.setMaster(edgeDest);
//								}
//								else{
//									ADestMesh.setSlave(edgeDest,(*AVarEdgeLID)[(*itEdge)->getID()],(*AVarEdgeOwner)[(*itEdge)->getID()]);
//								}
//							}
//						}
//					}
//				}
//			}
//		}
//
//		/* connection N->E*/
//		if(MeshDescriptor<M>::hasNodes && MeshDescriptor<M>::NKnowE){
//			typename IGMesh::nodes_iterator itNode  = mesh_.nodes_begin();
//			typename IGMesh::nodes_iterator itNodee = mesh_.nodes_end();
//
//			for(;itNode!=itNodee;itNode++){
//				if(mesh_.isMarked(*itNode,AMark)){
//
//					Node* nodeDest = ADestMesh.getLNode(mapID2Nodes[(*itNode)->getID()]);
//
//					std::vector<Edge*> edges = (*itNode)->getEdges();
//					std::vector<Edge*>::iterator itE  = edges.begin();
//					std::vector<Edge*>::iterator itEe = edges.end();
//					for(;itE!=itEe;itE++){
//						if(mesh_.isMarked(*itE,AMark)){
//							nodeDest->addEdge(ADestMesh.getLEdge(mapID2Edges[(*itE)->getID()]));
//						}
//					}
//				}
//			}
//		}
//
//		/* connection E->E*/
//		if(MeshDescriptor<M>::hasEdges && MeshDescriptor<M>::EKnowE){
//			typename IGMesh::edges_iterator itEdge  = mesh_.edges_begin();
//			typename IGMesh::edges_iterator itEdgee = mesh_.edges_end();
//
//			for(;itEdge!=itEdgee;itEdge++){
//				if(mesh_.isMarked(*itEdge,AMark)){
//
//					Edge* edgeDest = ADestMesh.getLEdge(mapID2Edges[(*itEdge)->getID()]);
//
//					std::vector<Edge*> edges = (*itEdge)->getEdges();
//					std::vector<Edge*>::iterator itE  = edges.begin();
//					std::vector<Edge*>::iterator itEe = edges.end();
//					for(;itE!=itEe;itE++){
//						if(mesh_.isMarked(*itE,AMark)){
//							edgeDest->addEdge(ADestMesh.getLEdge(mapID2Edges[(*itE)->getID()]));
//						}
//					}
//				}
//			}
//		}
//
//	} // if(MeshDescriptor<M>::hasEdges)
//
//	if(MeshDescriptor<M>::hasFaces){
//
//		typename IGMesh::faces_iterator itFace  = mesh_.faces_begin();
//		typename IGMesh::faces_iterator itFacee = mesh_.faces_end();
//
//		for(;itFace!=itFacee;itFace++){
//			if(mesh_.isMarked(*itFace,AMark)){
//
//				/* connection F->N*/
//				if(MeshDescriptor<M>::FKnowN){
//					std::vector<TCellID> nodesID = (*itFace)->getNodeIDs();
//					std::vector<Node*> nodes;
//
//					std::vector<TCellID>::iterator itNodeID  = nodesID.begin();
//					std::vector<TCellID>::iterator itNodeIDe = nodesID.end();
//
//					for(;itNodeID!=itNodeIDe;itNodeID++){
//						nodes.push_back(ADestMesh.getLNode(mapID2Nodes[*itNodeID]));
//					}
//
//					Face* faceDest = ADestMesh.newFace(nodes);
//
//					mapID2Faces[(*itFace)->getID()] = faceDest->getID();
//				}
//
//				/* connection F->E*/
//				if(MeshDescriptor<M>::FKnowE){
//					std::vector<TCellID> edgesID = (*itFace)->getEdgeIDs();
//					std::vector<Edge*> edges;
//					edges.reserve(edgesID.size());
//					std::vector<TCellID>::iterator itEdgeID  = edgesID.begin();
//					std::vector<TCellID>::iterator itEdgeIDe = edgesID.end();
//					for(;itEdgeID!=itEdgeIDe;itEdgeID++){
//						edges.push_back(ADestMesh.getLEdge(mapID2Edges[*itEdgeID]));
//					}
//
//					if(MeshDescriptor<M>::FKnowN){
//						Face* faceDest = ADestMesh.getLFace(mapID2Faces[(*itFace)->getID()]);
//						for(std::vector<Edge*>::iterator itEdge=edges.begin(); itEdge!=edges.end(); itEdge++){
//							faceDest->addEdge(*itEdge);
//						}
//						//faceDest->setEdges(edges);
//					}
//					else{
//						Face* faceDest = ADestMesh.newFace((*itFace)->getType());
//						for(std::vector<Edge*>::iterator itEdge=edges.begin(); itEdge!=edges.end(); itEdge++){
//							faceDest->addEdge(*itEdge);
//						}
//						//faceDest->setEdges(edges);
//						mapID2Faces[(*itFace)->getID()] = faceDest->getID();
//					}
//				}
//
//				// setting master/slave for faces.
//				Face* faceDest = ADestMesh.getLFace(mapID2Faces[(*itFace)->getID()]);
//
//				if ((*AVarFaceOwner)[(*itFace)->getID()]==PartID){
//					(*AVarFaceLID)[(*itFace)->getID()] = faceDest->getID();
//				}
//
//				if(mesh_.isMaster(*itFace)){
//					if((*AVarFaceOwner)[(*itFace)->getID()]==PartID){
//						ADestMesh.setMaster(faceDest);
//					}
//					else{
//						ADestMesh.setSlave(faceDest,(*AVarFaceLID)[(*itFace)->getID()],(*AVarFaceOwner)[(*itFace)->getID()]);
//					}
//				}
//				else{
//					if(mesh_.isSlave(*itFace)){
//						TInt masterPart;
//						id masterLID;
//						mesh_.getMasterData(*itFace,masterPart,masterLID);
//						ADestMesh.setSlave(faceDest,masterLID,masterPart);
//					}
//					else{
//						if(mesh_.isMarked(*itFace,AMarkShared)){
//							if ((*AVarFaceOwner)[(*itFace)->getID()]==PartID){
//								ADestMesh.setMaster(faceDest);
//							}
//							else{
//								ADestMesh.setSlave(faceDest,(*AVarFaceLID)[(*itFace)->getID()],(*AVarFaceOwner)[(*itFace)->getID()]);
//							}
//						}
//					}
//				}
//			}
//		}
//
//		/* connection N->F*/
//		if(MeshDescriptor<M>::hasNodes && MeshDescriptor<M>::NKnowF){
//			typename IGMesh::nodes_iterator itNode  = mesh_.nodes_begin();
//			typename IGMesh::nodes_iterator itNodee = mesh_.nodes_end();
//
//			for(;itNode!=itNodee;itNode++){
//				if(mesh_.isMarked(*itNode,AMark)){
//
//					Node* nodeDest = ADestMesh.getLNode(mapID2Nodes[(*itNode)->getID()]);
//
//					std::vector<Face*> faces = (*itNode)->getFaces();
//					std::vector<Face*>::iterator itF  = faces.begin();
//					std::vector<Face*>::iterator itFe = faces.end();
//					for(;itF!=itFe;itF++){
//						if(mesh_.isMarked(*itF,AMark)){
//							nodeDest->addFace(ADestMesh.getLFace(mapID2Faces[(*itF)->getID()]));
//						}
//					}
//				}
//			}
//		}
//
//		/* connection E->F*/
//		if(MeshDescriptor<M>::hasEdges && MeshDescriptor<M>::EKnowF){
//			typename IGMesh::edges_iterator itEdge  = mesh_.edges_begin();
//			typename IGMesh::edges_iterator itEdgee = mesh_.edges_end();
//
//			for(;itEdge!=itEdgee;itEdge++){
//				if(mesh_.isMarked(*itEdge,AMark)){
//
//					Edge* edgeDest = ADestMesh.getLEdge(mapID2Edges[(*itEdge)->getID()]);
//
//					std::vector<Face*> faces = (*itEdge)->getFaces();
//					std::vector<Face*>::iterator itF  = faces.begin();
//					std::vector<Face*>::iterator itFe = faces.end();
//					for(;itF!=itFe;itF++){
//						if(mesh_.isMarked(*itF,AMark)){
//							edgeDest->addFace(ADestMesh.getLFace(mapID2Faces[(*itF)->getID()]));
//						}
//					}
//				}
//			}
//		}
//
//		/* connection F->F*/
//		if(MeshDescriptor<M>::hasFaces && MeshDescriptor<M>::FKnowF){
//			typename IGMesh::faces_iterator itface  = mesh_.faces_begin();
//			typename IGMesh::faces_iterator itfacee = mesh_.faces_end();
//
//			for(;itFace!=itFacee;itFace++){
//				if(mesh_.isMarked(*itFace,AMark)){
//
//					Face* faceDest = ADestMesh.getLFace(mapID2Faces[(*itFace)->getID()]);
//
//					std::vector<Face*> faces = (*itFace)->getFaces();
//					std::vector<Face*>::iterator itF  = faces.begin();
//					std::vector<Face*>::iterator itFe = faces.end();
//					for(;itF!=itFe;itF++){
//						if(mesh_.isMarked(*itF,AMark)){
//							faceDest->addFace(ADestMesh.getLFace(mapID2Faces[(*itF)->getID()]));
//						}
//					}
//				}
//			}
//		}
//
//	} // if(MeshDescriptor<M>::hasFaces)
//
//	if(MeshDescriptor<M>::hasRegions){
//
//		typename IGMesh::regions_iterator itRegion  = mesh_.regions_begin();
//		typename IGMesh::regions_iterator itRegione = mesh_.regions_end();
//
//		for(;itRegion!=itRegione;itRegion++){
//			if(mesh_.isMarked(*itRegion,AMark)){
//
//				/* connection R->N*/
//				if(MeshDescriptor<M>::RKnowN){
//					std::vector<TCellID> nodesID = (*itRegion)->getNodeIDs();
//					std::vector<Node*> nodes;
//
//					std::vector<TCellID>::iterator itNodeID  = nodesID.begin();
//					std::vector<TCellID>::iterator itNodeIDe = nodesID.end();
//
//					for(;itNodeID!=itNodeIDe;itNodeID++){
//						nodes.push_back(ADestMesh.getLNode(mapID2Nodes[*itNodeID]));
//					}
//
//					Region* regionDest = ADestMesh.newRegion((*itRegion)->getType(),nodes);
//
//					mapID2Regions[(*itRegion)->getID()] = regionDest->getID();
//				}
//
//				/* connection R->E*/
//				if(MeshDescriptor<M>::RKnowE){
//					std::vector<TCellID> edgesID = (*itRegion)->getEdgeIDs();
//					std::vector<Edge*> edges;
//					edges.reserve(edgesID.size());
//					std::vector<TCellID>::iterator itEdgeID  = edgesID.begin();
//					std::vector<TCellID>::iterator itEdgeIDe = edgesID.end();
//					for(;itEdgeID!=itEdgeIDe;itEdgeID++){
//						edges.push_back(ADestMesh.getLEdge(mapID2Edges[*itEdgeID]));
//					}
//
//					if(MeshDescriptor<M>::RKnowN){
//						Region* regionDest = ADestMesh.getLRegion(mapID2Regions[(*itRegion)->getID()]);
//						for(std::vector<Edge*>::iterator itEdge=edges.begin(); itEdge!=edges.end(); itEdge++){
//							regionDest->addEdge(*itEdge);
//						}
//						//regionDest->setEdges(edges);
//					}
//					else{
//						Region* regionDest = ADestMesh.newRegion((*itRegion)->getType());
//						for(std::vector<Edge*>::iterator itEdge=edges.begin(); itEdge!=edges.end(); itEdge++){
//							regionDest->addEdge(*itEdge);
//						}
//						//regionDest->setEdges(edges);
//						mapID2Regions[(*itRegion)->getID()] = regionDest->getID();
//					}
//				}
//
//				/* connection R->F*/
//				if(MeshDescriptor<M>::RKnowF){
//					std::vector<TCellID> facesID = (*itRegion)->getFaceIDs();
//					std::vector<Face*> faces;
//					faces.reserve(facesID.size());
//					std::vector<TCellID>::iterator itFaceID  = facesID.begin();
//					std::vector<TCellID>::iterator itFaceIDe = facesID.end();
//					for(;itFaceID!=itFaceIDe;itFaceID++){
//						faces.push_back(ADestMesh.getLFace(mapID2Faces[*itFaceID]));
//					}
//
//					if(MeshDescriptor<M>::RKnowN || MeshDescriptor<M>::RKnowE){
//						Region* regionDest = ADestMesh.getLRegion(mapID2Regions[(*itRegion)->getID()]);
//						for(std::vector<Face*>::iterator itFace=faces.begin(); itFace!=faces.end(); itFace++){
//							regionDest->addFace(*itFace);
//						}
//						//regionDest->setFaces(faces);
//					}
//					else{
//						Region* regionDest = ADestMesh.newRegion((*itRegion)->getType());
//						for(std::vector<Face*>::iterator itFace=faces.begin(); itFace!=faces.end(); itFace++){
//							regionDest->addFace(*itFace);
//						}
//						//regionDest->setFaces(faces);
//						mapID2Regions[(*itRegion)->getID()] = regionDest->getID();
//					}
//				}
//
//				// setting master/slave for faces.
//				Region* regionDest = ADestMesh.getLRegion(mapID2Regions[(*itRegion)->getID()]);
//
//				if ((*AVarRegionOwner)[(*itRegion)->getID()]==PartID){
//					(*AVarRegionLID)[(*itRegion)->getID()] = regionDest->getID();
//				}
//
//				if(mesh_.isMaster(*itRegion)){
//					if((*AVarRegionOwner)[(*itRegion)->getID()]==PartID){
//						ADestMesh.setMaster(regionDest);
//					}
//					else{
//						ADestMesh.setSlave(regionDest,(*AVarRegionLID)[(*itRegion)->getID()],(*AVarRegionOwner)[(*itRegion)->getID()]);
//					}
//				}
//				else{
//					if(mesh_.isSlave(*itRegion)){
//						TInt masterPart;
//						id masterLID;
//						mesh_.getMasterData(*itRegion,masterPart,masterLID);
//						ADestMesh.setSlave(regionDest,masterLID,masterPart);
//					}
//					else{
//						if(mesh_.isMarked(*itRegion,AMarkShared)){
//							if ((*AVarRegionOwner)[(*itRegion)->getID()]==PartID){
//								ADestMesh.setMaster(regionDest);
//							}
//							else{
//								ADestMesh.setSlave(regionDest,(*AVarRegionLID)[(*itRegion)->getID()],(*AVarRegionOwner)[(*itRegion)->getID()]);
//							}
//						}
//					}
//				}
//			}
//		}
//
//		/* connection N->R*/
//		if(MeshDescriptor<M>::hasNodes && MeshDescriptor<M>::NKnowR){
//			typename IGMesh::nodes_iterator itNode  = mesh_.nodes_begin();
//			typename IGMesh::nodes_iterator itNodee = mesh_.nodes_end();
//
//			for(;itNode!=itNodee;itNode++){
//				if(mesh_.isMarked(*itNode,AMark)){
//
//					Node* nodeDest = ADestMesh.getLNode(mapID2Nodes[(*itNode)->getID()]);
//
//					std::vector<Region*> regions = (*itNode)->getRegions();
//					std::vector<Region*>::iterator itR  = regions.begin();
//					std::vector<Region*>::iterator itRe = regions.end();
//					for(;itR!=itRe;itR++){
//						if(mesh_.isMarked(*itR,AMark)){
//							nodeDest->addRegion(ADestMesh.getLRegion(mapID2Regions[(*itR)->getID()]));
//						}
//					}
//				}
//			}
//		}
//
//		/* connection E->R*/
//		if(MeshDescriptor<M>::hasEdges && MeshDescriptor<M>::EKnowR){
//			typename IGMesh::edges_iterator itEdge  = mesh_.edges_begin();
//			typename IGMesh::edges_iterator itEdgee = mesh_.edges_end();
//
//			for(;itEdge!=itEdgee;itEdge++){
//				if(mesh_.isMarked(*itEdge,AMark)){
//
//					Edge* edgeDest = ADestMesh.getLEdge(mapID2Edges[(*itEdge)->getID()]);
//
//					std::vector<Region*> regions = (*itEdge)->getRegions();
//					std::vector<Region*>::iterator itR  = regions.begin();
//					std::vector<Region*>::iterator itRe = regions.end();
//					for(;itR!=itRe;itR++){
//						if(mesh_.isMarked(*itR,AMark)){
//							edgeDest->addRegion(ADestMesh.getLRegion(mapID2Regions[(*itR)->getID()]));
//						}
//					}
//				}
//			}
//		}
//
//		/* connection F->R*/
//		if(MeshDescriptor<M>::hasFaces && MeshDescriptor<M>::FKnowR){
//			typename IGMesh::faces_iterator itFace  = mesh_.faces_begin();
//			typename IGMesh::faces_iterator itFacee = mesh_.faces_end();
//
//			for(;itFace!=itFacee;itFace++){
//				if(mesh_.isMarked(*itFace,AMark)){
//
//					Face* faceDest = ADestMesh.getLFace(mapID2Faces[(*itFace)->getID()]);
//
//					std::vector<Region*> regions = (*itFace)->getRegions();
//					std::vector<Region*>::iterator itR  = regions.begin();
//					std::vector<Region*>::iterator itRe = regions.end();
//					for(;itR!=itRe;itR++){
//						if(mesh_.isMarked(*itR,AMark)){
//							faceDest->addRegion(ADestMesh.getLRegion(mapID2Regions[(*itR)->getID()]));
//						}
//					}
//				}
//			}
//		}
//
//		/* connection R->R*/
//		if(MeshDescriptor<M>::hasRegions && MeshDescriptor<M>::RKnowR){
//			typename IGMesh::regions_iterator itRegion  = mesh_.regions_begin();
//			typename IGMesh::regions_iterator itRegione = mesh_.regions_end();
//
//			for(;itRegion!=itRegione;itRegion++){
//				if(mesh_.isMarked(*itRegion,AMark)){
//
//					Region* regionDest = ADestMesh.getLRegion(mapID2Regions[(*itRegion)->getID()]);
//
//					std::vector<Region*> regions = (*itRegion)->getRegions();
//					std::vector<Region*>::iterator itR  = regions.begin();
//					std::vector<Region*>::iterator itRe = regions.end();
//					for(;itR!=itRe;itR++){
//						if(mesh_.isMarked(*itR,AMark)){
//							regionDest->addRegion(ADestMesh.getLRegion(mapID2Regions[(*itR)->getID()]));
//						}
//					}
//				}
//			}
//		}
//
//	} // if(MeshDescriptor<M>::hasRegions)
//
//	return;
//}
///*----------------------------------------------------------------------------*/
//
//void SubMeshExtractor::
//getSubMeshRepartitionExp(std::vector<IGMesh >& SubMeshes,	std::vector<TInt>& partsOrdering, std::map<TInt,TInt>& partsOrder,
//		Variable<std::vector<TCellID> >* AVarNodeLIDs, Variable<std::vector<TCellID> >* AVarEdgeLIDs, Variable<std::vector<TCellID> >* AVarFaceLIDs, Variable<std::vector<TCellID> >* AVarRegionLIDs,
//		Variable<std::vector<TInt> >* AVarNodeOwners, Variable<std::vector<TInt> >* AVarEdgeOwners, Variable<std::vector<TInt> >* AVarFaceOwners, Variable<std::vector<TInt> >* AVarRegionOwners,
//		const int AMarkShared)
//{
//	// clear the destination meshes.
//	for(int i=0; i<SubMeshes.size(); i++){
//		SubMeshes[i].clear();
//	}
//
//	/* ADestMesh is built starting from nodes to regions. During this process,
//	 * downward connections can be built too. A second process from regions to
//	 * nodes is done to build upward connections. */
//
//	if(MeshDescriptor<M>::hasNodes){
//		typename IGMesh::nodes_iterator itNode  = mesh_.nodes_begin();
//		typename IGMesh::nodes_iterator itNodee = mesh_.nodes_end();
//
//		for(;itNode!=itNodee;itNode++){
//
//			bool toKeep = false;
//
//			for(TInt iOwner=0; iOwner<((*AVarNodeOwners)[(*itNode)->getID()]).size(); iOwner++){
//
//				TInt currentOwner = ((*AVarNodeOwners)[(*itNode)->getID()])[iOwner];
//
//				if(currentOwner==mesh_.getPartID()){
//					toKeep = true;
//					((*AVarNodeLIDs)[(*itNode)->getID()]).push_back((*itNode)->getID());
//					if(mesh_.isMarked(*itNode,AMarkShared)){
//						if(currentOwner==((*AVarNodeOwners)[(*itNode)->getID()])[0]){
//							mesh_.setMaster((*itNode));
//						}
//						else{
//							mesh_.setSlave(*itNode,((*AVarNodeLIDs)[(*itNode)->getID()])[0],((*AVarNodeOwners)[(*itNode)->getID()])[0]);
//						}
//					}
//				}
//				else{
//					TInt iSubMesh = partsOrder[currentOwner]-1;
//
//					//std::cout << (*itNode)->getID() << " " << SubMeshes.size() << " " << iSubMesh << " " << (*itNode)->getX() << " " << (*itNode)->getY() << " " << (*itNode)->getZ() << std::endl;
//					Node* nodeDest = (SubMeshes[iSubMesh]).newNode((*itNode)->getX(),(*itNode)->getY(),(*itNode)->getZ());
//					((*AVarNodeLIDs)[(*itNode)->getID()]).push_back(nodeDest->getID());
//					if(mesh_.isMarked(*itNode,AMarkShared)){
//						if(currentOwner==((*AVarNodeOwners)[(*itNode)->getID()])[0]){
//							SubMeshes[iSubMesh].setMaster(nodeDest);
//						}
//						else{
//							SubMeshes[iSubMesh].setSlave(nodeDest,((*AVarNodeLIDs)[(*itNode)->getID()])[0],((*AVarNodeOwners)[(*itNode)->getID()])[0]);
//						}
//					}
//				}
//			}
//
//			if(!toKeep){
//				//mesh_.deleteNode((*itNode));
//			}
//
//		} // for(;itNode!=itNodee;itNode++)
//
//	} // if(MeshDescriptor<M>::hasNodes)
//
//
//
////	if (MeshDescriptor<M>::hasEdges){
////		if(MeshDescriptor<M>::EKnowN){
////			typename IGMesh::edges_iterator itEdge  = mesh_.edges_begin();
////			typename IGMesh::edges_iterator itEdgee = mesh_.edges_end();
////
////			for(;itEdge!=itEdgee;itEdge++){
////				if(mesh_.isMarked(*itEdge,AMark)){
////					std::vector<TCellID> nodesID = (*itEdge)->getNodeIDs();
////
////					Edge* edgeDest = ADestMesh.newEdge(mapID2Nodes[nodesID[0]],mapID2Nodes[nodesID[1]]);
////
////					mapID2Edges[(*itEdge)->getID()] = edgeDest->getID();
////
////					if ((*AVarEdgeOwner)[(*itEdge)->getID()]==PartID){
////						(*AVarEdgeLID)[(*itEdge)->getID()] = edgeDest->getID();
////					}
////
////					if(mesh_.isMaster(*itEdge)){
////						if((*AVarEdgeOwner)[(*itEdge)->getID()]==PartID){
////							ADestMesh.setMaster(edgeDest);
////						}
////						else{
////							ADestMesh.setSlave(edgeDest,(*AVarEdgeLID)[(*itEdge)->getID()],(*AVarEdgeOwner)[(*itEdge)->getID()]);
////						}
////					}
////					else{
////						if(mesh_.isSlave(*itEdge)){
////							TInt masterPart;
////							id masterLID;
////							mesh_.getMasterData(*itEdge,masterPart,masterLID);
////							ADestMesh.setSlave(edgeDest,masterLID,masterPart);
////						}
////						else{
////							if(mesh_.isMarked(*itEdge,AMarkShared)){
////								if ((*AVarEdgeOwner)[(*itEdge)->getID()]==PartID){
////									ADestMesh.setMaster(edgeDest);
////								}
////								else{
////									ADestMesh.setSlave(edgeDest,(*AVarEdgeLID)[(*itEdge)->getID()],(*AVarEdgeOwner)[(*itEdge)->getID()]);
////								}
////							}
////						}
////					}
////				}
////			}
////		}
////
////		/* connection N->E*/
////		if(MeshDescriptor<M>::hasNodes && MeshDescriptor<M>::NKnowE){
////			typename IGMesh::nodes_iterator itNode  = mesh_.nodes_begin();
////			typename IGMesh::nodes_iterator itNodee = mesh_.nodes_end();
////
////			for(;itNode!=itNodee;itNode++){
////				if(mesh_.isMarked(*itNode,AMark)){
////
////					Node* nodeDest = ADestMesh.getLNode(mapID2Nodes[(*itNode)->getID()]);
////
////					std::vector<Edge*> edges = (*itNode)->getEdges();
////					std::vector<Edge*>::iterator itE  = edges.begin();
////					std::vector<Edge*>::iterator itEe = edges.end();
////					for(;itE!=itEe;itE++){
////						if(mesh_.isMarked(*itE,AMark)){
////							nodeDest->addEdge(ADestMesh.getLEdge(mapID2Edges[(*itE)->getID()]));
////						}
////					}
////				}
////			}
////		}
////
////		/* connection E->E*/
////		if(MeshDescriptor<M>::hasEdges && MeshDescriptor<M>::EKnowE){
////			typename IGMesh::edges_iterator itEdge  = mesh_.edges_begin();
////			typename IGMesh::edges_iterator itEdgee = mesh_.edges_end();
////
////			for(;itEdge!=itEdgee;itEdge++){
////				if(mesh_.isMarked(*itEdge,AMark)){
////
////					Edge* edgeDest = ADestMesh.getLEdge(mapID2Edges[(*itEdge)->getID()]);
////
////					std::vector<Edge*> edges = (*itEdge)->getEdges();
////					std::vector<Edge*>::iterator itE  = edges.begin();
////					std::vector<Edge*>::iterator itEe = edges.end();
////					for(;itE!=itEe;itE++){
////						if(mesh_.isMarked(*itE,AMark)){
////							edgeDest->addEdge(ADestMesh.getLEdge(mapID2Edges[(*itE)->getID()]));
////						}
////					}
////				}
////			}
////		}
////
////	} // if(MeshDescriptor<M>::hasEdges)
//
//	if(MeshDescriptor<M>::hasFaces){
//
//		typename IGMesh::faces_iterator itFace  = mesh_.faces_begin();
//		typename IGMesh::faces_iterator itFacee = mesh_.faces_end();
//
//		for(;itFace!=itFacee;itFace++){
//
//			bool toKeep = false;
//
//			for(TInt iOwner=0; iOwner<((*AVarFaceOwners)[(*itFace)->getID()]).size(); iOwner++){
//
//				TInt currentOwner = ((*AVarFaceOwners)[(*itFace)->getID()])[iOwner];
//
//				if(currentOwner==mesh_.getPartID()){
//					toKeep = true;
//					((*AVarFaceLIDs)[(*itFace)->getID()]).push_back((*itFace)->getID());
//					if(mesh_.isMarked(*itFace,AMarkShared)){
//						if(currentOwner==((*AVarFaceOwners)[(*itFace)->getID()])[0]){
//							mesh_.setMaster((*itFace));
//						}
//						else{
//							mesh_.setSlave(*itFace,((*AVarFaceLIDs)[(*itFace)->getID()])[0],((*AVarFaceOwners)[(*itFace)->getID()])[0]);
//						}
//					}
//				}
//				else{
//					TInt iSubMesh = partsOrder[currentOwner]-1;
//
//					/* connection F->N*/
//					if(MeshDescriptor<M>::FKnowN){
//						std::vector<TCellID> nodesID = (*itFace)->getNodeIDs();
//						std::vector<Node*> nodes;
//
//						std::vector<TCellID>::iterator itNodeID  = nodesID.begin();
//						std::vector<TCellID>::iterator itNodeIDe = nodesID.end();
//
//						for(;itNodeID!=itNodeIDe;itNodeID++){
//
//							TInt iNodeOwner = 0;
//							for(; iNodeOwner<(*AVarNodeOwners)[(*itNodeID)].size();iNodeOwner++){
//								if((*AVarNodeOwners)[(*itNodeID)][iNodeOwner] == currentOwner){
//									break;
//								}
//							}
//
//							nodes.push_back(SubMeshes[iSubMesh].getLNode((*AVarNodeLIDs)[(*itNodeID)][iNodeOwner]));
//						}
//
//						Face* faceDest = SubMeshes[iSubMesh].newFace(nodes);
//						((*AVarFaceLIDs)[(*itFace)->getID()]).push_back(faceDest->getID());
//
//					}
//				}
//			}
//
//			if(!toKeep){
//				//mesh_.deleteFace((*itFace));
//			}
//
//		} // for(;itFace!=itFacee;itFace++)
//
//
//
//
////					Node* nodeDest = SubMeshes[partsOrder[currentOwner]].newNode((*itNode)->getX(),(*itNode)->getY(),(*itNode)->getZ());
////					((*AVarNodeLIDs)[(*itNode)->getID()]).push_back(nodeDest->getID());
////					if(mesh_.isMarked(*itNode,AMarkShared)){
////						if(currentOwner==((*AVarNodeOwners)[(*itNode)->getID()])[0]){
////							SubMeshes[partsOrder[currentOwner]].setMaster(nodeDest);
////						}
////						else{
////							SubMeshes[partsOrder[currentOwner]].setSlave(nodeDest,((*AVarNodeLIDs)[(*itNode)->getID()])[0],((*AVarNodeOwners)[(*itNode)->getID()])[0]);
////						}
////					}
////				}
////			}
////
////
////
////
////
////
////
////
////			if(mesh_.isMarked(*itFace,AMark)){
////
////				/* connection F->N*/
////				if(MeshDescriptor<M>::FKnowN){
////					std::vector<TCellID> nodesID = (*itFace)->getNodeIDs();
////					std::vector<Node*> nodes;
////
////					std::vector<TCellID>::iterator itNodeID  = nodesID.begin();
////					std::vector<TCellID>::iterator itNodeIDe = nodesID.end();
////
////					for(;itNodeID!=itNodeIDe;itNodeID++){
////						nodes.push_back(ADestMesh.getLNode(mapID2Nodes[*itNodeID]));
////					}
////
////					Face* faceDest = ADestMesh.newFace(nodes);
////
////					mapID2Faces[(*itFace)->getID()] = faceDest->getID();
////				}
////
////				/* connection F->E*/
////				if(MeshDescriptor<M>::FKnowE){
////					std::vector<TCellID> edgesID = (*itFace)->getEdgeIDs();
////					std::vector<Edge*> edges;
////					edges.reserve(edgesID.size());
////					std::vector<TCellID>::iterator itEdgeID  = edgesID.begin();
////					std::vector<TCellID>::iterator itEdgeIDe = edgesID.end();
////					for(;itEdgeID!=itEdgeIDe;itEdgeID++){
////						edges.push_back(ADestMesh.getLEdge(mapID2Edges[*itEdgeID]));
////					}
////
////					if(MeshDescriptor<M>::FKnowN){
////						Face* faceDest = ADestMesh.getLFace(mapID2Faces[(*itFace)->getID()]);
////						for(std::vector<Edge*>::iterator itEdge=edges.begin(); itEdge!=edges.end(); itEdge++){
////							faceDest->addEdge(*itEdge);
////						}
////						//faceDest->setEdges(edges);
////					}
////					else{
////						Face* faceDest = ADestMesh.newFace((*itFace)->getType());
////						for(std::vector<Edge*>::iterator itEdge=edges.begin(); itEdge!=edges.end(); itEdge++){
////							faceDest->addEdge(*itEdge);
////						}
////						//faceDest->setEdges(edges);
////						mapID2Faces[(*itFace)->getID()] = faceDest->getID();
////					}
////				}
////
////				// setting master/slave for faces.
////				Face* faceDest = ADestMesh.getLFace(mapID2Faces[(*itFace)->getID()]);
////
////				if ((*AVarFaceOwner)[(*itFace)->getID()]==PartID){
////					(*AVarFaceLID)[(*itFace)->getID()] = faceDest->getID();
////				}
////
////				if(mesh_.isMaster(*itFace)){
////					if((*AVarFaceOwner)[(*itFace)->getID()]==PartID){
////						ADestMesh.setMaster(faceDest);
////					}
////					else{
////						ADestMesh.setSlave(faceDest,(*AVarFaceLID)[(*itFace)->getID()],(*AVarFaceOwner)[(*itFace)->getID()]);
////					}
////				}
////				else{
////					if(mesh_.isSlave(*itFace)){
////						TInt masterPart;
////						id masterLID;
////						mesh_.getMasterData(*itFace,masterPart,masterLID);
////						ADestMesh.setSlave(faceDest,masterLID,masterPart);
////					}
////					else{
////						if(mesh_.isMarked(*itFace,AMarkShared)){
////							if ((*AVarFaceOwner)[(*itFace)->getID()]==PartID){
////								ADestMesh.setMaster(faceDest);
////							}
////							else{
////								ADestMesh.setSlave(faceDest,(*AVarFaceLID)[(*itFace)->getID()],(*AVarFaceOwner)[(*itFace)->getID()]);
////							}
////						}
////					}
////				}
////			}
////		}
////
////		/* connection N->F*/
////		if(MeshDescriptor<M>::hasNodes && MeshDescriptor<M>::NKnowF){
////			typename IGMesh::nodes_iterator itNode  = mesh_.nodes_begin();
////			typename IGMesh::nodes_iterator itNodee = mesh_.nodes_end();
////
////			for(;itNode!=itNodee;itNode++){
////				if(mesh_.isMarked(*itNode,AMark)){
////
////					Node* nodeDest = ADestMesh.getLNode(mapID2Nodes[(*itNode)->getID()]);
////
////					std::vector<Face*> faces = (*itNode)->getFaces();
////					std::vector<Face*>::iterator itF  = faces.begin();
////					std::vector<Face*>::iterator itFe = faces.end();
////					for(;itF!=itFe;itF++){
////						if(mesh_.isMarked(*itF,AMark)){
////							nodeDest->addFace(ADestMesh.getLFace(mapID2Faces[(*itF)->getID()]));
////						}
////					}
////				}
////			}
////		}
////
////		/* connection E->F*/
////		if(MeshDescriptor<M>::hasEdges && MeshDescriptor<M>::EKnowF){
////			typename IGMesh::edges_iterator itEdge  = mesh_.edges_begin();
////			typename IGMesh::edges_iterator itEdgee = mesh_.edges_end();
////
////			for(;itEdge!=itEdgee;itEdge++){
////				if(mesh_.isMarked(*itEdge,AMark)){
////
////					Edge* edgeDest = ADestMesh.getLEdge(mapID2Edges[(*itEdge)->getID()]);
////
////					std::vector<Face*> faces = (*itEdge)->getFaces();
////					std::vector<Face*>::iterator itF  = faces.begin();
////					std::vector<Face*>::iterator itFe = faces.end();
////					for(;itF!=itFe;itF++){
////						if(mesh_.isMarked(*itF,AMark)){
////							edgeDest->addFace(ADestMesh.getLFace(mapID2Faces[(*itF)->getID()]));
////						}
////					}
////				}
////			}
////		}
////
////		/* connection F->F*/
////		if(MeshDescriptor<M>::hasFaces && MeshDescriptor<M>::FKnowF){
////			typename IGMesh::faces_iterator itface  = mesh_.faces_begin();
////			typename IGMesh::faces_iterator itfacee = mesh_.faces_end();
////
////			for(;itFace!=itFacee;itFace++){
////				if(mesh_.isMarked(*itFace,AMark)){
////
////					Face* faceDest = ADestMesh.getLFace(mapID2Faces[(*itFace)->getID()]);
////
////					std::vector<Face*> faces = (*itFace)->getFaces();
////					std::vector<Face*>::iterator itF  = faces.begin();
////					std::vector<Face*>::iterator itFe = faces.end();
////					for(;itF!=itFe;itF++){
////						if(mesh_.isMarked(*itF,AMark)){
////							faceDest->addFace(ADestMesh.getLFace(mapID2Faces[(*itF)->getID()]));
////						}
////					}
////				}
////			}
////		}
//
//	} // if(MeshDescriptor<M>::hasFaces)
//
////	if(MeshDescriptor<M>::hasRegions){
////
////		typename IGMesh::regions_iterator itRegion  = mesh_.regions_begin();
////		typename IGMesh::regions_iterator itRegione = mesh_.regions_end();
////
////		for(;itRegion!=itRegione;itRegion++){
////			if(mesh_.isMarked(*itRegion,AMark)){
////
////				/* connection R->N*/
////				if(MeshDescriptor<M>::RKnowN){
////					std::vector<TCellID> nodesID = (*itRegion)->getNodeIDs();
////					std::vector<Node*> nodes;
////
////					std::vector<TCellID>::iterator itNodeID  = nodesID.begin();
////					std::vector<TCellID>::iterator itNodeIDe = nodesID.end();
////
////					for(;itNodeID!=itNodeIDe;itNodeID++){
////						nodes.push_back(ADestMesh.getLNode(mapID2Nodes[*itNodeID]));
////					}
////
////					Region* regionDest = ADestMesh.newRegion((*itRegion)->getType(),nodes);
////
////					mapID2Regions[(*itRegion)->getID()] = regionDest->getID();
////				}
////
////				/* connection R->E*/
////				if(MeshDescriptor<M>::RKnowE){
////					std::vector<TCellID> edgesID = (*itRegion)->getEdgeIDs();
////					std::vector<Edge*> edges;
////					edges.reserve(edgesID.size());
////					std::vector<TCellID>::iterator itEdgeID  = edgesID.begin();
////					std::vector<TCellID>::iterator itEdgeIDe = edgesID.end();
////					for(;itEdgeID!=itEdgeIDe;itEdgeID++){
////						edges.push_back(ADestMesh.getLEdge(mapID2Edges[*itEdgeID]));
////					}
////
////					if(MeshDescriptor<M>::RKnowN){
////						Region* regionDest = ADestMesh.getLRegion(mapID2Regions[(*itRegion)->getID()]);
////						for(std::vector<Edge*>::iterator itEdge=edges.begin(); itEdge!=edges.end(); itEdge++){
////							regionDest->addEdge(*itEdge);
////						}
////						//regionDest->setEdges(edges);
////					}
////					else{
////						Region* regionDest = ADestMesh.newRegion((*itRegion)->getType());
////						for(std::vector<Edge*>::iterator itEdge=edges.begin(); itEdge!=edges.end(); itEdge++){
////							regionDest->addEdge(*itEdge);
////						}
////						//regionDest->setEdges(edges);
////						mapID2Regions[(*itRegion)->getID()] = regionDest->getID();
////					}
////				}
////
////				/* connection R->F*/
////				if(MeshDescriptor<M>::RKnowF){
////					std::vector<TCellID> facesID = (*itRegion)->getFaceIDs();
////					std::vector<Face*> faces;
////					faces.reserve(facesID.size());
////					std::vector<TCellID>::iterator itFaceID  = facesID.begin();
////					std::vector<TCellID>::iterator itFaceIDe = facesID.end();
////					for(;itFaceID!=itFaceIDe;itFaceID++){
////						faces.push_back(ADestMesh.getLFace(mapID2Faces[*itFaceID]));
////					}
////
////					if(MeshDescriptor<M>::RKnowN || MeshDescriptor<M>::RKnowE){
////						Region* regionDest = ADestMesh.getLRegion(mapID2Regions[(*itRegion)->getID()]);
////						for(std::vector<Face*>::iterator itFace=faces.begin(); itFace!=faces.end(); itFace++){
////							regionDest->addFace(*itFace);
////						}
////						//regionDest->setFaces(faces);
////					}
////					else{
////						Region* regionDest = ADestMesh.newRegion((*itRegion)->getType());
////						for(std::vector<Face*>::iterator itFace=faces.begin(); itFace!=faces.end(); itFace++){
////							regionDest->addFace(*itFace);
////						}
////						//regionDest->setFaces(faces);
////						mapID2Regions[(*itRegion)->getID()] = regionDest->getID();
////					}
////				}
////
////				// setting master/slave for faces.
////				Region* regionDest = ADestMesh.getLRegion(mapID2Regions[(*itRegion)->getID()]);
////
////				if ((*AVarRegionOwner)[(*itRegion)->getID()]==PartID){
////					(*AVarRegionLID)[(*itRegion)->getID()] = regionDest->getID();
////				}
////
////				if(mesh_.isMaster(*itRegion)){
////					if((*AVarRegionOwner)[(*itRegion)->getID()]==PartID){
////						ADestMesh.setMaster(regionDest);
////					}
////					else{
////						ADestMesh.setSlave(regionDest,(*AVarRegionLID)[(*itRegion)->getID()],(*AVarRegionOwner)[(*itRegion)->getID()]);
////					}
////				}
////				else{
////					if(mesh_.isSlave(*itRegion)){
////						TInt masterPart;
////						id masterLID;
////						mesh_.getMasterData(*itRegion,masterPart,masterLID);
////						ADestMesh.setSlave(regionDest,masterLID,masterPart);
////					}
////					else{
////						if(mesh_.isMarked(*itRegion,AMarkShared)){
////							if ((*AVarRegionOwner)[(*itRegion)->getID()]==PartID){
////								ADestMesh.setMaster(regionDest);
////							}
////							else{
////								ADestMesh.setSlave(regionDest,(*AVarRegionLID)[(*itRegion)->getID()],(*AVarRegionOwner)[(*itRegion)->getID()]);
////							}
////						}
////					}
////				}
////			}
////		}
////
////		/* connection N->R*/
////		if(MeshDescriptor<M>::hasNodes && MeshDescriptor<M>::NKnowR){
////			typename IGMesh::nodes_iterator itNode  = mesh_.nodes_begin();
////			typename IGMesh::nodes_iterator itNodee = mesh_.nodes_end();
////
////			for(;itNode!=itNodee;itNode++){
////				if(mesh_.isMarked(*itNode,AMark)){
////
////					Node* nodeDest = ADestMesh.getLNode(mapID2Nodes[(*itNode)->getID()]);
////
////					std::vector<Region*> regions = (*itNode)->getRegions();
////					std::vector<Region*>::iterator itR  = regions.begin();
////					std::vector<Region*>::iterator itRe = regions.end();
////					for(;itR!=itRe;itR++){
////						if(mesh_.isMarked(*itR,AMark)){
////							nodeDest->addRegion(ADestMesh.getLRegion(mapID2Regions[(*itR)->getID()]));
////						}
////					}
////				}
////			}
////		}
////
////		/* connection E->R*/
////		if(MeshDescriptor<M>::hasEdges && MeshDescriptor<M>::EKnowR){
////			typename IGMesh::edges_iterator itEdge  = mesh_.edges_begin();
////			typename IGMesh::edges_iterator itEdgee = mesh_.edges_end();
////
////			for(;itEdge!=itEdgee;itEdge++){
////				if(mesh_.isMarked(*itEdge,AMark)){
////
////					Edge* edgeDest = ADestMesh.getLEdge(mapID2Edges[(*itEdge)->getID()]);
////
////					std::vector<Region*> regions = (*itEdge)->getRegions();
////					std::vector<Region*>::iterator itR  = regions.begin();
////					std::vector<Region*>::iterator itRe = regions.end();
////					for(;itR!=itRe;itR++){
////						if(mesh_.isMarked(*itR,AMark)){
////							edgeDest->addRegion(ADestMesh.getLRegion(mapID2Regions[(*itR)->getID()]));
////						}
////					}
////				}
////			}
////		}
////
////		/* connection F->R*/
////		if(MeshDescriptor<M>::hasFaces && MeshDescriptor<M>::FKnowR){
////			typename IGMesh::faces_iterator itFace  = mesh_.faces_begin();
////			typename IGMesh::faces_iterator itFacee = mesh_.faces_end();
////
////			for(;itFace!=itFacee;itFace++){
////				if(mesh_.isMarked(*itFace,AMark)){
////
////					Face* faceDest = ADestMesh.getLFace(mapID2Faces[(*itFace)->getID()]);
////
////					std::vector<Region*> regions = (*itFace)->getRegions();
////					std::vector<Region*>::iterator itR  = regions.begin();
////					std::vector<Region*>::iterator itRe = regions.end();
////					for(;itR!=itRe;itR++){
////						if(mesh_.isMarked(*itR,AMark)){
////							faceDest->addRegion(ADestMesh.getLRegion(mapID2Regions[(*itR)->getID()]));
////						}
////					}
////				}
////			}
////		}
////
////		/* connection R->R*/
////		if(MeshDescriptor<M>::hasRegions && MeshDescriptor<M>::RKnowR){
////			typename IGMesh::regions_iterator itRegion  = mesh_.regions_begin();
////			typename IGMesh::regions_iterator itRegione = mesh_.regions_end();
////
////			for(;itRegion!=itRegione;itRegion++){
////				if(mesh_.isMarked(*itRegion,AMark)){
////
////					Region* regionDest = ADestMesh.getLRegion(mapID2Regions[(*itRegion)->getID()]);
////
////					std::vector<Region*> regions = (*itRegion)->getRegions();
////					std::vector<Region*>::iterator itR  = regions.begin();
////					std::vector<Region*>::iterator itRe = regions.end();
////					for(;itR!=itRe;itR++){
////						if(mesh_.isMarked(*itR,AMark)){
////							regionDest->addRegion(ADestMesh.getLRegion(mapID2Regions[(*itR)->getID()]));
////						}
////					}
////				}
////			}
////		}
////
////	} // if(MeshDescriptor<M>::hasRegions)
//
//	if(MeshDescriptor<M>::hasNodes){
//		typename IGMesh::nodes_iterator itNode  = mesh_.nodes_begin();
//		typename IGMesh::nodes_iterator itNodee = mesh_.nodes_end();
//
//		for(;itNode!=itNodee;itNode++){
//
//			bool toKeep = false;
//
//			for(TInt iOwner=0; iOwner<((*AVarNodeOwners)[(*itNode)->getID()]).size(); iOwner++){
//
//				TInt currentOwner = ((*AVarNodeOwners)[(*itNode)->getID()])[iOwner];
//
//				if(currentOwner==mesh_.getPartID()){
//					toKeep = true;
//				}
//			}
//
//			if(!toKeep){
//				mesh_.deleteNode(*itNode);
//			}
//		}
//	}
//
//	if(MeshDescriptor<M>::hasFaces){
//		typename IGMesh::faces_iterator itFace  = mesh_.faces_begin();
//		typename IGMesh::faces_iterator itFacee = mesh_.faces_end();
//
//		for(;itFace!=itFacee;itFace++){
//
//			bool toKeep = false;
//
//			for(TInt iOwner=0; iOwner<((*AVarFaceOwners)[(*itFace)->getID()]).size(); iOwner++){
//
//				TInt currentOwner = ((*AVarFaceOwners)[(*itFace)->getID()])[iOwner];
//
//				if(currentOwner==mesh_.getPartID()){
//					toKeep = true;
//				}
//			}
//
//			if(!toKeep){
//				mesh_.deleteFace(*itFace);
//			}
//		}
//	}
//
//	return;
//}
///*----------------------------------------------------------------------------*/
//
//void SubMeshExtractor::
//getSubMeshRepartitionExp2(std::vector<IGMesh >& SubMeshes,	std::vector<TInt>& partsOrdering, std::map<TInt,TInt>& partsOrder,
//		bool* isSharedNode,bool* isSharedEdge,bool* isSharedFace,bool* isSharedRegion,
//		id* nodes_LID, id* edges_LID, id* faces_LID, id* regions_LID,
//		TInt* nodes_owner,TInt* edges_owner, TInt* faces_owner,TInt* regions_owner,
//		std::vector<std::map<TInt,TCellID> >& nodes_owners2LIDs, std::vector<std::map<TInt,TCellID> >& edges_owners2LIDs, std::vector<std::map<TInt,TCellID> >& faces_owners2LIDs, std::vector<std::map<TInt,TCellID> >& regions_owners2LIDs)
//{
//	// clear the destination meshes.
//	for(int i=0; i<SubMeshes.size(); i++){
//		SubMeshes[i].clear();
//	}
//
//	// mark entities that will be kept in mesh_.
//	int entityToDeleteMark = mesh_.getNewMark();
//
//	if(MeshDescriptor<M>::hasNodes){
//		typename IGMesh::nodes_iterator itNode  = mesh_.nodes_begin();
//		typename IGMesh::nodes_iterator itNodee = mesh_.nodes_end();
//
//		for(;itNode!=itNodee;itNode++){
//			bool toKeep = false;
//
//			if(!isSharedNode[(*itNode)->getID()]){
//				TInt currentOwner = nodes_owner[(*itNode)->getID()];
//
//				if(currentOwner == mesh_.getPartID()){
//					toKeep = true;
//
//					nodes_LID[(*itNode)->getID()] = (*itNode)->getID();
//				}
//				else{
//					TInt iSubMesh = partsOrder[currentOwner]-1;
//
//					Node* nodeDest = (SubMeshes[iSubMesh]).newNode((*itNode)->getX(),(*itNode)->getY(),(*itNode)->getZ());
//					nodes_LID[(*itNode)->getID()] = nodeDest->getID();
//				}
//			}
//			else{
//				std::map<TInt,TCellID>::iterator itNodeOwner2LID = nodes_owners2LIDs[nodes_LID[(*itNode)->getID()]].begin();
//				std::map<TInt,TCellID>::iterator itNodeOwner2LIDe = nodes_owners2LIDs[nodes_LID[(*itNode)->getID()]].end();
//
//				for(;itNodeOwner2LID!=itNodeOwner2LIDe;itNodeOwner2LID++){
//					TInt currentOwner = itNodeOwner2LID->first;
//
//					if(currentOwner == mesh_.getPartID()){
//						toKeep = true;
//						itNodeOwner2LID->second = (*itNode)->getID();
//
//						if(nodes_owner[(*itNode)->getID()] == mesh_.getPartID()){
//							mesh_.setMaster((*itNode));
//						}
//						else{
//							mesh_.setSlave(*itNode,nodes_owners2LIDs[nodes_LID[(*itNode)->getID()]][nodes_owner[(*itNode)->getID()]],nodes_owner[(*itNode)->getID()]);
//						}
//					}
//					else{
//						TInt iSubMesh = partsOrder[currentOwner]-1;
//
//						Node* nodeDest = (SubMeshes[iSubMesh]).newNode((*itNode)->getX(),(*itNode)->getY(),(*itNode)->getZ());
//						nodes_owners2LIDs[nodes_LID[(*itNode)->getID()]][currentOwner] = nodeDest->getID();
//
//						if(currentOwner == nodes_owner[(*itNode)->getID()]){
//							SubMeshes[iSubMesh].setMaster(nodeDest);
//						}
//						else{
//							SubMeshes[iSubMesh].setSlave(nodeDest,nodes_owners2LIDs[nodes_LID[(*itNode)->getID()]][nodes_owner[(*itNode)->getID()]],nodes_owner[(*itNode)->getID()]);
//						}
//					}
//				}
//			} // if(!isSharedNode[(*itNode)->getID()])
//
//			if(!toKeep){
//				//mesh_.deleteNode((*itNode));
//				mesh_.mark((*itNode),entityToDeleteMark);
//			}
//
//		} // for(;itNode!=itNodee;itNode++)
//
//	} // if(MeshDescriptor<M>::hasNodes)
//
//	if(MeshDescriptor<M>::hasEdges){
//		typename IGMesh::edges_iterator itEdge  = mesh_.edges_begin();
//		typename IGMesh::edges_iterator itEdgee = mesh_.edges_end();
//
//		for(;itEdge!=itEdgee;itEdge++){
//			bool toKeep = false;
//
//			if(!isSharedEdge[(*itEdge)->getID()]){
//				TInt currentOwner = edges_owner[(*itEdge)->getID()];
//
//				if(currentOwner == mesh_.getPartID()){
//					toKeep = true;
//
//					edges_LID[(*itEdge)->getID()] = (*itEdge)->getID();
//				}
//				else{
//					TInt iSubMesh = partsOrder[currentOwner]-1;
//
//					/* connection E->N*/
//					if(MeshDescriptor<M>::EKnowN){
//						std::vector<TCellID> nodesID = (*itEdge)->getNodeIDs();
//						std::vector<Node*> nodes;
//
//						std::vector<TCellID>::iterator itNodeID  = nodesID.begin();
//						std::vector<TCellID>::iterator itNodeIDe = nodesID.end();
//
//						for(;itNodeID!=itNodeIDe;itNodeID++){
//							if(!isSharedNode[*itNodeID]){
//								nodes.push_back(SubMeshes[iSubMesh].getLNode(nodes_LID[*itNodeID]));
//							}
//							else{
//								nodes.push_back(SubMeshes[iSubMesh].getLNode(nodes_owners2LIDs[nodes_LID[*itNodeID]][currentOwner]));
//							}
//						}
//
//						Edge* edgeDest = SubMeshes[iSubMesh].newEdge(nodes[0],nodes[1]);
//						edges_LID[(*itEdge)->getID()] = edgeDest->getID();
//					}
//				}
//			}
//			else{
//				std::map<TInt,TCellID>::iterator itEdgeOwner2LID = edges_owners2LIDs[edges_LID[(*itEdge)->getID()]].begin();
//				std::map<TInt,TCellID>::iterator itEdgeOwner2LIDe = edges_owners2LIDs[edges_LID[(*itEdge)->getID()]].end();
//
//				for(;itEdgeOwner2LID!=itEdgeOwner2LIDe;itEdgeOwner2LID++){
//					TInt currentOwner = itEdgeOwner2LID->first;
//
//					if(currentOwner == mesh_.getPartID()){
//						toKeep = true;
//
//						itEdgeOwner2LID->second = (*itEdge)->getID();
//
//						if(edges_owner[(*itEdge)->getID()] == mesh_.getPartID()){
//							mesh_.setMaster((*itEdge));
//						}
//						else{
//							mesh_.setSlave(*itEdge,edges_owners2LIDs[edges_LID[(*itEdge)->getID()]][edges_owner[(*itEdge)->getID()]],edges_owner[(*itEdge)->getID()]);
//						}
//					}
//					else{
//						TInt iSubMesh = partsOrder[currentOwner]-1;
//
//						/* connection E->N*/
//						if(MeshDescriptor<M>::EKnowN){
//							std::vector<TCellID> nodesID = (*itEdge)->getNodeIDs();
//							std::vector<Node*> nodes;
//
//							std::vector<TCellID>::iterator itNodeID  = nodesID.begin();
//							std::vector<TCellID>::iterator itNodeIDe = nodesID.end();
//
//							for(;itNodeID!=itNodeIDe;itNodeID++){
//								if(!isSharedNode[*itNodeID]){
//									nodes.push_back(SubMeshes[iSubMesh].getLNode(nodes_LID[*itNodeID]));
//								}
//								else{
//									nodes.push_back(SubMeshes[iSubMesh].getLNode(nodes_owners2LIDs[nodes_LID[*itNodeID]][currentOwner]));
//								}
//							}
//
//							Edge* edgeDest = SubMeshes[iSubMesh].newEdge(nodes[0],nodes[1]);
//							edges_owners2LIDs[edges_LID[(*itEdge)->getID()]][currentOwner] = edgeDest->getID();
//
//							if(currentOwner == edges_owner[(*itEdge)->getID()]){
//								SubMeshes[iSubMesh].setMaster((edgeDest));
//							}
//							else{
//								SubMeshes[iSubMesh].setSlave(edgeDest,edges_owners2LIDs[edges_LID[(*itEdge)->getID()]][edges_owner[(*itEdge)->getID()]],edges_owner[(*itEdge)->getID()]);
//							}
//						}
//					}
//				}
//			} // if(!isSharedEdge[(*itEdge)->getID()])
//
//			if(!toKeep){
//				//mesh_.deleteEdge((*itEdge));
//				mesh_.mark((*itEdge),entityToDeleteMark);
//			}
//
//		} // for(;itEdge!=itEdgee;itEdge++)
//
//		/* connection N->E*/
//		if(MeshDescriptor<M>::hasNodes && MeshDescriptor<M>::NKnowE){
//			typename IGMesh::nodes_iterator itNode  = mesh_.nodes_begin();
//			typename IGMesh::nodes_iterator itNodee = mesh_.nodes_end();
//
//			for(;itNode!=itNodee;itNode++){
//
//				std::vector<Edge*> edges2add;
//				std::vector<TCellID> edgesID = (*itNode)->getEdgeIDs();
//				std::vector<TCellID>::iterator itEdgeID;
//				std::vector<TCellID>::iterator itEdgeIDb = edgesID.begin();
//				std::vector<TCellID>::iterator itEdgeIDe = edgesID.end();
//
//				if(!isSharedNode[(*itNode)->getID()]){
//					TInt currentOwner = nodes_owner[(*itNode)->getID()];
//
//					if(currentOwner == mesh_.getPartID()){
//						Node* nodeDest = mesh_.getLNode((*itNode)->getID());
//
//						for(edges2add.empty(),itEdgeID=itEdgeIDb;itEdgeID!=itEdgeIDe;itEdgeID++){
//							if(!isSharedEdge[*itEdgeID]){
//								if(edges_owner[*itEdgeID] == mesh_.getPartID()){
//									edges2add.push_back(mesh_.getLEdge(*itEdgeID));
//								}
//								else{
//								}
//							}
//							else{
//								if(edges_owners2LIDs[edges_LID[*itEdgeID]].find(currentOwner) != edges_owners2LIDs[edges_LID[*itEdgeID]].end()){
//									edges2add.push_back(mesh_.getLEdge(*itEdgeID));
//								}
//								else{
//								}
//							}
//						}
//						nodeDest->setEdges(edges2add);
//					}
//					else{
//						TInt iSubMesh = partsOrder[currentOwner]-1;
//
//						Node* nodeDest = SubMeshes[iSubMesh].getLNode(nodes_LID[(*itNode)->getID()]);
//
//						for(edges2add.empty(),itEdgeID=itEdgeIDb;itEdgeID!=itEdgeIDe;itEdgeID++){
//							if(!isSharedEdge[*itEdgeID]){
//								if(edges_owner[*itEdgeID] == currentOwner){
//									edges2add.push_back(SubMeshes[iSubMesh].getLEdge(edges_LID[*itEdgeID]));
//								}
//								else{
//								}
//							}
//							else{
//								if(edges_owners2LIDs[edges_LID[*itEdgeID]].find(currentOwner) != edges_owners2LIDs[edges_LID[*itEdgeID]].end()){
//									edges2add.push_back(SubMeshes[iSubMesh].getLEdge(edges_owners2LIDs[edges_LID[*itEdgeID]][currentOwner]));
//								}
//								else{
//								}
//							}
//						}
//						nodeDest->setEdges(edges2add);
//					}
//				}
//				else{
//					std::map<TInt,TCellID>::iterator itNodeOwner2LID = nodes_owners2LIDs[nodes_LID[(*itNode)->getID()]].begin();
//					std::map<TInt,TCellID>::iterator itNodeOwner2LIDe = nodes_owners2LIDs[nodes_LID[(*itNode)->getID()]].end();
//
//					for(;itNodeOwner2LID!=itNodeOwner2LIDe;itNodeOwner2LID++){
//						TInt currentOwner = itNodeOwner2LID->first;
//
//						if(currentOwner == mesh_.getPartID()){
//
//							Node* nodeDest = mesh_.getLNode(itNodeOwner2LID->second);
//
//							for(edges2add.empty(),itEdgeID=itEdgeIDb;itEdgeID!=itEdgeIDe;itEdgeID++){
//								if(!isSharedEdge[*itEdgeID]){
//									if(edges_owner[*itEdgeID] == mesh_.getPartID()){
//										edges2add.push_back(mesh_.getLEdge(*itEdgeID));
//									}
//									else{
//									}
//								}
//								else{
//									if(edges_owners2LIDs[edges_LID[*itEdgeID]].find(currentOwner) != edges_owners2LIDs[edges_LID[*itEdgeID]].end()){
//										edges2add.push_back(mesh_.getLEdge(*itEdgeID));
//									}
//									else{
//									}
//								}
//							}
//							nodeDest->setEdges(edges2add);
//						}
//						else{
//							TInt iSubMesh = partsOrder[currentOwner]-1;
//
//							Node* nodeDest = SubMeshes[iSubMesh].getLNode(nodes_owners2LIDs[nodes_LID[(*itNode)->getID()]][currentOwner]);
//
//							for(edges2add.empty(),itEdgeID=itEdgeIDb;itEdgeID!=itEdgeIDe;itEdgeID++){
//								if(!isSharedEdge[*itEdgeID]){
//									if(edges_owner[*itEdgeID] == currentOwner){
//										edges2add.push_back(SubMeshes[iSubMesh].getLEdge(edges_LID[*itEdgeID]));
//									}
//									else{
//									}
//								}
//								else{
//									if(edges_owners2LIDs[edges_LID[*itEdgeID]].find(currentOwner) != edges_owners2LIDs[edges_LID[*itEdgeID]].end()){
//										edges2add.push_back(SubMeshes[iSubMesh].getLEdge(edges_owners2LIDs[edges_LID[*itEdgeID]][currentOwner]));
//									}
//									else{
//									}
//								}
//							}
//							nodeDest->setEdges(edges2add);
//						}
//					}
//				}
//			} // for(;itNode!=itNodee;itNode++)
//		} // if(MeshDescriptor<M>::hasNodes && MeshDescriptor<M>::NKnowE)
//
//		/* connection E->E*/
//		if(MeshDescriptor<M>::hasEdges && MeshDescriptor<M>::EKnowE){
//			typename IGMesh::edges_iterator itEdge  = mesh_.edges_begin();
//			typename IGMesh::edges_iterator itEdgee = mesh_.edges_end();
//
//			for(;itEdge!=itEdgee;itEdge++){
//
//				std::vector<Edge*> edges2add;
//				std::vector<TCellID> edgesID = (*itEdge)->getEdgeIDs();
//				std::vector<TCellID>::iterator itEdgeID;
//				std::vector<TCellID>::iterator itEdgeIDb = edgesID.begin();
//				std::vector<TCellID>::iterator itEdgeIDe = edgesID.end();
//
//				if(!isSharedEdge[(*itEdge)->getID()]){
//					TInt currentOwner = edges_owner[(*itEdge)->getID()];
//
//					if(currentOwner == mesh_.getPartID()){
//						Edge* edgeDest = mesh_.getLEdge((*itEdge)->getID());
//
//						for(edges2add.empty(),itEdgeID=itEdgeIDb;itEdgeID!=itEdgeIDe;itEdgeID++){
//							if(!isSharedEdge[*itEdgeID]){
//								if(edges_owner[*itEdgeID] == mesh_.getPartID()){
//									edges2add.push_back(mesh_.getLEdge(*itEdgeID));
//								}
//								else{
//								}
//							}
//							else{
//								if(edges_owners2LIDs[edges_LID[*itEdgeID]].find(currentOwner) != edges_owners2LIDs[edges_LID[*itEdgeID]].end()){
//									edges2add.push_back(mesh_.getLEdge(*itEdgeID));
//								}
//								else{
//								}
//							}
//						}
//						edgeDest->setEdges(edges2add);
//					}
//					else{
//						TInt iSubMesh = partsOrder[currentOwner]-1;
//
//						Edge* edgeDest = SubMeshes[iSubMesh].getLEdge(edges_LID[(*itEdge)->getID()]);
//
//						for(edges2add.empty(),itEdgeID=itEdgeIDb;itEdgeID!=itEdgeIDe;itEdgeID++){
//							if(!isSharedEdge[*itEdgeID]){
//								if(edges_owner[*itEdgeID] == currentOwner){
//									edges2add.push_back(SubMeshes[iSubMesh].getLEdge(edges_LID[*itEdgeID]));
//								}
//								else{
//								}
//							}
//							else{
//								if(edges_owners2LIDs[edges_LID[*itEdgeID]].find(currentOwner) != edges_owners2LIDs[edges_LID[*itEdgeID]].end()){
//									edges2add.push_back(SubMeshes[iSubMesh].getLEdge(edges_owners2LIDs[edges_LID[*itEdgeID]][currentOwner]));
//								}
//								else{
//								}
//							}
//						}
//						edgeDest->setEdges(edges2add);
//					}
//				}
//				else{
//					std::map<TInt,TCellID>::iterator itEdgeOwner2LID = edges_owners2LIDs[edges_LID[(*itEdge)->getID()]].begin();
//					std::map<TInt,TCellID>::iterator itEdgeOwner2LIDe = edges_owners2LIDs[edges_LID[(*itEdge)->getID()]].end();
//
//					for(;itEdgeOwner2LID!=itEdgeOwner2LIDe;itEdgeOwner2LID++){
//						TInt currentOwner = itEdgeOwner2LID->first;
//
//						if(currentOwner == mesh_.getPartID()){
//
//							Edge* edgeDest = mesh_.getLEdge(itEdgeOwner2LID->second);
//
//							for(edges2add.empty(),itEdgeID=itEdgeIDb;itEdgeID!=itEdgeIDe;itEdgeID++){
//								if(!isSharedEdge[*itEdgeID]){
//									if(edges_owner[*itEdgeID] == mesh_.getPartID()){
//										edges2add.push_back(mesh_.getLEdge(*itEdgeID));
//									}
//									else{
//									}
//								}
//								else{
//									if(edges_owners2LIDs[edges_LID[*itEdgeID]].find(currentOwner) != edges_owners2LIDs[edges_LID[*itEdgeID]].end()){
//										edges2add.push_back(mesh_.getLEdge(*itEdgeID));
//									}
//									else{
//									}
//								}
//							}
//							edgeDest->setEdges(edges2add);
//						}
//						else{
//							TInt iSubMesh = partsOrder[currentOwner]-1;
//
//							Edge* edgeDest = SubMeshes[iSubMesh].getLEdge(edges_owners2LIDs[edges_LID[(*itEdge)->getID()]][currentOwner]);
//
//							for(edges2add.empty(),itEdgeID=itEdgeIDb;itEdgeID!=itEdgeIDe;itEdgeID++){
//								if(!isSharedEdge[*itEdgeID]){
//									if(edges_owner[*itEdgeID] == currentOwner){
//										edges2add.push_back(SubMeshes[iSubMesh].getLEdge(edges_LID[*itEdgeID]));
//									}
//									else{
//									}
//								}
//								else{
//									if(edges_owners2LIDs[edges_LID[*itEdgeID]].find(currentOwner) != edges_owners2LIDs[edges_LID[*itEdgeID]].end()){
//										edges2add.push_back(SubMeshes[iSubMesh].getLEdge(edges_owners2LIDs[edges_LID[*itEdgeID]][currentOwner]));
//									}
//									else{
//									}
//								}
//							}
//							edgeDest->setEdges(edges2add);
//						}
//					}
//				}
//
//			} // for(;itEdge!=itEdgee;itEdge++)
//		} // if(MeshDescriptor<M>::hasEdges && MeshDescriptor<M>::EKnowE)
//
//	} // if(MeshDescriptor<M>::hasEdges)
//
//	if(MeshDescriptor<M>::hasFaces){
//		typename IGMesh::faces_iterator itFace  = mesh_.faces_begin();
//		typename IGMesh::faces_iterator itFacee = mesh_.faces_end();
//
//		for(;itFace!=itFacee;itFace++){
//			bool toKeep = false;
//
//			if(!isSharedFace[(*itFace)->getID()]){
//				TInt currentOwner = faces_owner[(*itFace)->getID()];
//
//				if(currentOwner == mesh_.getPartID()){
//					toKeep = true;
//
//					faces_LID[(*itFace)->getID()] = (*itFace)->getID();
//				}
//				else{
//					TInt iSubMesh = partsOrder[currentOwner]-1;
//
//					/* connection F->N*/
//					if(MeshDescriptor<M>::FKnowN){
//						std::vector<TCellID> nodesID = (*itFace)->getNodeIDs();
//						std::vector<Node*> nodes;
//
//						std::vector<TCellID>::iterator itNodeID  = nodesID.begin();
//						std::vector<TCellID>::iterator itNodeIDe = nodesID.end();
//
//						for(;itNodeID!=itNodeIDe;itNodeID++){
//							if(!isSharedNode[*itNodeID]){
//								nodes.push_back(SubMeshes[iSubMesh].getLNode(nodes_LID[*itNodeID]));
//							}
//							else{
//								nodes.push_back(SubMeshes[iSubMesh].getLNode(nodes_owners2LIDs[nodes_LID[*itNodeID]][currentOwner]));
//							}
//						}
//
//						Face* faceDest = SubMeshes[iSubMesh].newFace(nodes);
//						faces_LID[(*itFace)->getID()] = faceDest->getID();
//					}
//
//					/* connection F->E*/
//					if(MeshDescriptor<M>::FKnowE){
//						std::vector<TCellID> edgesID = (*itFace)->getEdgeIDs();
//						std::vector<Edge*> edges;
//						edges.reserve(edgesID.size());
//						std::vector<TCellID>::iterator itEdgeID  = edgesID.begin();
//						std::vector<TCellID>::iterator itEdgeIDe = edgesID.end();
//						for(;itEdgeID!=itEdgeIDe;itEdgeID++){
//							if(!isSharedEdge[*itEdgeID]){
//								edges.push_back(SubMeshes[iSubMesh].getLEdge(edges_LID[*itEdgeID]));
//							}
//							else{
//								edges.push_back(SubMeshes[iSubMesh].getLEdge(edges_owners2LIDs[edges_LID[*itEdgeID]][currentOwner]));
//							}
//						}
//
//						if(MeshDescriptor<M>::FKnowN){
//							Face* faceDest = SubMeshes[iSubMesh].getLFace(faces_LID[(*itFace)->getID()]);
//							for(std::vector<Edge*>::iterator itEdge=edges.begin(); itEdge!=edges.end(); itEdge++){
//								faceDest->addEdge(*itEdge);
//							}
//							//faceDest->setEdges(edges);
//						}
//						else{
//							Face* faceDest = SubMeshes[iSubMesh].newFace((*itFace)->getType());
//							for(std::vector<Edge*>::iterator itEdge=edges.begin(); itEdge!=edges.end(); itEdge++){
//								faceDest->addEdge(*itEdge);
//							}
//							//faceDest->setEdges(edges);
//							faces_LID[(*itFace)->getID()] = faceDest->getID();
//						}
//					}
//				}
//			}
//			else{
//				std::map<TInt,TCellID>::iterator itFaceOwner2LID = faces_owners2LIDs[faces_LID[(*itFace)->getID()]].begin();
//				std::map<TInt,TCellID>::iterator itFaceOwner2LIDe = faces_owners2LIDs[faces_LID[(*itFace)->getID()]].end();
//
//				for(;itFaceOwner2LID!=itFaceOwner2LIDe;itFaceOwner2LID++){
//					TInt currentOwner = itFaceOwner2LID->first;
//
//					if(currentOwner == mesh_.getPartID()){
//						toKeep = true;
//
//						itFaceOwner2LID->second = (*itFace)->getID();
//
//						if(faces_owner[(*itFace)->getID()] == mesh_.getPartID()){
//							mesh_.setMaster((*itFace));
//						}
//						else{
//							mesh_.setSlave(*itFace,faces_owners2LIDs[faces_LID[(*itFace)->getID()]][faces_owner[(*itFace)->getID()]],faces_owner[(*itFace)->getID()]);
//						}
//					}
//					else{
//						TInt iSubMesh = partsOrder[currentOwner]-1;
//
//						/* connection F->N*/
//						if(MeshDescriptor<M>::FKnowN){
//							std::vector<TCellID> nodesID = (*itFace)->getNodeIDs();
//							std::vector<Node*> nodes;
//
//							std::vector<TCellID>::iterator itNodeID  = nodesID.begin();
//							std::vector<TCellID>::iterator itNodeIDe = nodesID.end();
//
//							for(;itNodeID!=itNodeIDe;itNodeID++){
//								if(!isSharedNode[*itNodeID]){
//									nodes.push_back(SubMeshes[iSubMesh].getLNode(nodes_LID[*itNodeID]));
//								}
//								else{
//									nodes.push_back(SubMeshes[iSubMesh].getLNode(nodes_owners2LIDs[nodes_LID[*itNodeID]][currentOwner]));
//								}
//							}
//
//							Face* faceDest = SubMeshes[iSubMesh].newFace(nodes);
//							faces_owners2LIDs[faces_LID[(*itFace)->getID()]][currentOwner] = faceDest->getID();
//
//							if(currentOwner == faces_owner[(*itFace)->getID()]){
//								SubMeshes[iSubMesh].setMaster((faceDest));
//							}
//							else{
//								SubMeshes[iSubMesh].setSlave(faceDest,faces_owners2LIDs[faces_LID[(*itFace)->getID()]][faces_owner[(*itFace)->getID()]],faces_owner[(*itFace)->getID()]);
//							}
//						}
//
//						/* connection F->E*/
//						if(MeshDescriptor<M>::FKnowE){
//							std::vector<TCellID> edgesID = (*itFace)->getEdgeIDs();
//							std::vector<Edge*> edges;
//							edges.reserve(edgesID.size());
//							std::vector<TCellID>::iterator itEdgeID  = edgesID.begin();
//							std::vector<TCellID>::iterator itEdgeIDe = edgesID.end();
//							for(;itEdgeID!=itEdgeIDe;itEdgeID++){
//								if(!isSharedEdge[*itEdgeID]){
//									edges.push_back(SubMeshes[iSubMesh].getLEdge(edges_LID[*itEdgeID]));
//								}
//								else{
//									edges.push_back(SubMeshes[iSubMesh].getLEdge(edges_owners2LIDs[edges_LID[*itEdgeID]][currentOwner]));
//								}
//							}
//
//							if(MeshDescriptor<M>::FKnowN){
//								Face* faceDest = SubMeshes[iSubMesh].getLFace(faces_owners2LIDs[faces_LID[(*itFace)->getID()]][currentOwner]);
//								for(std::vector<Edge*>::iterator itEdge=edges.begin(); itEdge!=edges.end(); itEdge++){
//									faceDest->addEdge(*itEdge);
//								}
//								//faceDest->setEdges(edges);
//							}
//							else{
//								Face* faceDest = SubMeshes[iSubMesh].newFace((*itFace)->getType());
//								for(std::vector<Edge*>::iterator itEdge=edges.begin(); itEdge!=edges.end(); itEdge++){
//									faceDest->addEdge(*itEdge);
//								}
//								//faceDest->setEdges(edges);
//								faces_owners2LIDs[faces_LID[(*itFace)->getID()]][currentOwner] = faceDest->getID();
//
//								if(currentOwner == faces_owner[(*itFace)->getID()]){
//									SubMeshes[iSubMesh].setMaster((faceDest));
//								}
//								else{
//									SubMeshes[iSubMesh].setSlave(faceDest,faces_owners2LIDs[faces_LID[(*itFace)->getID()]][faces_owner[(*itFace)->getID()]],faces_owner[(*itFace)->getID()]);
//								}
//							}
//						}
//					}
//				}
//			} // if(!isSharedFace[(*itFace)->getID()])
//
//			if(!toKeep){
//				//mesh_.deleteFace((*itFace));
//				mesh_.mark((*itFace),entityToDeleteMark);
//			}
//
//		} // for(;itFace!=itFacee;itFace++)
//
//		/* connection N->F*/
//		if(MeshDescriptor<M>::hasNodes && MeshDescriptor<M>::NKnowF){
//			typename IGMesh::nodes_iterator itNode  = mesh_.nodes_begin();
//			typename IGMesh::nodes_iterator itNodee = mesh_.nodes_end();
//
//			for(;itNode!=itNodee;itNode++){
//
//				std::vector<Face*> faces2add;
//				std::vector<TCellID> facesID = (*itNode)->getFaceIDs();
//				std::vector<TCellID>::iterator itFaceID;
//				std::vector<TCellID>::iterator itFaceIDb  = facesID.begin();
//				std::vector<TCellID>::iterator itFaceIDe = facesID.end();
//
//				if(!isSharedNode[(*itNode)->getID()]){
//					TInt currentOwner = nodes_owner[(*itNode)->getID()];
//
//					if(currentOwner == mesh_.getPartID()){
//						Node* nodeDest = mesh_.getLNode((*itNode)->getID());
//
//						for(faces2add.empty(),itFaceID=itFaceIDb;itFaceID!=itFaceIDe;itFaceID++){
//							if(!isSharedFace[*itFaceID]){
//								if(faces_owner[*itFaceID] == mesh_.getPartID()){
//									faces2add.push_back(mesh_.getLFace(*itFaceID));
//								}
//								else{
//								}
//							}
//							else{
//								if(faces_owners2LIDs[faces_LID[*itFaceID]].find(currentOwner) != faces_owners2LIDs[faces_LID[*itFaceID]].end()){
//									faces2add.push_back(mesh_.getLFace(*itFaceID));
//								}
//								else{
//								}
//							}
//						}
//						nodeDest->setFaces(faces2add);
//					}
//					else{
//						TInt iSubMesh = partsOrder[currentOwner]-1;
//
//						Node* nodeDest = SubMeshes[iSubMesh].getLNode(nodes_LID[(*itNode)->getID()]);
//
//						for(faces2add.empty(),itFaceID=itFaceIDb;itFaceID!=itFaceIDe;itFaceID++){
//							if(!isSharedFace[*itFaceID]){
//								if(faces_owner[*itFaceID] == currentOwner){
//									faces2add.push_back(SubMeshes[iSubMesh].getLFace(faces_LID[*itFaceID]));
//								}
//								else{
//								}
//							}
//							else{
//								if(faces_owners2LIDs[faces_LID[*itFaceID]].find(currentOwner) != faces_owners2LIDs[faces_LID[*itFaceID]].end()){
//									faces2add.push_back(SubMeshes[iSubMesh].getLFace(faces_owners2LIDs[faces_LID[*itFaceID]][currentOwner]));
//								}
//								else{
//								}
//							}
//						}
//						nodeDest->setFaces(faces2add);
//					}
//				}
//				else{
//					std::map<TInt,TCellID>::iterator itNodeOwner2LID = nodes_owners2LIDs[nodes_LID[(*itNode)->getID()]].begin();
//					std::map<TInt,TCellID>::iterator itNodeOwner2LIDe = nodes_owners2LIDs[nodes_LID[(*itNode)->getID()]].end();
//
//					for(;itNodeOwner2LID!=itNodeOwner2LIDe;itNodeOwner2LID++){
//						TInt currentOwner = itNodeOwner2LID->first;
//
//						if(currentOwner == mesh_.getPartID()){
//
//							Node* nodeDest = mesh_.getLNode(itNodeOwner2LID->second);
//
//							for(faces2add.empty(),itFaceID=itFaceIDb;itFaceID!=itFaceIDe;itFaceID++){
//								if(!isSharedFace[*itFaceID]){
//									if(faces_owner[*itFaceID] == mesh_.getPartID()){
//										faces2add.push_back(mesh_.getLFace(*itFaceID));
//									}
//									else{
//									}
//								}
//								else{
//									if(faces_owners2LIDs[faces_LID[*itFaceID]].find(currentOwner) != faces_owners2LIDs[faces_LID[*itFaceID]].end()){
//										faces2add.push_back(mesh_.getLFace(*itFaceID));
//									}
//									else{
//									}
//								}
//							}
//							nodeDest->setFaces(faces2add);
//						}
//						else{
//							TInt iSubMesh = partsOrder[currentOwner]-1;
//
//							Node* nodeDest = SubMeshes[iSubMesh].getLNode(nodes_owners2LIDs[nodes_LID[(*itNode)->getID()]][currentOwner]);
//
//							for(faces2add.empty(),itFaceID=itFaceIDb;itFaceID!=itFaceIDe;itFaceID++){
//								if(!isSharedFace[*itFaceID]){
//									if(faces_owner[*itFaceID] == currentOwner){
//										faces2add.push_back(SubMeshes[iSubMesh].getLFace(faces_LID[*itFaceID]));
//									}
//									else{
//									}
//								}
//								else{
//									if(faces_owners2LIDs[faces_LID[*itFaceID]].find(currentOwner) != faces_owners2LIDs[faces_LID[*itFaceID]].end()){
//										faces2add.push_back(SubMeshes[iSubMesh].getLFace(faces_owners2LIDs[faces_LID[*itFaceID]][currentOwner]));
//									}
//									else{
//									}
//								}
//							}
//							nodeDest->setFaces(faces2add);
//						}
//					}
//				}
//
//			} // for(;itNode!=itNodee;itNode++)
//		} // if(MeshDescriptor<M>::hasNodes && MeshDescriptor<M>::NKnowF)
//
//		/* connection E->F*/
//		if(MeshDescriptor<M>::hasEdges && MeshDescriptor<M>::EKnowF){
//			typename IGMesh::edges_iterator itEdge  = mesh_.edges_begin();
//			typename IGMesh::edges_iterator itEdgee = mesh_.edges_end();
//
//			for(;itEdge!=itEdgee;itEdge++){
//
//				std::vector<Face*> faces2add;
//				std::vector<TCellID> facesID = (*itEdge)->getFaceIDs();
//				std::vector<TCellID>::iterator itFaceID;
//				std::vector<TCellID>::iterator itFaceIDb = facesID.begin();
//				std::vector<TCellID>::iterator itFaceIDe = facesID.end();
//
//				if(!isSharedEdge[(*itEdge)->getID()]){
//					TInt currentOwner = edges_owner[(*itEdge)->getID()];
//
//					if(currentOwner == mesh_.getPartID()){
//						Edge* edgeDest = mesh_.getLEdge((*itEdge)->getID());
//
//						for(faces2add.empty(),itFaceID=itFaceIDb;itFaceID!=itFaceIDe;itFaceID++){
//							if(!isSharedFace[*itFaceID]){
//								if(faces_owner[*itFaceID] == mesh_.getPartID()){
//									faces2add.push_back(mesh_.getLFace(*itFaceID));
//								}
//								else{
//								}
//							}
//							else{
//								if(faces_owners2LIDs[faces_LID[*itFaceID]].find(currentOwner) != faces_owners2LIDs[faces_LID[*itFaceID]].end()){
//									faces2add.push_back(mesh_.getLFace(*itFaceID));
//								}
//								else{
//								}
//							}
//						}
//						edgeDest->setFaces(faces2add);
//					}
//					else{
//						TInt iSubMesh = partsOrder[currentOwner]-1;
//
//						Edge* edgeDest = SubMeshes[iSubMesh].getLEdge(edges_LID[(*itEdge)->getID()]);
//
//						for(faces2add.empty(),itFaceID=itFaceIDb;itFaceID!=itFaceIDe;itFaceID++){
//							if(!isSharedFace[*itFaceID]){
//								if(faces_owner[*itFaceID] == currentOwner){
//									faces2add.push_back(SubMeshes[iSubMesh].getLFace(faces_LID[*itFaceID]));
//								}
//								else{
//								}
//							}
//							else{
//								if(faces_owners2LIDs[faces_LID[*itFaceID]].find(currentOwner) != faces_owners2LIDs[faces_LID[*itFaceID]].end()){
//									faces2add.push_back(SubMeshes[iSubMesh].getLFace(faces_owners2LIDs[faces_LID[*itFaceID]][currentOwner]));
//								}
//								else{
//								}
//							}
//						}
//						edgeDest->setFaces(faces2add);
//					}
//				}
//				else{
//					std::map<TInt,TCellID>::iterator itEdgeOwner2LID = edges_owners2LIDs[edges_LID[(*itEdge)->getID()]].begin();
//					std::map<TInt,TCellID>::iterator itEdgeOwner2LIDe = edges_owners2LIDs[edges_LID[(*itEdge)->getID()]].end();
//
//					for(;itEdgeOwner2LID!=itEdgeOwner2LIDe;itEdgeOwner2LID++){
//						TInt currentOwner = itEdgeOwner2LID->first;
//
//						if(currentOwner == mesh_.getPartID()){
//
//							Edge* edgeDest = mesh_.getLEdge(itEdgeOwner2LID->second);
//
//							for(faces2add.empty(),itFaceID=itFaceIDb;itFaceID!=itFaceIDe;itFaceID++){
//								if(!isSharedFace[*itFaceID]){
//									if(faces_owner[*itFaceID] == mesh_.getPartID()){
//										faces2add.push_back(mesh_.getLFace(*itFaceID));
//									}
//									else{
//									}
//								}
//								else{
//									if(faces_owners2LIDs[faces_LID[*itFaceID]].find(currentOwner) != faces_owners2LIDs[faces_LID[*itFaceID]].end()){
//										faces2add.push_back(mesh_.getLFace(*itFaceID));
//									}
//									else{
//									}
//								}
//							}
//							edgeDest->setFaces(faces2add);
//						}
//						else{
//							TInt iSubMesh = partsOrder[currentOwner]-1;
//
//							Edge* edgeDest = SubMeshes[iSubMesh].getLEdge(edges_owners2LIDs[edges_LID[(*itEdge)->getID()]][currentOwner]);
//
//							for(faces2add.empty(),itFaceID=itFaceIDb;itFaceID!=itFaceIDe;itFaceID++){
//								if(!isSharedFace[*itFaceID]){
//									if(faces_owner[*itFaceID] == currentOwner){
//										faces2add.push_back(SubMeshes[iSubMesh].getLFace(faces_LID[*itFaceID]));
//									}
//									else{
//									}
//								}
//								else{
//									if(faces_owners2LIDs[faces_LID[*itFaceID]].find(currentOwner) != faces_owners2LIDs[faces_LID[*itFaceID]].end()){
//										faces2add.push_back(SubMeshes[iSubMesh].getLFace(faces_owners2LIDs[faces_LID[*itFaceID]][currentOwner]));
//									}
//									else{
//									}
//								}
//							}
//							edgeDest->setFaces(faces2add);
//						}
//					}
//				}
//
//			} // for(;itEdge!=itEdgee;itEdge++)
//		} // if(MeshDescriptor<M>::hasEdges && MeshDescriptor<M>::EKnowF)
//
//		/* connection F->F*/
//		if(MeshDescriptor<M>::hasFaces && MeshDescriptor<M>::FKnowF){
//			typename IGMesh::faces_iterator itFace  = mesh_.faces_begin();
//			typename IGMesh::faces_iterator itFacee = mesh_.faces_end();
//
//			for(;itFace!=itFacee;itFace++){
//
//				std::vector<Face*> faces2add;
//				std::vector<TCellID> facesID = (*itFace)->getFaceIDs();
//				std::vector<TCellID>::iterator itFaceID;
//				std::vector<TCellID>::iterator itFaceIDb = facesID.begin();
//				std::vector<TCellID>::iterator itFaceIDe = facesID.end();
//
//				if(!isSharedFace[(*itFace)->getID()]){
//					TInt currentOwner = faces_owner[(*itFace)->getID()];
//
//					if(currentOwner == mesh_.getPartID()){
//						Face* faceDest = mesh_.getLFace((*itFace)->getID());
//
//						for(faces2add.empty(),itFaceID=itFaceIDb;itFaceID!=itFaceIDe;itFaceID++){
//							if(!isSharedFace[*itFaceID]){
//								if(faces_owner[*itFaceID] == mesh_.getPartID()){
//									faces2add.push_back(mesh_.getLFace(*itFaceID));
//								}
//								else{
//								}
//							}
//							else{
//								if(faces_owners2LIDs[faces_LID[*itFaceID]].find(currentOwner) != faces_owners2LIDs[faces_LID[*itFaceID]].end()){
//									faces2add.push_back(mesh_.getLFace(*itFaceID));
//								}
//								else{
//								}
//							}
//						}
//						faceDest->setFaces(faces2add);
//					}
//					else{
//						TInt iSubMesh = partsOrder[currentOwner]-1;
//
//						Face* faceDest = SubMeshes[iSubMesh].getLFace(faces_LID[(*itFace)->getID()]);
//
//						for(faces2add.empty(),itFaceID=itFaceIDb;itFaceID!=itFaceIDe;itFaceID++){
//							if(!isSharedFace[*itFaceID]){
//								if(faces_owner[*itFaceID] == currentOwner){
//									faces2add.push_back(SubMeshes[iSubMesh].getLFace(faces_LID[*itFaceID]));
//								}
//								else{
//								}
//							}
//							else{
//								if(faces_owners2LIDs[faces_LID[*itFaceID]].find(currentOwner) != faces_owners2LIDs[faces_LID[*itFaceID]].end()){
//									faces2add.push_back(SubMeshes[iSubMesh].getLFace(faces_owners2LIDs[faces_LID[*itFaceID]][currentOwner]));
//								}
//								else{
//								}
//							}
//						}
//						faceDest->setFaces(faces2add);
//					}
//				}
//				else{
//					std::map<TInt,TCellID>::iterator itFaceOwner2LID = faces_owners2LIDs[faces_LID[(*itFace)->getID()]].begin();
//					std::map<TInt,TCellID>::iterator itFaceOwner2LIDe = faces_owners2LIDs[faces_LID[(*itFace)->getID()]].end();
//
//					for(;itFaceOwner2LID!=itFaceOwner2LIDe;itFaceOwner2LID++){
//						TInt currentOwner = itFaceOwner2LID->first;
//
//						if(currentOwner == mesh_.getPartID()){
//
//							Face* faceDest = mesh_.getLFace(itFaceOwner2LID->second);
//
//							for(faces2add.empty(),itFaceID=itFaceIDb;itFaceID!=itFaceIDe;itFaceID++){
//								if(!isSharedFace[*itFaceID]){
//									if(faces_owner[*itFaceID] == mesh_.getPartID()){
//										faces2add.push_back(mesh_.getLFace(*itFaceID));
//									}
//									else{
//									}
//								}
//								else{
//									if(faces_owners2LIDs[faces_LID[*itFaceID]].find(currentOwner) != faces_owners2LIDs[faces_LID[*itFaceID]].end()){
//										faces2add.push_back(mesh_.getLFace(*itFaceID));
//									}
//									else{
//									}
//								}
//							}
//							faceDest->setFaces(faces2add);
//						}
//						else{
//							TInt iSubMesh = partsOrder[currentOwner]-1;
//
//							Face* faceDest = SubMeshes[iSubMesh].getLFace(faces_owners2LIDs[faces_LID[(*itFace)->getID()]][currentOwner]);
//
//							for(faces2add.empty(),itFaceID=itFaceIDb;itFaceID!=itFaceIDe;itFaceID++){
//								if(!isSharedFace[*itFaceID]){
//									if(faces_owner[*itFaceID] == currentOwner){
//										faces2add.push_back(SubMeshes[iSubMesh].getLFace(faces_LID[*itFaceID]));
//									}
//									else{
//									}
//								}
//								else{
//									if(faces_owners2LIDs[faces_LID[*itFaceID]].find(currentOwner) != faces_owners2LIDs[faces_LID[*itFaceID]].end()){
//										faces2add.push_back(SubMeshes[iSubMesh].getLFace(faces_owners2LIDs[faces_LID[*itFaceID]][currentOwner]));
//									}
//									else{
//									}
//								}
//							}
//							faceDest->setFaces(faces2add);
//						}
//					}
//				}
//
//			} // for(;itFace!=itFacee;itFace++)
//		} // if(MeshDescriptor<M>::hasFaces && MeshDescriptor<M>::FKnowF)
//
//	} // if(MeshDescriptor<M>::hasFaces)
//
//	if(MeshDescriptor<M>::hasRegions){
//		typename IGMesh::regions_iterator itRegion  = mesh_.regions_begin();
//		typename IGMesh::regions_iterator itRegione = mesh_.regions_end();
//
//		for(;itRegion!=itRegione;itRegion++){
//			bool toKeep = false;
//
//			if(!isSharedRegion[(*itRegion)->getID()]){
//				TInt currentOwner = regions_owner[(*itRegion)->getID()];
//
//				if(currentOwner == mesh_.getPartID()){
//					toKeep = true;
//
//					regions_LID[(*itRegion)->getID()] = (*itRegion)->getID();
//				}
//				else{
//					TInt iSubMesh = partsOrder[currentOwner]-1;
//
//					/* connection R->N*/
//					if(MeshDescriptor<M>::RKnowN){
//						std::vector<TCellID> nodesID = (*itRegion)->getNodeIDs();
//						std::vector<Node*> nodes;
//
//						std::vector<TCellID>::iterator itNodeID  = nodesID.begin();
//						std::vector<TCellID>::iterator itNodeIDe = nodesID.end();
//
//						for(;itNodeID!=itNodeIDe;itNodeID++){
//							if(!isSharedNode[*itNodeID]){
//								nodes.push_back(SubMeshes[iSubMesh].getLNode(nodes_LID[*itNodeID]));
//							}
//							else{
//								nodes.push_back(SubMeshes[iSubMesh].getLNode(nodes_owners2LIDs[nodes_LID[*itNodeID]][currentOwner]));
//							}
//						}
//
//						Region* regionDest = SubMeshes[iSubMesh].newRegion((*itRegion)->getType(),nodes);
//						regions_LID[(*itRegion)->getID()] = regionDest->getID();
//					}
//
//					/* connection R->E*/
//					if(MeshDescriptor<M>::RKnowE){
//						std::vector<TCellID> edgesID = (*itRegion)->getEdgeIDs();
//						std::vector<Edge*> edges;
//						edges.reserve(edgesID.size());
//						std::vector<TCellID>::iterator itEdgeID  = edgesID.begin();
//						std::vector<TCellID>::iterator itEdgeIDe = edgesID.end();
//						for(;itEdgeID!=itEdgeIDe;itEdgeID++){
//							if(!isSharedEdge[*itEdgeID]){
//								edges.push_back(SubMeshes[iSubMesh].getLEdge(edges_LID[*itEdgeID]));
//							}
//							else{
//								edges.push_back(SubMeshes[iSubMesh].getLEdge(edges_owners2LIDs[edges_LID[*itEdgeID]][currentOwner]));
//							}
//						}
//
//						if(MeshDescriptor<M>::RKnowN){
//							Region* regionDest = SubMeshes[iSubMesh].getLRegion(regions_LID[(*itRegion)->getID()]);
//							for(std::vector<Edge*>::iterator itEdge=edges.begin(); itEdge!=edges.end(); itEdge++){
//								regionDest->addEdge(*itEdge);
//							}
//							//faceDest->setEdges(edges);
//						}
//						else{
//							Region* regionDest = SubMeshes[iSubMesh].newRegion((*itRegion)->getType());
//							for(std::vector<Edge*>::iterator itEdge=edges.begin(); itEdge!=edges.end(); itEdge++){
//								regionDest->addEdge(*itEdge);
//							}
//							//faceDest->setEdges(edges);
//							regions_LID[(*itRegion)->getID()] = regionDest->getID();
//						}
//					}
//
//					/* connection R->F*/
//					if(MeshDescriptor<M>::RKnowF){
//						std::vector<TCellID> facesID = (*itRegion)->getFaceIDs();
//						std::vector<Face*> faces;
//						faces.reserve(facesID.size());
//						std::vector<TCellID>::iterator itFaceID  = facesID.begin();
//						std::vector<TCellID>::iterator itFaceIDe = facesID.end();
//						for(;itFaceID!=itFaceIDe;itFaceID++){
//							if(!isSharedFace[*itFaceID]){
//								faces.push_back(SubMeshes[iSubMesh].getLFace(faces_LID[*itFaceID]));
//							}
//							else{
//								faces.push_back(SubMeshes[iSubMesh].getLFace(faces_owners2LIDs[faces_LID[*itFaceID]][currentOwner]));
//							}
//						}
//
//						if(MeshDescriptor<M>::RKnowN || MeshDescriptor<M>::RKnowE){
//							Region* regionDest = SubMeshes[iSubMesh].getLRegion(regions_LID[(*itRegion)->getID()]);
//							for(std::vector<Face*>::iterator itFace=faces.begin(); itFace!=faces.end(); itFace++){
//								regionDest->addFace(*itFace);
//							}
//							//faceDest->setEdges(edges);
//						}
//						else{
//							Region* regionDest = SubMeshes[iSubMesh].newRegion((*itRegion)->getType());
//							for(std::vector<Face*>::iterator itFace=faces.begin(); itFace!=faces.end(); itFace++){
//								regionDest->addFace(*itFace);
//							}
//							//faceDest->setEdges(edges);
//							regions_LID[(*itRegion)->getID()] = regionDest->getID();
//						}
//					}
//				}
//			}
//			else{
//				std::map<TInt,TCellID>::iterator itRegionOwner2LID = regions_owners2LIDs[regions_LID[(*itRegion)->getID()]].begin();
//				std::map<TInt,TCellID>::iterator itRegionOwner2LIDe = regions_owners2LIDs[regions_LID[(*itRegion)->getID()]].end();
//
//				for(;itRegionOwner2LID!=itRegionOwner2LIDe;itRegionOwner2LID++){
//					TInt currentOwner = itRegionOwner2LID->first;
//
//					if(currentOwner == mesh_.getPartID()){
//						toKeep = true;
//
//						itRegionOwner2LID->second = (*itRegion)->getID();
//
//						if(regions_owner[(*itRegion)->getID()] == mesh_.getPartID()){
//							mesh_.setMaster((*itRegion));
//						}
//						else{
//							mesh_.setSlave(*itRegion,regions_owners2LIDs[regions_LID[(*itRegion)->getID()]][regions_owner[(*itRegion)->getID()]],regions_owner[(*itRegion)->getID()]);
//						}
//					}
//					else{
//						TInt iSubMesh = partsOrder[currentOwner]-1;
//
//						/* connection R->N*/
//						if(MeshDescriptor<M>::RKnowN){
//							std::vector<TCellID> nodesID = (*itRegion)->getNodeIDs();
//							std::vector<Node*> nodes;
//
//							std::vector<TCellID>::iterator itNodeID  = nodesID.begin();
//							std::vector<TCellID>::iterator itNodeIDe = nodesID.end();
//
//							for(;itNodeID!=itNodeIDe;itNodeID++){
//								if(!isSharedNode[*itNodeID]){
//									nodes.push_back(SubMeshes[iSubMesh].getLNode(nodes_LID[*itNodeID]));
//								}
//								else{
//									nodes.push_back(SubMeshes[iSubMesh].getLNode(nodes_owners2LIDs[nodes_LID[*itNodeID]][currentOwner]));
//								}
//							}
//
//							Region* regionDest = SubMeshes[iSubMesh].newRegion((*itRegion)->getType(),nodes);
//							regions_owners2LIDs[regions_LID[(*itRegion)->getID()]][currentOwner] = regionDest->getID();
//
//							if(currentOwner == regions_owner[(*itRegion)->getID()]){
//								SubMeshes[iSubMesh].setMaster((regionDest));
//							}
//							else{
//								SubMeshes[iSubMesh].setSlave(regionDest,regions_owners2LIDs[regions_LID[(*itRegion)->getID()]][regions_owner[(*itRegion)->getID()]],regions_owner[(*itRegion)->getID()]);
//							}
//						}
//
//						/* connection R->E*/
//						if(MeshDescriptor<M>::RKnowE){
//							std::vector<TCellID> edgesID = (*itRegion)->getEdgeIDs();
//							std::vector<Edge*> edges;
//							edges.reserve(edgesID.size());
//							std::vector<TCellID>::iterator itEdgeID  = edgesID.begin();
//							std::vector<TCellID>::iterator itEdgeIDe = edgesID.end();
//							for(;itEdgeID!=itEdgeIDe;itEdgeID++){
//								if(!isSharedEdge[*itEdgeID]){
//									edges.push_back(SubMeshes[iSubMesh].getLEdge(edges_LID[*itEdgeID]));
//								}
//								else{
//									edges.push_back(SubMeshes[iSubMesh].getLEdge(edges_owners2LIDs[edges_LID[*itEdgeID]][currentOwner]));
//								}
//							}
//
//							if(MeshDescriptor<M>::RKnowN){
//								Region* regionDest = SubMeshes[iSubMesh].getLRegion(regions_owners2LIDs[regions_LID[(*itRegion)->getID()]][currentOwner]);
//								for(std::vector<Edge*>::iterator itEdge=edges.begin(); itEdge!=edges.end(); itEdge++){
//									regionDest->addEdge(*itEdge);
//								}
//								//faceDest->setEdges(edges);
//							}
//							else{
//								Region* regionDest = SubMeshes[iSubMesh].newRegion((*itRegion)->getType());
//								for(std::vector<Edge*>::iterator itEdge=edges.begin(); itEdge!=edges.end(); itEdge++){
//									regionDest->addEdge(*itEdge);
//								}
//								//faceDest->setEdges(edges);
//								regions_owners2LIDs[regions_LID[(*itRegion)->getID()]][currentOwner] = regionDest->getID();
//
//								if(currentOwner == regions_owner[(*itRegion)->getID()]){
//									SubMeshes[iSubMesh].setMaster((regionDest));
//								}
//								else{
//									SubMeshes[iSubMesh].setSlave(regionDest,regions_owners2LIDs[regions_LID[(*itRegion)->getID()]][regions_owner[(*itRegion)->getID()]],regions_owner[(*itRegion)->getID()]);
//								}
//							}
//						}
//
//						/* connection R->F*/
//						if(MeshDescriptor<M>::RKnowF){
//							std::vector<TCellID> facesID = (*itRegion)->getFaceIDs();
//							std::vector<Face*> faces;
//							faces.reserve(facesID.size());
//							std::vector<TCellID>::iterator itFaceID  = facesID.begin();
//							std::vector<TCellID>::iterator itFaceIDe = facesID.end();
//							for(;itFaceID!=itFaceIDe;itFaceID++){
//								if(!isSharedFace[*itFaceID]){
//									faces.push_back(SubMeshes[iSubMesh].getLFace(faces_LID[*itFaceID]));
//								}
//								else{
//									faces.push_back(SubMeshes[iSubMesh].getLFace(faces_owners2LIDs[faces_LID[*itFaceID]][currentOwner]));
//								}
//							}
//
//							if(MeshDescriptor<M>::RKnowN || MeshDescriptor<M>::RKnowE){
//								Region* regionDest = SubMeshes[iSubMesh].getLRegion(regions_owners2LIDs[regions_LID[(*itRegion)->getID()]][currentOwner]);
//								for(std::vector<Face*>::iterator itFace=faces.begin(); itFace!=faces.end(); itFace++){
//									regionDest->addFace(*itFace);
//								}
//								//faceDest->setEdges(edges);
//							}
//							else{
//								Region* regionDest = SubMeshes[iSubMesh].newRegion((*itRegion)->getType());
//								for(std::vector<Face*>::iterator itFace=faces.begin(); itFace!=faces.end(); itFace++){
//									regionDest->addFace(*itFace);
//								}
//								//faceDest->setEdges(edges);
//								regions_owners2LIDs[regions_LID[(*itRegion)->getID()]][currentOwner] = regionDest->getID();
//
//								if(currentOwner == regions_owner[(*itRegion)->getID()]){
//									SubMeshes[iSubMesh].setMaster((regionDest));
//								}
//								else{
//									SubMeshes[iSubMesh].setSlave(regionDest,regions_owners2LIDs[regions_LID[(*itRegion)->getID()]][regions_owner[(*itRegion)->getID()]],regions_owner[(*itRegion)->getID()]);
//								}
//							}
//						}
//					}
//				}
//			} // if(!isSharedRegion[(*itRegion)->getID()])
//
//			if(!toKeep){
//				//mesh_.deleteRegion((*itRegion));
//				mesh_.mark((*itRegion),entityToDeleteMark);
//			}
//
//		} // for(;itRegion!=itRegione;itRegion++)
//
//		/* connection N->R*/
//		if(MeshDescriptor<M>::hasNodes && MeshDescriptor<M>::NKnowR){
//			typename IGMesh::nodes_iterator itNode  = mesh_.nodes_begin();
//			typename IGMesh::nodes_iterator itNodee = mesh_.nodes_end();
//
//			for(;itNode!=itNodee;itNode++){
//
//				std::vector<Region*> regions2add;
//				std::vector<TCellID> regionsID = (*itNode)->getRegionIDs();
//				std::vector<TCellID>::iterator itRegionID;
//				std::vector<TCellID>::iterator itRegionIDb = regionsID.begin();
//				std::vector<TCellID>::iterator itRegionIDe = regionsID.end();
//
//				if(!isSharedNode[(*itNode)->getID()]){
//					TInt currentOwner = nodes_owner[(*itNode)->getID()];
//
//					if(currentOwner == mesh_.getPartID()){
//						Node* nodeDest = mesh_.getLNode((*itNode)->getID());
//
//						for(regions2add.empty(),itRegionID=itRegionIDb;itRegionID!=itRegionIDe;itRegionID++){
//							if(!isSharedRegion[*itRegionID]){
//								if(regions_owner[*itRegionID] == mesh_.getPartID()){
//									regions2add.push_back(mesh_.getLRegion(*itRegionID));
//								}
//								else{
//								}
//							}
//							else{
//								if(regions_owners2LIDs[regions_LID[*itRegionID]].find(currentOwner) != regions_owners2LIDs[regions_LID[*itRegionID]].end()){
//									regions2add.push_back(mesh_.getLRegion(*itRegionID));
//								}
//								else{
//								}
//							}
//						}
//						nodeDest->setRegions(regions2add);
//					}
//					else{
//						TInt iSubMesh = partsOrder[currentOwner]-1;
//
//						Node* nodeDest = SubMeshes[iSubMesh].getLNode(nodes_LID[(*itNode)->getID()]);
//
//						for(regions2add.empty(),itRegionID=itRegionIDb;itRegionID!=itRegionIDe;itRegionID++){
//							if(!isSharedRegion[*itRegionID]){
//								if(regions_owner[*itRegionID] == currentOwner){
//									regions2add.push_back(SubMeshes[iSubMesh].getLRegion(regions_LID[*itRegionID]));
//								}
//								else{
//								}
//							}
//							else{
//								if(regions_owners2LIDs[regions_LID[*itRegionID]].find(currentOwner) != regions_owners2LIDs[regions_LID[*itRegionID]].end()){
//									regions2add.push_back(SubMeshes[iSubMesh].getLRegion(regions_owners2LIDs[regions_LID[*itRegionID]][currentOwner]));
//								}
//								else{
//								}
//							}
//						}
//						nodeDest->setRegions(regions2add);
//					}
//				}
//				else{
//					std::map<TInt,TCellID>::iterator itNodeOwner2LID = nodes_owners2LIDs[nodes_LID[(*itNode)->getID()]].begin();
//					std::map<TInt,TCellID>::iterator itNodeOwner2LIDe = nodes_owners2LIDs[nodes_LID[(*itNode)->getID()]].end();
//
//					for(;itNodeOwner2LID!=itNodeOwner2LIDe;itNodeOwner2LID++){
//						TInt currentOwner = itNodeOwner2LID->first;
//
//						if(currentOwner == mesh_.getPartID()){
//
//							Node* nodeDest = mesh_.getLNode(itNodeOwner2LID->second);
//
//							for(regions2add.empty(),itRegionID=itRegionIDb;itRegionID!=itRegionIDe;itRegionID++){
//								if(!isSharedRegion[*itRegionID]){
//									if(regions_owner[*itRegionID] == mesh_.getPartID()){
//										regions2add.push_back(mesh_.getLRegion(*itRegionID));
//									}
//									else{
//									}
//								}
//								else{
//									if(regions_owners2LIDs[regions_LID[*itRegionID]].find(currentOwner) != regions_owners2LIDs[regions_LID[*itRegionID]].end()){
//										regions2add.push_back(mesh_.getLRegion(*itRegionID));
//									}
//									else{
//									}
//								}
//							}
//							nodeDest->setRegions(regions2add);
//						}
//						else{
//							TInt iSubMesh = partsOrder[currentOwner]-1;
//
//							Node* nodeDest = SubMeshes[iSubMesh].getLNode(nodes_owners2LIDs[nodes_LID[(*itNode)->getID()]][currentOwner]);
//
//							for(regions2add.empty(),itRegionID=itRegionIDb;itRegionID!=itRegionIDe;itRegionID++){
//								if(!isSharedRegion[*itRegionID]){
//									if(regions_owner[*itRegionID] == currentOwner){
//										regions2add.push_back(SubMeshes[iSubMesh].getLRegion(regions_LID[*itRegionID]));
//									}
//									else{
//									}
//								}
//								else{
//									if(regions_owners2LIDs[regions_LID[*itRegionID]].find(currentOwner) != regions_owners2LIDs[regions_LID[*itRegionID]].end()){
//										regions2add.push_back(SubMeshes[iSubMesh].getLRegion(regions_owners2LIDs[regions_LID[*itRegionID]][currentOwner]));
//									}
//									else{
//									}
//								}
//							}
//							nodeDest->setRegions(regions2add);
//						}
//					}
//				}
//
//			} // for(;itNode!=itNodee;itNode++)
//		} // if(MeshDescriptor<M>::hasNodes && MeshDescriptor<M>::NKnowR)
//
//		/* connection E->R*/
//		if(MeshDescriptor<M>::hasEdges && MeshDescriptor<M>::EKnowR){
//			typename IGMesh::edges_iterator itEdge  = mesh_.edges_begin();
//			typename IGMesh::edges_iterator itEdgee = mesh_.edges_end();
//
//			for(;itEdge!=itEdgee;itEdge++){
//
//				std::vector<Region*> regions2add;
//				std::vector<TCellID> regionsID = (*itEdge)->getRegionIDs();
//				std::vector<TCellID>::iterator itRegionID;
//				std::vector<TCellID>::iterator itRegionIDb = regionsID.begin();
//				std::vector<TCellID>::iterator itRegionIDe = regionsID.end();
//
//				if(!isSharedEdge[(*itEdge)->getID()]){
//					TInt currentOwner = edges_owner[(*itEdge)->getID()];
//
//					if(currentOwner == mesh_.getPartID()){
//						Edge* edgeDest = mesh_.getLEdge((*itEdge)->getID());
//
//						for(regions2add.empty(),itRegionID=itRegionIDb;itRegionID!=itRegionIDe;itRegionID++){
//							if(!isSharedRegion[*itRegionID]){
//								if(regions_owner[*itRegionID] == mesh_.getPartID()){
//									regions2add.push_back(mesh_.getLRegion(*itRegionID));
//								}
//								else{
//								}
//							}
//							else{
//								if(regions_owners2LIDs[regions_LID[*itRegionID]].find(currentOwner) != regions_owners2LIDs[regions_LID[*itRegionID]].end()){
//									regions2add.push_back(mesh_.getLRegion(*itRegionID));
//								}
//								else{
//								}
//							}
//						}
//						edgeDest->setRegions(regions2add);
//					}
//					else{
//						TInt iSubMesh = partsOrder[currentOwner]-1;
//
//						Edge* edgeDest = SubMeshes[iSubMesh].getLEdge(edges_LID[(*itEdge)->getID()]);
//
//						for(regions2add.empty(),itRegionID=itRegionIDb;itRegionID!=itRegionIDe;itRegionID++){
//							if(!isSharedRegion[*itRegionID]){
//								if(regions_owner[*itRegionID] == currentOwner){
//									regions2add.push_back(SubMeshes[iSubMesh].getLRegion(regions_LID[*itRegionID]));
//								}
//								else{
//								}
//							}
//							else{
//								if(regions_owners2LIDs[regions_LID[*itRegionID]].find(currentOwner) != regions_owners2LIDs[regions_LID[*itRegionID]].end()){
//									regions2add.push_back(SubMeshes[iSubMesh].getLRegion(regions_owners2LIDs[regions_LID[*itRegionID]][currentOwner]));
//								}
//								else{
//								}
//							}
//						}
//						edgeDest->setRegions(regions2add);
//					}
//				}
//				else{
//					std::map<TInt,TCellID>::iterator itEdgeOwner2LID = edges_owners2LIDs[edges_LID[(*itEdge)->getID()]].begin();
//					std::map<TInt,TCellID>::iterator itEdgeOwner2LIDe = edges_owners2LIDs[edges_LID[(*itEdge)->getID()]].end();
//
//					for(;itEdgeOwner2LID!=itEdgeOwner2LIDe;itEdgeOwner2LID++){
//						TInt currentOwner = itEdgeOwner2LID->first;
//
//						if(currentOwner == mesh_.getPartID()){
//
//							Edge* edgeDest = mesh_.getLEdge(itEdgeOwner2LID->second);
//
//							for(regions2add.empty(),itRegionID=itRegionIDb;itRegionID!=itRegionIDe;itRegionID++){
//								if(!isSharedRegion[*itRegionID]){
//									if(regions_owner[*itRegionID] == mesh_.getPartID()){
//										regions2add.push_back(mesh_.getLRegion(*itRegionID));
//									}
//									else{
//									}
//								}
//								else{
//									if(regions_owners2LIDs[regions_LID[*itRegionID]].find(currentOwner) != regions_owners2LIDs[regions_LID[*itRegionID]].end()){
//										regions2add.push_back(mesh_.getLRegion(*itRegionID));
//									}
//									else{
//									}
//								}
//							}
//							edgeDest->setRegions(regions2add);
//						}
//						else{
//							TInt iSubMesh = partsOrder[currentOwner]-1;
//
//							Edge* edgeDest = SubMeshes[iSubMesh].getLEdge(edges_owners2LIDs[edges_LID[(*itEdge)->getID()]][currentOwner]);
//
//							for(regions2add.empty(),itRegionID=itRegionIDb;itRegionID!=itRegionIDe;itRegionID++){
//								if(!isSharedRegion[*itRegionID]){
//									if(regions_owner[*itRegionID] == currentOwner){
//										regions2add.push_back(SubMeshes[iSubMesh].getLRegion(regions_LID[*itRegionID]));
//									}
//									else{
//									}
//								}
//								else{
//									if(regions_owners2LIDs[regions_LID[*itRegionID]].find(currentOwner) != regions_owners2LIDs[regions_LID[*itRegionID]].end()){
//										regions2add.push_back(SubMeshes[iSubMesh].getLRegion(regions_owners2LIDs[regions_LID[*itRegionID]][currentOwner]));
//									}
//									else{
//									}
//								}
//							}
//							edgeDest->setRegions(regions2add);
//						}
//					}
//				}
//
//			} // for(;itEdge!=itEdgee;itEdge++)
//		} // if(MeshDescriptor<M>::hasEdges && MeshDescriptor<M>::EKnowR)
//
//		/* connection F->R*/
//		if(MeshDescriptor<M>::hasFaces && MeshDescriptor<M>::FKnowR){
//			typename IGMesh::faces_iterator itFace  = mesh_.faces_begin();
//			typename IGMesh::faces_iterator itFacee = mesh_.faces_end();
//
//			for(;itFace!=itFacee;itFace++){
//
//				std::vector<Region*> regions2add;
//				std::vector<TCellID> regionsID = (*itFace)->getRegionIDs();
//				std::vector<TCellID>::iterator itRegionID;
//				std::vector<TCellID>::iterator itRegionIDb = regionsID.begin();
//				std::vector<TCellID>::iterator itRegionIDe = regionsID.end();
//
//				if(!isSharedFace[(*itFace)->getID()]){
//					TInt currentOwner = faces_owner[(*itFace)->getID()];
//
//					if(currentOwner == mesh_.getPartID()){
//						Face* faceDest = mesh_.getLFace((*itFace)->getID());
//
//						for(regions2add.empty(),itRegionID=itRegionIDb;itRegionID!=itRegionIDe;itRegionID++){
//							if(!isSharedRegion[*itRegionID]){
//								if(regions_owner[*itRegionID] == mesh_.getPartID()){
//									regions2add.push_back(mesh_.getLRegion(*itRegionID));
//								}
//								else{
//								}
//							}
//							else{
//								if(regions_owners2LIDs[regions_LID[*itRegionID]].find(currentOwner) != regions_owners2LIDs[regions_LID[*itRegionID]].end()){
//									regions2add.push_back(mesh_.getLRegion(*itRegionID));
//								}
//								else{
//								}
//							}
//						}
//						faceDest->setRegions(regions2add);
//					}
//					else{
//						TInt iSubMesh = partsOrder[currentOwner]-1;
//
//						Face* faceDest = SubMeshes[iSubMesh].getLFace(faces_LID[(*itFace)->getID()]);
//
//						for(regions2add.empty(),itRegionID=itRegionIDb;itRegionID!=itRegionIDe;itRegionID++){
//							if(!isSharedRegion[*itRegionID]){
//								if(regions_owner[*itRegionID] == currentOwner){
//									regions2add.push_back(SubMeshes[iSubMesh].getLRegion(regions_LID[*itRegionID]));
//								}
//								else{
//								}
//							}
//							else{
//								if(regions_owners2LIDs[regions_LID[*itRegionID]].find(currentOwner) != regions_owners2LIDs[regions_LID[*itRegionID]].end()){
//									regions2add.push_back(SubMeshes[iSubMesh].getLRegion(regions_owners2LIDs[regions_LID[*itRegionID]][currentOwner]));
//								}
//								else{
//								}
//							}
//						}
//						faceDest->setRegions(regions2add);
//					}
//				}
//				else{
//					std::map<TInt,TCellID>::iterator itFaceOwner2LID = faces_owners2LIDs[faces_LID[(*itFace)->getID()]].begin();
//					std::map<TInt,TCellID>::iterator itFaceOwner2LIDe = faces_owners2LIDs[faces_LID[(*itFace)->getID()]].end();
//
//					for(;itFaceOwner2LID!=itFaceOwner2LIDe;itFaceOwner2LID++){
//						TInt currentOwner = itFaceOwner2LID->first;
//
//						if(currentOwner == mesh_.getPartID()){
//
//							Face* faceDest = mesh_.getLFace(itFaceOwner2LID->second);
//
//							for(regions2add.empty(),itRegionID=itRegionIDb;itRegionID!=itRegionIDe;itRegionID++){
//								if(!isSharedRegion[*itRegionID]){
//									if(regions_owner[*itRegionID] == mesh_.getPartID()){
//										regions2add.push_back(mesh_.getLRegion(*itRegionID));
//									}
//									else{
//									}
//								}
//								else{
//									if(regions_owners2LIDs[regions_LID[*itRegionID]].find(currentOwner) != regions_owners2LIDs[regions_LID[*itRegionID]].end()){
//										regions2add.push_back(mesh_.getLRegion(*itRegionID));
//									}
//									else{
//									}
//								}
//							}
//							faceDest->setRegions(regions2add);
//						}
//						else{
//							TInt iSubMesh = partsOrder[currentOwner]-1;
//
//							Face* faceDest = SubMeshes[iSubMesh].getLFace(faces_owners2LIDs[faces_LID[(*itFace)->getID()]][currentOwner]);
//
//							for(regions2add.empty(),itRegionID=itRegionIDb;itRegionID!=itRegionIDe;itRegionID++){
//								if(!isSharedRegion[*itRegionID]){
//									if(regions_owner[*itRegionID] == currentOwner){
//										regions2add.push_back(SubMeshes[iSubMesh].getLRegion(regions_LID[*itRegionID]));
//									}
//									else{
//									}
//								}
//								else{
//									if(regions_owners2LIDs[regions_LID[*itRegionID]].find(currentOwner) != regions_owners2LIDs[regions_LID[*itRegionID]].end()){
//										regions2add.push_back(SubMeshes[iSubMesh].getLRegion(regions_owners2LIDs[regions_LID[*itRegionID]][currentOwner]));
//									}
//									else{
//									}
//								}
//							}
//							faceDest->setRegions(regions2add);
//						}
//					}
//				}
//
//			} // for(;itFace!=itFacee;itFace++)
//		} // if(MeshDescriptor<M>::hasFaces && MeshDescriptor<M>::FKnowR)
//
//		/* connection R->R*/
//		if(MeshDescriptor<M>::hasRegions && MeshDescriptor<M>::RKnowR){
//			typename IGMesh::regions_iterator itRegion  = mesh_.regions_begin();
//			typename IGMesh::regions_iterator itRegione = mesh_.regions_end();
//
//			for(;itRegion!=itRegione;itRegion++){
//
//				std::vector<Region*> regions2add;
//				std::vector<TCellID> regionsID = (*itRegion)->getRegionIDs();
//				std::vector<TCellID>::iterator itRegionID;
//				std::vector<TCellID>::iterator itRegionIDb = regionsID.begin();
//				std::vector<TCellID>::iterator itRegionIDe = regionsID.end();
//
//				if(!isSharedRegion[(*itRegion)->getID()]){
//					TInt currentOwner = regions_owner[(*itRegion)->getID()];
//
//					if(currentOwner == mesh_.getPartID()){
//						Region* regionDest = mesh_.getLRegion((*itRegion)->getID());
//
//						for(regions2add.empty(),itRegionID=itRegionIDb;itRegionID!=itRegionIDe;itRegionID++){
//							if(!isSharedRegion[*itRegionID]){
//								if(regions_owner[*itRegionID] == mesh_.getPartID()){
//									regions2add.push_back(mesh_.getLRegion(*itRegionID));
//								}
//								else{
//								}
//							}
//							else{
//								if(regions_owners2LIDs[regions_LID[*itRegionID]].find(currentOwner) != regions_owners2LIDs[regions_LID[*itRegionID]].end()){
//									regions2add.push_back(mesh_.getLRegion(*itRegionID));
//								}
//								else{
//								}
//							}
//						}
//						regionDest->setRegions(regions2add);
//					}
//					else{
//						TInt iSubMesh = partsOrder[currentOwner]-1;
//
//						Region* regionDest = SubMeshes[iSubMesh].getLRegion(regions_LID[(*itRegion)->getID()]);
//
//						for(regions2add.empty(),itRegionID=itRegionIDb;itRegionID!=itRegionIDe;itRegionID++){
//							if(!isSharedRegion[*itRegionID]){
//								if(regions_owner[*itRegionID] == currentOwner){
//									regions2add.push_back(SubMeshes[iSubMesh].getLRegion(regions_LID[*itRegionID]));
//								}
//								else{
//								}
//							}
//							else{
//								if(regions_owners2LIDs[regions_LID[*itRegionID]].find(currentOwner) != regions_owners2LIDs[regions_LID[*itRegionID]].end()){
//									regions2add.push_back(SubMeshes[iSubMesh].getLRegion(regions_owners2LIDs[regions_LID[*itRegionID]][currentOwner]));
//								}
//								else{
//								}
//							}
//						}
//						regionDest->setRegions(regions2add);
//					}
//				}
//				else{
//					std::map<TInt,TCellID>::iterator itRegionOwner2LID = regions_owners2LIDs[regions_LID[(*itRegion)->getID()]].begin();
//					std::map<TInt,TCellID>::iterator itRegionOwner2LIDe = regions_owners2LIDs[regions_LID[(*itRegion)->getID()]].end();
//
//					for(;itRegionOwner2LID!=itRegionOwner2LIDe;itRegionOwner2LID++){
//						TInt currentOwner = itRegionOwner2LID->first;
//
//						if(currentOwner == mesh_.getPartID()){
//
//							Region* regionDest = mesh_.getLRegion(itRegionOwner2LID->second);
//
//							for(regions2add.empty(),itRegionID=itRegionIDb;itRegionID!=itRegionIDe;itRegionID++){
//								if(!isSharedRegion[*itRegionID]){
//									if(regions_owner[*itRegionID] == mesh_.getPartID()){
//										regions2add.push_back(mesh_.getLRegion(*itRegionID));
//									}
//									else{
//									}
//								}
//								else{
//									if(regions_owners2LIDs[regions_LID[*itRegionID]].find(currentOwner) != regions_owners2LIDs[regions_LID[*itRegionID]].end()){
//										regions2add.push_back(mesh_.getLRegion(*itRegionID));
//									}
//									else{
//									}
//								}
//							}
//							regionDest->setRegions(regions2add);
//						}
//						else{
//							TInt iSubMesh = partsOrder[currentOwner]-1;
//
//							Region* regionDest = SubMeshes[iSubMesh].getLRegion(regions_owners2LIDs[regions_LID[(*itRegion)->getID()]][currentOwner]);
//
//							for(regions2add.empty(),itRegionID=itRegionIDb;itRegionID!=itRegionIDe;itRegionID++){
//								if(!isSharedRegion[*itRegionID]){
//									if(regions_owner[*itRegionID] == currentOwner){
//										regions2add.push_back(SubMeshes[iSubMesh].getLRegion(regions_LID[*itRegionID]));
//									}
//									else{
//									}
//								}
//								else{
//									if(regions_owners2LIDs[regions_LID[*itRegionID]].find(currentOwner) != regions_owners2LIDs[regions_LID[*itRegionID]].end()){
//										regions2add.push_back(SubMeshes[iSubMesh].getLRegion(regions_owners2LIDs[regions_LID[*itRegionID]][currentOwner]));
//									}
//									else{
//									}
//								}
//							}
//							regionDest->setRegions(regions2add);
//						}
//					}
//				}
//
//			} // for(;itRegion!=itRegione;itRegion++)
//		} // if(MeshDescriptor<M>::hasRegions && MeshDescriptor<M>::RKnowR)
//
//	} // if(MeshDescriptor<M>::hasRegions)
//
//	// deletion phase
//	if(MeshDescriptor<M>::hasNodes){
//		typename IGMesh::nodes_iterator it  = mesh_.nodes_begin();
//		typename IGMesh::nodes_iterator ite = mesh_.nodes_end();
//
//		for(;it!=ite;it++){
//			if(mesh_.isMarked((*it),entityToDeleteMark)){
//				mesh_.deleteNode(*it);
//			}
//		}
//	}
//	if(MeshDescriptor<M>::hasEdges){
//		typename IGMesh::edges_iterator it  = mesh_.edges_begin();
//		typename IGMesh::edges_iterator ite = mesh_.edges_end();
//
//		for(;it!=ite;it++){
//			if(mesh_.isMarked((*it),entityToDeleteMark)){
//				mesh_.deleteEdge(*it);
//			}
//		}
//	}
//	if(MeshDescriptor<M>::hasFaces){
//		typename IGMesh::faces_iterator it  = mesh_.faces_begin();
//		typename IGMesh::faces_iterator ite = mesh_.faces_end();
//
//		for(;it!=ite;it++){
//			if(mesh_.isMarked((*it),entityToDeleteMark)){
//				mesh_.deleteFace(*it);
//			}
//		}
//	}
//	if(MeshDescriptor<M>::hasRegions){
//		typename IGMesh::regions_iterator it  = mesh_.regions_begin();
//		typename IGMesh::regions_iterator ite = mesh_.regions_end();
//
//		for(;it!=ite;it++){
//			if(mesh_.isMarked((*it),entityToDeleteMark)){
//				mesh_.deleteRegion(*it);
//			}
//		}
//	}
//
//	mesh_.unmarkAll(entityToDeleteMark);
//	mesh_.freeMark(entityToDeleteMark);
//
//	return;
//}
///*----------------------------------------------------------------------------*/
//
//void SubMeshExtractor::extractNoGhostWithNativeSharedInfoTopLevelRepartition(const ECellType ATypePartAlong, const std::string& AName,
//												std::vector<IGMesh >& SubMeshes)
//{
//	// check that the model is coherent with this function.
//	if(MeshDescriptor<M>::NKnowN
//	|| MeshDescriptor<M>::NKnowE
//	|| MeshDescriptor<M>::NKnowF
//	|| MeshDescriptor<M>::NKnowR
//	|| MeshDescriptor<M>::EKnowE
//	|| MeshDescriptor<M>::EKnowF
//	|| MeshDescriptor<M>::EKnowR
//	|| MeshDescriptor<M>::FKnowF
//	|| MeshDescriptor<M>::FKnowR
//	|| MeshDescriptor<M>::RKnowR)
//	{
//		std::cout << "WARNING : Not a valid model, there are same-level/upwards connectivities!" << std::endl;
//		//throw GMDSException();
//	}
//
//
//	Variable<TInt> *var_dest = mesh_.template getVariable<TInt>(ATypePartAlong,AName);
//
//	// build parts order.
//	std::vector<TInt> partsOrdering;
//	std::map<TInt,TInt> partsOrder;
//
//	// the first one will be the local partition : that way we might reduce
//	// masters migrations since ownership is decided on a first-coming first-served principle.
//	partsOrdering.push_back(mesh_.getPartID());
//	partsOrder[mesh_.getPartID()] = 0;
//
//	{
//		typename IGMesh::cells_iterator it  = mesh_.cells_begin();
//		typename IGMesh::cells_iterator ite = mesh_.cells_end();
//
//		for(;it!=ite;it++){
//			if(partsOrder.find((*var_dest)[(*it)->getID()]) == partsOrder.end()){
//				partsOrdering.push_back((*var_dest)[(*it)->getID()]);
//				partsOrder[(*var_dest)[(*it)->getID()]] = partsOrdering.size();
//			}
//		}
//	}
//
//	// submeshes initialization. Of course, there needs to be at least one extracted
//	// submesh (current mesh does not count).
//	if(partsOrdering.size()>1){
//		SubMeshes.resize(partsOrdering.size()-1);
//	}
//
//	// we create two variables for each entity in the mesh to keep track of the owner and
//	// the local Id.
//	// the attribution of ownership is on a "first-come first-served" basis, except for the
//	// partitioned cell type.
//	Variable<TInt> *var_nodes_owner;
//	Variable<TCellID> *var_nodes_lId;
//	Variable<TInt> *var_faces_owner;
//	Variable<TCellID> *var_faces_lId;
//	Variable<TInt> *var_edges_owner;
//	Variable<TCellID> *var_edges_lId;
//	Variable<TInt> *var_regions_owner;
//	Variable<TCellID> *var_regions_lId;
//
//	// we create a mark that will characterize shared cells.
//	int cellSharedMark = mesh_.getNewMark();
//
//	if (MeshDescriptor<M>::hasNodes){
//		var_nodes_owner = mesh_.template newVariable<TInt>(GMDS_NODE,"var_nodes_owner");
//		var_nodes_lId = mesh_.template newVariable<TCellID>(GMDS_NODE,"var_nodes_lId");
//
//		typename IGMesh::nodes_iterator it  = mesh_.nodes_begin();
//		typename IGMesh::nodes_iterator ite = mesh_.nodes_end();
//		for(;it!=ite;it++){
//			(*var_nodes_owner)[(*it)->getID()] = NullTInt;
//			(*var_nodes_lId)[(*it)->getID()] = NullID;
//		}
//	}
//	if (MeshDescriptor<M>::hasEdges){
//		var_edges_owner = mesh_.template newVariable<TInt>(GMDS_EDGE,"var_edges_owner");
//		var_edges_lId = mesh_.template newVariable<TCellID>(GMDS_EDGE,"var_edges_lId");
//
//		typename IGMesh::edges_iterator it  = mesh_.edges_begin();
//		typename IGMesh::edges_iterator ite = mesh_.edges_end();
//		for(;it!=ite;it++){
//			(*var_edges_owner)[(*it)->getID()] = NullTInt;
//			(*var_edges_lId)[(*it)->getID()] = NullID;
//		}
//	}
//	if (MeshDescriptor<M>::hasFaces){
//		var_faces_owner = mesh_.template newVariable<TInt>(GMDS_FACE,"var_faces_owner");
//		var_faces_lId = mesh_.template newVariable<TCellID>(GMDS_FACE,"var_faces_lId");
//
//		typename IGMesh::faces_iterator it  = mesh_.faces_begin();
//		typename IGMesh::faces_iterator ite = mesh_.faces_end();
//		for(;it!=ite;it++){
//			(*var_faces_owner)[(*it)->getID()] = NullTInt;
//			(*var_faces_lId)[(*it)->getID()] = NullID;
//		}
//	}
//	if (MeshDescriptor<M>::hasRegions){
//		var_regions_owner = mesh_.template newVariable<TInt>(GMDS_REGION,"var_regions_owner");
//		var_regions_lId = mesh_.template newVariable<TCellID>(GMDS_REGION,"var_regions_lId");
//
//
//		typename IGMesh::regions_iterator it  = mesh_.regions_begin();
//		typename IGMesh::regions_iterator ite = mesh_.regions_end();
//		for(;it!=ite;it++){
//			(*var_regions_owner)[(*it)->getID()] = NullTInt;
//			(*var_regions_lId)[(*it)->getID()] = NullID;
//		}
//	}
//
//	// now we mark which entities will be kept in the submeshes, while
//	// updating the owner/local id variables.
//	// this is the first pass; it will determine which cells are shared.
//	for(TInt iPart=0;iPart<partsOrdering.size();iPart++){
//
//		TInt PartID = partsOrdering[iPart];
//
//		// local ids in the local mesh will be preserve; there is preemption for the local mesh.
//		bool setLID = (PartID == mesh_.getPartID());
//
//		int cellKeep = mesh_.getNewMark();
//
//		// first we mark the partitioned cells (depends on the type) as cellKeep.
//		// we also recursively mark the descendants of the partitioned cells as cellKeep.
//		{
//			typename IGMesh::cells_iterator it  = mesh_.cells_begin();
//			typename IGMesh::cells_iterator ite = mesh_.cells_end();
//
//			for(;it!=ite;it++){
//				if((*var_dest)[(*it)->getID()]==PartID){
//					recursiveMarkWithSharedRepartition(*it,cellKeep,PartID,
//							var_nodes_lId,var_edges_lId,var_faces_lId,var_regions_lId,
//							var_nodes_owner,var_edges_owner,var_faces_owner,var_regions_owner,
//							cellSharedMark,true,setLID);
//				}
//			}
//		}
//
//		mesh_.unmarkAll(cellKeep);
//		mesh_.freeMark(cellKeep);
//	}
//
//	// now we mark which entities will be kept in the submeshes, while
//	// updating the owner/local id variables.
//	// this is the second pass; it will extract the submeshes.
//	// We ignore the local mesh, as excess entities will removed from it later.
//	for(int iPart=1;iPart<partsOrdering.size();iPart++){
//
//		TInt PartID = partsOrdering[iPart];
//
//		int cellKeep = mesh_.getNewMark();
//
//		// first we mark the partitioned cells (depends on the type) as cellKeep.
//		// we also recursively mark the descendants of the partitioned cells as cellKeep.
//
//		typename IGMesh::cells_iterator it  = mesh_.cells_begin();
//		typename IGMesh::cells_iterator ite = mesh_.cells_end();
//
//		for(;it!=ite;it++){
//			if((*var_dest)[(*it)->getID()]==PartID){
//				recursiveMarkWithSharedRepartition(*it,cellKeep,PartID,
//						var_nodes_lId,var_edges_lId,var_faces_lId,var_regions_lId,
//						var_nodes_owner,var_edges_owner,var_faces_owner,var_regions_owner,
//						cellSharedMark,false,false);
//			}
//		}
//
//		getSubMeshRepartition(SubMeshes[iPart-1],cellKeep,partsOrdering[iPart],
//				var_nodes_lId,var_edges_lId,var_faces_lId,var_regions_lId,
//				var_nodes_owner,var_edges_owner,var_faces_owner,var_regions_owner,
//				cellSharedMark);
//
//		mesh_.unmarkAll(cellKeep);
//		mesh_.freeMark(cellKeep);
//
//	}
//
//	// we have to delete all the migrated entities in the local mesh.
//	{
//		TInt PartID = mesh_.getPartID();
//
//		int cellKeep = mesh_.getNewMark();
//
//		// first we mark the partitioned cells (depends on the type) as cellKeep.
//		// we also recursively mark the descendants of the partitioned cells as cellKeep.
//
//		typename IGMesh::cells_iterator it  = mesh_.cells_begin();
//		typename IGMesh::cells_iterator ite = mesh_.cells_end();
//
//		for(;it!=ite;it++){
//			if((*var_dest)[(*it)->getID()]==PartID){
//				recursiveMarkWithSharedRepartition(*it,cellKeep,PartID,
//						var_nodes_lId,var_edges_lId,var_faces_lId,var_regions_lId,
//						var_nodes_owner,var_edges_owner,var_faces_owner,var_regions_owner,
//						cellSharedMark,false,false);
//			}
//		}
//
//		getSubMeshRepartitionRemove(cellKeep,PartID,
//				var_nodes_lId,var_edges_lId,var_faces_lId,var_regions_lId,
//				var_nodes_owner,var_edges_owner,var_faces_owner,var_regions_owner,
//				cellSharedMark);
//
//		mesh_.unmarkAll(cellKeep);
//		mesh_.freeMark(cellKeep);
//	}
//
//
//	mesh_.unmarkAll(cellSharedMark);
//	mesh_.freeMark(cellSharedMark);
//
//	if (MeshDescriptor<M>::hasNodes){
//		mesh_.deleteVariable(GMDS_NODE,"var_nodes_owner");
//		mesh_.deleteVariable(GMDS_NODE,"var_nodes_lId");
//	}
//	if (MeshDescriptor<M>::hasEdges){
//		mesh_.deleteVariable(GMDS_EDGE,"var_edges_owner");
//		mesh_.deleteVariable(GMDS_EDGE,"var_edges_lId");
//	}
//	if (MeshDescriptor<M>::hasFaces){
//		mesh_.deleteVariable(GMDS_FACE,"var_faces_owner");
//		mesh_.deleteVariable(GMDS_FACE,"var_faces_lId");
//	}
//	if (MeshDescriptor<M>::hasRegions){
//		mesh_.deleteVariable(GMDS_REGION,"var_regions_owner");
//		mesh_.deleteVariable(GMDS_REGION,"var_regions_lId");
//	}
//
//	return;
//}
///*----------------------------------------------------------------------------*/
//
//void SubMeshExtractor::extractNoGhostWithNativeSharedInfoTopLevelRepartitionExp(const ECellType ATypePartAlong, const std::string& AName,
//												std::vector<IGMesh >& SubMeshes)
//{
//
//	// check that the model is coherent with this function.
//	if(MeshDescriptor<M>::NKnowN
//	|| MeshDescriptor<M>::NKnowE
//	|| MeshDescriptor<M>::NKnowF
//	|| MeshDescriptor<M>::NKnowR
//	|| MeshDescriptor<M>::EKnowE
//	|| MeshDescriptor<M>::EKnowF
//	|| MeshDescriptor<M>::EKnowR
//	|| MeshDescriptor<M>::FKnowF
//	|| MeshDescriptor<M>::FKnowR
//	|| MeshDescriptor<M>::RKnowR)
//	{
//		std::cout << "WARNING : Not a valid model, there are same-level/upwards connectivities!" << std::endl;
//		//throw GMDSException();
//	}
//
//	Variable<TInt> *var_dest = mesh_.template getVariable<TInt>(ATypePartAlong,AName);
//
//	// build parts order.
//	std::vector<TInt> partsOrdering;
//	std::map<TInt,TInt> partsOrder;
//
//	// the first one will be the local partition : that way we might reduce
//	// masters migrations since ownership is decided on a first-come first-served basis.
//	partsOrdering.push_back(mesh_.getPartID());
//	partsOrder[mesh_.getPartID()] = 0;
//
//	{
//		typename IGMesh::cells_iterator it  = mesh_.cells_begin();
//		typename IGMesh::cells_iterator ite = mesh_.cells_end();
//
//		for(;it!=ite;it++){
//			if(partsOrder.find((*var_dest)[(*it)->getID()]) == partsOrder.end()){
//				partsOrder[(*var_dest)[(*it)->getID()]] = partsOrdering.size();
//				partsOrdering.push_back((*var_dest)[(*it)->getID()]);
//			}
//		}
//	}
//
//	// submeshes initialization. Of course, there needs to be at least one extracted
//	// submesh (current mesh does not count).
//	if(partsOrdering.size()>1){
//		SubMeshes.resize(partsOrdering.size()-1);
//	}
//
//	// we create two variables for each entity in the mesh to keep track of the owner and
//	// the local Id.
//	// the attribution of ownership is on a "first-come first-served" basis, except for the
//	// partitioned cell type.
////	Variable<TInt> *var_nodes_owner;
////	Variable<TCellID> *var_nodes_lId;
////	Variable<TInt> *var_faces_owner;
////	Variable<TCellID> *var_faces_lId;
////	Variable<TInt> *var_edges_owner;
////	Variable<TCellID> *var_edges_lId;
////	Variable<TInt> *var_regions_owner;
////	Variable<TCellID> *var_regions_lId;
//	Variable<std::vector<TInt> > *var_nodes_owners;
//	Variable<std::vector<TCellID> > *var_nodes_lIds;
//	Variable<std::vector<TInt> > *var_faces_owners;
//	Variable<std::vector<TCellID> > *var_faces_lIds;
//	Variable<std::vector<TInt> > *var_edges_owners;
//	Variable<std::vector<TCellID> > *var_edges_lIds;
//	Variable<std::vector<TInt> > *var_regions_owners;
//	Variable<std::vector<TCellID> > *var_regions_lIds;
//
//	// we create a mark that will characterize shared cells.
//	int cellSharedMark = mesh_.getNewMark();
//
//	if (MeshDescriptor<M>::hasNodes){
////		var_nodes_owner = mesh_.template newVariable<TInt>(GMDS_NODE,"var_nodes_owner");
////		var_nodes_lId = mesh_.template newVariable<TCellID>(GMDS_NODE,"var_nodes_lId");
//		var_nodes_owners = mesh_.template newVariable<std::vector<TInt> >(GMDS_NODE,"var_nodes_owners");
//		var_nodes_lIds = mesh_.template newVariable<std::vector<TCellID> >(GMDS_NODE,"var_nodes_lIds");
//
////		typename IGMesh::nodes_iterator it  = mesh_.nodes_begin();
////		typename IGMesh::nodes_iterator ite = mesh_.nodes_end();
////		for(;it!=ite;it++){
////			(*var_nodes_owner)[(*it)->getID()] = NullTInt;
////			(*var_nodes_lId)[(*it)->getID()] = NullID;
////		}
//	}
//	if (MeshDescriptor<M>::hasEdges){
////		var_edges_owner = mesh_.template newVariable<TInt>(GMDS_EDGE,"var_edges_owner");
////		var_edges_lId = mesh_.template newVariable<TCellID>(GMDS_EDGE,"var_edges_lId");
//		var_edges_owners = mesh_.template newVariable<std::vector<TInt> >(GMDS_EDGE,"var_edges_owners");
//		var_edges_lIds = mesh_.template newVariable<std::vector<TCellID> >(GMDS_EDGE,"var_edges_lIds");
//
////		typename IGMesh::edges_iterator it  = mesh_.edges_begin();
////		typename IGMesh::edges_iterator ite = mesh_.edges_end();
////		for(;it!=ite;it++){
////			(*var_edges_owner)[(*it)->getID()] = NullTInt;
////			(*var_edges_lId)[(*it)->getID()] = NullID;
////		}
//	}
//	if (MeshDescriptor<M>::hasFaces){
////		var_faces_owner = mesh_.template newVariable<TInt>(GMDS_FACE,"var_faces_owner");
////		var_faces_lId = mesh_.template newVariable<TCellID>(GMDS_FACE,"var_faces_lId");
//		var_faces_owners = mesh_.template newVariable<std::vector<TInt> >(GMDS_FACE,"var_faces_owners");
//		var_faces_lIds = mesh_.template newVariable<std::vector<TCellID> >(GMDS_FACE,"var_faces_lIds");
//
////		typename IGMesh::faces_iterator it  = mesh_.faces_begin();
////		typename IGMesh::faces_iterator ite = mesh_.faces_end();
////		for(;it!=ite;it++){
////			(*var_faces_owner)[(*it)->getID()] = NullTInt;
////			(*var_faces_lId)[(*it)->getID()] = NullID;
////		}
//	}
//	if (MeshDescriptor<M>::hasRegions){
////		var_regions_owner = mesh_.template newVariable<TInt>(GMDS_REGION,"var_regions_owner");
////		var_regions_lId = mesh_.template newVariable<TCellID>(GMDS_REGION,"var_regions_lId");
//		var_regions_owners = mesh_.template newVariable<std::vector<TInt> >(GMDS_REGION,"var_regions_owners");
//		var_regions_lIds = mesh_.template newVariable<std::vector<TCellID> >(GMDS_REGION,"var_regions_lIds");
//
////		typename IGMesh::regions_iterator it  = mesh_.regions_begin();
////		typename IGMesh::regions_iterator ite = mesh_.regions_end();
////		for(;it!=ite;it++){
////			(*var_regions_owner)[(*it)->getID()] = NullTInt;
////			(*var_regions_lId)[(*it)->getID()] = NullID;
////		}
//	}
//
//	// now we mark which entities will be kept in the submeshes, while
//	// updating the owner/local id variables.
//	// this is the first pass; it will determine which cells are shared.
//	{
//		typename IGMesh::cells_iterator it  = mesh_.cells_begin();
//		typename IGMesh::cells_iterator ite = mesh_.cells_end();
//
//		for(;it!=ite;it++){
//			recursiveMarkWithSharedRepartitionExp(*it,(*var_dest)[(*it)->getID()],
//					var_nodes_owners,var_edges_owners,var_faces_owners,var_regions_owners,
//					cellSharedMark);
//		}
//	}
//
//	getSubMeshRepartitionExp(SubMeshes,partsOrdering,partsOrder,
//			var_nodes_lIds,var_edges_lIds,var_faces_lIds,var_regions_lIds,
//			var_nodes_owners,var_edges_owners,var_faces_owners,var_regions_owners,
//			cellSharedMark);
//
////
////	// now we mark which entities will be kept in the submeshes, while
////	// updating the owner/local id variables.
////	// this is the second pass; it will extract the submeshes.
////	// We ignore the local mesh, as excess entities will removed from it later.
////	for(int iPart=1;iPart<partsOrdering.size();iPart++){
////
////		TInt PartID = partsOrdering[iPart];
////
////		int cellKeep = mesh_.getNewMark();
////
////		// first we mark the partitioned cells (depends on the type) as cellKeep.
////		// we also recursively mark the descendants of the partitioned cells as cellKeep.
////
////		typename IGMesh::cells_iterator it  = mesh_.cells_begin();
////		typename IGMesh::cells_iterator ite = mesh_.cells_end();
////
////		for(;it!=ite;it++){
////			if((*var_dest)[(*it)->getID()]==PartID){
////				recursiveMarkWithSharedRepartition(*it,cellKeep,PartID,
////						var_nodes_lId,var_edges_lId,var_faces_lId,var_regions_lId,
////						var_nodes_owner,var_edges_owner,var_faces_owner,var_regions_owner,
////						cellSharedMark,false,false);
////			}
////		}
////
////		getSubMeshRepartition(SubMeshes[iPart-1],cellKeep,partsOrdering[iPart],
////				var_nodes_lId,var_edges_lId,var_faces_lId,var_regions_lId,
////				var_nodes_owner,var_edges_owner,var_faces_owner,var_regions_owner,
////				cellSharedMark);
////
////		mesh_.unmarkAll(cellKeep);
////		mesh_.freeMark(cellKeep);
////
////	}
////
////	// we have to delete all the migrated entities in the local mesh.
////	{
////		TInt PartID = mesh_.getPartitionID();
////
////		int cellKeep = mesh_.getNewMark();
////
////		// first we mark the partitioned cells (depends on the type) as cellKeep.
////		// we also recursively mark the descendants of the partitioned cells as cellKeep.
////
////		typename IGMesh::cells_iterator it  = mesh_.cells_begin();
////		typename IGMesh::cells_iterator ite = mesh_.cells_end();
////
////		for(;it!=ite;it++){
////			if((*var_dest)[(*it)->getID()]==PartID){
////				recursiveMarkWithSharedRepartition(*it,cellKeep,PartID,
////						var_nodes_lId,var_edges_lId,var_faces_lId,var_regions_lId,
////						var_nodes_owner,var_edges_owner,var_faces_owner,var_regions_owner,
////						cellSharedMark,false,false);
////			}
////		}
////
////		getSubMeshRepartitionRemove(cellKeep,PartID,
////				var_nodes_lId,var_edges_lId,var_faces_lId,var_regions_lId,
////				var_nodes_owner,var_edges_owner,var_faces_owner,var_regions_owner,
////				cellSharedMark);
////
////		mesh_.unmarkAll(cellKeep);
////		mesh_.freeMark(cellKeep);
////	}
////
////
//	mesh_.unmarkAll(cellSharedMark);
//	mesh_.freeMark(cellSharedMark);
//
//	if (MeshDescriptor<M>::hasNodes){
////		mesh_.deleteVariable(GMDS_NODE,"var_nodes_owner");
////		mesh_.deleteVariable(GMDS_NODE,"var_nodes_lId");
//		mesh_.deleteVariable(GMDS_NODE,"var_nodes_owners");
//		mesh_.deleteVariable(GMDS_NODE,"var_nodes_lIds");
//	}
//	if (MeshDescriptor<M>::hasEdges){
////		mesh_.deleteVariable(GMDS_EDGE,"var_edges_owner");
////		mesh_.deleteVariable(GMDS_EDGE,"var_edges_lId");
//		mesh_.deleteVariable(GMDS_EDGE,"var_edges_owners");
//		mesh_.deleteVariable(GMDS_EDGE,"var_edges_lIds");
//	}
//	if (MeshDescriptor<M>::hasFaces){
////		mesh_.deleteVariable(GMDS_FACE,"var_faces_owner");
////		mesh_.deleteVariable(GMDS_FACE,"var_faces_lId");
//		mesh_.deleteVariable(GMDS_FACE,"var_faces_owners");
//		mesh_.deleteVariable(GMDS_FACE,"var_faces_lIds");
//	}
//	if (MeshDescriptor<M>::hasRegions){
////		mesh_.deleteVariable(GMDS_REGION,"var_regions_owner");
////		mesh_.deleteVariable(GMDS_REGION,"var_regions_lId");
//		mesh_.deleteVariable(GMDS_REGION,"var_regions_owners");
//		mesh_.deleteVariable(GMDS_REGION,"var_regions_lIds");
//	}
//
//	return;
//}
///*----------------------------------------------------------------------------*/
//
//void SubMeshExtractor::extractNoGhostWithNativeSharedInfoTopLevelRepartitionExp2(const ECellType ATypePartAlong, const std::string& AName,
//												std::vector<IGMesh >& SubMeshes)
//{
//
//	// check that the model is coherent with this function.
//	// This function does not generate a ghost layer, so only
//	// downwards connections should be present to keep a coherent model.
//	if(MeshDescriptor<M>::NKnowN
//	|| MeshDescriptor<M>::NKnowE
//	|| MeshDescriptor<M>::NKnowF
//	|| MeshDescriptor<M>::NKnowR
//	|| MeshDescriptor<M>::EKnowE
//	|| MeshDescriptor<M>::EKnowF
//	|| MeshDescriptor<M>::EKnowR
//	|| MeshDescriptor<M>::FKnowF
//	|| MeshDescriptor<M>::FKnowR
//	|| MeshDescriptor<M>::RKnowR)
//	{
//		std::cout << "WARNING : Not a valid model, there are same-level/upwards connectivities!" << std::endl;
//		//throw GMDSException();
//	}
//
//	Variable<TInt> *var_dest = mesh_.template getVariable<TInt>(ATypePartAlong,AName);
//
//	// build parts order.
//	std::vector<TInt> partsOrdering;
//	std::map<TInt,TInt> partsOrder;
//
//	// the first one will be the local partition : that way we might reduce
//	// masters migrations since ownership is decided on a first-come first-served basis.
//	partsOrdering.push_back(mesh_.getPartID());
//	partsOrder[mesh_.getPartID()] = 0;
//
//	{
//		typename IGMesh::cells_iterator it  = mesh_.cells_begin();
//		typename IGMesh::cells_iterator ite = mesh_.cells_end();
//
//		for(;it!=ite;it++){
//			if(partsOrder.find((*var_dest)[(*it)->getID()]) == partsOrder.end()){
//				partsOrder[(*var_dest)[(*it)->getID()]] = partsOrdering.size();
//				partsOrdering.push_back((*var_dest)[(*it)->getID()]);
//			}
//		}
//	}
//
//	// submeshes initialization. Of course, there needs to be at least one extracted
//	// submesh (current mesh does not count).
//	if(partsOrdering.size()>1){
//		SubMeshes.resize(partsOrdering.size()-1);
//	}
//
//	// we create two variables for each entity in the mesh to keep track of the owner and
//	// the local Id.
//	// The attribution of ownership is on a "first-come first-served" basis, except for the
//	// partitioned cell type.
//	// The owners2LIDs containers will manage entities present on multiple parts.
//	std::cout << "begin initialize allocate" << std::endl;
//	bool isSharedNode[mesh_.getMaxLocalID(0)+1];
//	bool isSharedEdge[mesh_.getMaxLocalID(1)+1];
//	bool isSharedFace[mesh_.getMaxLocalID(2)+1];
//	bool isSharedRegion[mesh_.getMaxLocalID(3)+1];
//
//	TInt nodes_owner[mesh_.getMaxLocalID(0)+1];
//	TInt edges_owner[mesh_.getMaxLocalID(1)+1];
//	TInt faces_owner[mesh_.getMaxLocalID(2)+1];
//	TInt regions_owner[mesh_.getMaxLocalID(3)+1];
//
//	TInt nodes_LID[mesh_.getMaxLocalID(0)+1];
//	TInt edges_LID[mesh_.getMaxLocalID(1)+1];
//	TInt faces_LID[mesh_.getMaxLocalID(2)+1];
//	TInt regions_LID[mesh_.getMaxLocalID(3)+1];
//
//	std::vector<std::map<TInt,TCellID> > nodes_owners2LIDs;
//	std::vector<std::map<TInt,TCellID> > edges_owners2LIDs;
//	std::vector<std::map<TInt,TCellID> > faces_owners2LIDs;
//	std::vector<std::map<TInt,TCellID> > regions_owners2LIDs;
//
//	std::cout << "begin initialize values" << std::endl;
//	if (MeshDescriptor<M>::hasNodes){
//		id maxID = mesh_.getMaxLocalID(0)+1;
//
//		for(id i=0; i<maxID; i++){
//			isSharedNode[i] = false;
//			nodes_owner[i] = NullTInt;
//			nodes_LID[i] = NullID;
//		}
//	}
//
//	if (MeshDescriptor<M>::hasEdges){
//		id maxID = mesh_.getMaxLocalID(1)+1;
//
//		for(id i=0; i<maxID; i++){
//			isSharedEdge[i] = false;
//			edges_owner[i] = NullTInt;
//			edges_LID[i] = NullID;
//		}
//	}
//
//	if (MeshDescriptor<M>::hasFaces){
//		id maxID = mesh_.getMaxLocalID(2)+1;
//
//		for(id i=0; i<maxID; i++){
//			isSharedFace[i] = false;
//			faces_owner[i] = NullTInt;
//			faces_LID[i] = NullID;
//		}
//	}
//
//	if (MeshDescriptor<M>::hasRegions){
//		id maxID = mesh_.getMaxLocalID(3)+1;
//
//		for(id i=0; i<maxID; i++){
//			isSharedRegion[i] = false;
//			regions_owner[i] = NullTInt;
//			regions_LID[i] = NullID;
//		}
//	}
//
//	std::cout << "begin recursiveMarkWithSharedRepartitionExp2" << std::endl;
//	// now we mark which entities will be kept in the submeshes, while
//	// updating the owner variables.
//	{
//		typename IGMesh::cells_iterator it  = mesh_.cells_begin();
//		typename IGMesh::cells_iterator ite = mesh_.cells_end();
//
//		for(;it!=ite;it++){
//			recursiveMarkWithSharedRepartitionExp2(*it,(*var_dest)[(*it)->getID()],
//					isSharedNode,isSharedEdge,isSharedFace,isSharedRegion,
//					nodes_LID,edges_LID,faces_LID,regions_LID,
//					nodes_owner,edges_owner,faces_owner,regions_owner,
//					nodes_owners2LIDs,edges_owners2LIDs,faces_owners2LIDs,regions_owners2LIDs
//					);
//		}
//	}
//	std::cout << "end recursiveMarkWithSharedRepartitionExp2" << std::endl;
//	std::cout << "begin getSubMeshRepartitionExp2" << std::endl;
//
//	// submeshes extraction and deletion of entities no longer in mesh_
//	getSubMeshRepartitionExp2(SubMeshes,partsOrdering,partsOrder,
//			isSharedNode,isSharedEdge,isSharedFace,isSharedRegion,
//			nodes_LID,edges_LID,faces_LID,regions_LID,
//			nodes_owner,edges_owner,faces_owner,regions_owner,
//			nodes_owners2LIDs,edges_owners2LIDs,faces_owners2LIDs,regions_owners2LIDs);
//	std::cout << "end getSubMeshRepartitionExp2" << std::endl;
//
//	return;
//}
///*----------------------------------------------------------------------------*/
//
//void SubMeshExtractor::extractWithGhostWithNativeSharedInfoTopLevelRepartitionExp2(const ECellType ATypePartAlong, const std::string& AName,
//												std::vector<IGMesh >& SubMeshes)
//{
//	std::cout << "DO NOT CALL THAT METHOD FOR THE TIME BEING : extractWithGhostWithNativeSharedInfoTopLevelRepartitionExp2" << std::endl;
//	//throw GMDSException();
//
//	Variable<TInt> *var_dest = mesh_.template getVariable<TInt>(ATypePartAlong,AName);
//
//	// build parts order.
//	std::vector<TInt> partsOrdering;
//	std::map<TInt,TInt> partsOrder;
//
//	// the first one will be the local partition : that way we might reduce
//	// masters migrations since ownership is decided on a first-come first-served basis.
//	partsOrdering.push_back(mesh_.getPartID());
//	partsOrder[mesh_.getPartID()] = 0;
//
//	{
//		typename IGMesh::cells_iterator it  = mesh_.cells_begin();
//		typename IGMesh::cells_iterator ite = mesh_.cells_end();
//
//		for(;it!=ite;it++){
//			if(partsOrder.find((*var_dest)[(*it)->getID()]) == partsOrder.end()){
//				partsOrder[(*var_dest)[(*it)->getID()]] = partsOrdering.size();
//				partsOrdering.push_back((*var_dest)[(*it)->getID()]);
//			}
//		}
//	}
//
//	// submeshes initialization. Of course, there needs to be at least one extracted
//	// submesh (current mesh does not count).
//	if(partsOrdering.size()>1){
//		SubMeshes.resize(partsOrdering.size()-1);
//	}
//
//	// we create two variables for each entity in the mesh to keep track of the owner and
//	// the local Id.
//	// The attribution of ownership is on a "first-come first-served" basis, except for the
//	// partitioned cell type.
//	// The owners2LIDs containers will manage entities present on multiple parts.
//	std::cout << "begin initialize allocate" << std::endl;
//	bool isSharedNode[mesh_.getMaxLocalID(0)+1];
//	bool isSharedEdge[mesh_.getMaxLocalID(1)+1];
//	bool isSharedFace[mesh_.getMaxLocalID(2)+1];
//	bool isSharedRegion[mesh_.getMaxLocalID(3)+1];
//
//	TInt nodes_owner[mesh_.getMaxLocalID(0)+1];
//	TInt edges_owner[mesh_.getMaxLocalID(1)+1];
//	TInt faces_owner[mesh_.getMaxLocalID(2)+1];
//	TInt regions_owner[mesh_.getMaxLocalID(3)+1];
//
//	TInt nodes_LID[mesh_.getMaxLocalID(0)+1];
//	TInt edges_LID[mesh_.getMaxLocalID(1)+1];
//	TInt faces_LID[mesh_.getMaxLocalID(2)+1];
//	TInt regions_LID[mesh_.getMaxLocalID(3)+1];
//
//	std::vector<std::map<TInt,TCellID> > nodes_owners2LIDs;
//	std::vector<std::map<TInt,TCellID> > edges_owners2LIDs;
//	std::vector<std::map<TInt,TCellID> > faces_owners2LIDs;
//	std::vector<std::map<TInt,TCellID> > regions_owners2LIDs;
//
//	std::cout << "begin initialize values" << std::endl;
//	if (MeshDescriptor<M>::hasNodes){
//		id maxID = mesh_.getMaxLocalID(0)+1;
//
//		for(id i=0; i<maxID; i++){
//			isSharedNode[i] = false;
//			nodes_owner[i] = NullTInt;
//			nodes_LID[i] = NullID;
//		}
//	}
//
//	if (MeshDescriptor<M>::hasEdges){
//		id maxID = mesh_.getMaxLocalID(1)+1;
//
//		for(id i=0; i<maxID; i++){
//			isSharedEdge[i] = false;
//			edges_owner[i] = NullTInt;
//			edges_LID[i] = NullID;
//		}
//	}
//
//	if (MeshDescriptor<M>::hasFaces){
//		id maxID = mesh_.getMaxLocalID(2)+1;
//
//		for(id i=0; i<maxID; i++){
//			isSharedFace[i] = false;
//			faces_owner[i] = NullTInt;
//			faces_LID[i] = NullID;
//		}
//	}
//
//	if (MeshDescriptor<M>::hasRegions){
//		id maxID = mesh_.getMaxLocalID(3)+1;
//
//		for(id i=0; i<maxID; i++){
//			isSharedRegion[i] = false;
//			regions_owner[i] = NullTInt;
//			regions_LID[i] = NullID;
//		}
//	}
//
//	std::cout << "begin recursiveMarkWithSharedRepartitionExp2" << std::endl;
//	// now we mark which entities will be kept in the submeshes, while
//	// updating the owner variables.
//	{
//		typename IGMesh::cells_iterator it  = mesh_.cells_begin();
//		typename IGMesh::cells_iterator ite = mesh_.cells_end();
//
//		for(;it!=ite;it++){
//			recursiveMarkWithSharedRepartitionExp2(*it,(*var_dest)[(*it)->getID()],
//					isSharedNode,isSharedEdge,isSharedFace,isSharedRegion,
//					nodes_LID,edges_LID,faces_LID,regions_LID,
//					nodes_owner,edges_owner,faces_owner,regions_owner,
//					nodes_owners2LIDs,edges_owners2LIDs,faces_owners2LIDs,regions_owners2LIDs
//					);
//		}
//	}
//	std::cout << "end recursiveMarkWithSharedRepartitionExp2" << std::endl;
//
//	// extend the marking to the ghost layer.
//	{
//		Variable<std::vector<TInt> >* cells2Mark;
//
//		if(MeshDescriptor<M>::hasRegions){
//			cells2Mark = mesh_.template newVariable<std::vector<TInt> >(GMDS_REGION,"cells2Mark");
//		}
//		else{
//			if(MeshDescriptor<M>::hasFaces){
//				cells2Mark = mesh_.template newVariable<std::vector<TInt> >(GMDS_FACE,"cells2Mark");
//			}
//			else{
//				if(MeshDescriptor<M>::hasEdges){
//					cells2Mark = mesh_.template newVariable<std::vector<TInt> >(GMDS_FACE,"cells2Mark");
//				}
//				else{
//					std::cout << "no point in trying to get ghosts when there are only nodes." << std::endl;
//					throw GMDSException();
//				}
//			}
//		}
//
//		{
//			typename IGMesh::cells_iterator it  = mesh_.cells_begin();
//			typename IGMesh::cells_iterator ite = mesh_.cells_end();
//
//			for(;it!=ite;it++){
//				recursiveGhostMarkingExp2(*it,(*cells2Mark)[(*it)->getID()],
//						isSharedNode,isSharedEdge,isSharedFace,isSharedRegion,
//						nodes_LID,edges_LID,faces_LID,regions_LID,
//						nodes_owner,edges_owner,faces_owner,regions_owner,
//						nodes_owners2LIDs,edges_owners2LIDs,faces_owners2LIDs,regions_owners2LIDs
//						);
//			}
//		}
//
//		{
//			typename IGMesh::cells_iterator it  = mesh_.cells_begin();
//			typename IGMesh::cells_iterator ite = mesh_.cells_end();
//
//			for(;it!=ite;it++){
//				for(id iPart2add=0; iPart2add<(*cells2Mark)[(*it)->getID()].size(); iPart2add++){
//					TInt iPart = (*cells2Mark)[(*it)->getID()][iPart2add];
//
//					recursiveMarkWithSharedRepartitionExp2(*it,iPart,
//							isSharedNode,isSharedEdge,isSharedFace,isSharedRegion,
//							nodes_LID,edges_LID,faces_LID,regions_LID,
//							nodes_owner,edges_owner,faces_owner,regions_owner,
//							nodes_owners2LIDs,edges_owners2LIDs,faces_owners2LIDs,regions_owners2LIDs
//							);
//				}
//			}
//		}
//
//		if(MeshDescriptor<M>::hasRegions){
//			mesh_.deleteVariable(GMDS_REGION,"cells2Mark");
//		}
//		else{
//			if(MeshDescriptor<M>::hasFaces){
//				mesh_.deleteVariable(GMDS_FACE,"cells2Mark");
//			}
//			else{
//				if(MeshDescriptor<M>::hasEdges){
//					mesh_.deleteVariable(GMDS_EDGE,"cells2Mark");
//				}
//				else{
//					std::cout << "no point in trying to get ghosts when there are only nodes." << std::endl;
//					throw GMDSException();
//				}
//			}
//		}
//
//
//	}
//
//	std::cout << "begin getSubMeshRepartitionExp2" << std::endl;
//
//	// submeshes extraction and deletion of entities no longer in mesh_
//	getSubMeshRepartitionExp2(SubMeshes,partsOrdering,partsOrder,
//			isSharedNode,isSharedEdge,isSharedFace,isSharedRegion,
//			nodes_LID,edges_LID,faces_LID,regions_LID,
//			nodes_owner,edges_owner,faces_owner,regions_owner,
//			nodes_owners2LIDs,edges_owners2LIDs,faces_owners2LIDs,regions_owners2LIDs);
//	std::cout << "end getSubMeshRepartitionExp2" << std::endl;
//
//	return;
//}
///*----------------------------------------------------------------------------*/
//
//void SubMeshExtractor::extractWithGhostWithNativeSharedInfoTopLevelRepartition(const ECellType ATypePartAlong, const std::string& AName,
//												std::vector<IGMesh >& SubMeshes)
//{
//	// TODO check that the model is coherent with this function.
//	// all models are OK?
//
//	Variable<TInt> *var_dest = mesh_.template getVariable<TInt>(ATypePartAlong,AName);
//
//	// build parts order.
//	std::vector<TInt> partsOrdering;
//	std::map<TInt,TInt> partsOrder;
//
//	// the first one will be the local partition : that way we might reduce
//	// masters migrations since ownership is decided on a first-come first-served basis.
//	partsOrdering.push_back(mesh_.getPartID());
//	partsOrder[mesh_.getPartID()] = 0;
//
//	{
//		typename IGMesh::cells_iterator it  = mesh_.cells_begin();
//		typename IGMesh::cells_iterator ite = mesh_.cells_end();
//
//		for(;it!=ite;it++){
//			if(partsOrder.find((*var_dest)[(*it)->getID()]) == partsOrder.end()){
//				partsOrdering.push_back((*var_dest)[(*it)->getID()]);
//				partsOrder[(*var_dest)[(*it)->getID()]] = partsOrdering.size();
//			}
//		}
//	}
//
//	// submeshes initialization. Of course, there needs to be at least one extracted
//	// submesh (current mesh does not count).
//	if(partsOrdering.size()>1){
//		SubMeshes.resize(partsOrdering.size()-1);
//	}
//
//	// we create two variables for each entity in the mesh to keep track of the owner and
//	// the local Id.
//	// the attribution of ownership is on a "first-come first-served" basis, except for the
//	// partitioned cell type.
//	Variable<TInt> *var_nodes_owner;
//	Variable<TCellID> *var_nodes_lId;
//	Variable<TInt> *var_faces_owner;
//	Variable<TCellID> *var_faces_lId;
//	Variable<TInt> *var_edges_owner;
//	Variable<TCellID> *var_edges_lId;
//	Variable<TInt> *var_regions_owner;
//	Variable<TCellID> *var_regions_lId;
//
//	// we create a mark that will characterize shared cells.
//	int cellSharedMark = mesh_.getNewMark();
//
//	if(MeshDescriptor<M>::hasNodes){
//		var_nodes_owner = mesh_.template newVariable<TInt>(GMDS_NODE,"var_nodes_owner");
//		var_nodes_lId = mesh_.template newVariable<TCellID>(GMDS_NODE,"var_nodes_lId");
//
//		typename IGMesh::nodes_iterator it  = mesh_.nodes_begin();
//		typename IGMesh::nodes_iterator ite = mesh_.nodes_end();
//		for(;it!=ite;it++){
//			(*var_nodes_owner)[(*it)->getID()] = NullTInt;
//			(*var_nodes_lId)[(*it)->getID()] = NullID;
//		}
//	}
//	if(MeshDescriptor<M>::hasEdges){
//		var_edges_owner = mesh_.template newVariable<TInt>(GMDS_EDGE,"var_edges_owner");
//		var_edges_lId = mesh_.template newVariable<TCellID>(GMDS_EDGE,"var_edges_lId");
//
//		typename IGMesh::edges_iterator it  = mesh_.edges_begin();
//		typename IGMesh::edges_iterator ite = mesh_.edges_end();
//		for(;it!=ite;it++){
//			(*var_edges_owner)[(*it)->getID()] = NullTInt;
//			(*var_edges_lId)[(*it)->getID()] = NullID;
//		}
//	}
//	if (MeshDescriptor<M>::hasFaces){
//		var_faces_owner = mesh_.template newVariable<TInt>(GMDS_FACE,"var_faces_owner");
//		var_faces_lId = mesh_.template newVariable<TCellID>(GMDS_FACE,"var_faces_lId");
//
//		typename IGMesh::faces_iterator it  = mesh_.faces_begin();
//		typename IGMesh::faces_iterator ite = mesh_.faces_end();
//		for(;it!=ite;it++){
//			(*var_faces_owner)[(*it)->getID()] = NullTInt;
//			(*var_faces_lId)[(*it)->getID()] = NullID;
//		}
//	}
//	if (MeshDescriptor<M>::hasRegions){
//		var_regions_owner = mesh_.template newVariable<TInt>(GMDS_REGION,"var_regions_owner");
//		var_regions_lId = mesh_.template newVariable<TCellID>(GMDS_REGION,"var_regions_lId");
//
//		typename IGMesh::regions_iterator it  = mesh_.regions_begin();
//		typename IGMesh::regions_iterator ite = mesh_.regions_end();
//		for(;it!=ite;it++){
//			(*var_regions_owner)[(*it)->getID()] = NullTInt;
//			(*var_regions_lId)[(*it)->getID()] = NullID;
//		}
//	}
//
//	// now we mark which entities will be kept in the submeshes, while
//	// updating the owner/local id variables.
//	// this is the first pass; it will determine which cells are shared.
//	for(TInt iPart=0;iPart<partsOrdering.size();iPart++){
//
//		TInt PartID = partsOrdering[iPart];
//
//		// local ids in the local mesh will be preserve; there is preemption for the local mesh.
//		bool setLID = (PartID == mesh_.getPartID());
//
//		int cellKeep = mesh_.getNewMark();
//
//		// first we mark the partitioned cells (depends on the type) as cellKeep.
//		// we also recursively mark the descendants of the partitioned cells as cellKeep.
//		{
//			typename IGMesh::cells_iterator it  = mesh_.cells_begin();
//			typename IGMesh::cells_iterator ite = mesh_.cells_end();
//
//			for(;it!=ite;it++){
//				if((*var_dest)[(*it)->getID()]==PartID){
//					recursiveMarkWithSharedRepartition(*it,cellKeep,PartID,
//							var_nodes_lId,var_edges_lId,var_faces_lId,var_regions_lId,
//							var_nodes_owner,var_edges_owner,var_faces_owner,var_regions_owner,
//							cellSharedMark,true,setLID);
//				}
//			}
//		}
//
//		mesh_.unmarkAll(cellKeep);
//		mesh_.freeMark(cellKeep);
//	}
//
//	for(TInt iPart=0;iPart<partsOrdering.size();iPart++){
//
//		TInt PartID = partsOrdering[iPart];
//
//		// local ids in the local mesh will be preserve; there is preemption for the local mesh.
//		//bool setLID = (PartID == mesh_.getPartitionID());
//
//		int cellKeep = mesh_.getNewMark();
//		int cellKeepImplied = mesh_.getNewMark();
//
//		// first we mark the partitioned cells (depends on the type) as cellKeep.
//		// we also recursively mark the descendants of the partitioned cells as cellKeep.
//		{
//			typename IGMesh::cells_iterator it  = mesh_.cells_begin();
//			typename IGMesh::cells_iterator ite = mesh_.cells_end();
//
//			for(;it!=ite;it++){
//				if((*var_dest)[(*it)->getID()]==PartID){
//					//std::cout << "Keep " << PartID << " " << (*it)->getID() << std::endl;
//					recursiveMarkWithSharedRepartition(*it,cellKeep,PartID,
//							var_nodes_lId,var_edges_lId,var_faces_lId,var_regions_lId,
//							var_nodes_owner,var_edges_owner,var_faces_owner,var_regions_owner,
//							cellSharedMark,false,false);
//							//cellSharedMark,false,setLID);
//				}
//			}
//		}
//
//		// mark the entities in the ghost layer if there is need of one.
//		if(MeshDescriptor<M>::NKnowN
//		|| MeshDescriptor<M>::NKnowE
//		|| MeshDescriptor<M>::NKnowF
//		|| MeshDescriptor<M>::NKnowR
//		|| MeshDescriptor<M>::EKnowE
//		|| MeshDescriptor<M>::EKnowF
//		|| MeshDescriptor<M>::EKnowR
//		|| MeshDescriptor<M>::FKnowF
//		|| MeshDescriptor<M>::FKnowR
//		|| MeshDescriptor<M>::RKnowR)
//		{
//			typename IGMesh::cells_iterator it  = mesh_.cells_begin();
//			typename IGMesh::cells_iterator ite = mesh_.cells_end();
//
//			for(;it!=ite;it++){
//				if(recursiveIsMarked(*it,cellKeep)){
//					//std::cout << "Implied " << PartID << " " << (*it)->getID() << std::endl;
//					recursiveMarkWithSharedRepartition(*it,cellKeepImplied,PartID,
//							var_nodes_lId,var_edges_lId,var_faces_lId,var_regions_lId,
//							var_nodes_owner,var_edges_owner,var_faces_owner,var_regions_owner,
//							cellSharedMark,false,false);
//				}
//			}
//		}
//
//		mesh_.unmarkAll(cellKeep);
//		mesh_.freeMark(cellKeep);
//		mesh_.unmarkAll(cellKeepImplied);
//		mesh_.freeMark(cellKeepImplied);
//	}
//
//	// now we mark which entities will be kept in the submeshes, while
//	// updating the owner/local id variables.
//	// this is the second pass; it will extract the submeshes.
//	// We ignore the local mesh, as excess entities will be removed from it later.
//	for(int iPart=1;iPart<partsOrdering.size();iPart++){
//
//		TInt PartID = partsOrdering[iPart];
//
//		int cellKeep = mesh_.getNewMark();
//		int cellKeepImplied = mesh_.getNewMark();
//
//		// first we mark the partitioned cells (depends on the type) as cellKeep.
//		// we also recursively mark the descendants of the partitioned cells as cellKeep.
//
//		{
//		typename IGMesh::cells_iterator it  = mesh_.cells_begin();
//		typename IGMesh::cells_iterator ite = mesh_.cells_end();
//
//		for(;it!=ite;it++){
//			if((*var_dest)[(*it)->getID()]==PartID){
//				recursiveMarkWithSharedRepartition(*it,cellKeep,PartID,
//						var_nodes_lId,var_edges_lId,var_faces_lId,var_regions_lId,
//						var_nodes_owner,var_edges_owner,var_faces_owner,var_regions_owner,
//						cellSharedMark,false,false);
//			}
//		}
//		}
//
//		// mark the entities in the ghost layer if there is need of one.
//		if(MeshDescriptor<M>::NKnowN
//		|| MeshDescriptor<M>::NKnowE
//		|| MeshDescriptor<M>::NKnowF
//		|| MeshDescriptor<M>::NKnowR
//		|| MeshDescriptor<M>::EKnowE
//		|| MeshDescriptor<M>::EKnowF
//		|| MeshDescriptor<M>::EKnowR
//		|| MeshDescriptor<M>::FKnowF
//		|| MeshDescriptor<M>::FKnowR
//		|| MeshDescriptor<M>::RKnowR)
//		{
//			typename IGMesh::cells_iterator it  = mesh_.cells_begin();
//			typename IGMesh::cells_iterator ite = mesh_.cells_end();
//
//			for(;it!=ite;it++){
//				if(recursiveIsMarked(*it,cellKeep)){
//					mesh_.mark(*it,cellKeepImplied);
////					recursiveMarkWithSharedRepartition(*it,cellKeepImplied,PartID,
////							var_nodes_lId,var_edges_lId,var_faces_lId,var_regions_lId,
////							var_nodes_owner,var_edges_owner,var_faces_owner,var_regions_owner,
////							cellSharedMark,false,false);
//				}
//			}
//		}
//
//		{
//			typename IGMesh::cells_iterator it  = mesh_.cells_begin();
//			typename IGMesh::cells_iterator ite = mesh_.cells_end();
//
//			for(;it!=ite;it++){
//				if(mesh_.isMarked(*it,cellKeepImplied)){
//					recursiveMarkWithSharedRepartition(*it,cellKeep,PartID,
//							var_nodes_lId,var_edges_lId,var_faces_lId,var_regions_lId,
//							var_nodes_owner,var_edges_owner,var_faces_owner,var_regions_owner,
//							cellSharedMark,false,false);
//				}
//			}
//		}
//
//		getSubMeshRepartition(SubMeshes[iPart-1],cellKeep,partsOrdering[iPart],
//				var_nodes_lId,var_edges_lId,var_faces_lId,var_regions_lId,
//				var_nodes_owner,var_edges_owner,var_faces_owner,var_regions_owner,
//				cellSharedMark);
//
//		mesh_.unmarkAll(cellKeep);
//		mesh_.freeMark(cellKeep);
//		mesh_.unmarkAll(cellKeepImplied);
//		mesh_.freeMark(cellKeepImplied);
//
//	}
//	// TODO a second time TESTING
//	for(int iPart=1;iPart<partsOrdering.size();iPart++){
//
//		TInt PartID = partsOrdering[iPart];
//
//		int cellKeep = mesh_.getNewMark();
//		int cellKeepImplied = mesh_.getNewMark();
//
//		// first we mark the partitioned cells (depends on the type) as cellKeep.
//		// we also recursively mark the descendants of the partitioned cells as cellKeep.
//
//		{
//		typename IGMesh::cells_iterator it  = mesh_.cells_begin();
//		typename IGMesh::cells_iterator ite = mesh_.cells_end();
//
//		for(;it!=ite;it++){
//			if((*var_dest)[(*it)->getID()]==PartID){
//				recursiveMarkWithSharedRepartition(*it,cellKeep,PartID,
//						var_nodes_lId,var_edges_lId,var_faces_lId,var_regions_lId,
//						var_nodes_owner,var_edges_owner,var_faces_owner,var_regions_owner,
//						cellSharedMark,false,false);
//			}
//		}
//		}
//
//		// mark the entities in the ghost layer if there is need of one.
//		if(MeshDescriptor<M>::NKnowN
//		|| MeshDescriptor<M>::NKnowE
//		|| MeshDescriptor<M>::NKnowF
//		|| MeshDescriptor<M>::NKnowR
//		|| MeshDescriptor<M>::EKnowE
//		|| MeshDescriptor<M>::EKnowF
//		|| MeshDescriptor<M>::EKnowR
//		|| MeshDescriptor<M>::FKnowF
//		|| MeshDescriptor<M>::FKnowR
//		|| MeshDescriptor<M>::RKnowR)
//		{
//			typename IGMesh::cells_iterator it  = mesh_.cells_begin();
//			typename IGMesh::cells_iterator ite = mesh_.cells_end();
//
//			for(;it!=ite;it++){
//				if(recursiveIsMarked(*it,cellKeep)){
//					mesh_.mark(*it,cellKeepImplied);
////					recursiveMarkWithSharedRepartition(*it,cellKeepImplied,PartID,
////							var_nodes_lId,var_edges_lId,var_faces_lId,var_regions_lId,
////							var_nodes_owner,var_edges_owner,var_faces_owner,var_regions_owner,
////							cellSharedMark,false,false);
//				}
//			}
//		}
//
//		{
//			typename IGMesh::cells_iterator it  = mesh_.cells_begin();
//			typename IGMesh::cells_iterator ite = mesh_.cells_end();
//
//			for(;it!=ite;it++){
//				if(mesh_.isMarked(*it,cellKeepImplied)){
//					recursiveMarkWithSharedRepartition(*it,cellKeep,PartID,
//							var_nodes_lId,var_edges_lId,var_faces_lId,var_regions_lId,
//							var_nodes_owner,var_edges_owner,var_faces_owner,var_regions_owner,
//							cellSharedMark,false,false);
//				}
//			}
//		}
//
//		getSubMeshRepartition(SubMeshes[iPart-1],cellKeep,partsOrdering[iPart],
//				var_nodes_lId,var_edges_lId,var_faces_lId,var_regions_lId,
//				var_nodes_owner,var_edges_owner,var_faces_owner,var_regions_owner,
//				cellSharedMark);
//
//		mesh_.unmarkAll(cellKeep);
//		mesh_.freeMark(cellKeep);
//		mesh_.unmarkAll(cellKeepImplied);
//		mesh_.freeMark(cellKeepImplied);
//
//	}
//
//	// we have to delete all the migrated entities in the local mesh.
//	{
//		TInt PartID = mesh_.getPartID();
//
//		int cellKeep = mesh_.getNewMark();
//		int cellKeepImplied = mesh_.getNewMark();
//
//		// first we mark the partitioned cells (depends on the type) as cellKeep.
//		// we also recursively mark the descendants of the partitioned cells as cellKeep.
//
//		typename IGMesh::cells_iterator it  = mesh_.cells_begin();
//		typename IGMesh::cells_iterator ite = mesh_.cells_end();
//
//		for(;it!=ite;it++){
//			if((*var_dest)[(*it)->getID()]==PartID){
//				recursiveMarkWithSharedRepartition(*it,cellKeep,PartID,
//						var_nodes_lId,var_edges_lId,var_faces_lId,var_regions_lId,
//						var_nodes_owner,var_edges_owner,var_faces_owner,var_regions_owner,
//						cellSharedMark,false,false);
//			}
//		}
//
//		// mark the entities in the ghost layer if there is need of one.
//		if(MeshDescriptor<M>::NKnowN
//		|| MeshDescriptor<M>::NKnowE
//		|| MeshDescriptor<M>::NKnowF
//		|| MeshDescriptor<M>::NKnowR
//		|| MeshDescriptor<M>::EKnowE
//		|| MeshDescriptor<M>::EKnowF
//		|| MeshDescriptor<M>::EKnowR
//		|| MeshDescriptor<M>::FKnowF
//		|| MeshDescriptor<M>::FKnowR
//		|| MeshDescriptor<M>::RKnowR)
//		{
//			typename IGMesh::cells_iterator it  = mesh_.cells_begin();
//			typename IGMesh::cells_iterator ite = mesh_.cells_end();
//
//			for(;it!=ite;it++){
//				// we mark every entity that has a cellKeep descendant as
//				// cellKeepImplied. Thus, we will not miss any sibling/upward
//				// connectivity.
//				// Problem is, we could restrict further the selection.
//				if(recursiveIsMarked(*it,cellKeep)){
//					recursiveMarkWithSharedRepartition(*it,cellKeepImplied,PartID,
//							var_nodes_lId,var_edges_lId,var_faces_lId,var_regions_lId,
//							var_nodes_owner,var_edges_owner,var_faces_owner,var_regions_owner,
//							cellSharedMark,false,false);
//				}
//			}
//		}
//
//		{
//			typename IGMesh::cells_iterator it  = mesh_.cells_begin();
//			typename IGMesh::cells_iterator ite = mesh_.cells_end();
//
//			for(;it!=ite;it++){
//				if(mesh_.isMarked(*it,cellKeepImplied)){
//					recursiveMarkWithSharedRepartition(*it,cellKeep,PartID,
//							var_nodes_lId,var_edges_lId,var_faces_lId,var_regions_lId,
//							var_nodes_owner,var_edges_owner,var_faces_owner,var_regions_owner,
//							cellSharedMark,false,false);
//				}
//			}
//		}
//
//		getSubMeshRepartitionRemove(cellKeep,PartID,
//				var_nodes_lId,var_edges_lId,var_faces_lId,var_regions_lId,
//				var_nodes_owner,var_edges_owner,var_faces_owner,var_regions_owner,
//				cellSharedMark);
//
//		mesh_.unmarkAll(cellKeep);
//		mesh_.freeMark(cellKeep);
//		mesh_.unmarkAll(cellKeepImplied);
//		mesh_.freeMark(cellKeepImplied);
//	}
//
//	mesh_.unmarkAll(cellSharedMark);
//	mesh_.freeMark(cellSharedMark);
//
//	if (MeshDescriptor<M>::hasNodes){
//		mesh_.deleteVariable(GMDS_NODE,"var_nodes_owner");
//		mesh_.deleteVariable(GMDS_NODE,"var_nodes_lId");
//	}
//	if (MeshDescriptor<M>::hasEdges){
//		mesh_.deleteVariable(GMDS_EDGE,"var_edges_owner");
//		mesh_.deleteVariable(GMDS_EDGE,"var_edges_lId");
//	}
//	if (MeshDescriptor<M>::hasFaces){
//		mesh_.deleteVariable(GMDS_FACE,"var_faces_owner");
//		mesh_.deleteVariable(GMDS_FACE,"var_faces_lId");
//	}
//	if (MeshDescriptor<M>::hasRegions){
//		mesh_.deleteVariable(GMDS_REGION,"var_regions_owner");
//		mesh_.deleteVariable(GMDS_REGION,"var_regions_lId");
//	}
//
//	return;
//}
///*----------------------------------------------------------------------------*/
//
//void SubMeshExtractor::partitioning(const ECellType ATypePartAlong, const TInt nParts)
//{
//	Variable<TInt> *var_dest = mesh_.template newVariable<TInt>(ATypePartAlong,"part");
//
//	id nCells = mesh_.getNbCells();
//	id nCellsParPart = ceil((float)nCells/(float)nParts);
//
//	typename IGMesh::cells_iterator it  = mesh_.cells_begin();
//	typename IGMesh::cells_iterator ite = mesh_.cells_end();
//
//	for(;it!=ite;it++){
//		(*var_dest)[(*it)->getID()] = (*it)->getID()/nCellsParPart;
//	}
//
//	return;
//}
///*----------------------------------------------------------------------------*/
//
//void SubMeshExtractor::
//recursiveMarkWithSharedRepartition(AbstractCell* ACell, const int AMark, const TInt PartID,
//		Variable<TCellID>* AVarNodeLID, Variable<TCellID>* AVarEdgeLID, Variable<TCellID>* AVarFaceLID, Variable<TCellID>* AVarRegionLID,
//		Variable<TInt>* AVarNodeOwner, Variable<TInt>* AVarEdgeOwner, Variable<TInt>* AVarFaceOwner, Variable<TInt>* AVarRegionOwner,
//		const int AMarkShared, const bool setPartID, const bool setLID)
//{
//	if(mesh_.isMarked(ACell,AMark)){
//		return;
//	}
//
//	mesh_.mark(ACell,AMark);
//
//	switch(ACell->getDim())
//	{
//	case 0 :
//	{
//		if ((*AVarNodeOwner)[ACell->getID()] == NullTInt){
//			if(setPartID){
//				(*AVarNodeOwner)[ACell->getID()] = PartID;
//			}
//			if(setLID){
//				(*AVarNodeLID)[ACell->getID()] = ACell->getID();
//			}
//		}
//		if ((*AVarNodeOwner)[ACell->getID()] != PartID){
//			mesh_.mark(ACell,AMarkShared);
//		}
//		return;
//		break;
//	}
//	case 1 :
//	{
//		if ((*AVarEdgeOwner)[ACell->getID()] == NullTInt){
//			if(setPartID){
//				(*AVarEdgeOwner)[ACell->getID()] = PartID;
//			}
//			if(setLID){
//				(*AVarEdgeLID)[ACell->getID()] = ACell->getID();
//			}
//		}
//		if ((*AVarEdgeOwner)[ACell->getID()] != PartID){
//			mesh_.mark(ACell,AMarkShared);
//		}
//
//		if(MeshDescriptor<M>::hasNodes && MeshDescriptor<M>::EKnowN){
//			std::vector<Node*> nodes = ACell->getNodes();
//
//			std::vector<Node*>::iterator it  = nodes.begin();
//			std::vector<Node*>::iterator ite = nodes.end();
//
//			for(;it!=ite;it++){
//				recursiveMarkWithSharedRepartition(*it,AMark,PartID,
//						AVarNodeLID,AVarEdgeLID,AVarFaceLID,AVarRegionLID,
//						AVarNodeOwner, AVarEdgeOwner, AVarFaceOwner, AVarRegionOwner,
//						AMarkShared,setPartID,setLID);
//			}
//			return;
//		}
//		break;
//	}
//	case 2 :
//	{
//		if ((*AVarFaceOwner)[ACell->getID()] == NullTInt){
//			if(setPartID){
//				(*AVarFaceOwner)[ACell->getID()] = PartID;
//			}
//			if(setLID){
//				(*AVarFaceLID)[ACell->getID()] = ACell->getID();
//			}
//		}
//		if ((*AVarFaceOwner)[ACell->getID()] != PartID){
//			mesh_.mark(ACell,AMarkShared);
//		}
//
//		if(MeshDescriptor<M>::hasEdges && MeshDescriptor<M>::FKnowE){
//			std::vector<Edge*> edges = ACell->getEdges();
//
//			std::vector<Edge*>::iterator it  = edges.begin();
//			std::vector<Edge*>::iterator ite = edges.end();
//
//			for(;it!=ite;it++){
//				recursiveMarkWithSharedRepartition(*it,AMark,PartID,
//						AVarNodeLID,AVarEdgeLID,AVarFaceLID,AVarRegionLID,
//						AVarNodeOwner, AVarEdgeOwner, AVarFaceOwner, AVarRegionOwner,
//						AMarkShared,setPartID,setLID);
//			}
//			//return;
//		}
//
//		if(MeshDescriptor<M>::hasNodes && MeshDescriptor<M>::FKnowN){
//			std::vector<Node*> nodes = ACell->getNodes();
//
//			std::vector<Node*>::iterator it  = nodes.begin();
//			std::vector<Node*>::iterator ite = nodes.end();
//
//			for(;it!=ite;it++){
//				recursiveMarkWithSharedRepartition(*it,AMark,PartID,
//						AVarNodeLID,AVarEdgeLID,AVarFaceLID,AVarRegionLID,
//						AVarNodeOwner, AVarEdgeOwner, AVarFaceOwner, AVarRegionOwner,
//						AMarkShared,setPartID,setLID);
//			}
//			//return;
//		}
//
//		break;
//	}
//	case 3 :
//	{
//		if ((*AVarRegionOwner)[ACell->getID()] == NullTInt){
//			if(setPartID){
//				(*AVarRegionOwner)[ACell->getID()] = PartID;
//			}
//			if(setLID){
//				(*AVarRegionLID)[ACell->getID()] = ACell->getID();
//			}
//		}
//		if ((*AVarRegionOwner)[ACell->getID()] != PartID){
//			mesh_.mark(ACell,AMarkShared);
//		}
//
//		if(MeshDescriptor<M>::hasFaces && MeshDescriptor<M>::RKnowF){
//			std::vector<Face*> faces = ACell->getFaces();
//
//			std::vector<Face*>::iterator it  = faces.begin();
//			std::vector<Face*>::iterator ite = faces.end();
//
//			for(;it!=ite;it++){
//				recursiveMarkWithSharedRepartition(*it,AMark,PartID,
//						AVarNodeLID,AVarEdgeLID,AVarFaceLID,AVarRegionLID,
//						AVarNodeOwner, AVarEdgeOwner, AVarFaceOwner, AVarRegionOwner,
//						AMarkShared,setPartID,setLID);
//			}
//			//return;
//		}
//
//		if(MeshDescriptor<M>::hasEdges && MeshDescriptor<M>::RKnowE){
//			std::vector<Edge*> edges = ACell->getEdges();
//
//			std::vector<Edge*>::iterator it  = edges.begin();
//			std::vector<Edge*>::iterator ite = edges.end();
//
//			for(;it!=ite;it++){
//				recursiveMarkWithSharedRepartition(*it,AMark,PartID,
//						AVarNodeLID,AVarEdgeLID,AVarFaceLID,AVarRegionLID,
//						AVarNodeOwner, AVarEdgeOwner, AVarFaceOwner, AVarRegionOwner,
//						AMarkShared,setPartID,setLID);
//			}
//			//return;
//		}
//
//		if(MeshDescriptor<M>::hasNodes && MeshDescriptor<M>::RKnowN){
//			std::vector<Node*> nodes = ACell->getNodes();
//
//			std::vector<Node*>::iterator it  = nodes.begin();
//			std::vector<Node*>::iterator ite = nodes.end();
//
//			for(;it!=ite;it++){
//				recursiveMarkWithSharedRepartition(*it,AMark,PartID,
//						AVarNodeLID,AVarEdgeLID,AVarFaceLID,AVarRegionLID,
//						AVarNodeOwner, AVarEdgeOwner, AVarFaceOwner, AVarRegionOwner,
//						AMarkShared,setPartID,setLID);
//			}
//
//			//return;
//		}
//
//		break;
//	}
//	default :
//		std::cout << "Not a valid part type!!" << std::endl;
//		throw GMDSException();
//		return;
//	}
//
//	return;
//}
///*----------------------------------------------------------------------------*/
//
//void SubMeshExtractor::
//recursiveMarkWithSharedRepartitionExp(AbstractCell* ACell, TInt PartID,
//		Variable<std::vector<TInt> >* AVarNodeOwners, Variable<std::vector<TInt> >* AVarEdgeOwners, Variable<std::vector<TInt> >* AVarFaceOwners, Variable<std::vector<TInt> >* AVarRegionOwners,
//		const int AMarkShared)
//{
//	switch(ACell->getDim())
//	{
//	case 0 :
//	{
//		for(int i=0; i<(*AVarNodeOwners)[ACell->getID()].size(); i++){
//			if(((*AVarNodeOwners)[ACell->getID()])[i]==PartID){
//				return;
//			}
//		}
//
//		(*AVarNodeOwners)[ACell->getID()].push_back(PartID);
//
//		if((*AVarNodeOwners)[ACell->getID()].size()>1){
//			mesh_.mark(ACell,AMarkShared);
//		}
//
//		return;
//		break;
//	}
//	case 1 :
//	{
//		for(int i=0; i<(*AVarEdgeOwners)[ACell->getID()].size(); i++){
//			if(((*AVarEdgeOwners)[ACell->getID()])[i]==PartID){
//				return;
//			}
//		}
//
//		(*AVarEdgeOwners)[ACell->getID()].push_back(PartID);
//
//		if((*AVarEdgeOwners)[ACell->getID()].size()>1){
//			mesh_.mark(ACell,AMarkShared);
//		}
//
//		if(MeshDescriptor<M>::hasNodes && MeshDescriptor<M>::EKnowN){
//			std::vector<Node*> nodes = ACell->getNodes();
//
//			std::vector<Node*>::iterator it  = nodes.begin();
//			std::vector<Node*>::iterator ite = nodes.end();
//
//			for(;it!=ite;it++){
//				recursiveMarkWithSharedRepartitionExp(*it,PartID,
//						AVarNodeOwners, AVarEdgeOwners, AVarFaceOwners, AVarRegionOwners,
//						AMarkShared);
//			}
//			return;
//		}
//		break;
//	}
//	case 2 :
//	{
//		for(int i=0; i<(*AVarFaceOwners)[ACell->getID()].size(); i++){
//			if(((*AVarFaceOwners)[ACell->getID()])[i]==PartID){
//				return;
//			}
//		}
//
//		(*AVarFaceOwners)[ACell->getID()].push_back(PartID);
//
//		if((*AVarFaceOwners)[ACell->getID()].size()>1){
//			mesh_.mark(ACell,AMarkShared);
//		}
//
//		if(MeshDescriptor<M>::hasEdges && MeshDescriptor<M>::FKnowE){
//			std::vector<Edge*> edges = ACell->getEdges();
//
//			std::vector<Edge*>::iterator it  = edges.begin();
//			std::vector<Edge*>::iterator ite = edges.end();
//
//			for(;it!=ite;it++){
//				recursiveMarkWithSharedRepartitionExp(*it,PartID,
//						AVarNodeOwners, AVarEdgeOwners, AVarFaceOwners, AVarRegionOwners,
//						AMarkShared);
//			}
//			//return;
//		}
//
//		if(MeshDescriptor<M>::hasNodes && MeshDescriptor<M>::FKnowN){
//			std::vector<Node*> nodes = ACell->getNodes();
//
//			std::vector<Node*>::iterator it  = nodes.begin();
//			std::vector<Node*>::iterator ite = nodes.end();
//
//			for(;it!=ite;it++){
//				recursiveMarkWithSharedRepartitionExp(*it,PartID,
//						AVarNodeOwners, AVarEdgeOwners, AVarFaceOwners, AVarRegionOwners,
//						AMarkShared);
//			}
//			//return;
//		}
//
//		break;
//	}
//	case 3 :
//	{
//		for(int i=0; i<(*AVarRegionOwners)[ACell->getID()].size(); i++){
//			if(((*AVarRegionOwners)[ACell->getID()])[i]==PartID){
//				return;
//			}
//		}
//
//		(*AVarRegionOwners)[ACell->getID()].push_back(PartID);
//
//		if((*AVarRegionOwners)[ACell->getID()].size()>1){
//			mesh_.mark(ACell,AMarkShared);
//		}
//
//		if(MeshDescriptor<M>::hasFaces && MeshDescriptor<M>::RKnowF){
//			std::vector<Face*> faces = ACell->getFaces();
//
//			std::vector<Face*>::iterator it  = faces.begin();
//			std::vector<Face*>::iterator ite = faces.end();
//
//			for(;it!=ite;it++){
//				recursiveMarkWithSharedRepartitionExp(*it,PartID,
//						AVarNodeOwners, AVarEdgeOwners, AVarFaceOwners, AVarRegionOwners,
//						AMarkShared);
//			}
//			//return;
//		}
//
//		if(MeshDescriptor<M>::hasEdges && MeshDescriptor<M>::RKnowE){
//			std::vector<Edge*> edges = ACell->getEdges();
//
//			std::vector<Edge*>::iterator it  = edges.begin();
//			std::vector<Edge*>::iterator ite = edges.end();
//
//			for(;it!=ite;it++){
//				recursiveMarkWithSharedRepartitionExp(*it,PartID,
//						AVarNodeOwners, AVarEdgeOwners, AVarFaceOwners, AVarRegionOwners,
//						AMarkShared);
//			}
//			//return;
//		}
//
//		if(MeshDescriptor<M>::hasNodes && MeshDescriptor<M>::RKnowN){
//			std::vector<Node*> nodes = ACell->getNodes();
//
//			std::vector<Node*>::iterator it  = nodes.begin();
//			std::vector<Node*>::iterator ite = nodes.end();
//
//			for(;it!=ite;it++){
//				recursiveMarkWithSharedRepartitionExp(*it,PartID,
//						AVarNodeOwners, AVarEdgeOwners, AVarFaceOwners, AVarRegionOwners,
//						AMarkShared);
//			}
//
//			//return;
//		}
//
//		break;
//	}
//	default :
//		std::cout << "Not a valid part type!!" << std::endl;
//		throw GMDSException();
//		return;
//	}
//
//	return;
//}
///*----------------------------------------------------------------------------*/
//
//void SubMeshExtractor::
//recursiveMarkWithSharedRepartitionExp2(AbstractCell* ACell, TInt PartID,
//		bool* isSharedNode,bool* isSharedEdge,bool* isSharedFace,bool* isSharedRegion,
//		id* nodes_LID, id* edges_LID, id* faces_LID, id* regions_LID,
//		TInt* nodes_owner,TInt* edges_owner, TInt* faces_owner,TInt* regions_owner,
//		std::vector<std::map<TInt,TCellID> >& nodes_owners2LIDs, std::vector<std::map<TInt,TCellID> >& edges_owners2LIDs, std::vector<std::map<TInt,TCellID> >& faces_owners2LIDs, std::vector<std::map<TInt,TCellID> >& regions_owners2LIDs)
//{
//	switch(ACell->getDim())
//	{
//	case 0 :
//	{
//		// entity encountered for the first time
//		if(nodes_owner[ACell->getID()] == NullTInt){
//			nodes_owner[ACell->getID()] = PartID;
//		}
//		else{
//			// entity already marked by another part
//			if(nodes_owner[ACell->getID()] != PartID){
//				// multiple occurrences of this entity are known
//				if(isSharedNode[ACell->getID()]){
//					if(nodes_owners2LIDs[nodes_LID[ACell->getID()]].find(PartID) == nodes_owners2LIDs[nodes_LID[ACell->getID()]].end()){
//						nodes_owners2LIDs[nodes_LID[ACell->getID()]][PartID] = NullID;
//					}
//					// entity already encountered by this part
//					else{
//					}
//				}
//				// only one occurrence is known
//				else{
//					isSharedNode[ACell->getID()] = true;
//					nodes_LID[ACell->getID()] = nodes_owners2LIDs.size();
//					nodes_owners2LIDs.resize(nodes_owners2LIDs.size()+1);
//					nodes_owners2LIDs[nodes_LID[ACell->getID()]][nodes_owner[ACell->getID()]] = NullID;
//					nodes_owners2LIDs[nodes_LID[ACell->getID()]][PartID] = NullID;
//				}
//			}
//			// entity already encountered by this part
//			else{
//			}
//		}
//
//		return;
//		break;
//	}
//	case 1 :
//	{
//		// entity encountered for the first time
//		if(edges_owner[ACell->getID()] == NullTInt){
//			edges_owner[ACell->getID()] = PartID;
//		}
//		else{
//			// entity already marked by another part
//			if(edges_owner[ACell->getID()] != PartID){
//				// multiple occurrences of this entity are known
//				if(isSharedEdge[ACell->getID()]){
//					if(edges_owners2LIDs[edges_LID[ACell->getID()]].find(PartID) == edges_owners2LIDs[edges_LID[ACell->getID()]].end()){
//						edges_owners2LIDs[edges_LID[ACell->getID()]][PartID] = NullID;
//					}
//				}
//				// only one occurrence is known
//				else{
//					isSharedEdge[ACell->getID()] = true;
//					edges_LID[ACell->getID()] = edges_owners2LIDs.size();
//					edges_owners2LIDs.resize(edges_owners2LIDs.size()+1);
//					edges_owners2LIDs[edges_LID[ACell->getID()]][edges_owner[ACell->getID()]] = NullID;
//					edges_owners2LIDs[edges_LID[ACell->getID()]][PartID] = NullID;
//				}
//			}
//			// entity already encountered by this part
//			else{
//			}
//		}
//
//		if(MeshDescriptor<M>::hasNodes && MeshDescriptor<M>::EKnowN){
//			std::vector<Node*> nodes = ACell->getNodes();
//
//			std::vector<Node*>::iterator it  = nodes.begin();
//			std::vector<Node*>::iterator ite = nodes.end();
//
//			for(;it!=ite;it++){
//				recursiveMarkWithSharedRepartitionExp2(*it,PartID,
//						isSharedNode,isSharedEdge,isSharedFace,isSharedRegion,
//						nodes_LID,edges_LID,faces_LID,regions_LID,
//						nodes_owner,edges_owner,faces_owner,regions_owner,
//						nodes_owners2LIDs,edges_owners2LIDs,faces_owners2LIDs,regions_owners2LIDs);
//			}
//		}
//
//		return;
//		break;
//	}
//	case 2 :
//	{
//		// entity encountered for the first time
//		if(faces_owner[ACell->getID()] == NullTInt){
//			faces_owner[ACell->getID()] = PartID;
//		}
//		else{
//			// entity already marked by another part
//			if(faces_owner[ACell->getID()] != PartID){
//				// multiple occurrences of this entity are known
//				if(isSharedFace[ACell->getID()]){
//					if(faces_owners2LIDs[faces_LID[ACell->getID()]].find(PartID) == faces_owners2LIDs[faces_LID[ACell->getID()]].end()){
//						faces_owners2LIDs[faces_LID[ACell->getID()]][PartID] = NullID;
//					}
//				}
//				// only one occurence is known
//				else{
//					isSharedFace[ACell->getID()] = true;
//					faces_LID[ACell->getID()] = faces_owners2LIDs.size();
//					faces_owners2LIDs.resize(faces_owners2LIDs.size()+1);
//					faces_owners2LIDs[faces_LID[ACell->getID()]][faces_owner[ACell->getID()]] = NullID;
//					faces_owners2LIDs[faces_LID[ACell->getID()]][PartID] = NullID;
//				}
//			}
//			// entity already encountered by this part
//			else{
//			}
//		}
//
//		if(MeshDescriptor<M>::hasEdges && MeshDescriptor<M>::FKnowE){
//			std::vector<Edge*> edges = ACell->getEdges();
//
//			std::vector<Edge*>::iterator it  = edges.begin();
//			std::vector<Edge*>::iterator ite = edges.end();
//
//			for(;it!=ite;it++){
//				recursiveMarkWithSharedRepartitionExp2(*it,PartID,
//						isSharedNode,isSharedEdge,isSharedFace,isSharedRegion,
//						nodes_LID,edges_LID,faces_LID,regions_LID,
//						nodes_owner,edges_owner,faces_owner,regions_owner,
//						nodes_owners2LIDs,edges_owners2LIDs,faces_owners2LIDs,regions_owners2LIDs);
//			}
//		}
//
//		if(MeshDescriptor<M>::hasNodes && MeshDescriptor<M>::FKnowN){
//			std::vector<Node*> nodes = ACell->getNodes();
//
//			std::vector<Node*>::iterator it  = nodes.begin();
//			std::vector<Node*>::iterator ite = nodes.end();
//
//			for(;it!=ite;it++){
//				recursiveMarkWithSharedRepartitionExp2(*it,PartID,
//						isSharedNode,isSharedEdge,isSharedFace,isSharedRegion,
//						nodes_LID,edges_LID,faces_LID,regions_LID,
//						nodes_owner,edges_owner,faces_owner,regions_owner,
//						nodes_owners2LIDs,edges_owners2LIDs,faces_owners2LIDs,regions_owners2LIDs);
//			}
//		}
//
//		return;
//		break;
//	}
//	case 3 :
//	{
//		// entity encountered for the first time
//		if(regions_owner[ACell->getID()] == NullTInt){
//			regions_owner[ACell->getID()] = PartID;
//		}
//		else{
//			// entity already marked by another part
//			if(regions_owner[ACell->getID()] != PartID){
//				// multiple occurrences of this entity are known
//				if(isSharedRegion[ACell->getID()]){
//					if(regions_owners2LIDs[regions_LID[ACell->getID()]].find(PartID) == regions_owners2LIDs[regions_LID[ACell->getID()]].end()){
//						regions_owners2LIDs[regions_LID[ACell->getID()]][PartID] = NullID;
//					}
//				}
//				// only one occurrence is known
//				else{
//					isSharedRegion[ACell->getID()] = true;
//					regions_LID[ACell->getID()] = regions_owners2LIDs.size();
//					regions_owners2LIDs.resize(regions_owners2LIDs.size()+1);
//					regions_owners2LIDs[regions_LID[ACell->getID()]][regions_owner[ACell->getID()]] = NullID;
//					regions_owners2LIDs[regions_LID[ACell->getID()]][PartID] = NullID;
//				}
//			}
//			// entity already encountered by this part
//			else{
//			}
//		}
//
//		if(MeshDescriptor<M>::hasFaces && MeshDescriptor<M>::RKnowF){
//			std::vector<Face*> faces = ACell->getFaces();
//
//			std::vector<Face*>::iterator it  = faces.begin();
//			std::vector<Face*>::iterator ite = faces.end();
//
//			for(;it!=ite;it++){
//				recursiveMarkWithSharedRepartitionExp2(*it,PartID,
//						isSharedNode,isSharedEdge,isSharedFace,isSharedRegion,
//						nodes_LID,edges_LID,faces_LID,regions_LID,
//						nodes_owner,edges_owner,faces_owner,regions_owner,
//						nodes_owners2LIDs,edges_owners2LIDs,faces_owners2LIDs,regions_owners2LIDs);
//			}
//		}
//
//		if(MeshDescriptor<M>::hasEdges && MeshDescriptor<M>::RKnowE){
//			std::vector<Edge*> edges = ACell->getEdges();
//
//			std::vector<Edge*>::iterator it  = edges.begin();
//			std::vector<Edge*>::iterator ite = edges.end();
//
//			for(;it!=ite;it++){
//				recursiveMarkWithSharedRepartitionExp2(*it,PartID,
//						isSharedNode,isSharedEdge,isSharedFace,isSharedRegion,
//						nodes_LID,edges_LID,faces_LID,regions_LID,
//						nodes_owner,edges_owner,faces_owner,regions_owner,
//						nodes_owners2LIDs,edges_owners2LIDs,faces_owners2LIDs,regions_owners2LIDs);
//			}
//		}
//
//		if(MeshDescriptor<M>::hasNodes && MeshDescriptor<M>::RKnowN){
//			std::vector<Node*> nodes = ACell->getNodes();
//
//			std::vector<Node*>::iterator it  = nodes.begin();
//			std::vector<Node*>::iterator ite = nodes.end();
//
//			for(;it!=ite;it++){
//				recursiveMarkWithSharedRepartitionExp2(*it,PartID,
//						isSharedNode,isSharedEdge,isSharedFace,isSharedRegion,
//						nodes_LID,edges_LID,faces_LID,regions_LID,
//						nodes_owner,edges_owner,faces_owner,regions_owner,
//						nodes_owners2LIDs,edges_owners2LIDs,faces_owners2LIDs,regions_owners2LIDs);
//			}
//		}
//
//		return;
//		break;
//	}
//	default :
//		std::cout << "Not a valid part type!!" << std::endl;
//		throw GMDSException();
//		return;
//	}
//
//	return;
//}
///*----------------------------------------------------------------------------*/
//
//bool SubMeshExtractor::
//recursiveIsMarked(AbstractCell* ACell, const int AMark)
//{
//	if(mesh_.isMarked(ACell,AMark)){
//		return true;
//	}
//
//	switch(ACell->getDim())
//	{
//	case 0 :
//	{
//		break;
//	}
//	case 1 :
//	{
//		if(MeshDescriptor<M>::hasNodes && MeshDescriptor<M>::EKnowN){
//			std::vector<Node*> nodes = ACell->getNodes();
//
//			std::vector<Node*>::iterator it  = nodes.begin();
//			std::vector<Node*>::iterator ite = nodes.end();
//
//			for(;it!=ite;it++){
//				if(recursiveIsMarked(*it,AMark)){
//					return true;
//				}
//			}
//		}
//		break;
//	}
//	case 2 :
//	{
//		if(MeshDescriptor<M>::hasEdges && MeshDescriptor<M>::FKnowE){
//			std::vector<Edge*> edges = ACell->getEdges();
//
//			std::vector<Edge*>::iterator it  = edges.begin();
//			std::vector<Edge*>::iterator ite = edges.end();
//
//			for(;it!=ite;it++){
//				if(recursiveIsMarked(*it,AMark)){
//					return true;
//				}
//			}
//		}
//
//		if(MeshDescriptor<M>::hasNodes && MeshDescriptor<M>::FKnowN){
//			std::vector<Node*> nodes = ACell->getNodes();
//
//			std::vector<Node*>::iterator it  = nodes.begin();
//			std::vector<Node*>::iterator ite = nodes.end();
//
//			for(;it!=ite;it++){
//				if(recursiveIsMarked(*it,AMark)){
//					return true;
//				}
//			}
//		}
//
//		break;
//	}
//	case 3 :
//	{
//		if(MeshDescriptor<M>::hasFaces && MeshDescriptor<M>::RKnowF){
//			std::vector<Face*> faces = ACell->getFaces();
//
//			std::vector<Face*>::iterator it  = faces.begin();
//			std::vector<Face*>::iterator ite = faces.end();
//
//			for(;it!=ite;it++){
//				if(recursiveIsMarked(*it,AMark)){
//					return true;
//				}
//			}
//		}
//
//		if(MeshDescriptor<M>::hasEdges && MeshDescriptor<M>::RKnowE){
//			std::vector<Edge*> edges = ACell->getEdges();
//
//			std::vector<Edge*>::iterator it  = edges.begin();
//			std::vector<Edge*>::iterator ite = edges.end();
//
//			for(;it!=ite;it++){
//				if(recursiveIsMarked(*it,AMark)){
//					return true;
//				}
//			}
//		}
//
//		if(MeshDescriptor<M>::hasNodes && MeshDescriptor<M>::RKnowN){
//			std::vector<Node*> nodes = ACell->getNodes();
//
//			std::vector<Node*>::iterator it  = nodes.begin();
//			std::vector<Node*>::iterator ite = nodes.end();
//
//			for(;it!=ite;it++){
//				if(recursiveIsMarked(*it,AMark)){
//					return true;
//				}
//			}
//		}
//
//		break;
//	}
//	default :
//		std::cout << "Not a valid part type!!" << std::endl;
//		throw GMDSException();
//		return false;
//	}
//
//	return false;
//
//}
///*----------------------------------------------------------------------------*/
//
//void SubMeshExtractor::
//recursiveGhostMarkingExp2(AbstractCell* ACell,std::vector<TInt>& parts2add,
//		bool* isSharedNode,bool* isSharedEdge,bool* isSharedFace,bool* isSharedRegion,
//		id* nodes_LID, id* edges_LID, id* faces_LID, id* regions_LID,
//		TInt* nodes_owner,TInt* edges_owner, TInt* faces_owner,TInt* regions_owner,
//		std::vector<std::map<TInt,TCellID> >& nodes_owners2LIDs, std::vector<std::map<TInt,TCellID> >& edges_owners2LIDs, std::vector<std::map<TInt,TCellID> >& faces_owners2LIDs, std::vector<std::map<TInt,TCellID> >& regions_owners2LIDs)
//{
//	if(!MeshDescriptor<M>::CKnowN){
//		std::cout << "CELLS SHOULD KNOW NODES : recursiveGhostMarkingExp2" << std::endl;
//		throw GMDSException();
//	}
//
//	std::list<TInt> allParts;
//
//	std::vector<TCellID> nodesIDs = ACell->getNodeIDs();
//
//	std::vector<TCellID>::iterator nodeID = nodesIDs.begin();
//	std::vector<TCellID>::iterator nodeIDe = nodesIDs.end();
//
//	for(;nodeID!=nodeIDe;nodeID++){
//		if(isSharedNode[*nodeID]){
//
//			std::map<TInt,TCellID>::iterator itNodeOwner2LID = nodes_owners2LIDs[nodes_LID[*nodeID]].begin();
//			std::map<TInt,TCellID>::iterator itNodeOwner2LIDe = nodes_owners2LIDs[nodes_LID[*nodeID]].end();
//
//			for(;itNodeOwner2LID!=itNodeOwner2LIDe;itNodeOwner2LID++){
//				allParts.push_back(itNodeOwner2LID->first);
//			}
//		}
//	}
//
//	allParts.sort();
//	allParts.unique();
//
//	std::list<TInt>::iterator allPartsID = allParts.begin();
//	std::list<TInt>::iterator allPartsIDe = allParts.end();
//
//	for(;allPartsID!=allPartsIDe;allPartsID++){
//		parts2add.push_back(*allPartsID);
//	}
//
//	return;
//}

/*----------------------------------------------------------------------------*/
