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
/** \file    GETMe.cpp
 *  \author  legoff
 *  \date    10/16/2015
 */
/*----------------------------------------------------------------------------*/
#include "GMDS/Algo/GETMe.h"
#include "GMDS/CAD/GeomEntity.h"
#include "GMDS/Math/Segment.h"
#include "GMDS/Math/Triangle.h"
/*----------------------------------------------------------------------------*/
#include <sstream>
#include <set>
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/
GETMe::GETMe(IGMesh& m_mesh)
  :m_mesh(m_mesh),m_isVolumeTargetOn(false),m_volumeTargetVariableName("")
{}
/*----------------------------------------------------------------------------*/
GETMe::~GETMe()
{}
/*----------------------------------------------------------------------------*/
void
GETMe::setVolumeTarget(bool ABool)
{
  m_isVolumeTargetOn = ABool;
}
/*----------------------------------------------------------------------------*/
void
GETMe::setVolumeTargetVariableName(std::string AVariableName)
{
  m_volumeTargetVariableName = AVariableName;
}
/*----------------------------------------------------------------------------*/
void 
GETMe::exec(
		int ANbMaxIterSimult,
		int ANbMaxIterSeq,
		double AQualityThreshold,
		bool AAvoidTangledElements,
		int AMarkFixedNodes)
{
	// simultaneous GETMe
	execSimult(
			ANbMaxIterSimult,
			AQualityThreshold,
			0.5,
			1.0,
			1.0,
			1.0,
			false,
			false,
			AAvoidTangledElements,
			AMarkFixedNodes
	);

	execSeq(
		ANbMaxIterSeq,
		AQualityThreshold,
		0.01,
		1.0,
		0.0002,
		0.0004,
		0.0002,
		false,
		false,
		false,
		AMarkFixedNodes
	  );
}
/*----------------------------------------------------------------------------*/
void
GETMe::execSimult(
		int ANbMaxIter,
		double AQualityThreshold,
		double ARelaxFactor,
		double AQualityDependantFactorMin,
		double AQualityDependantFactorMax,
		double AWeightExponent,
		bool AIsGeomAssociated,
		bool AIsGeomAssocForced,
		bool AAvoidTangledElements,
		int AMarkFixedNodes)
{
	// check for parameters consistency
	if(AAvoidTangledElements) {
		if(computeMinNormalizedScaledJacobian() < 0.) {
			throw GMDSException("GETMe::execSimult AAvoidTangledElements can not be set to true when min scaled jacobian is negative");
		}
	}

	// execution control variables
	int numIter = 0;	
	
	for(; numIter<ANbMaxIter && computeMinNormalizedScaledJacobian()<AQualityThreshold; numIter++) {

std::cout<<"GETMe::execSimult numIter/ANbMaxIter "<<numIter<<" "<<ANbMaxIter<<std::endl;
std::cout<<"computeMinNormalizedScaledJacobian "<<computeMinNormalizedScaledJacobian()<<std::endl;
	// 
        std::map<TCellID, gmds::math::Point> nodesNewPos;
	std::map<TCellID, gmds::math::Point> nodesOldPos;
        std::map<TCellID, double> nodesNewPosWeight;

	// transformation of all elements
        {
	gmds::IGMesh::region_iterator itr  = this->m_mesh.regions_begin();
        for(;!itr.isDone();itr.next()) {
                gmds::Region current_region = itr.value();
		
		// cellQuality as used here should range between [0.,1.]
		double cellQuality = current_region.computeNormalizedScaledJacobian();
		cellQuality = (cellQuality + 1.)/2.;

                // compute a quality dependant factor
                double qualityDependantFactor = AQualityDependantFactorMin + (AQualityDependantFactorMax - AQualityDependantFactorMin) * (1. - cellQuality);

		// compute new Positions
		std::vector<gmds::Node> nodes = current_region.get<gmds::Node>();
                std::vector<gmds::math::Point> newPos(nodes.size());
		computeGETMePoint(&current_region,qualityDependantFactor,newPos);

		// apply scaling 
		applyScaling(&current_region,newPos,1.);

		// apply relaxation
		for(int iPoint=0; iPoint<nodes.size(); iPoint++) {
                        newPos[iPoint] = (1. - ARelaxFactor) * nodes[iPoint].getPoint() + ARelaxFactor*newPos[iPoint]; 
                }

		// compute weight
		double newPosWeight = std::pow((1. - cellQuality),AWeightExponent);	

		// store weight and position
		for(int iPoint=0; iPoint<nodes.size(); iPoint++) {
			if(nodesNewPosWeight.find(nodes[iPoint].getID()) != nodesNewPosWeight.end()) {
				nodesNewPosWeight[nodes[iPoint].getID()] += newPosWeight;
				nodesNewPos[nodes[iPoint].getID()] = nodesNewPos[nodes[iPoint].getID()] + newPosWeight*newPos[iPoint];
			} else {
				nodesNewPosWeight[nodes[iPoint].getID()] = newPosWeight;
				nodesNewPos[nodes[iPoint].getID()] = newPosWeight*newPos[iPoint];
			}
		}
	
        } // for(;!itr.isDone();itr.next())

	// divide by the weight
	gmds::IGMesh::node_iterator itn  = this->m_mesh.nodes_begin();
        for(;!itn.isDone();itn.next()) {
                gmds::Node current_node = itn.value();

		if(nodesNewPosWeight[current_node.getID()] != 0.) { 
	                nodesNewPos[current_node.getID()] = (1./nodesNewPosWeight[current_node.getID()]) * nodesNewPos[current_node.getID()];
		} else {
			nodesNewPos[current_node.getID()] = current_node.getPoint();
		}
         
        }
	
	}

	// apply geometric association if any
        {
	if(AIsGeomAssociated) {
		double geomAssocFactor = (numIter+1)/ANbMaxIter;

		if(AIsGeomAssocForced) {
			geomAssocFactor = 1;
		}

        	gmds::IGMesh::node_iterator itn  = this->m_mesh.nodes_begin();
	        for(;!itn.isDone();itn.next()) {
        	        gmds::Node current_node = itn.value();

			if(AMarkFixedNodes == -1 || (AMarkFixedNodes != -1 && !this->m_mesh.isMarked(current_node,AMarkFixedNodes))) {
			  geom::GeomEntity* refGeometry = m_mesh.getGeometricClassification(current_node);
				if(refGeometry != NULL && refGeometry->getDim() != 3) {
					gmds::math::Point pt = nodesNewPos[current_node.getID()];
					refGeometry->project(pt);

					nodesNewPos[current_node.getID()] = (1. - geomAssocFactor) * nodesNewPos[current_node.getID()] + geomAssocFactor * pt;
				}
			}
		}	
        }
        }


	// update nodes positions
	{
	gmds::IGMesh::node_iterator itn  = this->m_mesh.nodes_begin();
        for(;!itn.isDone();itn.next()) {
                gmds::Node current_node = itn.value();

		nodesOldPos[current_node.getID()] = current_node.getPoint();

		if(AMarkFixedNodes == -1 || (AMarkFixedNodes != -1 && !this->m_mesh.isMarked(current_node,AMarkFixedNodes))) {
			current_node.setPoint(nodesNewPos[current_node.getID()]);
		}
	}
	}

	// invalid element handling
	if (AAvoidTangledElements) {
 
	bool tangledCellDetected = false;
	do {
		tangledCellDetected = false;		
		std::set<gmds::TCellID> nodes2RollBack;

		// detect invalid cells
		gmds::IGMesh::region_iterator itr  = this->m_mesh.regions_begin();
        	for(;!itr.isDone();itr.next()) {
                	gmds::Region current_region = itr.value();

			double scaledJacobian = current_region.computeNormalizedScaledJacobian();
			if(scaledJacobian <= 0.) {
				tangledCellDetected = true;
				std::vector<gmds::Node> nodes = current_region.get<gmds::Node>();
				for(int iNode=0; iNode<nodes.size(); iNode++) {
					nodes2RollBack.insert(nodes[iNode].getID());
				}	
			}
		}
	
		// rollback nodes
		std::set<gmds::TCellID>::iterator itn = nodes2RollBack.begin();
		for(; itn!=nodes2RollBack.end(); itn++) {
			(this->m_mesh.get<gmds::Node>(*itn)).setPoint(nodesOldPos[(*itn)]);
		}

	} while(tangledCellDetected);

	} // if (AAvoidTangledElements)	

std::cout<<"computeMinNormalizedScaledJacobian "<<computeMinNormalizedScaledJacobian()<<std::endl;
	} // for(; numIter<ANbMaxIter; numIter++)

}
/*----------------------------------------------------------------------------*/
void
GETMe::execSimult2D(
		int ANbMaxIter,
		double AQualityThreshold,
		double ARelaxFactor,
		double AQualityDependantFactorMin,
		double AQualityDependantFactorMax,
		double AWeightExponent,
		bool AIsGeomAssociated,
		bool AIsGeomAssocForced,
		bool AAvoidTangledElements,
		int AMarkFixedNodes)
{
	// check for parameters consistency
        if(AAvoidTangledElements) {
                if(computeMinScaledJacobian2D() < 0.) {
                        throw GMDSException("GETMe::execSimult2D AAvoidTangledElements can not be set to true when min scaled jacobian is negative");
                }
        }

	// execution control variables
	int numIter = 0;	
	
	for(; numIter<ANbMaxIter && computeMinScaledJacobian2D()<AQualityThreshold; numIter++) {

std::cout<<"GETMe::execSimult2D numIter/ANbMaxIter "<<numIter<<" "<<ANbMaxIter<<std::endl;
std::cout<<"computeMinScaledJacobian2D "<<computeMinScaledJacobian2D()<<std::endl;
	// 
        std::map<TCellID, gmds::math::Point> nodesNewPos;
	std::map<TCellID, gmds::math::Point> nodesOldPos;
        std::map<TCellID, double> nodesNewPosWeight;

	// transformation of all elements
        {
	gmds::IGMesh::face_iterator itf  = this->m_mesh.faces_begin();
        for(;!itf.isDone();itf.next()) {
                gmds::Face current_cell = itf.value();
		
		// cellQuality should range between [0.,1.]
		double cellQuality = current_cell.computeScaledJacobian2D();
		cellQuality = (cellQuality + 1.)/2.;

                // compute a quality dependant factor
                double qualityDependantFactor = AQualityDependantFactorMin + (AQualityDependantFactorMax - AQualityDependantFactorMin) * (1. - cellQuality);

		// compute new Positions
		std::vector<gmds::Node> nodes = current_cell.get<gmds::Node>();
                std::vector<gmds::math::Point> newPos(nodes.size());
		computeGETMePoint(&current_cell,qualityDependantFactor,newPos);

		// apply scaling 
		// TODO

		// apply relaxation
		for(int iPoint=0; iPoint<nodes.size(); iPoint++) {
                        newPos[iPoint] = (1. - ARelaxFactor) * nodes[iPoint].getPoint() + ARelaxFactor*newPos[iPoint]; 
                }

		// compute weight
		double newPosWeight = std::pow((1. - cellQuality),AWeightExponent);	

		// store weight and position
		for(int iPoint=0; iPoint<nodes.size(); iPoint++) {
			if(nodesNewPosWeight.find(nodes[iPoint].getID()) != nodesNewPosWeight.end()) {
				nodesNewPosWeight[nodes[iPoint].getID()] += newPosWeight;
				nodesNewPos[nodes[iPoint].getID()] = nodesNewPos[nodes[iPoint].getID()] + newPosWeight*newPos[iPoint];
			} else {
				nodesNewPosWeight[nodes[iPoint].getID()] = newPosWeight;
				nodesNewPos[nodes[iPoint].getID()] = newPosWeight*newPos[iPoint];
			}
		}
	
        } // for(;!itr.isDone();itr.next())

	// divide by the weight
	gmds::IGMesh::node_iterator itn  = this->m_mesh.nodes_begin();
        for(;!itn.isDone();itn.next()) {
                gmds::Node current_node = itn.value();

		if(nodesNewPosWeight[current_node.getID()] != 0.) { 
	                nodesNewPos[current_node.getID()] = (1./nodesNewPosWeight[current_node.getID()]) * nodesNewPos[current_node.getID()];
		} else {
			nodesNewPos[current_node.getID()] = current_node.getPoint();
		}
         
        }
	
	}

	// apply geometric association if any
        {
	if(AIsGeomAssociated) {

		double geomAssocFactor = (numIter+1)/ANbMaxIter;

		if(AIsGeomAssocForced) {
			geomAssocFactor = 1;
		}

        	gmds::IGMesh::node_iterator itn  = this->m_mesh.nodes_begin();
	        for(;!itn.isDone();itn.next()) {
        	        gmds::Node current_node = itn.value();

			if(AMarkFixedNodes == -1 || (AMarkFixedNodes != -1 && !this->m_mesh.isMarked(current_node,AMarkFixedNodes))) {
					
				geom::GeomEntity* refGeometry = m_mesh.getGeometricClassification(current_node);
				if(refGeometry != NULL && refGeometry->getDim() != 3) {
		
					gmds::math::Point pt = nodesNewPos[current_node.getID()];
					refGeometry->project(pt);

					nodesNewPos[current_node.getID()] = (1. - geomAssocFactor) * nodesNewPos[current_node.getID()] + geomAssocFactor * pt;
				}
			}
		}	
        }
        }


	// update nodes positions
	{
	gmds::IGMesh::node_iterator itn  = this->m_mesh.nodes_begin();
        for(;!itn.isDone();itn.next()) {
                gmds::Node current_node = itn.value();

		nodesOldPos[current_node.getID()] = current_node.getPoint();

		if(AMarkFixedNodes == -1 || (AMarkFixedNodes != -1 && !this->m_mesh.isMarked(current_node,AMarkFixedNodes))) {
			  current_node.setPoint(nodesNewPos[current_node.getID()]);
		}
	}
	}

	// invalid element handling
        if (AAvoidTangledElements) {

        bool tangledCellDetected = false;
        do {
                tangledCellDetected = false;
                std::set<gmds::TCellID> nodes2RollBack;

                // detect invalid cells
                gmds::IGMesh::region_iterator itr  = this->m_mesh.regions_begin();
                for(;!itr.isDone();itr.next()) {
                        gmds::Region current_region = itr.value();

                        double scaledJacobian = current_region.computeScaledJacobian();
                        if(scaledJacobian <= 0.) {
                                tangledCellDetected = true;
                                std::vector<gmds::Node> nodes = current_region.get<gmds::Node>();
                                for(int iNode=0; iNode<nodes.size(); iNode++) {
                                        nodes2RollBack.insert(nodes[iNode].getID());
                                }
                        }
                }

                // rollback nodes
                std::set<gmds::TCellID>::iterator itn = nodes2RollBack.begin();
                for(; itn!=nodes2RollBack.end(); itn++) {
                        (this->m_mesh.get<gmds::Node>(*itn)).setPoint(nodesOldPos[(*itn)]);
                }

        } while(tangledCellDetected);

        } // if (AAvoidTangledElements)

//std::cout<<"computeMinScaledJacobian "<<computeMinScaledJacobian()<<std::endl;
	} // for(; numIter<ANbMaxIter; numIter++)

}
/*----------------------------------------------------------------------------*/
void
GETMe::execSeq(
		int ANbMaxIter,
                double AQualityThreshold,
                double ARelaxFactor,
                double AQualityDependantFactor,
                double ADeltaInvalid,
                double ADeltaReselect,
                double ADeltaSuccess,
                bool AIsGeomAssociated,
                bool AIsGeomAssocForced,
                bool AAvoidTangledElements,
                int AMarkFixedNodes)
{
	// check for parameters consistency
        if(AAvoidTangledElements) {
                if(computeMinNormalizedScaledJacobian() < 0.) {
                        throw GMDSException("GETMe::execSeq AAvoidTangledElements can not be set to true when min scaled jacobian is negative");
                }
        }
	if(ADeltaInvalid<=0.) {
	  throw GMDSException("GETMe::execSeq ADeltaInvalid must be > 0.");
	}
	if(ADeltaReselect<=0.) {
	  throw GMDSException("GETMe::execSeq ADeltaReselect must be > 0.");
	}
	if(ADeltaSuccess<=0.) {
	  throw GMDSException("GETMe::execSeq ADeltaSuccess must be > 0.");
	}

        // execution control variables
	const int nbIterPerRound = 100;

	// build ordered set of cells quality
	// we rely on the ordered property of the multimap
	// We also build the N2R adjacency
	std::multimap<double, gmds::Region> orderedByQualityRegions;
	std::map<gmds::Region, double> regions2qualityPenalty;
	std::map<gmds::Region, std::multimap<double, gmds::Region>::iterator> regions2Iterators;
	std::map<gmds::Node, std::set<gmds::Region> > node2Regions;
	gmds::TCellID previousRegionID = gmds::NullID;

	gmds::IGMesh::region_iterator itr  = this->m_mesh.regions_begin();
	for(;!itr.isDone();itr.next()) {

                gmds::Region current_region = itr.value();

		std::vector<gmds::Node> nodes = current_region.get<gmds::Node> ();
		for(size_t i=0; i<nodes.size(); i++) {
		  node2Regions[nodes[i]].insert(current_region);
		}
	}

        for(int numIter=0; (numIter<ANbMaxIter) && (computeMinNormalizedScaledJacobian()<AQualityThreshold); numIter++) {
	 
	  // reset quality penalty; this makes reordering the regions necessary
	  // Must be executed on the first pass
	  if((numIter % nbIterPerRound) == 0) {

	    std::cout<<"numIter "<<numIter<<" of "<<ANbMaxIter<<" "<<computeMinNormalizedScaledJacobian()<<std::endl;

	    orderedByQualityRegions.clear();
	    regions2qualityPenalty.clear();
	    regions2Iterators.clear();
	    gmds::IGMesh::region_iterator itr  = this->m_mesh.regions_begin();
	    for(;!itr.isDone();itr.next()) {
	      gmds::Region current_region = itr.value();

	      double scaledJacobian = current_region.computeNormalizedScaledJacobian();
	      std::multimap<double, gmds::Region>::iterator newIt = orderedByQualityRegions.insert(std::pair<double, gmds::Region> (scaledJacobian,current_region));
	      regions2Iterators.insert(std::pair<gmds::Region, std::multimap<double, gmds::Region>::iterator > (current_region, newIt));
	      regions2qualityPenalty.insert(std::pair<gmds::Region, double> (current_region, 0.));
	    }
	  }

	  std::multimap<double, gmds::Region>::iterator it = orderedByQualityRegions.begin();
	  	  
	  double quality = it->first;
	  gmds::Region current_region = it->second;
	  orderedByQualityRegions.erase(it);

	  // compute new positions
	  std::vector<gmds::Node> nodes = current_region.get<gmds::Node>();
	  std::vector<gmds::math::Point> newPos(nodes.size());
	  computeGETMePoint(&current_region,AQualityDependantFactor,newPos);
	  
	  // apply scaling
	  applyScaling(&current_region,newPos);

	  // apply relaxation
	  for(int iPoint=0; iPoint<nodes.size(); iPoint++) {
	    newPos[iPoint] = (1. - ARelaxFactor) * nodes[iPoint].getPoint() + ARelaxFactor*newPos[iPoint];
	  }

	  // apply geometric association if any
	  if(AIsGeomAssociated) {

	    double geomAssocFactor = (numIter+1)/ANbMaxIter;

	    if(AIsGeomAssocForced) {
	      geomAssocFactor = 1;
	    }
	    
	    for(size_t i=0; i<nodes.size(); i++) {
	      gmds::Node current_node = nodes[i];

	      if(AMarkFixedNodes == -1 || (AMarkFixedNodes != -1 && !this->m_mesh.isMarked(current_node,AMarkFixedNodes))) {
		geom::GeomEntity* refGeometry = m_mesh.getGeometricClassification(current_node);
		if(refGeometry != NULL && refGeometry->getDim() != 3) {
		  
		  gmds::math::Point pt = newPos[i];
		  refGeometry->project(pt);
		  newPos[i] = (1. - geomAssocFactor) * newPos[i] + geomAssocFactor * pt;
		}
	      }
	    }
	  }

	  // do not modify position if node is fixed
	  for(size_t i=0; i<nodes.size(); i++) {
	    if(AMarkFixedNodes != -1 && this->m_mesh.isMarked(nodes[i],AMarkFixedNodes)) {
	      newPos[i] = nodes[i].getPoint();
	    }
	  }
	  
	  // invalid element handling
	  bool elementsInverted = false;
	  if (AAvoidTangledElements) {
	    
	    std::map<gmds::Node, gmds::math::Point> nodes2Pts;
	    for(size_t i=0; i<nodes.size(); i++) {
	      nodes2Pts[nodes[i]] = nodes[i].getPoint();
	    }

	    std::set<gmds::Region> adjacentRegions;
	    for(size_t i=0; i<nodes.size(); i++) {
	      adjacentRegions.insert(node2Regions[nodes[i]].begin(),node2Regions[nodes[i]].end());
	    }

	    std::set<gmds::Region>::iterator itr = adjacentRegions.begin();
	    for(; itr != adjacentRegions.end(); itr++) {
	      
	      std::vector<gmds::Node> nodesTmp = itr->get<gmds::Node> ();
	      std::vector<gmds::math::Point> pts(nodesTmp.size());
	      for(size_t i=0; i<nodesTmp.size(); i++) {
		if(nodes2Pts.find(nodesTmp[i]) != nodes2Pts.end()) {
		  pts[i] = nodesTmp[i].getPoint();
		} else {
		  pts[i] = nodes2Pts[nodesTmp[i]];
		}
		
	      }

	      // check the would-be new cell quality
	      double scaledJacobian = 0.;
	      switch(itr->getType()) {
	      case GMDS_TETRA : 
		{
		  gmds::math::Tetrahedron tet (pts[0],pts[1],pts[2],pts[3]);
		  scaledJacobian = tet.computeScaledJacobian();
		}
		break;
	      case GMDS_HEX :
		{
		  gmds::math::Hexahedron hex (pts[0],pts[1],pts[2],pts[3],pts[4],pts[5],pts[6],pts[7]);
		  scaledJacobian = hex.computeScaledJacobian();
		}
		break;
	      case GMDS_PYRAMID :
		{
		  gmds::math::Pyramid pyr (pts[0],pts[1],pts[2],pts[3],pts[4]);
		  scaledJacobian = pyr.computeScaledJacobian();
		}
		break;
	      case GMDS_PRISM3 :
		{
		  gmds::math::Prism3 prism (pts[0],pts[1],pts[2],pts[3],pts[4],pts[5]);
		  scaledJacobian = prism.computeScaledJacobian();
		}
		break;
	      default :
		throw GMDSException("GETMe::execSeq could not roll-back a cell of unknown type.");
		break;
	      }
	
	      if(scaledJacobian <= 0.) {
		elementsInverted = true;
		break;
	      }
	    }

	  }

	  // update nodes positions                                                                                                                                            
	  {
	    if(!elementsInverted) {
	      for(size_t i=0; i<nodes.size(); i++) {
		
		gmds::Node current_node = nodes[i];
	
		if(AMarkFixedNodes == -1 || (AMarkFixedNodes != -1 && !this->m_mesh.isMarked(current_node,AMarkFixedNodes))) {
		  current_node.setPoint(newPos[i]);
		}
	      }
	    }
	    
	    // update the ordered container of regions
	    if(!elementsInverted) {
	      regions2qualityPenalty[current_region] -= ADeltaSuccess;
	    } else {
	      regions2qualityPenalty[current_region] += ADeltaInvalid;
	    }
	    if(previousRegionID == current_region.getID()) {
	      regions2qualityPenalty[current_region] += ADeltaReselect;
	    }
	    previousRegionID = current_region.getID();
	    
	    double newQualityPenalty = current_region.computeNormalizedScaledJacobian() + regions2qualityPenalty[current_region];
	    
	    regions2Iterators[current_region] = orderedByQualityRegions.insert(std::pair<double, gmds::Region> (newQualityPenalty,current_region));

	    // also update the adjacent regions
	    if(!elementsInverted) {
	      std::set<gmds::Region> adjacentRegions;
	      for(size_t i=0; i<nodes.size(); i++) {
		adjacentRegions.insert(node2Regions[nodes[i]].begin(),node2Regions[nodes[i]].end());
	      }

	      std::set<gmds::Region>::iterator itr = adjacentRegions.begin();
	      for(; itr != adjacentRegions.end(); itr++) {
		if(*itr != current_region) {
		  gmds::Region neighborRegion = *itr;
		  double newQualityPenalty = neighborRegion.computeNormalizedScaledJacobian() + regions2qualityPenalty[neighborRegion];
		  orderedByQualityRegions.erase(regions2Iterators[neighborRegion]);
		  regions2Iterators[neighborRegion] = orderedByQualityRegions.insert(std::pair<double, gmds::Region> (newQualityPenalty, neighborRegion));
		}
	      }
	    }

	  }
	  
	}

}
/*----------------------------------------------------------------------------*/
double
GETMe::computeMinMeanRatio()
{
	double minMeanRatio = HUGE_VALF;

	gmds::IGMesh::region_iterator itr  = this->m_mesh.regions_begin();
        for(;!itr.isDone();itr.next()) {
                gmds::Region current_region = itr.value();
	
		double meanRatio = current_region.computeMeanRatio();

		if(meanRatio < minMeanRatio) {
			minMeanRatio = meanRatio;
		}
	}

	return minMeanRatio;
}
/*----------------------------------------------------------------------------*/
double
GETMe::computeMinNormalizedScaledJacobian()
{
        double minValue = HUGE_VALF;

        gmds::IGMesh::region_iterator itr  = this->m_mesh.regions_begin();
        for(;!itr.isDone();itr.next()) {
                gmds::Region current_region = itr.value();

                double value = current_region.computeNormalizedScaledJacobian();

                if(value < minValue) {
                        minValue = value;
                }
        }

        return minValue;
}
/*----------------------------------------------------------------------------*/
double
GETMe::computeMinScaledJacobian2D()
{       
        double minValue = HUGE_VALF;
        
        gmds::IGMesh::face_iterator itf  = this->m_mesh.faces_begin();
        for(;!itf.isDone();itf.next()) {
                gmds::Face current_face = itf.value();
                
                double value = current_face.computeScaledJacobian2D();
                
                if(value < minValue) {
                        minValue = value;
                }
        }
        
        return minValue;
}
/*----------------------------------------------------------------------------*/
double
GETMe::computeVolumeTargetDiscrepancy()
{       
        if(!m_isVolumeTargetOn) {
    	        throw GMDSException("GETMe::computeVolumeTargetDiscrepancy can not be called when volume target is not active.");
        }

	gmds::Variable<double>* volumetarget = m_mesh.getVariable<double>(GMDS_REGION,m_volumeTargetVariableName);

        double maxValue = -HUGE_VALF;

        gmds::IGMesh::region_iterator itr  = this->m_mesh.regions_begin();
        for(;!itr.isDone();itr.next()) {
                gmds::Region current_region = itr.value();

                double volume = current_region.volume();
		double value = std::fabs(volume - (*volumetarget)[current_region.getID()]);

                if(value > maxValue) {
                        maxValue = value;
                }
        }

        return maxValue;
}
/*----------------------------------------------------------------------------*/
void
GETMe::computeGETMePoint(
                        const gmds::Cell* ACell,
                        double AScalingFactor, 
                        std::vector<gmds::math::Point>& APoints) const
{
	std::vector<gmds::Node> nodes = ACell->get<gmds::Node>();

	switch(ACell->getType()) {
		case GMDS_HEX :
			{
			gmds::math::Hexahedron hex(
				nodes[0].getPoint(),
	                        nodes[1].getPoint(),
        	                nodes[2].getPoint(),
                	        nodes[3].getPoint(),
                        	nodes[4].getPoint(),
	                        nodes[5].getPoint(),
        	                nodes[6].getPoint(),
                	        nodes[7].getPoint()
			);
			computeGETMePoint(hex,AScalingFactor,APoints);
			}
			break;
		case GMDS_TETRA :
			{
			gmds::math::Tetrahedron tet(
                                nodes[0].getPoint(),
                                nodes[1].getPoint(),
                                nodes[2].getPoint(),
                                nodes[3].getPoint()
			);
			computeGETMePoint(tet,AScalingFactor,APoints);
			}
			break;
		case GMDS_PYRAMID :
                        {
                        gmds::math::Pyramid pyr(
                                nodes[0].getPoint(),
                                nodes[1].getPoint(),
                                nodes[2].getPoint(),
                                nodes[3].getPoint(),
                                nodes[4].getPoint()
                        );
                        computeGETMePoint(pyr,AScalingFactor,APoints);
                        }
                        break;
		case GMDS_PRISM3 :
                        {
                        gmds::math::Prism3 prism(
                                nodes[0].getPoint(),
                                nodes[1].getPoint(),
                                nodes[2].getPoint(),
                                nodes[3].getPoint(),
				nodes[4].getPoint(),
                                nodes[5].getPoint()
                        );
                        computeGETMePoint(prism,AScalingFactor,APoints);
                        }
                        break;
		case GMDS_QUAD :
			{
			gmds::math::Quadrilateral quad(
				nodes[0].getPoint(),
                                nodes[1].getPoint(),
                                nodes[2].getPoint(),
                                nodes[3].getPoint()
			);
			computeGETMePoint(quad,AScalingFactor,APoints);	
			}
			break;		
		default :
			throw GMDSException("GETMe::computeGETMePoint not implemented for this cell type.");
			break;
	}
}
/*----------------------------------------------------------------------------*/
void
GETMe::computeGETMePoint(
			const gmds::math::Hexahedron AHex, 
			double AScalingFactor, 
			std::vector<gmds::math::Point>& APoints) const
{       
        // nodes of each face
        const int facesNodes[6][4] = {
                {0,3,2,1}, // bottom
                {0,1,5,4}, // front
                {1,2,6,5}, // right
                {2,3,7,6}, // back
                {3,0,4,7}, // left
                {4,5,6,7}  // top
        };
        
        // centroids of which faces is taken to build the "dual" triangle to the node
        // the triangle will be outward oriented
        const int node2faces[8][3] = {
                {0,1,4},
                {0,2,1},
                {0,3,2},
                {0,4,3},
                {4,1,5},
                {1,2,5},
                {2,3,5},
                {3,4,5}
        };

	gmds::math::Point pnts[8];
	for(int iPoint=0; iPoint<8; iPoint++) {
		pnts[iPoint] = AHex.getPoint(iPoint);
	}

	for(int iPoint=0; iPoint<8; iPoint++) {

        	// first build the "dual" triangle
	        gmds::math::Point pt0 = pnts[facesNodes[node2faces[iPoint][0]][0]] + pnts[facesNodes[node2faces[iPoint][0]][1]] +
        	                        pnts[facesNodes[node2faces[iPoint][0]][2]] + pnts[facesNodes[node2faces[iPoint][0]][3]];
	        pt0 = 1./4. * pt0;
        	gmds::math::Point pt1 = pnts[facesNodes[node2faces[iPoint][1]][0]] + pnts[facesNodes[node2faces[iPoint][1]][1]] +
                	                pnts[facesNodes[node2faces[iPoint][1]][2]] + pnts[facesNodes[node2faces[iPoint][1]][3]];
	        pt1 = 1./4. * pt1;
        	gmds::math::Point pt2 = pnts[facesNodes[node2faces[iPoint][2]][0]] + pnts[facesNodes[node2faces[iPoint][2]][1]] +
                	                pnts[facesNodes[node2faces[iPoint][2]][2]] + pnts[facesNodes[node2faces[iPoint][2]][3]];
	        pt2 = 1./4. * pt2;
        
        	gmds::math::Triangle tri(pt0,pt1,pt2);
	        gmds::math::Point center = tri.getCenter();
        
        	gmds::math::Vector normal = tri.getNormal();
        
	        gmds::math::Point newPos = center + (AScalingFactor / sqrt(normal.norm())) * normal;

		APoints[iPoint] = newPos;
	}
}
/*----------------------------------------------------------------------------*/
void
GETMe::computeGETMePoint(
                        const gmds::math::Tetrahedron ATet,
                        double AScalingFactor,
                        std::vector<gmds::math::Point>& APoints) const
{
        // nodes of each face
        const int facesNodes[4][3] = {
                {0,2,1}, // bottom
                {0,1,3}, // front
                {1,2,3}, // right
                {2,0,3}  // back
        };

        // centroids of which faces is taken to build the "dual" triangle to the node
        // the triangle will be outward oriented
        const int node2faces[4][3] = {
                {0,1,3},
                {0,2,1},
                {0,3,2},
                {1,2,3}
        };

        gmds::math::Point pnts[4];
        for(int iPoint=0; iPoint<4; iPoint++) {
                pnts[iPoint] = ATet.getPoint(iPoint);
        }

        for(int iPoint=0; iPoint<4; iPoint++) {

                // first build the "dual" triangle
                gmds::math::Point pt0 = pnts[facesNodes[node2faces[iPoint][0]][0]] + pnts[facesNodes[node2faces[iPoint][0]][1]] +
                                        pnts[facesNodes[node2faces[iPoint][0]][2]];
                pt0 = 1./3. * pt0;
                gmds::math::Point pt1 = pnts[facesNodes[node2faces[iPoint][1]][0]] + pnts[facesNodes[node2faces[iPoint][1]][1]] +
                                        pnts[facesNodes[node2faces[iPoint][1]][2]];
                pt1 = 1./3. * pt1;
                gmds::math::Point pt2 = pnts[facesNodes[node2faces[iPoint][2]][0]] + pnts[facesNodes[node2faces[iPoint][2]][1]] +
                                        pnts[facesNodes[node2faces[iPoint][2]][2]];
                pt2 = 1./3. * pt2;

                gmds::math::Triangle tri(pt0,pt1,pt2);
                gmds::math::Point center = tri.getCenter();

                gmds::math::Vector normal = tri.getNormal();

                gmds::math::Point newPos = center + (AScalingFactor / sqrt(normal.norm())) * normal;

                APoints[iPoint] = newPos;
        }

}
/*----------------------------------------------------------------------------*/
void
GETMe::computeGETMePoint(
                        const gmds::math::Pyramid APyr,
                        double AScalingFactor,
                        std::vector<gmds::math::Point>& APoints) const
{
        // nodes of each face
        const int facesNodes[5][4] = {
                {0,3,2, 1}, // bottom
                {0,1,4,-1}, // front
                {1,2,4,-1}, // right
                {2,3,4,-1}, // back
		{3,0,4,-1}  // left
        };

	const int facesNodesSize[5] = {
		4,  // bottom
		3,  // front
		3,  // right
		3,  // back
		3   // left
	};

        // centroids of which faces is taken to build the "dual" triangle to the node
        // the triangle will be outward oriented
        const int node2faces[5][4] = {
                {0,1,4,-1},
                {0,2,1,-1},
                {0,3,2,-1},
                {0,4,3,-1},
		{1,2,3, 4}
        };

	const int node2facesSize[5] = {
		3,
		3,
		3,
		3,
		4
	};

	const int nbPoints = 5;

        gmds::math::Point pnts[nbPoints];
        for(int iPoint=0; iPoint<nbPoints; iPoint++) {
                pnts[iPoint] = APyr.getPoint(iPoint);
        }

        for(int iPoint=0; iPoint<nbPoints; iPoint++) {

		std::vector<gmds::math::Point> pnts_tmp(node2facesSize[iPoint]);

		gmds::math::Point center;
		for(int i=0; i<node2facesSize[iPoint]; i++) {
			for(int j=0; j<facesNodesSize[node2faces[iPoint][i]]; j++) {
				pnts_tmp[i] = pnts_tmp[i] + pnts[facesNodes[node2faces[iPoint][i]][j]];
			}

			pnts_tmp[i] = 1./facesNodesSize[node2faces[iPoint][i]] * pnts_tmp[i];
			center = center + pnts_tmp[i];
		}
		center = 1./node2facesSize[iPoint] * center;

		gmds::math::Vector normal;
		if(node2facesSize[iPoint] == 3) {
			normal = gmds::math::Vector(pnts_tmp[1] - pnts_tmp[0]).cross(gmds::math::Vector(pnts_tmp[2] - pnts_tmp[0]));
		} else {
			normal = 1./2. * gmds::math::Vector(pnts_tmp[2] - pnts_tmp[0]).cross(gmds::math::Vector(pnts_tmp[3] - pnts_tmp[1]));
		}

                gmds::math::Point newPos = center + (AScalingFactor / sqrt(normal.norm())) * normal;

                APoints[iPoint] = newPos;
		//std::cout<<"poyop "<<pnts[iPoint]<<" "<<newPos<<std::endl;
        }

}
/*----------------------------------------------------------------------------*/
void
GETMe::computeGETMePoint(
                        const gmds::math::Prism3 APrism,
                        double AScalingFactor,
                        std::vector<gmds::math::Point>& APoints) const
{
        // nodes of each face
        const int facesNodes[5][4] = {
                {2,0,1,-1}, // bottom
                {0,1,4, 3}, // front
                {1,2,5, 4}, // right
                {2,1,3, 5}, // back
                {3,5,4,-1}  // top
        };

        const int facesNodesSize[5] = {
                3,  // bottom
                4,  // front
                4,  // right
                4,  // back
                3   // top
        };

        // centroids of which faces is taken to build the "dual" triangle to the node
        // the triangle will be outward oriented
        const int node2faces[6][3] = {
                {0,1,3},
                {0,2,1},
                {0,3,2},
                {4,3,1},
                {4,1,2},
		{4,2,3}
        };

	const int node2facesSize[6] = {
		3,3,3,3,3,3
	};

        const int nbPoints = 6;

        gmds::math::Point pnts[nbPoints];
        for(int iPoint=0; iPoint<nbPoints; iPoint++) {
                pnts[iPoint] = APrism.getPoint(iPoint);
        }

        for(int iPoint=0; iPoint<nbPoints; iPoint++) {

                std::vector<gmds::math::Point> pnts_tmp(node2facesSize[iPoint]);

                gmds::math::Point center;
                for(int i=0; i<node2facesSize[iPoint]; i++) {
                        for(int j=0; j<facesNodesSize[node2faces[iPoint][i]]; j++) {
                                pnts_tmp[i] = pnts_tmp[i] + pnts[facesNodes[node2faces[iPoint][i]][j]];
                        }

                        pnts_tmp[i] = 1./facesNodesSize[node2faces[iPoint][i]] * pnts_tmp[i];
                        center = center + pnts_tmp[i];
                }
		center = 1./node2facesSize[iPoint] * center;

                gmds::math::Vector normal;
                if(node2facesSize[iPoint] == 3) {
                        normal = gmds::math::Vector(pnts_tmp[1] - pnts_tmp[0]).cross(gmds::math::Vector(pnts_tmp[2] - pnts_tmp[0]));
                } else {
                        normal = 1./2. * gmds::math::Vector(pnts_tmp[2] - pnts_tmp[0]).cross(gmds::math::Vector(pnts_tmp[3] - pnts_tmp[1]));
                }

                gmds::math::Point newPos = center + (AScalingFactor / sqrt(normal.norm())) * normal;

                APoints[iPoint] = newPos;
        }

}
/*----------------------------------------------------------------------------*/
void
GETMe::computeGETMePoint(
                        const gmds::math::Quadrilateral AQuad,
                        double AScalingFactor, 
                        std::vector<gmds::math::Point>& APoints) const
{
        // nodes of each side
        int facesNodes[4][2] = {
                {0,1}, // bottom
                {1,2}, // right
                {2,3}, // top
                {3,0}  // left
        };

        // centroids of which faces is taken to build the "dual" triangle to the node
        // the triangle will be outward oriented
        int node2faces[4][2] = {
                {0,3},
                {1,0},
                {2,1},
                {3,2}
        };

	gmds::math::Point pnts[4];
        for(int iPoint=0; iPoint<4; iPoint++) {
                pnts[iPoint] = AQuad.getPoint(iPoint);
        }

        for(int iPoint=0; iPoint<4; iPoint++) {

                // first build the "dual" segment
                gmds::math::Point pt0 = pnts[facesNodes[node2faces[iPoint][0]][0]] + pnts[facesNodes[node2faces[iPoint][0]][1]];
                pt0 = 1./2. * pt0;
                gmds::math::Point pt1 = pnts[facesNodes[node2faces[iPoint][1]][0]] + pnts[facesNodes[node2faces[iPoint][1]][1]];
                pt1 = 1./2. * pt1;

                gmds::math::Segment seg(pt0,pt1);
                gmds::math::Point center = seg.computeCenter();
			
		// we are in 2D, the "normal" here is simply the normal to the segment counter-clockwise 
		gmds::math::Vector vec = seg.getUnitVector();
                gmds::math::Vector normal(-vec.Y(),vec.X());

                gmds::math::Point newPos = center + (AScalingFactor / sqrt(normal.norm())) * normal;

                APoints[iPoint] = newPos;
        }
}
/*----------------------------------------------------------------------------*/
void
GETMe::applyScaling(
		    const gmds::Cell* ACell,
		    std::vector<gmds::math::Point>& APoints,
		    double AModifySizeFactor) const
{  
        gmds::math::Point center = ACell->center();

	std::vector<gmds::Node> nodes = ACell->get<gmds::Node>();

	double factor_old;

	switch(ACell->getType()) {
		case GMDS_HEX :
			{
			gmds::math::Hexahedron hex(
				nodes[0].getPoint(),
	                        nodes[1].getPoint(),
        	                nodes[2].getPoint(),
                	        nodes[3].getPoint(),
                        	nodes[4].getPoint(),
	                        nodes[5].getPoint(),
        	                nodes[6].getPoint(),
                	        nodes[7].getPoint()
			);
			factor_old = hex.computeMeanEdgeLength();
			}
			break;
		case GMDS_TETRA :
			{
			gmds::math::Tetrahedron tet(
                                nodes[0].getPoint(),
                                nodes[1].getPoint(),
                                nodes[2].getPoint(),
                                nodes[3].getPoint()
			);
			factor_old = tet.computeMeanEdgeLength();
			}
			break;
		case GMDS_PYRAMID :
                        {
                        gmds::math::Pyramid pyr(
                                nodes[0].getPoint(),
                                nodes[1].getPoint(),
                                nodes[2].getPoint(),
                                nodes[3].getPoint(),
                                nodes[4].getPoint()
                        );
			factor_old = pyr.computeMeanEdgeLength();
                        }
                        break;
		case GMDS_PRISM3 :
                        {
                        gmds::math::Prism3 prism(
                                nodes[0].getPoint(),
                                nodes[1].getPoint(),
                                nodes[2].getPoint(),
                                nodes[3].getPoint(),
				nodes[4].getPoint(),
                                nodes[5].getPoint()
                        );
			factor_old = prism.computeMeanEdgeLength();
                        }
                        break;
		case GMDS_TRIANGLE :
			{
			gmds::math::Triangle tri(
				nodes[0].getPoint(),
                                nodes[1].getPoint(),
                                nodes[2].getPoint()
			);
			factor_old = tri.computeMeanEdgeLength();
			}
			break;		
		case GMDS_QUAD :
			{
			gmds::math::Quadrilateral quad(
				nodes[0].getPoint(),
                                nodes[1].getPoint(),
                                nodes[2].getPoint(),
                                nodes[3].getPoint()
			);
			factor_old = quad.computeMeanEdgeLength();
			}
			break;		
		default :
			throw GMDSException("GETMe::applyScaling not implemented for this cell type.");
			break;
	}

	double factor_new;

	switch(ACell->getType()) {
		case GMDS_HEX :
			{
			  gmds::math::Hexahedron hex(
				APoints[0],
				APoints[1],
				APoints[2],
				APoints[3],
				APoints[4],
				APoints[5],
				APoints[6],
				APoints[7]);
			  
			  factor_new = hex.computeMeanEdgeLength();
			}
			break;
		case GMDS_TETRA :
			{
			gmds::math::Tetrahedron tet(
				APoints[0],
				APoints[1],
				APoints[2],
				APoints[3]);

			  
			  factor_new = tet.computeMeanEdgeLength();
			}
			break;
		case GMDS_PYRAMID :
                        {
			gmds::math::Pyramid pyr(
				APoints[0],
				APoints[1],
				APoints[2],
				APoints[3],
				APoints[4]);

			factor_new = pyr.computeMeanEdgeLength();
                        }
                        break;
		case GMDS_PRISM3 :
                        {
			gmds::math::Prism3 prism(
				APoints[0],
				APoints[1],
				APoints[2],
				APoints[3],
				APoints[4],
				APoints[5]);

			factor_new = prism.computeMeanEdgeLength();
                        }
                        break;
	        case GMDS_TRIANGLE :
			{
			gmds::math::Triangle tri(
				APoints[0],
				APoints[1],
				APoints[2]);
			  
			  factor_new = tri.computeMeanEdgeLength();
			}
			break;		
		case GMDS_QUAD :
			{
			gmds::math::Quadrilateral quad(
				APoints[0],
				APoints[1],
				APoints[2],
				APoints[3]);
			  
			  factor_new = quad.computeMeanEdgeLength();
			}
			break;		
		default :
			throw GMDSException("GETMe::applyScaling not implemented for this cell type.");
			break;
	}

	// the factor used is the mean edge length, so a simple ratio enables the size preservation
	double factor = factor_old/factor_new;
	factor *= AModifySizeFactor;
	
	for(int i=0; i<APoints.size(); i++) {
	  APoints[i] = center + factor * (APoints[i] - center);
	}
}
/*----------------------------------------------------------------------------*/
