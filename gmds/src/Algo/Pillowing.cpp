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
/** \file    Pillowing.cpp
 *  \author  legoff
 *  \date    22/07/2016
 */
/*----------------------------------------------------------------------------*/
#include "GMDS/Algo/Pillowing.h"
/*----------------------------------------------------------------------------*/
#include <sstream>
#include <set>
#include <map>
#include <vector>
/*----------------------------------------------------------------------------*/
#include "GMDS/CAD/GeomEntity.h"
#include "GMDS/Math/Triangle.h"
/*----------------------------------------------------------------------------*/
namespace gmds {
/*----------------------------------------------------------------------------*/
Pillowing::Pillowing(IGMesh& m_mesh)
:m_mesh(m_mesh)
 {}
/*----------------------------------------------------------------------------*/
Pillowing::~Pillowing()
{}
/*----------------------------------------------------------------------------*/
void
Pillowing::pillow(
	const std::vector<gmds::Face>& AFaces,
	const std::vector<gmds::Region>& ARegions,
	std::vector<gmds::Region>& ANewRegions,
	bool ADisplacement)
{
	ANewRegions.clear();

	m_displacement = ADisplacement;

	// check whether the model meets the requirements for this algorithm
	if(!this->checkModel()) {
		throw GMDSException("Pillowing::pillow can not run on this mesh model.");
	}
	
	// create the fake faces and identify those on the boundary of the shrink sets
	this->retrieveBoundary(ARegions);

	// differentiate between the disjoined shrink sets
	this->identifyShrinkSets(ARegions);

	// create the new nodes
	this->createNewNodes();

	// create the new cells
	this->createNewCells(ANewRegions);
	
	// set the new nodes for the affected regions
	this->setNewNodes();

}
/*----------------------------------------------------------------------------*/
bool
Pillowing::checkModel()
{
	gmds::MeshModel mod = m_mesh.getModel();

	if(
	!mod.has(gmds::R) ||
	!mod.has(gmds::N) ||
	!mod.has(gmds::R2N)) {
		return false;
	}

	return true;
}
/*----------------------------------------------------------------------------*/
void
Pillowing::identifyShrinkSets(const std::vector<gmds::Region>& ARegions)
{
	std::set<gmds::Region> unTreatedRegions;
	std::set<gmds::Region> treatedRegions;

	for(size_t i=0; i<ARegions.size(); i++) {
		unTreatedRegions.insert(ARegions[i]);
	}

	uint current_set = 0;

	while (!unTreatedRegions.empty()) {
		// find an untreated cell
		gmds::Region current_region = *(unTreatedRegions.begin());
		this->propagateRecurse(current_region,current_set,unTreatedRegions,treatedRegions);
		current_set++;
	}
}
/*----------------------------------------------------------------------------*/
void
Pillowing::propagateRecurse(gmds::Region ARegion, uint ASetIndex, std::set<gmds::Region>& AUnTreatedRegions, std::set<gmds::Region>& ATreatedRegions)
{       
        uint current_set = 0;

	if(AUnTreatedRegions.find(ARegion) != AUnTreatedRegions.end()) {
		// mark this region
		m_regionsSets[ARegion] = ASetIndex;
		AUnTreatedRegions.erase(ARegion);
		ATreatedRegions.insert(ARegion);

		// propagate to neighbors
		for(size_t i=0; i<m_regions2FakeFaces[ARegion].size(); i++) {
			if(m_fakeFaces2Regions[m_regions2FakeFaces[ARegion][i]].size() == 2) {
				if(m_fakeFaces2Regions[m_regions2FakeFaces[ARegion][i]][0] != ARegion) {
					propagateRecurse(m_fakeFaces2Regions[m_regions2FakeFaces[ARegion][i]][0],ASetIndex,AUnTreatedRegions,ATreatedRegions);
				} else {
					propagateRecurse(m_fakeFaces2Regions[m_regions2FakeFaces[ARegion][i]][1],ASetIndex,AUnTreatedRegions,ATreatedRegions);
				}
			}
		}
	}
}
/*----------------------------------------------------------------------------*/
void
Pillowing::retrieveBoundary(const std::vector<gmds::Region>& ARegions)
{
	// build the fake faces
        for(size_t iRegion=0; iRegion<ARegions.size(); iRegion++) {

                gmds::Region current_region = ARegions[iRegion];
                std::vector<gmds::Node> nodes = current_region.get<gmds::Node>();

                std::vector<gmds::FakeFace> fakeFaces = current_region.getFakeFaces();

                for(size_t iFace=0; iFace<fakeFaces.size(); iFace++) {
                        m_fakeFaces2Regions[fakeFaces[iFace]].push_back(current_region);
                }
                m_regions2FakeFaces[current_region] = fakeFaces;
        }

	// get the boundary fake faces and 
	// determines whether the fake face is ordered outwardly from the cell
        {
                std::map<gmds::FakeFace, std::vector<gmds::Region> >::iterator it = m_fakeFaces2Regions.begin();
                for(; it!=m_fakeFaces2Regions.end(); it++) {
                        if(it->second.size() == 1) {
                                m_pillowedFakeFaces[it->first] = it->second[0];

                                std::vector<TCellID> ids = it->first.node_ids();
                                bool isOutward = it->second[0].isFaceOrientedOutward(it->first.node_ids());
                                m_pillowedFakeFacesIsOutward[it->first] = isOutward;
                        }
                }
        }
	// mark the nodes of the boundary of the shrink set
        {
                std::map<gmds::FakeFace, gmds::Region>::iterator it = m_pillowedFakeFaces.begin();
                for(; it!=m_pillowedFakeFaces.end(); it++) {
                        std::vector<TCellID> ids = it->first.node_ids();

                        for(size_t i=0; i<ids.size(); i++) {
                                m_node2FakeFaces[ids[i]].push_back(it->first);
                        }
                }
        }

}
/*----------------------------------------------------------------------------*/
void
Pillowing::createNewNodes()
	{
		std::map<gmds::TCellID, std::vector<gmds::FakeFace> >::iterator it = m_node2FakeFaces.begin();
                for(; it!=m_node2FakeFaces.end(); it++) {

			gmds::Node oldNode = m_mesh.get<gmds::Node>(it->first);

			gmds::TCellID oldNodeID = it->first;

				// use one vector per shrink set this node is a part of
				std::map<uint, gmds::math::Vector> v;
				std::map<uint, double> length;

				for(size_t iFace=0; iFace<it->second.size(); iFace++) {

					// get the index of the shrink set
					uint shrinkSet = m_regionsSets[m_pillowedFakeFaces[it->second[iFace]]];
					
					if(v.find(shrinkSet) == v.end()) {
						v[shrinkSet] = gmds::math::Vector (0,0,0);
						length[shrinkSet] = HUGE_VALF;
					}


					// locate the node in the face and build the triangle
					// in order to compute the normal
					std::vector<gmds::TCellID> ids = it->second[iFace].node_ids();

					int index = 0;
					for(; index<ids.size(); index++) {
						if(oldNodeID == ids[index]) {
							break;
						}
					}

					int indexNext = (index == ids.size()-1)? 0 : index+1;
					int indexPrev = (index == 0)? ids.size()-1 : index-1;

					gmds::math::Point pt0 = m_mesh.get<gmds::Node>(ids[indexPrev]).getPoint();
					gmds::math::Point pt1 = m_mesh.get<gmds::Node>(ids[index]).getPoint();
					gmds::math::Point pt2 = m_mesh.get<gmds::Node>(ids[indexNext]).getPoint();

					gmds::math::Triangle tri (pt0,pt1,pt2);

					gmds::math::Vector normal = tri.getNormal();
					normal.normalize();

					if(!m_pillowedFakeFacesIsOutward[it->second[iFace]]) {
						normal = (-1) * normal;
					}

					v[shrinkSet] = v[shrinkSet] + normal;

					if(length[shrinkSet] > pt1.distance(pt0)) {
						length[shrinkSet] = pt1.distance(pt0);
					}
					if(length[shrinkSet] > pt1.distance(pt2)) {
                                                length[shrinkSet] = pt1.distance(pt2);
                                        }
				}

				std::map<uint, gmds::math::Vector>::iterator it = v.begin();
				for(; it != v.end(); it++) {

					uint shrinkSet = it->first;
					gmds::math::Point newPos(oldNode.getPoint());
					if(m_displacement) {
					  newPos = newPos - (1./10.) * length[shrinkSet] * v[shrinkSet];
					} 
					gmds::Node newNode = m_mesh.newNode(newPos);

					m_nodesOld2New[shrinkSet][oldNode] = newNode;
					m_newNodes.insert(newNode);
				}
			
		}
	}
/*----------------------------------------------------------------------------*/
void
Pillowing::createNewCells(std::vector<gmds::Region>& ANewRegions)
{
		std::map<gmds::FakeFace, gmds::Region>::iterator it = m_pillowedFakeFaces.begin();
                for(; it!=m_pillowedFakeFaces.end(); it++) {
			std::vector<TCellID> ids = it->first.node_ids();

			uint shrinkSet = m_regionsSets[it->second];

			if(ids.size() == 3) {

				gmds::Node n0 = m_mesh.get<gmds::Node>(ids[0]);
                                gmds::Node n1 = m_mesh.get<gmds::Node>(ids[1]);
                                gmds::Node n2 = m_mesh.get<gmds::Node>(ids[2]);
                                gmds::Node n3 = m_nodesOld2New[shrinkSet][n0];
                                gmds::Node n4 = m_nodesOld2New[shrinkSet][n1];
                                gmds::Node n5 = m_nodesOld2New[shrinkSet][n2];

                                if(!m_pillowedFakeFacesIsOutward[it->first]) {
                                        gmds::Region region = m_mesh.newPrism3(n0,n1,n2,n3,n4,n5);
                                        ANewRegions.push_back(region);
                                } else {
                                        gmds::Region region = m_mesh.newPrism3(n3,n4,n5,n0,n1,n2);
                                        ANewRegions.push_back(region);
                                }
			} else if(ids.size() == 4) {

				gmds::Node n0 = m_mesh.get<gmds::Node>(ids[0]);
				gmds::Node n1 = m_mesh.get<gmds::Node>(ids[1]);
				gmds::Node n2 = m_mesh.get<gmds::Node>(ids[2]);
				gmds::Node n3 = m_mesh.get<gmds::Node>(ids[3]);
				gmds::Node n4 = m_nodesOld2New[shrinkSet][n0];
				gmds::Node n5 = m_nodesOld2New[shrinkSet][n1];
				gmds::Node n6 = m_nodesOld2New[shrinkSet][n2];
				gmds::Node n7 = m_nodesOld2New[shrinkSet][n3];
			
				if(!m_pillowedFakeFacesIsOutward[it->first]) {
					gmds::Region region = m_mesh.newHex(n0,n1,n2,n3,n4,n5,n6,n7);
					ANewRegions.push_back(region);
				} else {
					gmds::Region region = m_mesh.newHex(n4,n5,n6,n7,n0,n1,n2,n3);
                                	ANewRegions.push_back(region);
				}
			} else {
				throw GMDSException("Pillowing::pillow not implemented for fake faces other than triangles and quads.");
			}
		}
}
/*----------------------------------------------------------------------------*/
void
Pillowing::setNewNodes()
{
                std::map<gmds::FakeFace, gmds::Region>::iterator it = m_pillowedFakeFaces.begin();
                for(; it!=m_pillowedFakeFaces.end(); it++) {
                        gmds::Region current_region = it->second;
                        uint shrinkSet = m_regionsSets[current_region];

                        std::vector<gmds::Node> nodes = current_region.get<gmds::Node>();
                        for(size_t i=0; i<nodes.size(); i++) {
                                if(m_nodesOld2New.find(shrinkSet) != m_nodesOld2New.end()) {
                                        if(m_nodesOld2New[shrinkSet].find(nodes[i]) != m_nodesOld2New[shrinkSet].end()) {
                                                nodes[i] = m_nodesOld2New[shrinkSet][nodes[i]];
                                        }
                                }
                        }
                        current_region.set<gmds::Node>(nodes);
                }
}
/*----------------------------------------------------------------------------*/
std::set<Node>
Pillowing::getNewInternalNodes()
{
	return m_newNodes;
}
/*----------------------------------------------------------------------------*/
}  // namespace gmds
/*----------------------------------------------------------------------------*/
