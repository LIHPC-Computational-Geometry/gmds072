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
 * Edge.cpp
 *
 *  Created on: 19 may 2014
 *      Author: ledouxf
 */
/*----------------------------------------------------------------------------*/
#include <GMDS/IG/Edge.h>
#include <GMDS/IG/IGMesh.h>
#include <GMDS/IG/EdgeContainer.h>
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/
Edge::
Edge()
: Cell(0,GMDS_EDGE,NullID)
{
	m_edges_container = 0;
}
/*----------------------------------------------------------------------------*/
Edge::
Edge(IGMesh* AMesh, const TCellID& AID)
: Cell(AMesh,GMDS_EDGE,AID)
{
	if(AMesh!=0)
		m_edges_container  = AMesh->m_edges_container;
	else
		m_edges_container  = 0;
}
/*----------------------------------------------------------------------------*/
Edge::
Edge(const Edge& AEdge)
: Cell(AEdge.m_owner,GMDS_EDGE,AEdge.m_id)
{
	if(m_owner!=0)
		m_edges_container = m_owner->m_edges_container;
	else
		m_edges_container = 0;
}
/*----------------------------------------------------------------------------*/
Edge::~Edge()
{}
/*----------------------------------------------------------------------------*/
bool Edge::operator==(const Edge& AEdge) const
{
        return (m_owner == AEdge.m_owner && m_id==AEdge.m_id);
}
/*----------------------------------------------------------------------------*/
bool Edge::operator!=(const Edge& AEdge) const
{
        return (!(*this == AEdge));
}
/*----------------------------------------------------------------------------*/
TInt Edge::getNbNodes() const
{
	return 2;
}
/*----------------------------------------------------------------------------*/
TInt Edge::getNbEdges() const
{
	TabCellID<size_undef> cells = (*m_edges_container->m_E2E)[m_id];
	return cells.size();
}
/*----------------------------------------------------------------------------*/
TInt Edge::getNbFaces() const
{
	TabCellID<size_undef> cells = (*m_edges_container->m_E2F)[m_id];
	return cells.size();
}
/*----------------------------------------------------------------------------*/
TInt Edge::getNbRegions() const
{
	TabCellID<size_undef> cells = (*m_edges_container->m_E2R)[m_id];
	return cells.size();
}
/*----------------------------------------------------------------------------*/
TCoord Edge::length() const
{
	std::vector<Node> nodes;
	get<Node>(nodes);
	return math::Vector(nodes[0].getPoint(),nodes[1].getPoint()).norm();
}
/*----------------------------------------------------------------------------*/
math::Point Edge::center() const
{
	std::vector<Node> nodes;
	get<Node>(nodes);
	return 0.5*nodes[0].getPoint()+0.5*nodes[1].getPoint();
}
/*----------------------------------------------------------------------------*/
void Edge::delegateGet(std::vector<Node>&   ACells) const
{
	if(!m_owner->m_model.has(E2N))
		throw GMDSException("E2N adjacency is not supported by the mesh model");

	std::vector<TCellID> cellIDs;
	(*m_edges_container->m_E2N)[m_id].values(cellIDs);
	ACells.resize(cellIDs.size());
	for(unsigned int i=0;i<cellIDs.size();i++){
		ACells[i] =m_owner->m_nodes_container->buildNode(cellIDs[i]);
	}
}
/*----------------------------------------------------------------------------*/
void Edge::delegateGet(std::vector<Edge>&   ACells) const
{
	if(!m_owner->m_model.has(E2E))
		throw GMDSException("E2E adjacency is not supported by the mesh model");
	std::vector<TCellID> cellIDs;
	(*m_edges_container->m_E2E)[m_id].values(cellIDs);
	ACells.resize(cellIDs.size());
	for(unsigned int i=0;i<cellIDs.size();i++){
		ACells[i] =m_edges_container->buildEdge(cellIDs[i]);
	}
}
/*----------------------------------------------------------------------------*/
void Edge::delegateGet(std::vector<Face>&   ACells) const
{
	if(!m_owner->m_model.has(E2F))
		throw GMDSException("E2F adjacency is not supported by the mesh model");
	std::vector<TCellID> cellIDs;
	(*m_edges_container->m_E2F)[m_id].values(cellIDs);
	ACells.resize(cellIDs.size());
	for(unsigned int i=0;i<cellIDs.size();i++){
		ACells[i] =m_owner->m_faces_container->buildFace(cellIDs[i]);
	}
}
/*----------------------------------------------------------------------------*/
void Edge::delegateGet(std::vector<Region>& ACells) const
{
	if(!m_owner->m_model.has(E2R))
		throw GMDSException("E2R adjacency is not supported by the mesh model");
	std::vector<TCellID> cellIDs;
	(*m_edges_container->m_E2R)[m_id].values(cellIDs);
	ACells.resize(cellIDs.size());
	for(unsigned int i=0;i<cellIDs.size();i++){
		ACells[i] =m_owner->m_regions_container->buildRegion(cellIDs[i]);
	}
}
/*----------------------------------------------------------------------------*/
void Edge::delegateGetNodeIDs(std::vector<TCellID>&   ACells) const
{
	if(!m_owner->m_model.has(E2N))
		throw GMDSException("E2N adjacency is not supported by the mesh model");
	(*m_edges_container->m_E2N)[m_id].values(ACells);
}
/*----------------------------------------------------------------------------*/
void Edge::delegateGetEdgeIDs(std::vector<TCellID>&   ACells) const
{
	if(!m_owner->m_model.has(E2E))
		throw GMDSException("E2E adjacency is not supported by the mesh model");
	(*m_edges_container->m_E2E)[m_id].values(ACells);
}
/*----------------------------------------------------------------------------*/
void Edge::delegateGetFaceIDs(std::vector<TCellID>&   ACells) const
{
	if(!m_owner->m_model.has(E2F))
		throw GMDSException("E2F adjacency is not supported by the mesh model");
	(*m_edges_container->m_E2F)[m_id].values(ACells);
}
/*----------------------------------------------------------------------------*/
void Edge::delegateGetRegionIDs(std::vector<TCellID>&   ACells) const
{
	if(!m_owner->m_model.has(E2R))
		throw GMDSException("E2R adjacency is not supported by the mesh model");
	(*m_edges_container->m_E2R)[m_id].values(ACells);
}
/*----------------------------------------------------------------------------*/
void Edge::delegateGetAll(std::vector<Node>&   ACells) const
{
	if(!m_owner->m_model.has(E2N))
		throw GMDSException("E2N adjacency is not supported by the mesh model");
	std::vector<TCellID> cellIDs;
	(*m_edges_container->m_E2N)[m_id].allValues(cellIDs);
	ACells.resize(cellIDs.size());
	for(unsigned int i=0;i<cellIDs.size();i++){
		ACells[i] =m_owner->m_nodes_container->buildNode(cellIDs[i]);
	}
}
/*----------------------------------------------------------------------------*/
void Edge::delegateGetAll(std::vector<Edge>&   ACells) const
{
	if(!m_owner->m_model.has(E2E))
		throw GMDSException("E2E adjacency is not supported by the mesh model");
	std::vector<TCellID> cellIDs;
	(*m_edges_container->m_E2N)[m_id].allValues(cellIDs);
	ACells.resize(cellIDs.size());
	for(unsigned int i=0;i<cellIDs.size();i++){
		ACells[i] =m_edges_container->buildEdge(cellIDs[i]);
	}
}
/*----------------------------------------------------------------------------*/
void Edge::delegateGetAll(std::vector<Face>&   ACells) const
{
	if(!m_owner->m_model.has(E2F))
		throw GMDSException("E2F adjacency is not supported by the mesh model");
	std::vector<TCellID> cellIDs;
	(*m_edges_container->m_E2F)[m_id].allValues(cellIDs);
	ACells.resize(cellIDs.size());
	for(unsigned int i=0;i<cellIDs.size();i++){
		ACells[i] =m_owner->m_faces_container->buildFace(cellIDs[i]);
	}
}
/*----------------------------------------------------------------------------*/
void Edge::delegateGetAll(std::vector<Region>& ACells) const
{
	if(!m_owner->m_model.has(E2R))
		throw GMDSException("E2R adjacency is not supported by the mesh model");
	std::vector<TCellID> cellIDs;
	(*m_edges_container->m_E2R)[m_id].allValues(cellIDs);
	ACells.resize(cellIDs.size());
	for(unsigned int i=0;i<cellIDs.size();i++){
		ACells[i] =m_owner->m_regions_container->buildRegion(cellIDs[i]);
	}
}
/*----------------------------------------------------------------------------*/
void Edge::delegateGetAllNodeIDs(std::vector<TCellID>&   ACells) const
{
        if(!m_owner->m_model.has(E2N))
                throw GMDSException("E2N adjacency is not supported by the mesh model");
                
        (*m_edges_container->m_E2N)[m_id].allValues(ACells);
}
/*----------------------------------------------------------------------------*/
void Edge::delegateGetAllEdgeIDs(std::vector<TCellID>&   ACells) const
{
        if(!m_owner->m_model.has(E2E))
                throw GMDSException("E2E adjacency is not supported by the mesh model");
        (*m_edges_container->m_E2E)[m_id].allValues(ACells);
}
/*----------------------------------------------------------------------------*/
void Edge::delegateGetAllFaceIDs(std::vector<TCellID>&   ACells) const
{
        if(!m_owner->m_model.has(E2F))
                throw GMDSException("E2F adjacency is not supported by the mesh model");
        (*m_edges_container->m_E2F)[m_id].allValues(ACells);
}
/*----------------------------------------------------------------------------*/
void Edge::delegateGetAllRegionIDs(std::vector<TCellID>&   ACells) const
{
        if(!m_owner->m_model.has(E2R))
                throw GMDSException("E2R adjacency is not supported by the mesh model");
        (*m_edges_container->m_E2R)[m_id].allValues(ACells);
}

/*----------------------------------------------------------------------------*/
void Edge::delegateSetNodeIDs(const std::vector<TCellID>&   ACells)
{
	if(!m_owner->m_model.has(E2N))
		throw GMDSException("E2N adjacency is not supported by the mesh model");
	(*m_edges_container->m_E2N)[m_id]=ACells;
}
/*----------------------------------------------------------------------------*/
void Edge::delegateSetEdgeIDs(const std::vector<TCellID>&   ACells)
{
	if(!m_owner->m_model.has(E2E))
		throw GMDSException("E2E adjacency is not supported by the mesh model");
	(*m_edges_container->m_E2E)[m_id]=ACells;
}
/*----------------------------------------------------------------------------*/
void Edge::delegateSetFaceIDs(const std::vector<TCellID>&   ACells)
{
	if(!m_owner->m_model.has(E2F))
		throw GMDSException("E2F adjacency is not supported by the mesh model");
	(*m_edges_container->m_E2F)[m_id]=ACells;
}
/*----------------------------------------------------------------------------*/
void Edge::delegateSetRegionIDs(const std::vector<TCellID>&   ACells)
{
	if(!m_owner->m_model.has(E2R))
		throw GMDSException("E2R adjacency is not supported by the mesh model");
	(*m_edges_container->m_E2R)[m_id]=ACells;
}
/*----------------------------------------------------------------------------*/
void Edge::delegateNodeAdd(TCellID AID)
{
	if(!m_owner->m_model.has(E2N))
		throw GMDSException("E2N adjacency is not supported by the mesh model");
	(*m_edges_container->m_E2N)[m_id].add(AID);
}
/*----------------------------------------------------------------------------*/
void Edge::delegateEdgeAdd(TCellID AID)
{
	if(!m_owner->m_model.has(E2E))
		throw GMDSException("E2E adjacency is not supported by the mesh model");

	(*m_edges_container->m_E2E)[m_id].add(AID);
}
/*----------------------------------------------------------------------------*/
void Edge::delegateFaceAdd(TCellID AID)
{
	if(!m_owner->m_model.has(E2F))
		throw GMDSException("E2F adjacency is not supported by the mesh model");
	(*m_edges_container->m_E2F)[m_id].add(AID);
}
/*----------------------------------------------------------------------------*/
void Edge::delegateRegionAdd(TCellID AID)
{
	if(!m_owner->m_model.has(E2R))
		throw GMDSException("E2R adjacency is not supported by the mesh model");
	(*m_edges_container->m_E2R)[m_id].add(AID);
}
/*----------------------------------------------------------------------------*/
void Edge::delegateNodeRemove(TCellID AID)
{
	if(!m_owner->m_model.has(E2N))
		throw GMDSException("E2N adjacency is not supported by the mesh model");
	(*m_edges_container->m_E2N)[m_id].del(AID);
}
/*----------------------------------------------------------------------------*/
void Edge::delegateEdgeRemove(TCellID AID)
{
	if(!m_owner->m_model.has(E2E))
		throw GMDSException("E2E adjacency is not supported by the mesh model");
	(*m_edges_container->m_E2E)[m_id].del(AID);

}
/*----------------------------------------------------------------------------*/
void Edge::delegateFaceRemove(TCellID AID)
{
	if(!m_owner->m_model.has(E2F))
		throw GMDSException("E2F adjacency is not supported by the mesh model");

	(*m_edges_container->m_E2F)[m_id].del(AID);

}
/*----------------------------------------------------------------------------*/
void Edge::delegateRegionRemove(TCellID AID)
{
	if(!m_owner->m_model.has(E2R))
		throw GMDSException("E2R adjacency is not supported by the mesh model");
	(*m_edges_container->m_E2R)[m_id].del(AID);
}
/*----------------------------------------------------------------------------*/
void Edge::delegateNodeReplace(TCellID AID1, TCellID AID2)
{
	if(!m_owner->m_model.has(E2N))
		throw GMDSException("E2N adjacency is not supported by the mesh model");

	(*m_edges_container->m_E2N)[m_id].replace(AID1, AID2);
}
/*----------------------------------------------------------------------------*/
void Edge::delegateEdgeReplace(TCellID AID1, TCellID AID2)
{
	if(!m_owner->m_model.has(E2E))
		throw GMDSException("E2E adjacency is not supported by the mesh model");

	(*m_edges_container->m_E2E)[m_id].replace(AID1, AID2);
}
/*----------------------------------------------------------------------------*/
void Edge::delegateFaceReplace(TCellID AID1, TCellID AID2)
{
	if(!m_owner->m_model.has(E2F))
		throw GMDSException("E2F adjacency is not supported by the mesh model");

	(*m_edges_container->m_E2F)[m_id].replace(AID1, AID2);
}
/*----------------------------------------------------------------------------*/
void Edge::delegateRegionReplace(TCellID AID1, TCellID AID2)
{
	if(!m_owner->m_model.has(E2R))
		throw GMDSException("E2R adjacency is not supported by the mesh model");

	(*m_edges_container->m_E2R)[m_id].replace(AID1, AID2);
}
/*----------------------------------------------------------------------------*/
std::ostream & operator << (std::ostream & AStream, const Edge & AN)
{
	AStream<<"Edge "<<AN.getID();
	return AStream;
}
/*----------------------------------------------------------------------------*/
