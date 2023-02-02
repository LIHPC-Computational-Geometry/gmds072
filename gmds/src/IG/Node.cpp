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
 * Node.cpp
 *
 *  Created on: 5 févr. 2014
 *      Author: ledouxf
 */
/*----------------------------------------------------------------------------*/
#include <GMDS/IG/Node.h>
#include <GMDS/IG/IGMesh.h>
#include <GMDS/IG/EdgeContainer.h>
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
Node::
Node()
: Cell(0,GMDS_NODE,NullID)
{
	m_nodes_container = 0;
}
/*----------------------------------------------------------------------------*/
Node::
Node(IGMesh* AMesh,  const TCellID& AID,
		const TCoord& AX, const TCoord& AY, const TCoord& AZ)
: Cell(AMesh,GMDS_NODE,AID)
{
	if(AMesh!=0){
		m_nodes_container  = AMesh->m_nodes_container;
		m_nodes_container->m_node_coords[m_id].setXYZ(AX,AY,AZ);
	}
	else
		m_nodes_container=0;
}
/*----------------------------------------------------------------------------*/
Node::
Node(IGMesh* AMesh, const TCellID& AID,
		const math::Point& APt)
: Cell(AMesh,GMDS_NODE,AID)
{
	if(AMesh!=0){
		m_nodes_container  = AMesh->m_nodes_container;
		m_nodes_container->m_node_coords[m_id]=APt;
	}
	else
		m_nodes_container=0;
}
/*----------------------------------------------------------------------------*/
Node::
Node(const Node& ANode)
: Cell(ANode.m_owner,GMDS_NODE,ANode.m_id)
{
	if(m_owner!=0)
		m_nodes_container  = m_owner->m_nodes_container;
	else
		m_nodes_container = 0;
}
/*----------------------------------------------------------------------------*/
void Node::operator=(const Node& ANode)
{
	m_owner = ANode.m_owner;
	m_id = ANode.m_id;
	if(m_owner!=0)
		m_nodes_container  = m_owner->m_nodes_container;
	else
		m_nodes_container = 0;
}
/*----------------------------------------------------------------------------*/
bool Node::operator==(const Node& ANode) const
{
	return (m_owner == ANode.m_owner && m_id==ANode.m_id);
}
/*----------------------------------------------------------------------------*/
bool Node::operator!=(const Node& ANode) const
{
        return (!(*this == ANode));
}
/*----------------------------------------------------------------------------*/
Node::~Node()
{}
/*----------------------------------------------------------------------------*/
TInt Node::getNbNodes() const
{
	TabCellID<size_undef> cells = (*m_nodes_container->m_N2N)[m_id];
	return cells.size();
}
/*----------------------------------------------------------------------------*/
TInt Node::getNbEdges() const
{
	TabCellID<size_undef> cells = (*m_nodes_container->m_N2E)[m_id];
	return cells.size();
}
/*----------------------------------------------------------------------------*/
TInt Node::getNbFaces() const
{
	TabCellID<size_undef> cells = (*m_nodes_container->m_N2F)[m_id];
	return cells.size();
}
/*----------------------------------------------------------------------------*/
TInt Node::getNbRegions() const
{
	TabCellID<size_undef> cells = (*m_nodes_container->m_N2R)[m_id];
	return cells.size();
}
/*----------------------------------------------------------------------------*/
gmds::math::Point Node::center() const
{
        return m_nodes_container->m_node_coords[m_id];
}
/*----------------------------------------------------------------------------*/
void Node::delegateGet(std::vector<Node>&   ACells) const
{
	if(!m_owner->m_model.has(N2N))
		throw GMDSException("N2N adjacency is not supported by the mesh model");

	std::vector<TCellID> cellIDs;
	(*m_nodes_container->m_N2N)[m_id].values(cellIDs);
	ACells.resize(cellIDs.size());
	for(unsigned int i=0;i<cellIDs.size();i++){
		ACells[i] =m_nodes_container->buildNode(cellIDs[i]);
	}
}
/*----------------------------------------------------------------------------*/
void Node::delegateGet(std::vector<Edge>&   ACells) const
{
	if(!m_owner->m_model.has(N2E))
		throw GMDSException("N2E adjacency is not supported by the mesh model");
	std::vector<TCellID> cellIDs;
	(*m_nodes_container->m_N2E)[m_id].values(cellIDs);
	ACells.resize(cellIDs.size());
	for(unsigned int i=0;i<cellIDs.size();i++){
		ACells[i] =m_owner->m_edges_container->buildEdge(cellIDs[i]);
	}
}
/*----------------------------------------------------------------------------*/
void Node::delegateGet(std::vector<Face>&   ACells) const
{
	if(!m_owner->m_model.has(N2F))
		throw GMDSException("N2F adjacency is not supported by the mesh model");
	std::vector<TCellID> cellIDs;
	(*m_nodes_container->m_N2F)[m_id].values(cellIDs);
	ACells.resize(cellIDs.size());
	for(unsigned int i=0;i<cellIDs.size();i++){
		ACells[i] =m_owner->m_faces_container->buildFace(cellIDs[i]);
	}
}
/*----------------------------------------------------------------------------*/
void Node::delegateGet(std::vector<Region>& ACells) const
{
	if(!m_owner->m_model.has(N2R))
		throw GMDSException("N2R adjacency is not supported by the mesh model");
	std::vector<TCellID> cellIDs;
	(*m_nodes_container->m_N2R)[m_id].values(cellIDs);
	ACells.resize(cellIDs.size());
	for(unsigned int i=0;i<cellIDs.size();i++){
		ACells[i] =m_owner->m_regions_container->buildRegion(cellIDs[i]);
	}
}
/*----------------------------------------------------------------------------*/
void Node::delegateGetNodeIDs(std::vector<TCellID>&   ACells) const
{
	if(!m_owner->m_model.has(N2N))
		throw GMDSException("N2N adjacency is not supported by the mesh model");

	(*m_nodes_container->m_N2N)[m_id].values(ACells);
}
/*----------------------------------------------------------------------------*/
void Node::delegateGetEdgeIDs(std::vector<TCellID>&   ACells) const
{
	if(!m_owner->m_model.has(N2E))
		throw GMDSException("N2E adjacency is not supported by the mesh model");
	(*m_nodes_container->m_N2E)[m_id].values(ACells);
}
/*----------------------------------------------------------------------------*/
void Node::delegateGetFaceIDs(std::vector<TCellID>&   ACells) const
{
	if(!m_owner->m_model.has(N2F))
		throw GMDSException("N2F adjacency is not supported by the mesh model");
	(*m_nodes_container->m_N2F)[m_id].values(ACells);
}
/*----------------------------------------------------------------------------*/
void Node::delegateGetRegionIDs(std::vector<TCellID>&   ACells) const
{
	if(!m_owner->m_model.has(N2R))
		throw GMDSException("N2R adjacency is not supported by the mesh model");
	(*m_nodes_container->m_N2R)[m_id].values(ACells);
}
/*----------------------------------------------------------------------------*/
void Node::delegateGetAll(std::vector<Node>&   ACells) const
{
	if(!m_owner->m_model.has(N2N))
		throw GMDSException("N2N adjacency is not supported by the mesh model");
	std::vector<TCellID> cellIDs;
	(*m_nodes_container->m_N2N)[m_id].allValues(cellIDs);
	ACells.resize(cellIDs.size());
	for(unsigned int i=0;i<cellIDs.size();i++){
		ACells[i] =m_nodes_container->buildNode(cellIDs[i]);
	}
}
/*----------------------------------------------------------------------------*/
void Node::delegateGetAll(std::vector<Edge>&   ACells) const
{
	if(!m_owner->m_model.has(N2E))
		throw GMDSException("N2E adjacency is not supported by the mesh model");
	std::vector<TCellID> cellIDs;
	(*m_nodes_container->m_N2E)[m_id].allValues(cellIDs);
	ACells.resize(cellIDs.size());
	for(unsigned int i=0;i<cellIDs.size();i++){
		ACells[i] =m_owner->m_edges_container->buildEdge(cellIDs[i]);
	}
}
/*----------------------------------------------------------------------------*/
void Node::delegateGetAll(std::vector<Face>&   ACells) const
{
	if(!m_owner->m_model.has(N2F))
		throw GMDSException("N2F adjacency is not supported by the mesh model");
	std::vector<TCellID> cellIDs;
	(*m_nodes_container->m_N2F)[m_id].allValues(cellIDs);
	ACells.resize(cellIDs.size());
	for(unsigned int i=0;i<cellIDs.size();i++){
		ACells[i] =m_owner->m_faces_container->buildFace(cellIDs[i]);
	}
}
/*----------------------------------------------------------------------------*/
void Node::delegateGetAll(std::vector<Region>& ACells) const
{
	if(!m_owner->m_model.has(N2R))
		throw GMDSException("N2R adjacency is not supported by the mesh model");
	std::vector<TCellID> cellIDs;
	(*m_nodes_container->m_N2R)[m_id].allValues(cellIDs);
	ACells.resize(cellIDs.size());
	for(unsigned int i=0;i<cellIDs.size();i++){
		ACells[i] =m_owner->m_regions_container->buildRegion(cellIDs[i]);
	}
}
/*----------------------------------------------------------------------------*/
void Node::delegateGetAllNodeIDs(std::vector<TCellID>&   ACells) const
{
        if(!m_owner->m_model.has(N2N))
                throw GMDSException("N2N adjacency is not supported by the mesh model");

        (*m_nodes_container->m_N2N)[m_id].allValues(ACells);
}
/*----------------------------------------------------------------------------*/
void Node::delegateGetAllEdgeIDs(std::vector<TCellID>&   ACells) const
{
        if(!m_owner->m_model.has(N2E))
                throw GMDSException("N2E adjacency is not supported by the mesh model");
        (*m_nodes_container->m_N2E)[m_id].allValues(ACells);
}
/*----------------------------------------------------------------------------*/
void Node::delegateGetAllFaceIDs(std::vector<TCellID>&   ACells) const
{
        if(!m_owner->m_model.has(N2F))
                throw GMDSException("N2F adjacency is not supported by the mesh model");
        (*m_nodes_container->m_N2F)[m_id].allValues(ACells);
}
/*----------------------------------------------------------------------------*/
void Node::delegateGetAllRegionIDs(std::vector<TCellID>&   ACells) const
{
        if(!m_owner->m_model.has(N2R))
                throw GMDSException("N2R adjacency is not supported by the mesh model");
        (*m_nodes_container->m_N2R)[m_id].allValues(ACells);
}
/*----------------------------------------------------------------------------*/
void Node::delegateSetNodeIDs(const std::vector<TCellID>&   ACells)
{
	if(!m_owner->m_model.has(N2N))
		throw GMDSException("N2N adjacency is not supported by the mesh model");
	(*m_nodes_container->m_N2N)[m_id]=ACells;
}
/*----------------------------------------------------------------------------*/
void Node::delegateSetEdgeIDs(const std::vector<TCellID>&   ACells)
{
	if(!m_owner->m_model.has(N2E))
		throw GMDSException("N2E adjacency is not supported by the mesh model");
	(*m_nodes_container->m_N2E)[m_id]=ACells;
}
/*----------------------------------------------------------------------------*/
void Node::delegateSetFaceIDs(const std::vector<TCellID>&   ACells)
{
	if(!m_owner->m_model.has(N2F))
			throw GMDSException("N2F adjacency is not supported by the mesh model");
	(*m_nodes_container->m_N2F)[m_id]=ACells;
}
/*----------------------------------------------------------------------------*/
void Node::delegateSetRegionIDs(const std::vector<TCellID>& ACells)
{
	if(!m_owner->m_model.has(N2R))
		throw GMDSException("N2R adjacency is not supported by the mesh model");
	(*m_nodes_container->m_N2R)[m_id]=ACells;
}
/*----------------------------------------------------------------------------*/
void Node::delegateNodeAdd(TCellID AID)
{
	if(!m_owner->m_model.has(N2N))
		throw GMDSException("N2N adjacency is not supported by the mesh model");
	(*m_nodes_container->m_N2N)[m_id].add(AID);
}
/*----------------------------------------------------------------------------*/
void Node::delegateEdgeAdd(TCellID AID)
{
	if(!m_owner->m_model.has(N2E))
		throw GMDSException("N2E adjacency is not supported by the mesh model");
	(*m_nodes_container->m_N2E)[m_id].add(AID);
}
/*----------------------------------------------------------------------------*/
void Node::delegateFaceAdd(TCellID AID)
{
	if(!m_owner->m_model.has(N2F))
		throw GMDSException("N2F adjacency is not supported by the mesh model");
	(*m_nodes_container->m_N2F)[m_id].add(AID);
}
/*----------------------------------------------------------------------------*/
void Node::delegateRegionAdd(TCellID AID)
{
	if(!m_owner->m_model.has(N2R))
		throw GMDSException("N2R adjacency is not supported by the mesh model");
	(*m_nodes_container->m_N2R)[m_id].add(AID);
}
/*----------------------------------------------------------------------------*/
void Node::delegateNodeRemove(TCellID AID)
{
	if(!m_owner->m_model.has(N2N))
		throw GMDSException("N2N adjacency is not supported by the mesh model");
	(*m_nodes_container->m_N2N)[m_id].del(AID);
}
/*----------------------------------------------------------------------------*/
void Node::delegateEdgeRemove(TCellID AID)
{
	if(!m_owner->m_model.has(N2E))
		throw GMDSException("N2E adjacency is not supported by the mesh model");
	(*m_nodes_container->m_N2E)[m_id].del(AID);
}
/*----------------------------------------------------------------------------*/
void Node::delegateFaceRemove(TCellID AID)
{
	if(!m_owner->m_model.has(N2F))
		throw GMDSException("N2F adjacency is not supported by the mesh model");
	(*m_nodes_container->m_N2F)[m_id].del(AID);
}
/*----------------------------------------------------------------------------*/
void Node::delegateRegionRemove(TCellID AID)
{
	if(!m_owner->m_model.has(N2R))
		throw GMDSException("N2R adjacency is not supported by the mesh model");
	(*m_nodes_container->m_N2R)[m_id].del(AID);
}
/*----------------------------------------------------------------------------*/
void Node::delegateNodeReplace(TCellID AID1, TCellID AID2)
{
	if(!m_owner->m_model.has(N2N))
		throw GMDSException("N2N adjacency is not supported by the mesh model");
	(*m_nodes_container->m_N2N)[m_id].replace(AID1, AID2);
}
/*----------------------------------------------------------------------------*/
void Node::delegateEdgeReplace(TCellID AID1, TCellID AID2)
{
	if(!m_owner->m_model.has(N2E))
		throw GMDSException("N2E adjacency is not supported by the mesh model");
	(*m_nodes_container->m_N2E)[m_id].replace(AID1, AID2);
}
/*----------------------------------------------------------------------------*/
void Node::delegateFaceReplace(TCellID AID1, TCellID AID2)
{
	if(!m_owner->m_model.has(N2F))
		throw GMDSException("N2F adjacency is not supported by the mesh model");
	(*m_nodes_container->m_N2F)[m_id].replace(AID1, AID2);
}
/*----------------------------------------------------------------------------*/
void Node::delegateRegionReplace(TCellID AID1, TCellID AID2)
{
	if(!m_owner->m_model.has(N2R))
		throw GMDSException("N2R adjacency is not supported by the mesh model");
	(*m_nodes_container->m_N2R)[m_id].replace(AID1, AID2);
}
/*----------------------------------------------------------------------------*/
math::Point Node::getPoint() const
{
	return  m_nodes_container->m_node_coords[m_id];
}
/*----------------------------------------------------------------------------*/
void Node::setPoint(const math::Point& APnt)
{
	m_nodes_container->m_node_coords[m_id] = APnt;
}
/*----------------------------------------------------------------------------*/
TCoord Node::X() const
{
	return m_nodes_container->m_node_coords[m_id].X();
}
/*----------------------------------------------------------------------------*/
TCoord Node::Y() const
{
	return m_nodes_container->m_node_coords[m_id].Y();
}
/*----------------------------------------------------------------------------*/
TCoord Node::Z() const
{
	return m_nodes_container->m_node_coords[m_id].Z();
}
/*----------------------------------------------------------------------------*/
TCoord& Node::X()
{
	return m_nodes_container->m_node_coords[m_id].X();
}
/*----------------------------------------------------------------------------*/
TCoord& Node::Y()
{
	return m_nodes_container->m_node_coords[m_id].Y();
}
/*----------------------------------------------------------------------------*/
TCoord& Node::Z()
{
	return m_nodes_container->m_node_coords[m_id].Z();
}
/*----------------------------------------------------------------------------*/
void Node::setX(const TCoord AVal)
{
	m_nodes_container->m_node_coords[m_id].setX(AVal);
}
/*----------------------------------------------------------------------------*/
void Node::setY(const TCoord AVal)
{
	m_nodes_container->m_node_coords[m_id].setY(AVal);
}
/*----------------------------------------------------------------------------*/
void Node::setZ(const TCoord AVal)
{
	m_nodes_container->m_node_coords[m_id].setZ(AVal);
}
/*----------------------------------------------------------------------------*/
void Node::setXYZ(const TCoord AX, const TCoord AY,const TCoord AZ)
{
	m_nodes_container->m_node_coords[m_id].setXYZ(AX,AY,AZ);
}
/*----------------------------------------------------------------------------*/
std::ostream& operator<<(std::ostream& AStream, const Node& AN)
{
        AStream<<"Node "<<AN.getID()<<" ("
                        <<AN.getPoint().X()<<", "
                        <<AN.getPoint().Y()<<", "
                        <<AN.getPoint().Z()<<")";
        return AStream;
}
/*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/
