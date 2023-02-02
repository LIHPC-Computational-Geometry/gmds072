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
 * NodeContainer.cpp
 *
 *  Created on: 26 mars 2014
 *      Author: ledouxf
 */
/*----------------------------------------------------------------------------*/
#include <GMDS/IG/NodeContainer.h>
#include <GMDS/IG/IGMesh.h>
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
NodeContainer::NodeContainer( IGMesh* AMesh)
:m_mesh(AMesh),m_model(AMesh->getModel())
{
	if(m_model.has(N2N))
		m_N2N = new SmartVector<TabCellID<size_undef> >();
	else
		m_N2N = 0;
	if(m_model.has(N2E))
		m_N2E = new SmartVector<TabCellID<size_undef> >();
	else
		m_N2E = 0;
	if(m_model.has(N2F))
		m_N2F = new SmartVector<TabCellID<size_undef> >();
	else
		m_N2F = 0;
	if(m_model.has(N2R))
		m_N2R = new SmartVector<TabCellID<size_undef> >();
	else
		m_N2R = 0;
}
/*----------------------------------------------------------------------------*/
NodeContainer::~NodeContainer()
{
	if(m_N2N)
		delete m_N2N;
	if(m_N2E)
		delete m_N2E;
	if(m_N2F)
		delete m_N2F;
	if(m_N2R)
		delete m_N2R;
}
/*----------------------------------------------------------------------------*/
Node NodeContainer::add(const TCoord& AX, const TCoord& AY, const TCoord& AZ)
{
	TInt index = m_node_ids.getFreeIndex();
	return add(AX,AY,AZ,index);
}
/*----------------------------------------------------------------------------*/
Node NodeContainer::add(const TCoord& AX, const TCoord& AY, const TCoord& AZ,
		const TCellID& AGID)
{
	TCellID index = AGID;
	m_node_ids.assign(index);

	math::Point p(AX,AY,AZ);
	m_node_coords.assign(p, index);

	if(m_model.has(N2N)){
		//here only nodes are updated
		TabCellID<size_undef> t;
		TInt index_N2N = m_N2N->selectNewIndex();
		if(index_N2N!=index)
			throw GMDSException("Consistency trouble in the NodeContainer class");
		m_N2N->assign(t,index);
	}
	if(m_model.has(N2E)){
		//here only nodes are updated
		TabCellID<size_undef> t;
		TInt index_N2E = m_N2E->selectNewIndex();
		if(index_N2E!=index)
			throw GMDSException("Consistency trouble in the NodeContainer class");
		m_N2E->assign(t,index);
	}
	if(m_model.has(N2F)){
		//here only nodes are updated
		TabCellID<size_undef> t;
		TInt index_N2F = m_N2F->selectNewIndex();
		if(index_N2F!=index)
			throw GMDSException("Consistency trouble in the NodeContainer class");
		m_N2F->assign(t,index);
	}
	if(m_model.has(N2R)){
		//here only nodes are updated
		TabCellID<size_undef> t;
		TInt index_N2R = m_N2R->selectNewIndex();
		if(index_N2R!=index)
			throw GMDSException("Consistency trouble in the NodeContainer class");
		m_N2R->assign(t,index);
	}

	Node n(m_mesh,index,p);
	return n;
}
/*----------------------------------------------------------------------------*/
void NodeContainer::addConnectivityContainers(const TInt ADim)
{
	/** WARNING
	 * When wa add a new adjacency container, we must define its size. If the
	 * cell was already used, then X2N exists and is filled, so we resize to
	 * X2N->size()
	 */
	if(ADim==0){
		if(m_N2N==0){
			m_N2N = new SmartVector<TabCellID<size_undef> >(
				m_node_ids.getBits());
			m_node_ids.update();
		}
	}
	else if (ADim==1){
		if (m_N2E==0){
                        m_N2E = new SmartVector<TabCellID<size_undef> >(
                                m_node_ids.getBits());
			m_node_ids.update();
		}
	}
	else if (ADim==2){
		if (m_N2F==0){
                        m_N2F = new SmartVector<TabCellID<size_undef> >(
                                m_node_ids.getBits());
			m_node_ids.update();
		}
	}
	else if (ADim==3){
		if(m_N2R==0){
                        m_N2R = new SmartVector<TabCellID<size_undef> >(
                                m_node_ids.getBits());
			m_node_ids.update();
		}
	}
}
/*----------------------------------------------------------------------------*/
void NodeContainer::removeConnectivityContainers(const TInt ADim)
{
	if(ADim==0){
		if(m_N2N!=0){
			delete m_N2N;
			m_N2N=0;
		}
	}
	else if (ADim==1){
		if(m_N2E!=0){
			delete m_N2E;
			m_N2E=0;
		}
	}
	else if (ADim==2){
		if(m_N2F!=0){
			delete m_N2F;
			m_N2F=0;
		}
	}
	else if (ADim==3)
	{
		if(m_N2R!=0){
			delete m_N2R;
			m_N2R=0;
		}
	}
}
/*----------------------------------------------------------------------------*/
void NodeContainer::getNodesData(const TCellID& AID, int& ANbNodes)
{
	if(m_N2N==0)
		ANbNodes = 0;
	else
		ANbNodes = (*m_N2N)[AID].size();
}
/*----------------------------------------------------------------------------*/
void NodeContainer::getEdgesData(const TCellID& AID, int& ANbEdges)
{
	if(m_N2E==0)
		ANbEdges = 0;
	else
		ANbEdges = (*m_N2E)[AID].size();
}
/*----------------------------------------------------------------------------*/
void NodeContainer::getFacesData(const TCellID& AID, int& ANbFaces)
{
	if(m_N2F==0)
		ANbFaces = 0;
	else
		ANbFaces = (*m_N2F)[AID].size();
}
/*----------------------------------------------------------------------------*/
void NodeContainer::getRegionsData(const TCellID& AID, int& ANbRegions)
{
	if(m_N2R==0)
		ANbRegions = 0;
	else
		ANbRegions = (*m_N2R)[AID].size();
}
/*------------------------------------------------------------------------*/
math::Point NodeContainer::getPoint(const TCellID& AID)
{
	return m_node_coords[AID];
}
/*----------------------------------------------------------------------------*/
Node NodeContainer::buildNode(const TInt index)
{
	math::Point p = m_node_coords[index];
	Node n(m_mesh,index,p);
	return n;
}
/*----------------------------------------------------------------------------*/
bool NodeContainer::has(const TCellID& AID)
{
	return m_node_ids[AID];
}

/*----------------------------------------------------------------------------*/
void NodeContainer::update()
{
	m_node_ids.update();
}
/*----------------------------------------------------------------------------*/
void NodeContainer::clear()
{
	m_node_ids.clear();
	m_node_coords.clear();
	if(m_N2N)
		m_N2N->clear();
	if(m_N2E)
		m_N2E->clear();
	if(m_N2F)
		m_N2F->clear();
	if(m_N2R)
		m_N2R->clear();
}
/*----------------------------------------------------------------------------*/
void NodeContainer::resize(const TInt ASize)
{
	m_node_ids.resize(ASize);
	m_node_coords.resize(ASize);
}
/*----------------------------------------------------------------------------*/
NodeContainer::iterator NodeContainer::getIterator()
{
	return iterator(this);
}
/*----------------------------------------------------------------------------*/
void NodeContainer::serialize(std::ostream& AStr)
{
	m_node_ids.serialize(AStr);
	m_node_coords.serialize(AStr);
#ifdef GMDS_PARALLEL
	m_distributed_data.serialize(AStr);
#endif //GMDS_PARALLEL

	if(m_N2N)
		m_N2N->serialize(AStr);
	if(m_N2E)
		m_N2E->serialize(AStr);
	if(m_N2F)
		m_N2F->serialize(AStr);
	if(m_N2R)
		m_N2R->serialize(AStr);
}
/*----------------------------------------------------------------------------*/
void NodeContainer::unserialize(std::istream& AStr)
{
	clear();
	m_node_ids.unserialize(AStr);
	m_node_coords.unserialize(AStr);
#ifdef GMDS_PARALLEL
	m_distributed_data.unserialize(AStr);
#endif //GMDS_PARALLEL

	if(m_N2N)
		m_N2N->unserialize(AStr);
	if(m_N2E)
		m_N2E->unserialize(AStr);
	if(m_N2F)
		m_N2F->unserialize(AStr);
	if(m_N2R)
		m_N2R->unserialize(AStr);
}
/*----------------------------------------------------------------------------*/
NodeContainer::iterator::iterator():m_container(0){;}
/*----------------------------------------------------------------------------*/
NodeContainer::iterator::iterator(NodeContainer* nc):m_container(nc)
{
	m_iterator = m_container->m_node_ids.begin();
}
/*----------------------------------------------------------------------------*/
Node NodeContainer::iterator::value() const
{
	TInt node_id = m_iterator.value();
	return m_container->buildNode(node_id);
}
/*----------------------------------------------------------------------------*/
void NodeContainer::iterator::next()
{
	m_iterator.next();
}
/*----------------------------------------------------------------------------*/
bool NodeContainer::iterator::isDone() const
{
	return m_iterator.isDone();
}
/*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/
