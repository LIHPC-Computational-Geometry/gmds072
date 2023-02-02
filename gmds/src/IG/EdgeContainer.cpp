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
 * EdgeContainer.cpp
 *
 *  Created on: 26 mars 2014
 *      Author: ledouxf
 */
/*----------------------------------------------------------------------------*/
#include <GMDS/IG/EdgeContainer.h>
#include <GMDS/IG/IGMesh.h>
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
EdgeContainer::EdgeContainer( IGMesh* AMesh)
:m_mesh(AMesh),m_model(AMesh->getModel()),
 m_E2N(0),m_E2E(0),m_E2F(0),m_E2R(0)
{
	if(m_model.has(E2N))
	{
		m_E2N = new SmartVector<TabCellID<2> >();
	}
	if(m_model.has(E2E)){
		m_E2E = new SmartVector<TabCellID<size_undef> >();
	}
	if(m_model.has(E2F)){
		m_E2F = new SmartVector<TabCellID<size_undef> >();
	}
	if(m_model.has(E2R)){
		m_E2R = new SmartVector<TabCellID<size_undef> >();
	}
}
/*----------------------------------------------------------------------------*/
EdgeContainer::~EdgeContainer()
{
	if(m_E2N)
		delete m_E2N;
	if(m_E2E)
		delete m_E2E;
	if(m_E2F)
		delete m_E2F;
	if(m_E2R)
		delete m_E2R;
}
/*----------------------------------------------------------------------------*/
void EdgeContainer::addConnectivityContainers(const TInt ADim)
{
	if(ADim==0){
		if(m_E2N==0){
			m_E2N = new SmartVector<TabCellID<2> >(
				m_edge_ids.getBits());
			m_edge_ids.update();
		}
	}
	else if (ADim==1){
		if (m_E2E==0){
			m_E2E = new SmartVector<TabCellID<size_undef> >(
				m_edge_ids.getBits());
			m_edge_ids.update();
		}
	}
	else if (ADim==2){
		if (m_E2F==0){
			m_E2F = new SmartVector<TabCellID<size_undef> >(
				m_edge_ids.getBits());
			m_edge_ids.update();			
		}
	}
	else if (ADim==3){
		if(m_E2R==0){
			m_E2R = new SmartVector<TabCellID<size_undef> >(
				m_edge_ids.getBits());
			m_edge_ids.update();
		}
	}
}
/*----------------------------------------------------------------------------*/
void EdgeContainer::removeConnectivityContainers(const TInt ADim)
{
	if(ADim==0){
		if(m_E2N!=0){
			delete m_E2N;
			m_E2N=0;
		}
	}
	else if (ADim==1){
		if(m_E2E!=0){
			delete m_E2E;
			m_E2E=0;
		}
	}
	else if (ADim==2){
		if(m_E2F!=0){
			delete m_E2F;
			m_E2F=0;
		}
	}
	else if (ADim==3)
	{
		if(m_E2R!=0){
			delete m_E2R;
			m_E2R=0;
		}
	}
}
/*----------------------------------------------------------------------------*/
void EdgeContainer::getNodesData(const TCellID& AID, int& ANbNodes)
{
	if(m_E2N==0)
		ANbNodes = 0;
	else
		ANbNodes = 2;
}
/*----------------------------------------------------------------------------*/
void EdgeContainer::getEdgesData(const TCellID& AID, int& ANbEdges)
{
	if(m_E2E==0)
		ANbEdges = 0;
	else
		ANbEdges = (*m_E2E)[AID].size();
}
/*----------------------------------------------------------------------------*/
void EdgeContainer::getFacesData(const TCellID& AID, int& ANbFaces)
{
	if(m_E2F==0)
		ANbFaces = 0;
	else
		ANbFaces = (*m_E2F)[AID].size();
}
/*----------------------------------------------------------------------------*/
void EdgeContainer::getRegionsData(const TCellID& AID, int& ANbRegions)
{
	if(m_E2R==0)
		ANbRegions = 0;
	else
		ANbRegions = (*m_E2R)[AID].size();
}

/*----------------------------------------------------------------------------*/
Edge EdgeContainer::add(const TCellID& AN1,const TCellID& AN2)
{
	//===========================================
	// STEP 1 - we get the edgeindex
	//===========================================
	TCellID index = m_edge_ids.getFreeIndex();
	return add(AN1,AN2,index);
}
/*----------------------------------------------------------------------------*/
bool EdgeContainer::has(const TCellID& AID)
{
	return m_edge_ids[AID];
}
/*----------------------------------------------------------------------------*/
Edge EdgeContainer::
add(const TCellID& AN1,const TCellID& AN2, const TCellID& AID)
{
	TCellID index = AID;
	m_edge_ids.assign(index);
	//============================================
	// STEP 2 - we fill in F2N
	//============================================
	if(m_model.has(E2N)){
		//here only nodes are updated
		TabCellID<2> t;
		t.add(AN1);
		t.add(AN2);
		TInt index_E2N = m_E2N->selectNewIndex();
		if(index_E2N!=index)
			throw GMDSException("Consistency trouble in the FaceContainer class");
		m_E2N->assign(t,index);
	}
	if(m_model.has(E2E)){
		//here only nodes are updated
		TabCellID<size_undef> t;
		TInt index_E2E = m_E2E->selectNewIndex();
		if(index_E2E!=index)
			throw GMDSException("Consistency trouble in the EdgeContainer class");
		m_E2E->assign(t,index);
	}
	if(m_model.has(E2F)){
		//here only nodes are updated
		TabCellID<size_undef> t;
		TInt index_E2F = m_E2F->selectNewIndex();
		if(index_E2F!=index)
			throw GMDSException("Consistency trouble in the EdgeContainer class");
		m_E2F->assign(t,index);
	}
	if(m_model.has(E2R)){
		//here only nodes are updated
		TabCellID<size_undef> t;
		TInt index_E2R = m_E2R->selectNewIndex();
		if(index_E2R!=index)
			throw GMDSException("Consistency trouble in the EdgeContainer class");
		m_E2R->assign(t,index);
	}

	//============================================
	// STEP 3 - we create an edge object
	//============================================
	Edge e(m_mesh,index);
	return e;
}
/*----------------------------------------------------------------------------*/
Edge EdgeContainer::buildEdge(const TInt index)
{
	Edge e(m_mesh,index);
	return e;
}
/*----------------------------------------------------------------------------*/
EdgeContainer::iterator EdgeContainer::getIterator()
{
	return iterator(this);
}
/*----------------------------------------------------------------------------*/
EdgeContainer::iterator::iterator():m_container(0){;}
/*----------------------------------------------------------------------------*/
EdgeContainer::iterator::iterator(EdgeContainer* nc):m_container(nc)
{
	m_iterator = m_container->m_edge_ids.begin();
}
/*----------------------------------------------------------------------------*/
Edge EdgeContainer::iterator::value() const
{
	TInt node_id = m_iterator.value();
	return m_container->buildEdge(node_id);
}
/*----------------------------------------------------------------------------*/
void EdgeContainer::iterator::next()
{
	m_iterator.next();
}
/*----------------------------------------------------------------------------*/
bool EdgeContainer::iterator::isDone() const
{
	return m_iterator.isDone();
}
/*----------------------------------------------------------------------------*/
void EdgeContainer::clear()
{
	m_edge_ids.clear();
	if(m_E2N)
		m_E2N->clear();
	if(m_E2E)
		m_E2E->clear();
	if(m_E2F)
		m_E2F->clear();
	if(m_E2R)
		m_E2R->clear();
}
/*----------------------------------------------------------------------------*/
void EdgeContainer::resize(const TInt ASize)
{
	m_edge_ids.resize(ASize);
}

/*----------------------------------------------------------------------------*/
void EdgeContainer::update()
{
	m_edge_ids.update();
}
/*----------------------------------------------------------------------------*/
void EdgeContainer::serialize(std::ostream& AStr)
{
	m_edge_ids.serialize(AStr);
#ifdef GMDS_PARALLEL
	m_distributed_data.serialize(AStr);
#endif //GMDS_PARALLEL

	if(m_E2N)
		m_E2N->serialize(AStr);
	if(m_E2E)
		m_E2E->serialize(AStr);
	if(m_E2F)
		m_E2F->serialize(AStr);
	if(m_E2R)
		m_E2R->serialize(AStr);
}
/*----------------------------------------------------------------------------*/
void EdgeContainer::unserialize(std::istream& AStr)
{
	m_edge_ids.unserialize(AStr);
#ifdef GMDS_PARALLEL
	m_distributed_data.unserialize(AStr);
#endif //GMDS_PARALLEL

	if(m_E2N)
		m_E2N->unserialize(AStr);
	if(m_E2E)
		m_E2E->unserialize(AStr);
	if(m_E2F)
		m_E2F->unserialize(AStr);
	if(m_E2R)
		m_E2R->unserialize(AStr);

}
/*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/
