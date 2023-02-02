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
 * IGMesh.cpp
 *
 *  Created on: 5 févr. 2014
 *      Author: ledouxf
 */
/*----------------------------------------------------------------------------*/
#include <GMDS/IG/IGMesh.h>
#include <GMDS/IG/IGMeshDoctor.h>
#include <GMDS/CAD/GeomEntity.h>
/*----------------------------------------------------------------------------*/
#include <sstream>
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
template<> TInt IGMesh::getNewMark<Node>()
{
	if ( m_nbUsedMarks_nodes==31 )
		// not enough marks
		throw GMDSException("Limit of Boolean marks reached");
#ifdef _DEBUG

	if (m_nbUsedMarks_nodes==m_maxNbUsedMarks_nodes)
		m_maxNbUsedMarks_nodes = m_nbUsedMarks_nodes + 1;
#endif // _DEBUG

	TInt mark = m_marks_nodes[m_nbUsedMarks_nodes++];
	m_usedMarks_nodes.set(mark, true);
	return mark;
}
/*----------------------------------------------------------------------------*/
template<> TInt IGMesh::getNewMark<Edge>()
{
	if ( m_nbUsedMarks_edges==31 )
		// not enough marks
		throw GMDSException("Limit of Boolean marks reached");
#ifdef _DEBUG

	if (m_nbUsedMarks_edges==m_maxNbUsedMarks_edges)
		m_maxNbUsedMarks_edges = m_nbUsedMarks_edges + 1;
#endif // _DEBUG

	TInt mark = m_marks_edges[m_nbUsedMarks_edges++];
	m_usedMarks_edges.set(mark, true);
	return mark;
}
/*----------------------------------------------------------------------------*/
template<> TInt IGMesh::getNewMark<Face>()
{
	if ( m_nbUsedMarks_faces==31 )
		// not enough marks
		throw GMDSException("Limit of Boolean marks reached");
#ifdef _DEBUG

	if (m_nbUsedMarks_faces==m_maxNbUsedMarks_faces)
		m_maxNbUsedMarks_faces = m_nbUsedMarks_faces + 1;
#endif // _DEBUG

	TInt mark = m_marks_faces[m_nbUsedMarks_faces++];
	m_usedMarks_faces.set(mark, true);
	return mark;
}
/*----------------------------------------------------------------------------*/
template<> TInt IGMesh::getNewMark<Region>()
{
	if ( m_nbUsedMarks_regions==31 )
		// not enough marks
		throw GMDSException("Limit of Boolean marks reached");
#ifdef _DEBUG

	if (m_nbUsedMarks_regions==m_maxNbUsedMarks_regions)
		m_maxNbUsedMarks_regions = m_nbUsedMarks_regions + 1;
#endif // _DEBUG

	TInt mark = m_marks_regions[m_nbUsedMarks_regions++];
	m_usedMarks_regions.set(mark, true);
	return mark;
}
/*----------------------------------------------------------------------------*/
template<> void IGMesh::negateMaskMark<Node>(const TInt AMarkNumber)
{
	m_maskMarks_nodes.flip(AMarkNumber);
}
/*----------------------------------------------------------------------------*/
template<> void IGMesh::negateMaskMark<Edge>(const TInt AMarkNumber)
{
	m_maskMarks_edges.flip(AMarkNumber);
}
/*----------------------------------------------------------------------------*/
template<> void IGMesh::negateMaskMark<Face>(const TInt AMarkNumber)
{
	m_maskMarks_faces.flip(AMarkNumber);
}
/*----------------------------------------------------------------------------*/
template<> void IGMesh::negateMaskMark<Region>(const TInt AMarkNumber)
{
	m_maskMarks_regions.flip(AMarkNumber);
}
/*----------------------------------------------------------------------------*/
template<> bool IGMesh::isMarked<Node>(const TCellID& ACellID, int AMarkNumber) const
{
	return (*m_marks[0])[ACellID][AMarkNumber]!=m_maskMarks_nodes[AMarkNumber];
}
/*----------------------------------------------------------------------------*/
bool IGMesh::isMarked(const Node& ACell, int AMarkNumber) const
{
	return isMarked<Node>(ACell.getID(),AMarkNumber);
}
/*----------------------------------------------------------------------------*/
template<> bool IGMesh::
isMarked<Edge>(const TCellID& ACellID, int AMarkNumber) const
{
	return (*m_marks[1])[ACellID][AMarkNumber]!=m_maskMarks_edges[AMarkNumber];
}
/*----------------------------------------------------------------------------*/
bool IGMesh::isMarked(const Edge& ACell, int AMarkNumber) const
{
	return isMarked<Edge>(ACell.getID(),AMarkNumber);
}
/*----------------------------------------------------------------------------*/
template<> bool IGMesh::
isMarked<Face>(const TCellID& ACellID, int AMarkNumber) const
{
	return (*m_marks[2])[ACellID][AMarkNumber]!=m_maskMarks_faces[AMarkNumber];
}
/*----------------------------------------------------------------------------*/
bool IGMesh::isMarked(const Face& ACell, int AMarkNumber) const
{
	return isMarked<Face>(ACell.getID(),AMarkNumber);
}
/*----------------------------------------------------------------------------*/
template<> bool IGMesh::
isMarked<Region>(const TCellID& ACellID, int AMarkNumber) const
{
	return (*m_marks[3])[ACellID][AMarkNumber]!=m_maskMarks_regions[AMarkNumber];
}
/*----------------------------------------------------------------------------*/
bool IGMesh::isMarked(const Region& ACell, int AMarkNumber) const
{
	return isMarked<Region>(ACell.getID(),AMarkNumber);
}
/*----------------------------------------------------------------------------*/
template<> void IGMesh::
markTo<Node>(const TCellID& ACellID, int AMarkNumber, bool AState)
{
	(*m_marks[0])[ACellID].set(AMarkNumber, AState ^ m_maskMarks_nodes[AMarkNumber]);
}
/*----------------------------------------------------------------------------*/
void IGMesh::markTo(const Node& ACell, int AMarkNumber, bool AState)
{
	markTo<Node>(ACell.getID(), AMarkNumber, AState);
}
/*----------------------------------------------------------------------------*/
template<> void IGMesh::
markTo<Edge>(const TCellID& ACellID, int AMarkNumber, bool AState)
{
	(*m_marks[1])[ACellID].set(AMarkNumber, AState ^ m_maskMarks_edges[AMarkNumber]);
}
/*----------------------------------------------------------------------------*/
void IGMesh::markTo(const Edge& ACell, int AMarkNumber, bool AState)
{
	markTo<Edge>(ACell.getID(), AMarkNumber, AState);
}
/*----------------------------------------------------------------------------*/
template<> void IGMesh::
markTo<Face>(const TCellID& ACellID, int AMarkNumber, bool AState)
{
	(*m_marks[2])[ACellID].set(AMarkNumber, AState ^ m_maskMarks_faces[AMarkNumber]);
}
/*----------------------------------------------------------------------------*/
void IGMesh::markTo(const Face& ACell, int AMarkNumber, bool AState)
{
	markTo<Face>(ACell.getID(), AMarkNumber, AState);
}
/*----------------------------------------------------------------------------*/
template<> void IGMesh::
markTo<Region>(const TCellID& ACellID, int AMarkNumber, bool AState)
{
	(*m_marks[3])[ACellID].set(AMarkNumber, AState ^ m_maskMarks_regions[AMarkNumber]);
}
/*----------------------------------------------------------------------------*/
void IGMesh::markTo(const Region& ACell, int AMarkNumber, bool AState)
{
	markTo<Region>(ACell.getID(), AMarkNumber, AState);
}
/*----------------------------------------------------------------------------*/
void IGMesh::mark(const Node& ACell, int AMarkNumber)
{
	markTo(ACell,AMarkNumber,true);
}
/*----------------------------------------------------------------------------*/
template<> void IGMesh::
mark<Node>(const TCellID& ACellID, int AMarkNumber)
{
	markTo<Node>(ACellID,AMarkNumber,true);
}
/*----------------------------------------------------------------------------*/
void IGMesh::mark(const Edge& ACell, int AMarkNumber)
{
	markTo(ACell,AMarkNumber,true);
}
/*----------------------------------------------------------------------------*/
template<> void IGMesh::
mark<Edge>(const TCellID& ACellID, int AMarkNumber)
{
	markTo<Edge>(ACellID,AMarkNumber,true);
}
/*----------------------------------------------------------------------------*/
void IGMesh::mark(const Face& ACell, int AMarkNumber)
{
	markTo(ACell,AMarkNumber,true);
}
/*----------------------------------------------------------------------------*/
template<> void IGMesh::
mark<Face>(const TCellID& ACellID, int AMarkNumber)
{
	markTo<Face>(ACellID,AMarkNumber,true);
}
/*----------------------------------------------------------------------------*/
void IGMesh::mark(const Region& ACell, int AMarkNumber)
{
	markTo(ACell,AMarkNumber,true);
}
/*----------------------------------------------------------------------------*/
template<> void IGMesh::
mark<Region>(const TCellID& ACellID, int AMarkNumber)
{
	markTo<Region>(ACellID,AMarkNumber,true);
}
/*----------------------------------------------------------------------------*/
void IGMesh::unmark(const Node& ACell, int AMarkNumber)
{
	markTo(ACell,AMarkNumber,false);
}
/*----------------------------------------------------------------------------*/
template<> void IGMesh::
unmark<Node>(const TCellID& ACellID, int AMarkNumber)
{
	markTo<Node>(ACellID,AMarkNumber,false);
}
/*----------------------------------------------------------------------------*/
void IGMesh::unmark(const Edge& ACell, int AMarkNumber)
{
	markTo(ACell,AMarkNumber,false);
}
/*----------------------------------------------------------------------------*/
template<> void IGMesh::
unmark<Edge>(const TCellID& ACellID, int AMarkNumber)
{
	markTo<Edge>(ACellID,AMarkNumber,false);
}
/*----------------------------------------------------------------------------*/
void IGMesh::unmark(const Face& ACell, int AMarkNumber)
{
	markTo(ACell,AMarkNumber,false);
}
/*----------------------------------------------------------------------------*/
template<> void IGMesh::
unmark<Face>(const TCellID& ACellID, int AMarkNumber)
{
	markTo<Face>(ACellID,AMarkNumber,false);
}
/*----------------------------------------------------------------------------*/
void IGMesh::unmark(const Region& ACell, int AMarkNumber)
{
	markTo(ACell,AMarkNumber,false);
}
/*----------------------------------------------------------------------------*/
template<> void IGMesh::
unmark<Region>(const TCellID& ACellID, int AMarkNumber)
{
	markTo<Region>(ACellID,AMarkNumber,false);
}
/*----------------------------------------------------------------------------*/
template<> void IGMesh::unmarkAll<Node>(const TInt AMarkNumber)
{
	node_iterator it = nodes_begin();
	for(;!it.isDone();it.next())
		unmark(it.value(),AMarkNumber);
}
/*----------------------------------------------------------------------------*/
template<> void IGMesh::unmarkAll<Edge>(const TInt AMarkNumber)
{
	edge_iterator it = edges_begin();
	for(;!it.isDone();it.next())
		unmark(it.value(),AMarkNumber);
}
/*----------------------------------------------------------------------------*/
template<> void IGMesh::unmarkAll<Face>(const TInt AMarkNumber)
{
	face_iterator it = faces_begin();
	for(;!it.isDone();it.next())
		unmark(it.value(),AMarkNumber);
}
/*----------------------------------------------------------------------------*/
template<> void IGMesh::unmarkAll<Region>(const TInt AMarkNumber)
{
	region_iterator it = regions_begin();
	for(;!it.isDone();it.next())
		unmark(it.value(),AMarkNumber);
}
/*----------------------------------------------------------------------------*/
template<> void IGMesh::freeMark<Node>(int AMarkNumber)
{
	m_usedMarks_nodes.set(AMarkNumber, false);
	m_marks_nodes[-- m_nbUsedMarks_nodes] = AMarkNumber;
}
/*----------------------------------------------------------------------------*/
template<> void IGMesh::freeMark<Edge>(int AMarkNumber)
{
	m_usedMarks_edges.set(AMarkNumber, false);
	m_marks_edges[-- m_nbUsedMarks_edges] = AMarkNumber;
}
/*----------------------------------------------------------------------------*/
template<> void IGMesh::freeMark<Face>(int AMarkNumber)
{
	m_usedMarks_faces.set(AMarkNumber, false);
	m_marks_faces[-- m_nbUsedMarks_faces] = AMarkNumber;
}
/*----------------------------------------------------------------------------*/
template<> void IGMesh::freeMark<Region>(int AMarkNumber)
{
	m_usedMarks_regions.set(AMarkNumber, false);
	m_marks_regions[-- m_nbUsedMarks_regions] = AMarkNumber;
}
/*----------------------------------------------------------------------------*/
template<> Node IGMesh::get<Node>(const TCellID& AID)
{
	return m_nodes_container->buildNode(AID);
}
/*----------------------------------------------------------------------------*/
template<> Edge IGMesh::get<Edge>(const TCellID& AID)
{
	return m_edges_container->buildEdge(AID);
}
/*----------------------------------------------------------------------------*/
template<> Face IGMesh::get<Face>(const TCellID& AID)
{
	return m_faces_container->buildFace(AID);
}
/*----------------------------------------------------------------------------*/
template<> Region IGMesh::get<Region>(const TCellID& AID)
{
	return m_regions_container->buildRegion(AID);
}
/*----------------------------------------------------------------------------*/
math::Point IGMesh::getPoint(const TCellID& AID)
{
        return m_nodes_container->getPoint(AID);
}
/*----------------------------------------------------------------------------*/
template<> bool IGMesh::has<Node>(const TCellID& AID)
{
	return m_nodes_container->has(AID);
}
/*----------------------------------------------------------------------------*/
template<> bool IGMesh::has<Edge>(const TCellID& AID)
{
	return m_edges_container->has(AID);
}
/*----------------------------------------------------------------------------*/
template<> bool IGMesh::has<Face>(const TCellID& AID)
{
	return m_faces_container->has(AID);
}
/*----------------------------------------------------------------------------*/
template<> bool IGMesh::has<Region>(const TCellID& AID)
{
	return m_regions_container->has(AID);
}
/*----------------------------------------------------------------------------*/
template<> void IGMesh::getAll<Node>(std::vector<Node>& AVec)
{
	TInt nb_cells = getNbNodes();
	AVec.clear();
	AVec.resize(nb_cells);
	int i=0;
	for(node_iterator it = nodes_begin();!it.isDone();it.next())
		AVec[i++] = it.value();

}
/*----------------------------------------------------------------------------*/
template<> void IGMesh::getAll<Edge>(std::vector<Edge>& AVec)
{
	TInt nb_cells = getNbEdges();
	AVec.clear();
	AVec.resize(nb_cells);
	int i=0;
	for(edge_iterator it = edges_begin();!it.isDone();it.next())
		AVec[i++] = it.value();

}
/*----------------------------------------------------------------------------*/
template<> void IGMesh::getAll<Face>(std::vector<Face>& AVec)
{
	TInt nb_cells = getNbFaces();
	AVec.clear();
	AVec.resize(nb_cells);
	int i=0;
	for(face_iterator it = faces_begin();!it.isDone();it.next())
		AVec[i++] = it.value();
}
/*----------------------------------------------------------------------------*/
template<> void IGMesh::getAll<Region>(std::vector<Region>& AVec)
{
	TInt nb_cells = getNbRegions();
	AVec.clear();
	AVec.resize(nb_cells);
	int i=0;
	for(region_iterator it = regions_begin();!it.isDone();it.next())
		AVec[i++] = it.value();
}
/*----------------------------------------------------------------------------*/
template<> Marks32 IGMesh::getMarks<Node>(const Node& ACell)
{
  return (*m_marks[0])[ACell.getID()];
}
/*----------------------------------------------------------------------------*/
template<> Marks32 IGMesh::getMarks<Edge>(const Edge& ACell)
{
  return (*m_marks[1])[ACell.getID()];
}
/*----------------------------------------------------------------------------*/
template<> Marks32 IGMesh::getMarks<Face>(const Face& ACell)
{
  return (*m_marks[2])[ACell.getID()];
}
/*----------------------------------------------------------------------------*/
template<> Marks32 IGMesh::getMarks<Region>(const Region& ACell)
{
  return (*m_marks[3])[ACell.getID()];
}
/*----------------------------------------------------------------------------*/
IGMesh::IGMesh(MeshModel model)
: m_model(model)
{
	m_nodes_container   = new NodeContainer  (this);
	m_edges_container   = new EdgeContainer  (this);
	m_faces_container   = new FaceContainer  (this);
	m_regions_container = new RegionContainer(this);


	/* init all the bits to false*/
	m_maskMarks_nodes.reset();
	m_maskMarks_edges.reset();
	m_maskMarks_faces.reset();
	m_maskMarks_regions.reset();
	m_usedMarks_nodes.reset();
	m_usedMarks_edges.reset();
	m_usedMarks_faces.reset();
	m_usedMarks_regions.reset();

	m_nbUsedMarks_nodes = 0;
	m_nbUsedMarks_edges = 0;
	m_nbUsedMarks_faces = 0;
	m_nbUsedMarks_regions = 0;
#ifdef _DEBUG
	m_maxNbUsedMarks_nodes = 0;
	m_maxNbUsedMarks_edges = 0;
	m_maxNbUsedMarks_faces = 0;
	m_maxNbUsedMarks_regions = 0;
#endif // _DEBUG


	m_marks[0] = newVariable<Marks32>(GMDS_NODE  , "mark");
	m_marks[1] = newVariable<Marks32>(GMDS_EDGE  , "mark");
	m_marks[2] = newVariable<Marks32>(GMDS_FACE  , "mark");
	m_marks[3] = newVariable<Marks32>(GMDS_REGION, "mark");

	for (int i=0; i<32; ++i){
		m_marks_nodes  [i] = i;
		m_marks_edges  [i] = i;
		m_marks_faces  [i] = i;
		m_marks_regions[i] = i;
	}
}
/*----------------------------------------------------------------------------*/
IGMesh::~IGMesh()
{
	if(m_nodes_container)
		delete m_nodes_container;
	if(m_edges_container)
		delete m_edges_container;
	if(m_faces_container)
		delete m_faces_container;
	if(m_regions_container)
		delete m_regions_container;
}
/*----------------------------------------------------------------------------*/
MeshModel IGMesh::getModel()const
{
	return m_model;
}
/*----------------------------------------------------------------------------*/
void IGMesh::clear()
{
	if(m_nodes_container)
		m_nodes_container->clear();
	if(m_edges_container)
		m_edges_container->clear();
	if(m_faces_container)
		m_faces_container->clear();
	if(m_regions_container)
		m_regions_container->clear();

	m_node_variable_manager.clearVariables();
	m_edge_variable_manager.clearVariables();
	m_face_variable_manager.clearVariables();
	m_region_variable_manager.clearVariables();

	m_clouds.clear();
	m_lines.clear();
	m_surfaces.clear();
	m_volumes.clear();
}
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void IGMesh::changeModel(const MeshModel& AModel, const bool& ACallDoctor)
{
	if(ACallDoctor) {
		changeModelWithDoctor(AModel);
	} else {
		changeModelWithoutDoctor(AModel);	
	}
}
/*----------------------------------------------------------------------------*/
void IGMesh::changeModelWithoutDoctor(const MeshModel& AModel)
{
	MeshModel model_old = m_model;
	MeshModel model_new = AModel;

	/* the reference mesh model must be changed now since base cell classes
	 * (as the Node or Face ones for instance), will check some model
	 * properties before applying an operation.
	 */
	m_model = AModel;
	m_nodes_container->setModel(m_model);
	m_edges_container->setModel(m_model);
	m_faces_container->setModel(m_model);
	m_regions_container->setModel(m_model);
	// We look for the cells and relationships that differ between the 2 models
	MeshModel model_oldBUTnew = MeshModel::exclusion	(model_old, model_new);
	MeshModel model_newBUTold = MeshModel::exclusion	(model_new, model_old);

	/* Containers for cells of any dimension remain available all the time
	 * (see m_X_container), but adjacency and incidence relations X2Y must be
	 * added or removed
	*/

	/* We add the necesseray containers
	 */
	if(model_newBUTold.has(N2N))
		m_nodes_container->addConnectivityContainers(0);
	if(model_newBUTold.has(N2E))
		m_nodes_container->addConnectivityContainers(1);
	if(model_newBUTold.has(N2F))
		m_nodes_container->addConnectivityContainers(2);
	if(model_newBUTold.has(N2R))
		m_nodes_container->addConnectivityContainers(3);

	if(model_newBUTold.has(E2N))
		m_edges_container->addConnectivityContainers(0);
	if(model_newBUTold.has(E2E))
		m_edges_container->addConnectivityContainers(1);
	if(model_newBUTold.has(E2F))
		m_edges_container->addConnectivityContainers(2);
	if(model_newBUTold.has(E2R))
		m_edges_container->addConnectivityContainers(3);

	if(model_newBUTold.has(F2N))
		m_faces_container->addConnectivityContainers(0);
	if(model_newBUTold.has(F2E))
		m_faces_container->addConnectivityContainers(1);
	if(model_newBUTold.has(F2F))
		m_faces_container->addConnectivityContainers(2);
	if(model_newBUTold.has(F2R))
		m_faces_container->addConnectivityContainers(3);

	if(model_newBUTold.has(R2N))
		m_regions_container->addConnectivityContainers(0);
	if(model_newBUTold.has(R2E))
		m_regions_container->addConnectivityContainers(1);
	if(model_newBUTold.has(R2F))
		m_regions_container->addConnectivityContainers(2);
	if(model_newBUTold.has(R2R))
		m_regions_container->addConnectivityContainers(3);

	//REMOVING
	if(model_oldBUTnew.has(N)) {
		m_nodes_container->clear();
	}
	if(model_oldBUTnew.has(N2N))
		m_nodes_container->removeConnectivityContainers(0);
	if(model_oldBUTnew.has(N2E))
		m_nodes_container->removeConnectivityContainers(1);
	if(model_oldBUTnew.has(N2F))
		m_nodes_container->removeConnectivityContainers(2);
	if(model_oldBUTnew.has(N2R))
		m_nodes_container->removeConnectivityContainers(3);

	if(model_oldBUTnew.has(E)) {
		m_edges_container->clear();	
	}
	if(model_oldBUTnew.has(E2N))
		m_edges_container->removeConnectivityContainers(0);
	if(model_oldBUTnew.has(E2E))
		m_edges_container->removeConnectivityContainers(1);
	if(model_oldBUTnew.has(E2F))
		m_edges_container->removeConnectivityContainers(2);
	if(model_oldBUTnew.has(E2R))
		m_edges_container->removeConnectivityContainers(3);

	if(model_oldBUTnew.has(F)) {
		m_faces_container->clear();
	}
	if(model_oldBUTnew.has(F2N))
		m_faces_container->removeConnectivityContainers(0);
	if(model_oldBUTnew.has(F2E))
		m_faces_container->removeConnectivityContainers(1);
	if(model_oldBUTnew.has(F2F))
		m_faces_container->removeConnectivityContainers(2);
	if(model_oldBUTnew.has(F2R))
		m_faces_container->removeConnectivityContainers(3);

	if(model_oldBUTnew.has(R)) {
                m_regions_container->clear();
        }
	if(model_oldBUTnew.has(R2N))
		m_regions_container->removeConnectivityContainers(0);
	if(model_oldBUTnew.has(R2E))
		m_regions_container->removeConnectivityContainers(1);
	if(model_oldBUTnew.has(R2F))
		m_regions_container->removeConnectivityContainers(2);
	if(model_oldBUTnew.has(R2R))
		m_regions_container->removeConnectivityContainers(3);

}
/*----------------------------------------------------------------------------*/
void IGMesh::changeModelWithDoctor(const MeshModel& AModel)
{
	MeshModel model_old = m_model;
	MeshModel model_new = AModel;

	/* the reference mesh model must be changed now since base cell classes
	 * (as the Node or Face ones for instance), will check some model
	 * properties before applying an operation.
	 */
	m_model = AModel;
	m_nodes_container->setModel(m_model);
	m_edges_container->setModel(m_model);
	m_faces_container->setModel(m_model);
	m_regions_container->setModel(m_model);
	// We look for the cells and relationships that differ between the 2 models
	MeshModel model_oldBUTnew = MeshModel::exclusion	(model_old, model_new);
	MeshModel model_newBUTold = MeshModel::exclusion	(model_new, model_old);

	/* Containers for cells of any dimension remain available all the time
	 * (see m_X_container), but adjacency and incidence relations X2Y must be
	 * added or removed
	*/
	IGMeshDoctor doc(this);

	/* We add the necesseray containers
	 */
	if(model_newBUTold.has(N2N))
		m_nodes_container->addConnectivityContainers(0);
	if(model_newBUTold.has(N2E))
		m_nodes_container->addConnectivityContainers(1);
	if(model_newBUTold.has(N2F))
		m_nodes_container->addConnectivityContainers(2);
	if(model_newBUTold.has(N2R))
		m_nodes_container->addConnectivityContainers(3);

	if(model_newBUTold.has(E2N))
		m_edges_container->addConnectivityContainers(0);
	if(model_newBUTold.has(E2E))
		m_edges_container->addConnectivityContainers(1);
	if(model_newBUTold.has(E2F))
		m_edges_container->addConnectivityContainers(2);
	if(model_newBUTold.has(E2R))
		m_edges_container->addConnectivityContainers(3);

	if(model_newBUTold.has(F2N))
		m_faces_container->addConnectivityContainers(0);
	if(model_newBUTold.has(F2E))
		m_faces_container->addConnectivityContainers(1);
	if(model_newBUTold.has(F2F))
		m_faces_container->addConnectivityContainers(2);
	if(model_newBUTold.has(F2R))
		m_faces_container->addConnectivityContainers(3);

	if(model_newBUTold.has(R2N))
		m_regions_container->addConnectivityContainers(0);
	if(model_newBUTold.has(R2E))
		m_regions_container->addConnectivityContainers(1);
	if(model_newBUTold.has(R2F))
		m_regions_container->addConnectivityContainers(2);
	if(model_newBUTold.has(R2R))
		m_regions_container->addConnectivityContainers(3);

	// Now missing cells are created under the assumption that the X2N
	// adjacencies exists for X= R, F or E
	if(model_newBUTold.has(F))
		doc.buildF();
	if(model_newBUTold.has(E))
		doc.buildE();

	//Data in newBUTold must be added, then those in old_butNew will be removed
	//The addition is done beginning with downard ajacencies then upward.

	model_old = model_new;
	if(model_newBUTold.has(F2E))
		doc.buildF2E(model_old);

	if(model_newBUTold.has(R2N))
		doc.buildR2N(model_old);
	if(model_newBUTold.has(R2E))
		doc.buildR2E(model_old);
	if(model_newBUTold.has(R2F))
		doc.buildR2F(model_old);

	if(model_newBUTold.has(N2N))
		doc.buildN2N(model_old);
	if(model_newBUTold.has(N2E))
		doc.buildN2E(model_old);
	if(model_newBUTold.has(N2F))
		doc.buildN2F(model_old);
	if(model_newBUTold.has(N2R))
		doc.buildN2R(model_old);

	if(model_newBUTold.has(E2E))
		doc.buildE2E(model_old);
	if(model_newBUTold.has(E2F))
		doc.buildE2F(model_old);
	if(model_newBUTold.has(E2R))
		doc.buildE2R(model_old);

	if(model_newBUTold.has(F2F))
		doc.buildF2F(model_old);
	if(model_newBUTold.has(F2R))
		doc.buildF2R(model_old);

	if(model_newBUTold.has(R2R))
		doc.buildR2R(model_old);




	//REMOVING
	if(model_oldBUTnew.has(N)) {
		m_nodes_container->clear();
	}
	if(model_oldBUTnew.has(N2N))
		m_nodes_container->removeConnectivityContainers(0);
	if(model_oldBUTnew.has(N2E))
		m_nodes_container->removeConnectivityContainers(1);
	if(model_oldBUTnew.has(N2F))
		m_nodes_container->removeConnectivityContainers(2);
	if(model_oldBUTnew.has(N2R))
		m_nodes_container->removeConnectivityContainers(3);

	if(model_oldBUTnew.has(E)) {
                m_edges_container->clear();
        }
	if(model_oldBUTnew.has(E2N))
		m_edges_container->removeConnectivityContainers(0);
	if(model_oldBUTnew.has(E2E))
		m_edges_container->removeConnectivityContainers(1);
	if(model_oldBUTnew.has(E2F))
		m_edges_container->removeConnectivityContainers(2);
	if(model_oldBUTnew.has(E2R))
		m_edges_container->removeConnectivityContainers(3);

	if(model_oldBUTnew.has(F)) {
                m_faces_container->clear();
        }
	if(model_oldBUTnew.has(F2N))
		m_faces_container->removeConnectivityContainers(0);
	if(model_oldBUTnew.has(F2E))
		m_faces_container->removeConnectivityContainers(1);
	if(model_oldBUTnew.has(F2F))
		m_faces_container->removeConnectivityContainers(2);
	if(model_oldBUTnew.has(F2R))
		m_faces_container->removeConnectivityContainers(3);

	if(model_oldBUTnew.has(R)) {
                m_regions_container->clear();
        }
	if(model_oldBUTnew.has(R2N))
		m_regions_container->removeConnectivityContainers(0);
	if(model_oldBUTnew.has(R2E))
		m_regions_container->removeConnectivityContainers(1);
	if(model_oldBUTnew.has(R2F))
		m_regions_container->removeConnectivityContainers(2);
	if(model_oldBUTnew.has(R2R))
		m_regions_container->removeConnectivityContainers(3);

}
/*----------------------------------------------------------------------------*/
Node IGMesh::newNode(const TCoord& AX, const TCoord& AY, const TCoord AZ)
{
	Node n = m_nodes_container->add(AX,AY,AZ);
	m_node_variable_manager.addEntry(n.getID());
	return n;
}
/*----------------------------------------------------------------------------*/
Node IGMesh::newNode(const math::Point& APnt)
{
	return newNode(APnt.X(),APnt.Y(),APnt.Z());
}
/*----------------------------------------------------------------------------*/
Edge IGMesh::newEdge(const Node& AN1, const Node& AN2)
{
	return newEdge(AN1.getID(),AN2.getID());
}
/*----------------------------------------------------------------------------*/
Edge IGMesh::newEdge(const TCellID& AN1, const TCellID& AN2)
{
	Edge e = m_edges_container->add(AN1,AN2);
	m_edge_variable_manager.addEntry(e.getID());
	return e;
}
/*----------------------------------------------------------------------------*/
Face IGMesh::newTriangle(const Node& AN1, const Node& AN2, const Node& AN3)
{
	return newTriangle(AN1.getID(),AN2.getID(),AN3.getID());
}
/*----------------------------------------------------------------------------*/
Face IGMesh::
newTriangle(const TCellID& AN1, const TCellID& AN2, const TCellID& AN3)
{
	Face f = m_faces_container->addTriangle(AN1,AN2,AN3);
	m_face_variable_manager.addEntry(f.getID());
	return f;
}
/*----------------------------------------------------------------------------*/
Face IGMesh::
newQuad(const Node& AN1, const Node& AN2, const Node& AN3, const Node& AN4)
{
	return newQuad(AN1.getID(),AN2.getID(),AN3.getID(),AN4.getID());

}
/*----------------------------------------------------------------------------*/
Face IGMesh::
newQuad(const TCellID& AN1, const TCellID& AN2,
		const TCellID& AN3, const TCellID& AN4)
{
	Face f = m_faces_container->addQuad(AN1,AN2,AN3,AN4);
	m_face_variable_manager.addEntry(f.getID());
	return f;
}
/*----------------------------------------------------------------------------*/
Face IGMesh::
newPolygon(const std::vector<Node>& ANodes)
{
	std::vector<TCellID> ids;
	ids.resize(ANodes.size());
	for(unsigned int i=0;i<ANodes.size();i++)
		ids[i]=ANodes[i].getID();

	return newPolygon(ids);
}
/*----------------------------------------------------------------------------*/
Face IGMesh::
newPolygon(const std::vector<TCellID>& ANodes)
{
	Face f = m_faces_container->addPolygon(ANodes);
	m_face_variable_manager.addEntry(f.getID());
	return f;

}
/*----------------------------------------------------------------------------*/
Face IGMesh::
newFace(const std::vector<Node>& ANodes)
{
	if(ANodes.size()==3) {
		return newTriangle(ANodes[0].getID(),ANodes[1].getID(),ANodes[2].getID());
	} else if(ANodes.size()==4) {
		return newQuad(ANodes[0].getID(),ANodes[1].getID(),ANodes[2].getID(),ANodes[3].getID());
	} else {
		std::vector<TCellID> ids(ANodes.size());
		for(unsigned int iNode=0; iNode<ANodes.size(); iNode++) {
			ids[iNode] = ANodes[iNode].getID();	
		}
		return newPolygon(ids);
	}
}
/*----------------------------------------------------------------------------*/
Face IGMesh::
newFace(const std::vector<TCellID>& ANodes)
{
	Face f;
	if(ANodes.size()==3)
		f = newTriangle(ANodes[0],ANodes[1],ANodes[2]);
	else if(ANodes.size()==4)
		f = newQuad(ANodes[0],ANodes[1],ANodes[2],ANodes[3]);
	else
		f = newPolygon(ANodes);

	return f;
}
/*----------------------------------------------------------------------------*/
Region IGMesh::newTet(const Node& AN1, const Node& AN2, const Node& AN3,
		const Node& AN4)
{
	return newTet(AN1.getID(),AN2.getID(),AN3.getID(),AN4.getID());
}
/*----------------------------------------------------------------------------*/
Region IGMesh::newTet(const TCellID& AN1, const TCellID& AN2,
		const TCellID& AN3, const TCellID& AN4)
{
	Region r = m_regions_container->addTet(AN1,AN2,AN3,AN4);
	m_region_variable_manager.addEntry(r.getID());
	return r;
}
/*----------------------------------------------------------------------------*/
Region IGMesh::newPyramid(const Node& AN1, const Node& AN2, const Node& AN3,
		const Node& AN4, const Node& AN5)
{
	return newPyramid(AN1.getID(),AN2.getID(),AN3.getID(),AN4.getID(), AN5.getID());
}
/*----------------------------------------------------------------------------*/
Region IGMesh::newPyramid(const TCellID& AN1, const TCellID& AN2,
		const TCellID& AN3, const TCellID& AN4,
		const TCellID& AN5)
{
	Region r = m_regions_container->addPyramid(AN1,AN2,AN3,AN4,AN5);
	m_region_variable_manager.addEntry(r.getID());
	return r;
}
/*----------------------------------------------------------------------------*/
Region IGMesh::newPrism3(const Node& AN1, const Node& AN2, const Node& AN3,
		const Node& AN4, const Node& AN5, const Node& AN6)
{
	return newPrism3(
			AN1.getID(), AN2.getID(), AN3.getID(),
			AN4.getID(), AN5.getID(), AN6.getID());
}
/*----------------------------------------------------------------------------*/
Region IGMesh::newPrism3(const TCellID& AN1, const TCellID& AN2,
		const TCellID& AN3, const TCellID& AN4,
		const TCellID& AN5, const TCellID& AN6)
{
	Region r = m_regions_container->addPrism3(AN1,AN2,AN3,AN4,AN5,AN6);
	m_region_variable_manager.addEntry(r.getID());
	return r;
}
/*----------------------------------------------------------------------------*/
Region IGMesh::newHex(
		const Node& AN1, const Node& AN2, const Node& AN3, const Node& AN4,
		const Node& AN5, const Node& AN6, const Node& AN7, const Node& AN8)
{
	return newHex(AN1.getID(),AN2.getID(),AN3.getID(),AN4.getID(),
			AN5.getID(),AN6.getID(),AN7.getID(),AN8.getID());
}
/*----------------------------------------------------------------------------*/
Region IGMesh::newHex(
		const TCellID& AN1, const TCellID& AN2, const TCellID& AN3, const TCellID& AN4,
		const TCellID& AN5, const TCellID& AN6, const TCellID& AN7, const TCellID& AN8)
{
	Region r = m_regions_container->addHex(AN1,AN2,AN3,AN4,AN5,AN6,AN7,AN8);
	m_region_variable_manager.addEntry(r.getID());
	return r;
}
/*----------------------------------------------------------------------------*/
TCellID IGMesh::getMaxLocalID(const TInt& ADim) const {
	TInt max_id=0;
	if(ADim==0)
		max_id = m_nodes_container->getMaxID();
	else if(ADim==1)
		max_id = m_edges_container->getMaxID();
	else if(ADim==2)
		max_id = m_faces_container->getMaxID();
	else if(ADim==3)
		max_id = m_regions_container->getMaxID();

	return max_id;
}
/*----------------------------------------------------------------------------*/
void IGMesh::deleteVariable(ECellType AType, const std::string& AName)
{
	switch(AType){
	case GMDS_NODE:
		m_node_variable_manager.deleteVariable(AName);
		break;
	case GMDS_EDGE:
		m_edge_variable_manager.deleteVariable(AName);
		break;
	case GMDS_FACE:
		m_face_variable_manager.deleteVariable(AName);
		break;
	case GMDS_REGION:
		m_region_variable_manager.deleteVariable(AName);
		break;
	default:
		throw GMDSException("Unmanaged type of cell -> impossible to delete a variable");
	}
}
/*----------------------------------------------------------------------------*/
void IGMesh::deleteVariable(ECellType AType, VariableItf* AVar)
{
	switch(AType){
	case GMDS_NODE:
		m_node_variable_manager.deleteVariable(AVar);
		break;
	case GMDS_EDGE:
		m_edge_variable_manager.deleteVariable(AVar);
		break;
	case GMDS_FACE:
		m_face_variable_manager.deleteVariable(AVar);
		break;
	case GMDS_REGION:
		m_region_variable_manager.deleteVariable(AVar);
		break;
	default:
		throw GMDSException("Unmanaged type of cell -> impossible to delete a variable");
	}
}


/*----------------------------------------------------------------------------*/
std::vector<VariableItf*> IGMesh::getAllVariables(ECellType AType)
{
	switch (AType){
	case GMDS_NODE:{
					   return m_node_variable_manager.getAllVariables();
	}
		break;
	case GMDS_EDGE:{
					   return  m_edge_variable_manager.getAllVariables();
	}
		break;
	case GMDS_FACE:{
					   return m_face_variable_manager.getAllVariables();
	}
		break;
	case GMDS_REGION:{
						 return  m_region_variable_manager.getAllVariables();
	}
		break;
	default:
		throw GMDSException("Unmanaged type of cell -> impossible to access to a variable");
	}
}
/*----------------------------------------------------------------------------*/
bool IGMesh::doesVariableExist(ECellType AType, const std::string& AName)
{
	switch(AType){
	case GMDS_NODE:{
		return m_node_variable_manager.doesVariableExist(AName);
	}
	break;
	case GMDS_EDGE:{
		return m_edge_variable_manager.doesVariableExist(AName);
	}
	break;
	case GMDS_FACE:{
		return m_face_variable_manager.doesVariableExist(AName);
	}
	break;
	case GMDS_REGION:{
		return m_region_variable_manager.doesVariableExist(AName);
	}
	break;
	default:
		throw GMDSException("Unmanaged type of cell -> impossible to access to a variable");
	}

	return false;
}
/*----------------------------------------------------------------------------*/
void IGMesh::clearAndResizeNodeIDContainer(const TInt AMaxID)
{
	m_nodes_container->clear();
	m_nodes_container->resize(AMaxID);
	m_node_variable_manager.clearVariables();
}
/*----------------------------------------------------------------------------*/
void IGMesh::clearAndResizeEdgeIDContainer(const TInt AMaxID)
{
	m_edges_container->clear();
	m_edges_container->resize(AMaxID);
	m_edge_variable_manager.clearVariables();
}
/*----------------------------------------------------------------------------*/
void IGMesh::clearAndResizeFaceIDContainer(const TInt AMaxID)
{
	m_faces_container->clear();
	m_faces_container->resize(AMaxID);
	m_face_variable_manager.clearVariables();
}
/*----------------------------------------------------------------------------*/
void IGMesh::clearAndResizeRegionIDContainer(const TInt AMaxID)
{
	m_regions_container->clear();
	m_regions_container->resize(AMaxID);
	m_region_variable_manager.clearVariables();
}
/*----------------------------------------------------------------------------*/
void IGMesh::updateIDContainers()
{
	m_nodes_container->update();
	m_edges_container->update();
	m_faces_container->update();
	m_regions_container->update();
}
/*----------------------------------------------------------------------------*/
Node IGMesh::newNodeWithID(const TCoord& AX, const TCoord& AY,
		const TCoord& AZ,
		const TCellID& AGID)
{
	Node n = m_nodes_container->add(AX,AY,AZ,AGID);
	m_node_variable_manager.addEntry(AGID);
	return n;
}
/*----------------------------------------------------------------------------*/
Edge IGMesh::
newEdgeWithID(const TCellID& AV1, const TCellID& AV2, const TCellID& AGID)
{
	Edge e = m_edges_container->add(AV1,AV2,AGID);
	m_edge_variable_manager.addEntry(AGID);
	return e;
}
/*----------------------------------------------------------------------------*/
Face IGMesh::
newTriangleWithID(const TCellID& AN1,const TCellID& AN2,
		const TCellID& AN3, const TCellID& AGID)
{
	Face f = m_faces_container->addTriangle(AN1,AN2,AN3,AGID);
	m_face_variable_manager.addEntry(AGID);
	return f;
}
/*----------------------------------------------------------------------------*/
Face IGMesh::
newQuadWithID(const TCellID& AN1, const TCellID& AN2,
		const TCellID& AN3, const TCellID& AN4, const TCellID& AGID)
{
	Face f = m_faces_container->addQuad(AN1,AN2,AN3,AN4,AGID);
	m_face_variable_manager.addEntry(AGID);
	return f;
}
/*----------------------------------------------------------------------------*/
Face IGMesh::
newPolygonWithID(const std::vector<TCellID>& ANodes, const TCellID& AGID)
{
	Face f = m_faces_container->addPolygon(ANodes,AGID);
	m_face_variable_manager.addEntry(AGID);
	return f;
}
/*----------------------------------------------------------------------------*/
Face IGMesh::
newFaceWithID(const std::vector<TCellID>& ANodes, const TCellID& AGID)
{
	Face f;
	if(ANodes.size()==3)
		f = newTriangleWithID(ANodes[0],ANodes[1],ANodes[2],AGID);
	else if(ANodes.size()==4)
		f = newQuadWithID(ANodes[0],ANodes[1],ANodes[2],ANodes[3],AGID);
	else
		f = newPolygonWithID(ANodes,AGID);

	return f;
}
/*----------------------------------------------------------------------------*/
Region IGMesh::
newTetWithID(const TCellID& AN1, const TCellID& AN2,
		const TCellID& AN3, const TCellID& AN4, const TCellID& AGID)
{
	Region r = m_regions_container->addTet(AN1,AN2,AN3,AN4,AGID);
	m_region_variable_manager.addEntry(AGID);
	return r;
}
/*----------------------------------------------------------------------------*/
Region IGMesh::newPyramidWithID(const TCellID& AN1, const TCellID& AN2,
		   const TCellID& AN3, const TCellID& AN4,
		   const TCellID& AN5, const TCellID& AGID)
{
	Region r = m_regions_container->addPyramid(AN1,AN2,AN3,AN4,AN5,AGID);
	m_region_variable_manager.addEntry(AGID);
	return r;
}
/*----------------------------------------------------------------------------*/
Region IGMesh::
newHexWithID(const TCellID& AN1, const TCellID& AN2,
				   const TCellID& AN3, const TCellID& AN4,
				   const TCellID& AN5, const TCellID& AN6,
				   const TCellID& AN7, const TCellID& AN8,
				   const TCellID& AGID)
{
	Region r = m_regions_container->addHex(AN1,AN2,AN3,AN4,AN5,AN6,AN7,AN8,AGID);
	m_region_variable_manager.addEntry(AGID);
	return r;
}
/*----------------------------------------------------------------------------*/
void IGMesh::initializeGeometryClassification()
{
	classification[0] =
	  newVariable<geom::GeomEntity*>(GMDS_NODE,"NodeClassification");
	classification[1] =
	  newVariable<geom::GeomEntity*>(GMDS_EDGE,"EdgeClassification");
	classification[2] =
	  newVariable<geom::GeomEntity*>(GMDS_FACE,"FaceClassification");
	classification[3] =
	  newVariable<geom::GeomEntity*>(GMDS_REGION,"RegionClassification");

	for(int i=0;i<4;i++){
		(classification[i])->setValuesTo(0);
	}
}
/*----------------------------------------------------------------------------*/
bool IGMesh::doesGeometricClassificationExist(const int ADim){

	switch(ADim){
	case 0:{
		return m_node_variable_manager.doesVariableExist("NodeClassification");
	}
	break;
	case 1:{
		return m_edge_variable_manager.doesVariableExist("EdgeClassification");
	}
	break;
	case 2:{
		return m_face_variable_manager.doesVariableExist("FaceClassification");
	}
	break;
	case 3:{
		return m_region_variable_manager.doesVariableExist("RegionClassification");
	}
	break;
	default:
		throw GMDSException("Mesh::doesGeometricClassificationExist : bad ADim");
	}
}
/*----------------------------------------------------------------------------*/
Variable<geom::GeomEntity* >* IGMesh::
getGeometricClassification(const int ADim)
{
	return classification[ADim];
}
/*----------------------------------------------------------------------------*/
geom::GeomEntity* IGMesh::getGeometricClassification(const Node ACell)
{
	return (*(classification[0]))[ACell.getID()];
}
/*----------------------------------------------------------------------------*/
geom::GeomEntity* IGMesh::getGeometricClassification(const Edge ACell)
{
	return (*(classification[1]))[ACell.getID()];
}
/*----------------------------------------------------------------------------*/
geom::GeomEntity* IGMesh::getGeometricClassification(const Face ACell)
{
	return (*(classification[2]))[ACell.getID()];
}
/*----------------------------------------------------------------------------*/
geom::GeomEntity* IGMesh::getGeometricClassification(const Region ACell)
{
	return (*(classification[3]))[ACell.getID()];
}
/*----------------------------------------------------------------------------*/
void  IGMesh::
setGeometricClassification(const Node ACell, geom::GeomEntity* AGE){
	(classification[0])->set(ACell.getID(),AGE);
}
/*----------------------------------------------------------------------------*/
void  IGMesh::
setGeometricClassification(const Edge ACell, geom::GeomEntity* AGE){
	(classification[1])->set(ACell.getID(),AGE);
}
/*----------------------------------------------------------------------------*/
void  IGMesh::
setGeometricClassification(const Face ACell, geom::GeomEntity* AGE){
	(classification[2])->set(ACell.getID(),AGE);
}
/*----------------------------------------------------------------------------*/
void  IGMesh::
setGeometricClassification(const Region ACell, geom::GeomEntity* AGE){
	(classification[3])->set(ACell.getID(),AGE);
}
/*----------------------------------------------------------------------------*/
IGMesh::cloud& IGMesh::newCloud(const std::string& AName)
{
	std::list<cloud>::iterator it = m_clouds.begin();
	for(;it!=m_clouds.end();it++)
		if((*it).name()==AName)
		{
			const std::string mess = "A cloud named "+AName+" already exists !";
			throw GMDSException(mess);
		}

	cloud c(this,AName);
	m_clouds.push_back(c);
	return m_clouds.back();
}
/*----------------------------------------------------------------------------*/
void IGMesh::deleteCloud(cloud& ACloud)
{
	std::list<cloud>::iterator it = m_clouds.begin();
	for(;it!=m_clouds.end();it++)
		if(ACloud==*it){
			m_clouds.erase(it);
			return;
		}
}
/*----------------------------------------------------------------------------*/
IGMesh::cloud& IGMesh::getCloud(const std::string& AName)
{
	std::list<cloud>::iterator it = m_clouds.begin();
	for (; it != m_clouds.end(); it++)
	if ((*it).name() == AName)
		return *it;
    
    std::stringstream mess;
    mess<<"There is no cloud named "<<AName;
	throw GMDSException(mess.str());
	//useless but allows to remove compilation warnings
	return *it;
}
/*----------------------------------------------------------------------------*/
IGMesh::volume& IGMesh::getVolume(const unsigned int AIndex)
{
	if (AIndex >= m_volumes.size())
	{
        std::stringstream mess;
        mess<<"There is no volume indexed "<<AIndex;
        throw GMDSException(mess.str());
	}
	int index = 0;
	std::list<IGMesh::volume>::iterator it = m_volumes.begin();
	for (; it != m_volumes.end(); it++, index++)
	if (index == AIndex)
		return *it;

	//useless but allows to remove compilation warnings
	return *it;
}
/*----------------------------------------------------------------------------*/
IGMesh::cloud& IGMesh::getCloud(const unsigned int AIndex)
{
	if (AIndex >= m_clouds.size())
	{
        std::stringstream mess;
        mess<<"There is no cloud indexed "<<AIndex;
        throw GMDSException(mess.str());
	}
	int index = 0;
	std::list<IGMesh::cloud>::iterator it = m_clouds.begin();
	for (; it != m_clouds.end(); it++, index++)
	if (index == AIndex)
		return *it;

	//useless but allows to remove compilation warnings
	return *it;
}
/*----------------------------------------------------------------------------*/
IGMesh::line& IGMesh::getLine(const unsigned int AIndex)
{
	if (AIndex >= m_lines.size())
	{
        std::stringstream mess;
        mess<<"There is no line indexed "<<AIndex;
        throw GMDSException(mess.str());
	}
	int index = 0;
	std::list<IGMesh::line>::iterator it = m_lines.begin();
	for (; it != m_lines.end(); it++, index++)
	if (index == AIndex)
		return *it;

	//useless but allows to remove compilation warnings
	return *it;
}
/*----------------------------------------------------------------------------*/
IGMesh::surface& IGMesh::getSurface(const unsigned int AIndex)
{
	if (AIndex >= m_surfaces.size())
	{
        std::stringstream mess;
        mess<<"There is no surface indexed "<<AIndex;
        throw GMDSException(mess.str());
	}
	int index = 0;
	std::list<IGMesh::surface >::iterator it = m_surfaces.begin();
	for (; it != m_surfaces.end(); it++, index++)
	if (index == AIndex)
		return *it;

	//useless but allows to remove compilation warnings
	return *it;
}
/*----------------------------------------------------------------------------*/
IGMesh::clouds_iterator IGMesh::clouds_begin()
{
	return m_clouds.begin();
}
/*----------------------------------------------------------------------------*/
IGMesh::clouds_iterator IGMesh::clouds_end()
{
	return m_clouds.end();
}
/*----------------------------------------------------------------------------*/
unsigned int IGMesh::getNbClouds () const
{
	return m_clouds.size();
}
/*----------------------------------------------------------------------------*/
IGMesh::line& IGMesh::newLine(const std::string& AName)
{
	std::list<line>::iterator it = m_lines.begin();
	for(;it!=m_lines.end();it++)
		if((*it).name()==AName)
		{
			const std::string mess = "A line named "+AName+" already exists !";
			throw GMDSException(mess);
		}

	line l(this,AName);
	m_lines.push_back(l);
	return m_lines.back();
}
/*----------------------------------------------------------------------------*/
void IGMesh::deleteLine(line& ALine)
{
	std::list<line>::iterator it = m_lines.begin();
	for(;it!=m_lines.end();it++)
		if(*it==ALine)
		{
			m_lines.erase(it);
			return;
		}
}
/*----------------------------------------------------------------------------*/
IGMesh::line& IGMesh::getLine(const std::string& AName)
{
	std::list<line>::iterator it = m_lines.begin();
	for(;it!=m_lines.end();it++)
		if((*it).name()==AName)
			return *it;

	const std::string mess = "There is no line named "+AName;
	throw GMDSException(mess);
	//useless but allows to remove compilation warnings
	return *it;
}
/*----------------------------------------------------------------------------*/
unsigned int IGMesh::getNbLines() const
{
	return m_lines.size();
}
/*----------------------------------------------------------------------------*/
IGMesh::lines_iterator IGMesh::lines_begin()
{
	return m_lines.begin();
}
/*----------------------------------------------------------------------------*/
IGMesh::lines_iterator IGMesh::lines_end()
{
	return m_lines.end();
}

/*----------------------------------------------------------------------------*/
IGMesh::surface& IGMesh::newSurface(const std::string& AName)
{
	std::list<surface>::iterator it = m_surfaces.begin();
	for(;it!=m_surfaces.end();it++)
		if((*it).name()==AName)
		{
			const std::string mess = "A surface named "+AName+" already exists !";
			throw GMDSException(mess);
		}

	surface s(this,AName);
	m_surfaces.push_back(s);
	return m_surfaces.back();
}
/*----------------------------------------------------------------------------*/
void IGMesh::deleteSurface(surface& ASurf)
{
	std::list<surface>::iterator it = m_surfaces.begin();
	for(;it!=m_surfaces.end();it++)
		if(*it==ASurf)
		{
			m_surfaces.erase(it);
			return;
		}
}
/*----------------------------------------------------------------------------*/
IGMesh::surface& IGMesh::getSurface(const std::string& AName)
{
	std::list<surface>::iterator it = m_surfaces.begin();
	for(;it!=m_surfaces.end();it++)
		if((*it).name()==AName)
			return *it;

	const std::string mess = "There is no surface named "+AName;
	throw GMDSException(mess);
	//useless but allows to remove compilation warnings
	return *it;
}
/*----------------------------------------------------------------------------*/
unsigned int IGMesh::getNbSurfaces () const
{
	return m_surfaces.size();
}
/*----------------------------------------------------------------------------*/
IGMesh::surfaces_iterator IGMesh::surfaces_begin()
{
	return m_surfaces.begin();
}
/*----------------------------------------------------------------------------*/
IGMesh::surfaces_iterator IGMesh::surfaces_end()
{
	return m_surfaces.end();
}
/*----------------------------------------------------------------------------*/
IGMesh::volume& IGMesh::newVolume(const std::string& AName)
{
	std::list<volume>::iterator it = m_volumes.begin();
	for(;it!=m_volumes.end();it++)
		if((*it).name()==AName)
		{
			const std::string mess = "A volume named "+AName+" already exists !";
			throw GMDSException(mess);
		}

	volume v(this,AName);
	m_volumes.push_back(v);
	return m_volumes.back();
}
/*----------------------------------------------------------------------------*/
void IGMesh::deleteVolume(volume& AVol)
{
	std::list<volume>::iterator it = m_volumes.begin();
	for(;it!=m_volumes.end();it++)
		if(*it==AVol)
		{
			m_volumes.erase(it);
			return;
		}
}
/*----------------------------------------------------------------------------*/
IGMesh::volume& IGMesh::getVolume(const std::string& AName)
{
	std::list<volume>::iterator it = m_volumes.begin();
	for(;it!=m_volumes.end();it++)
		if((*it).name()==AName)
			return *it;

	const std::string mess = "There is no volume named "+AName;
	throw GMDSException(mess);
	//useless but allows to remove compilation warnings
	return *it;
}
/*----------------------------------------------------------------------------*/
unsigned int IGMesh::getNbVolumes () const
{
	return m_volumes.size();
}
/*----------------------------------------------------------------------------*/
IGMesh::volumes_iterator IGMesh::volumes_begin()
{
	return m_volumes.begin();
}
/*----------------------------------------------------------------------------*/
IGMesh::volumes_iterator IGMesh::volumes_end()
{
	return m_volumes.end();
}
/*------------------------------------------------------------------------*/
void IGMesh::serialize(std::ostream& AStr)
{
	int model = m_model.getDef();
	AStr.write((char*)&model,sizeof(int));
	m_nodes_container->serialize(AStr);
	m_edges_container->serialize(AStr);
	m_faces_container->serialize(AStr);
	m_regions_container->serialize(AStr);

	m_node_variable_manager.serialize(AStr);
	m_edge_variable_manager.serialize(AStr);
	m_face_variable_manager.serialize(AStr);
	m_region_variable_manager.serialize(AStr);

	for(std::list<cloud>::iterator it= m_clouds.begin(); it!=m_clouds.end();it++)
		it->serialize(AStr);
	for(std::list<line>::iterator it= m_lines.begin(); it!=m_lines.end();it++)
		it->serialize(AStr);
	for(std::list<surface>::iterator it= m_surfaces.begin(); it!=m_surfaces.end();it++)
		it->serialize(AStr);
	for(std::list<volume>::iterator it= m_volumes.begin(); it!=m_volumes.end();it++)
		it->serialize(AStr);

	/** classification has not to serialized, it is done during variable managers'
	 *  serializations. Idem for m_marks.
	 */
	AStr.write((char*)&m_usedMarks_nodes  , sizeof(int32_t));
	AStr.write((char*)&m_usedMarks_edges  , sizeof(int32_t));
		AStr.write((char*)&m_usedMarks_faces  , sizeof(int32_t));
	AStr.write((char*)&m_usedMarks_regions, sizeof(int32_t));

	AStr.write((char*)&m_maskMarks_nodes  , sizeof(int32_t));
	AStr.write((char*)&m_maskMarks_edges  , sizeof(int32_t));
	AStr.write((char*)&m_maskMarks_faces  , sizeof(int32_t));
	AStr.write((char*)&m_maskMarks_regions, sizeof(int32_t));

	AStr.write((char*)&m_marks_nodes[0]  , 32*sizeof(TInt));
	AStr.write((char*)&m_marks_edges[0]  , 32*sizeof(TInt));
	AStr.write((char*)&m_marks_faces[0]  , 32*sizeof(TInt));
	AStr.write((char*)&m_marks_regions[0], 32*sizeof(TInt));

	AStr.write((char*)&m_nbUsedMarks_nodes  , sizeof(TInt));
	AStr.write((char*)&m_nbUsedMarks_edges  , sizeof(TInt));
	AStr.write((char*)&m_nbUsedMarks_faces  , sizeof(TInt));
	AStr.write((char*)&m_nbUsedMarks_regions, sizeof(TInt));

}
/*------------------------------------------------------------------------*/
void IGMesh:: unserialize(std::istream& AStr) {

	int def_model;
	AStr.read((char*)&def_model, sizeof(int));
	m_model = MeshModel(def_model);
	m_nodes_container->unserialize(AStr);
	m_edges_container->unserialize(AStr);
	m_faces_container->unserialize(AStr);
	m_regions_container->unserialize(AStr);

	m_node_variable_manager.unserialize(AStr);
	m_edge_variable_manager.unserialize(AStr);
	m_face_variable_manager.unserialize(AStr);
	m_region_variable_manager.unserialize(AStr);

	for(std::list<cloud>::iterator it= m_clouds.begin(); it!=m_clouds.end();it++)
		it->unserialize(AStr);
	for(std::list<line>::iterator it= m_lines.begin(); it!=m_lines.end();it++)
		it->unserialize(AStr);
	for(std::list<surface>::iterator it= m_surfaces.begin(); it!=m_surfaces.end();it++)
		it->unserialize(AStr);
	for(std::list<volume>::iterator it= m_volumes.begin(); it!=m_volumes.end();it++)
		it->unserialize(AStr);

	/** classification has not to serialized, it is done during variable managers'
	 *  serializations. Idem for m_marks.
	 */
	AStr.read((char*)&m_usedMarks_nodes  , sizeof(int32_t));
	AStr.read((char*)&m_usedMarks_edges  , sizeof(int32_t));
	AStr.read((char*)&m_usedMarks_faces  , sizeof(int32_t));
	AStr.read((char*)&m_usedMarks_regions, sizeof(int32_t));

	AStr.read((char*)&m_maskMarks_nodes  , sizeof(int32_t));
	AStr.read((char*)&m_maskMarks_edges  , sizeof(int32_t));
	AStr.read((char*)&m_maskMarks_faces  , sizeof(int32_t));
	AStr.read((char*)&m_maskMarks_regions, sizeof(int32_t));

	AStr.read((char*)&m_marks_nodes[0]  , 32*sizeof(TInt));
	AStr.read((char*)&m_marks_edges[0]  , 32*sizeof(TInt));
	AStr.read((char*)&m_marks_faces[0]  , 32*sizeof(TInt));
	AStr.read((char*)&m_marks_regions[0], 32*sizeof(TInt));

	AStr.read((char*)&m_nbUsedMarks_nodes  , sizeof(TInt));
	AStr.read((char*)&m_nbUsedMarks_edges  , sizeof(TInt));
	AStr.read((char*)&m_nbUsedMarks_faces  , sizeof(TInt));
	AStr.read((char*)&m_nbUsedMarks_regions, sizeof(TInt));
}
/*----------------------------------------------------------------------------*/
void IGMesh::
getCommonNodes(const Face& AF1, const Face& AF2, std::vector<Node>& ANodes)
{
	ANodes.clear();
	std::vector<Node> nodes1 = AF1.get<Node>();
	std::vector<Node> nodes2 = AF2.get<Node>();


	bool found;
	for(unsigned int i1 = 0; i1<nodes1.size();i1++)
	{
		found=false;
		for(unsigned int i2 = 0; i2<nodes2.size() && !found;i2++)
			if(nodes1[i1]==nodes2[i2])
			{
				ANodes.push_back(nodes1[i1]);
				found=true;
			}
	}
}
/*----------------------------------------------------------------------------*/
std::vector<Node> IGMesh::
getCommonNodes(const Face& AF1, const Face& AF2)
{
	std::vector<Node> nodes;
	getCommonNodes(AF1,AF2,nodes);
	return nodes;
}
/*----------------------------------------------------------------------------*/
#ifdef GMDS_PARALLEL
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
TInt IGMesh::getPartID() const
{
	return m_part_id;
}
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
#endif //GMDS_PARALLEL
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/
