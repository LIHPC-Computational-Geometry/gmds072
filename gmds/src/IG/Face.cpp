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
 * Face.cpp
 *
 *  Created on: 20 févr. 2014
 *      Author: ledouxf
 */
/*----------------------------------------------------------------------------*/
#include <map>
/*----------------------------------------------------------------------------*/
#include <GMDS/IG/Face.h>
#include <GMDS/IG/FaceContainer.h>
#include <GMDS/IG/IGMesh.h>
#include <GMDS/Utils/Exception.h>
#include <GMDS/Math/Plane.h>
#include <GMDS/Math/Segment.h>
#include <GMDS/Math/Quadrilateral.h>
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
Face::
Face()
: Cell(0,GMDS_FACE,NullID), m_faces_container(0), m_type_id(NullID)
{}
/*----------------------------------------------------------------------------*/
TInt Face::getNbNodes() const
{
	TInt nb = 0;
	m_faces_container->getNodesData(m_id,nb);
	return nb;
}
/*----------------------------------------------------------------------------*/
TInt Face::getNbEdges() const
{
	TInt nb = 0;
	m_faces_container->getEdgesData(m_id,nb);
	return nb;
}
/*----------------------------------------------------------------------------*/
TInt Face::getNbFaces() const
{
	TInt nb = 0;
	m_faces_container->getFacesData(m_id,nb);
	return nb;
}
/*----------------------------------------------------------------------------*/
TInt Face::getNbRegions() const
{
	TInt nb = 0;
	m_faces_container->getRegionsData(m_id,nb);
	return nb;
}
/*----------------------------------------------------------------------------*/
math::Vector Face::normal() const
{
	math::Vector n;
	std::vector<Node> nodes = this->get<Node>();
    int nb_nodes = nodes.size();
    
    if(nb_nodes==3){
        math::Point p1 = nodes[0].getPoint();
        math::Point p2 = nodes[1].getPoint();
        math::Point p3 = nodes[2].getPoint();
        
        math::Vector v1(p1,p2);
        math::Vector v3(p1,p3);
        
        n = (v1.cross(v3));
        n.normalize();
    }
    else if (nb_nodes>3){
        math::Point p = center();
        math::Vector3d n_temp(0,0,0);
        for(auto i=0; i<nb_nodes; i++){
            math::Point pi = nodes[i].getPoint();
            math::Point pj = nodes[(i+1)%nb_nodes].getPoint();
            math::Vector3d vi(p,pi);
            math::Vector3d vj(p,pj);
            math::Vector3d nij = vi.cross(vj);
            nij.normalize();
            n_temp += nij;
        }
        n_temp /= nb_nodes;
        n_temp.normalize();
        n = math::Vector(n_temp.X(), n_temp.Y(), n_temp.Z());
    }
	return n;
}
/*----------------------------------------------------------------------------*/
math::Vector Face::normal(const Node& ANode) const
{
        math::Vector n;
        std::vector<Node> nodes = this->get<Node>();
	
	// first find ANode among the faces nodes
	unsigned int nodeIndex = NullID;

	for(unsigned int iNode=0; iNode<nodes.size(); iNode++) {
		if(ANode == nodes[iNode]) {
			nodeIndex = iNode;
		}
	}

	if(nodeIndex == NullID) {
		throw GMDSException("math::Vector Face::normal could not find ANode.");
	}

        Node n1 = nodes[nodeIndex];
        Node n2 = nodes[(nodeIndex+1)%nodes.size()];
        Node n3 = nodes[(nodeIndex+2)%nodes.size()];

        math::Vector v1(n2.X()-n1.X(), n2.Y()-n1.Y(), n2.Z()-n1.Z());
        math::Vector v3(n3.X()-n1.X(), n3.Y()-n1.Y(), n3.Z()-n1.Z());

        n = (v1.cross(v3));
        n.normalize();
        return n;
}
/*----------------------------------------------------------------------------*/
TCoord Face::area() const
{
	math::Vector n;
	std::vector<Node> nodes = this->get<Node>();
	int nb_nodes = nodes.size();

	if(nb_nodes!=3)
		throw GMDSException("Not yet implemented!");

	Node n1 = nodes[0];
	Node n2 = nodes[1];
	Node n3 = nodes[2];

	math::Vector v1(n2.X()-n1.X(), n2.Y()-n1.Y(), n2.Z()-n1.Z());
	math::Vector v3(n3.X()-n1.X(), n3.Y()-n1.Y(), n3.Z()-n1.Z());

	n = (v1.cross(v3));

	return 0.5 * n.norm();
}
/*----------------------------------------------------------------------------*/
math::Point Face::center() const
{
	TCoord p_coords[3] = {0.0,0.0,0.0};

	std::vector<Node> nodes = this->get<Node>();
	int nb_nodes = nodes.size();

	for(int i=0; i<nb_nodes; i++)
	{
		Node n = nodes[i];
		p_coords[0] += n.X();
		p_coords[1] += n.Y();
		p_coords[2] += n.Z();
	}

	p_coords[0] = p_coords[0] / nodes.size();
	p_coords[1] = p_coords[1] / nodes.size();
	p_coords[2] = p_coords[2] / nodes.size();

	math::Point p(p_coords[0],p_coords[1],p_coords[2]);

	return p;
}
/*----------------------------------------------------------------------------*/
double
Face::computeScaledJacobian2D() const
{
        std::vector<Node> nodes = this->get<Node>();

        switch(this->getType()) {
        case GMDS_TRIANGLE:
		throw GMDSException("Region::computeScaledJacobian2D not implemented yet for triangles.");
                break;
        case GMDS_QUAD:
                {
                        gmds::math::Quadrilateral quad(nodes[0].getPoint(),nodes[1].getPoint(),nodes[2].getPoint(),nodes[3].getPoint());
                        return quad.computeScaledJacobian2D();
                }
                break;
        default:
                throw GMDSException("Region::computeScaledJacobian2D not implemented yet for this cell type.");
                break;
        }
}
/*----------------------------------------------------------------------------*/
math::Point Face::project(const math::Point& AP) const
{
	if(m_type!=GMDS_TRIANGLE)
		throw GMDSException("Face::project only implemented for triangular faces");

	std::vector<Node> nodes = this->get<Node>();
	int nb_nodes = nodes.size();
	if(nb_nodes!=3)
		throw GMDSException("Face::project. Error a triangle has not 3 vertices");

	math::Point p0 = nodes[0].getPoint();
	math::Point p1 = nodes[1].getPoint();
	math::Point p2 = nodes[2].getPoint();

	//we get the projected point
	math::Point X = math::Plane(p0,p1,p2).project(AP);

	TCoord x,y,z;
	math::Point::computeBarycentric(p0,p1,p2,X,x,y,z);

	if(x<0.0)
	{
		if(y<0.0)
		{
			return p2;
		}
		else if(y<0.0)
		{
			return p1;
		}
		else
		{
			return math::Segment(p1,p2).project(X);
		}
	}
	else if(y<0.0)
	{
		if(z<0.0)
		{
			return p0;
		}
		else
		{
			return math::Segment(p0,p2).project(X);
		}
	}
	else if(z<0.0)
	{
		return math::Segment(p0,p1).project(X);
	}
	// we are in the triangle
	return X;
}
/*----------------------------------------------------------------------------*/
void Face::getAdjacentNodes(const Node& ANode1, Node& ANode2, Node& ANode3)
{
	std::vector<gmds::Node> nodes = get<Node>();
	if(m_type==GMDS_QUAD)
	{
		int index1=-1, index2, index3;
		for(int i=0;i<4;i++)
			if(nodes[i].getID()==ANode1.getID())
				index1=i;

		switch(index1) {
		case -1:
#ifdef _DEBUG
			std::cout<<"Index error"<<std::endl;
#endif //_DEBUG
			throw GMDSException("getAdjacentNodes: node 1 is not adjacent to the face");
			break;
		case 0:
			index2=3;
			index3=1;
			break;
		case 1:
			index2=0;
			index3=2;
			break;
		case 2:
			index2=1;
			index3=3;
			break;
		case 3:
			index2=2;
			index3=0;
			break;
		}

		ANode2 = nodes[index2];
		ANode3 = nodes[index3];
	}
	else if(m_type==GMDS_TRIANGLE)
	{
		int index1=-1, index2, index3;
		for(int i=0;i<3;i++)
			if(nodes[i].getID()==ANode1.getID())
				index1=i;

		switch(index1) {
		case -1:
#ifdef _DEBUG
			std::cout<<"Index error"<<std::endl;
			for(int i=0;i<3;i++)
				std::cout<<ATab[i]<<" ";
			std::cout<<"<- "<<ANode1<<std::endl;
#endif //_DEBUG
			throw GMDSException("getAdjacentNodes: node 1 is not adjacent to the face");
			break;
		case 0:
			index2=2;
			index3=1;
			break;
		case 1:
			index2=0;
			index3=2;
			break;
		case 2:
			index2=1;
			index3=0;
			break;
		default:
                        throw GMDSException("Face::getAdjacentNodes index incorrectly computed.");
                        break;
		}

		ANode2 = nodes[index2];
		ANode3 = nodes[index3];
	}
	else if(m_type==GMDS_POLYGON)
	{
		bool found=false;
		unsigned int index1, index2,index3;
		for(unsigned int i=0; !found; i++)
			if(nodes[i].getID()==ANode1.getID())
			{
				found=true;
				index1 = i;
			}

		const unsigned int last_index = nodes.size()-1;

		if(!found)
		{
			std::cout<<"Index error"<<std::endl;
			throw GMDSException("getNodeOfOppositeEdge: node 1 is not adjacent to the face");
		}
		else if(index1==0){

			index2=last_index;
			index3=1;
		}
		else if (index1==last_index){
			index2=last_index-1;
			index3=0;
		}
		else{
			index2=index1-1;
			index3=index1+1;
		}

		ANode2 = nodes[index2];
		ANode3 = nodes[index3];
	}
	else
		throw GMDSException("getAdjacent node not availabe");
}
/*----------------------------------------------------------------------------*/
void Face::getOrderedEdges(std::vector<Edge>& AEdges) const
{
	AEdges.clear();

	std::map<FakeEdge,Edge> fakeEdgeMap;

	std::vector<Edge> edges = this->get<Edge>();
	for(int iEdge=0; iEdge<edges.size(); iEdge++) {
		std::vector<Node> nodes = edges[iEdge].get<Node>();
		fakeEdgeMap[FakeEdge(nodes[0].getID(),nodes[1].getID())] = edges[iEdge];
	}

	std::vector<Node> nodes = this->get<Node>();
	for(int iNode=0; iNode<nodes.size(); iNode++) {
		if(fakeEdgeMap.find(FakeEdge(nodes[iNode].getID(),nodes[(iNode+1)%nodes.size()].getID())) == fakeEdgeMap.end()) {
			throw GMDSException("Face::getOrderedEdges could not find edge");
		}
		AEdges.push_back(fakeEdgeMap[FakeEdge(nodes[iNode].getID(),nodes[(iNode+1)%nodes.size()].getID())]);
        }
}
/*----------------------------------------------------------------------------*/
TCoord Face::distance(const math::Point& AP) const
{
	return project(AP).distance(AP);
}
/*----------------------------------------------------------------------------*/
Face::
Face(IGMesh* AMesh,  const ECellType AType, const TCellID& AID)
: Cell(AMesh,AType,AID) //mesh, type and id are filled in
{
	//============================================
	// we keep a reference on the face container
	if(AMesh!=0){
		m_faces_container = AMesh->m_faces_container;
		m_type_id = m_faces_container->getTypeID(AID);
	}
	else{
		m_faces_container = 0;
		m_type_id = NullID;
	}
}
/*----------------------------------------------------------------------------*/
Face::
Face(const Face& AF)
: Cell(AF.m_owner,AF.m_type,AF.m_id)
{
	if(m_owner!=0)
		m_faces_container = m_owner->m_faces_container;
	else
		m_faces_container = 0;

	m_type_id = AF.m_type_id;
}
/*----------------------------------------------------------------------------*/
void Face::operator=(const Face& AF)
{
	m_owner =AF.m_owner;
	m_type=AF.m_type;
	m_id = AF.m_id;
	if(m_owner!=0)
		m_faces_container = m_owner->m_faces_container;
	else
		m_faces_container = 0;
	m_type_id = AF.m_type_id;
}
/*----------------------------------------------------------------------------*/
bool Face::operator==(const Face& AFace) const
{
        return (m_owner == AFace.m_owner && m_id==AFace.m_id);
}
/*----------------------------------------------------------------------------*/
bool Face::operator!=(const Face& AFace) const
{
        return (!(*this == AFace));
}
/*----------------------------------------------------------------------------*/
Face::~Face(){}
/*----------------------------------------------------------------------------*/
void Face::delegateGet(std::vector<Node>&   ACells) const
{
	if(!m_owner->m_model.has(F2N))
		throw GMDSException("F2N adjacency is not supported by the mesh model");

	std::vector<TCellID> cellIDs;
	if(m_type==GMDS_TRIANGLE){
		(*m_faces_container->m_T2N)[m_type_id].values(cellIDs);
	}
	else if (m_type==GMDS_QUAD){
		(*m_faces_container->m_Q2N)[m_type_id].values(cellIDs);
	}
	else if (m_type==GMDS_POLYGON){
		(*m_faces_container->m_P2N)[m_type_id].values(cellIDs);
	}
	else
		throw GMDSException("Not yet implemented");

	ACells.resize(cellIDs.size());
	for(unsigned int i=0;i<cellIDs.size();i++){
		ACells[i] =m_owner->m_nodes_container->buildNode(cellIDs[i]);
	}
}
/*----------------------------------------------------------------------------*/
void Face::delegateGet(std::vector<Edge>&   ACells) const
{
	if(!m_owner->m_model.has(F2E))
		throw GMDSException("F2E adjacency is not supported by the mesh model");

	std::vector<TCellID> cellIDs;
	if(m_type==GMDS_TRIANGLE){
		(*m_faces_container->m_T2E)[m_type_id].values(cellIDs);
	}
	else if (m_type==GMDS_QUAD){
		(*m_faces_container->m_Q2E)[m_type_id].values(cellIDs);
	}
	else if (m_type==GMDS_POLYGON){
		(*m_faces_container->m_P2E)[m_type_id].values(cellIDs);
	}
	else
		throw GMDSException("Not yet implemented");

	ACells.resize(cellIDs.size());
	for(unsigned int i=0;i<cellIDs.size();i++){
		ACells[i] =m_owner->m_edges_container->buildEdge(cellIDs[i]);
	}
}
/*----------------------------------------------------------------------------*/
void Face::delegateGet(std::vector<Face>&   ACells) const
{
	if(!m_owner->m_model.has(F2F))
		throw GMDSException("F2F adjacency is not supported by the mesh model");

	std::vector<TCellID> cellIDs;
	if(m_type==GMDS_TRIANGLE){
		(*m_faces_container->m_T2F)[m_type_id].values(cellIDs);
	}
	else if (m_type==GMDS_QUAD){
		(*m_faces_container->m_Q2F)[m_type_id].values(cellIDs);
	}
	else if (m_type==GMDS_POLYGON){
		(*m_faces_container->m_P2F)[m_type_id].values(cellIDs);
	}
	else
		throw GMDSException("Not yet implemented");

	ACells.resize(cellIDs.size());
	for(unsigned int i=0;i<cellIDs.size();i++){
		ACells[i] =m_faces_container->buildFace(cellIDs[i]);
	}
}
/*----------------------------------------------------------------------------*/
void Face::delegateGet(std::vector<Region>& ACells) const
{
	if(!m_owner->m_model.has(F2R))
		throw GMDSException("F2R adjacency is not supported by the mesh model");

	std::vector<TCellID> cellIDs;
	if(m_type==GMDS_TRIANGLE){
		(*m_faces_container->m_T2R)[m_type_id].values(cellIDs);
	}
	else if (m_type==GMDS_QUAD){
		(*m_faces_container->m_Q2R)[m_type_id].values(cellIDs);
	}
	else if (m_type==GMDS_POLYGON){
		(*m_faces_container->m_P2R)[m_type_id].values(cellIDs);
	}
	else
		throw GMDSException("Not yet implemented");

	ACells.resize(cellIDs.size());
	for(unsigned int i=0;i<cellIDs.size();i++){
		ACells[i] =m_owner->m_regions_container->buildRegion(cellIDs[i]);
	}
}
/*----------------------------------------------------------------------------*/
void Face::delegateGetNodeIDs(std::vector<TCellID>&   ACells) const
{
	if(!m_owner->m_model.has(F2N))
		throw GMDSException("F2N adjacency is not supported by the mesh model");

	if(m_type==GMDS_TRIANGLE){
		(*m_faces_container->m_T2N)[m_type_id].values(ACells);
	}
	else if (m_type==GMDS_QUAD){
		(*m_faces_container->m_Q2N)[m_type_id].values(ACells);
	}
	else if (m_type==GMDS_POLYGON){
		(*m_faces_container->m_P2N)[m_type_id].values(ACells);
	}
}
/*----------------------------------------------------------------------------*/
void Face::delegateGetEdgeIDs(std::vector<TCellID>&   ACells) const
{
	if(!m_owner->m_model.has(F2E))
		throw GMDSException("F2E adjacency is not supported by the mesh model");

	if(m_type==GMDS_TRIANGLE){
		(*m_faces_container->m_T2E)[m_type_id].values(ACells);
	}
	else if (m_type==GMDS_QUAD){
		(*m_faces_container->m_Q2E)[m_type_id].values(ACells);
	}
	else if (m_type==GMDS_POLYGON){
		(*m_faces_container->m_P2E)[m_type_id].values(ACells);
	}
}
/*----------------------------------------------------------------------------*/
void Face::delegateGetFaceIDs(std::vector<TCellID>&   ACells) const
{
	if(!m_owner->m_model.has(F2F))
		throw GMDSException("F2F adjacency is not supported by the mesh model");

	if(m_type==GMDS_TRIANGLE){
		(*m_faces_container->m_T2F)[m_type_id].values(ACells);
	}
	else if (m_type==GMDS_QUAD){
		(*m_faces_container->m_Q2F)[m_type_id].values(ACells);
	}
	else if (m_type==GMDS_POLYGON){
		(*m_faces_container->m_P2F)[m_type_id].values(ACells);
	}
}
/*----------------------------------------------------------------------------*/
void Face::delegateGetRegionIDs(std::vector<TCellID>&   ACells) const
{
	if(!m_owner->m_model.has(F2R))
		throw GMDSException("F2R adjacency is not supported by the mesh model");

	if(m_type==GMDS_TRIANGLE){
		(*m_faces_container->m_T2R)[m_type_id].values(ACells);
	}
	else if (m_type==GMDS_QUAD){
		(*m_faces_container->m_Q2R)[m_type_id].values(ACells);
	}
	else if (m_type==GMDS_POLYGON){
		(*m_faces_container->m_P2R)[m_type_id].values(ACells);
	}
}

/*----------------------------------------------------------------------------*/
void Face::delegateGetAll(std::vector<Node>&   ACells) const
{
	if(!m_owner->m_model.has(F2N))
		throw GMDSException("F2N adjacency is not supported by the mesh model");

	std::vector<TCellID> cellIDs;
	if(m_type==GMDS_TRIANGLE){
		(*m_faces_container->m_T2N)[m_type_id].allValues(cellIDs);
	}
	else if (m_type==GMDS_QUAD){
		(*m_faces_container->m_Q2N)[m_type_id].allValues(cellIDs);
	}
	else if (m_type==GMDS_POLYGON){
		(*m_faces_container->m_P2N)[m_type_id].allValues(cellIDs);
	}
	else
		throw GMDSException("Not yet implemented");

	ACells.resize(cellIDs.size());
	for(unsigned int i=0;i<cellIDs.size();i++){
		ACells[i] =m_owner->m_nodes_container->buildNode(cellIDs[i]);
	}
}
/*----------------------------------------------------------------------------*/
void Face::delegateGetAll(std::vector<Edge>&   ACells) const
{
	if(!m_owner->m_model.has(F2E))
		throw GMDSException("F2E adjacency is not supported by the mesh model");

	std::vector<TCellID> cellIDs;
	if(m_type==GMDS_TRIANGLE){
		(*m_faces_container->m_T2E)[m_type_id].allValues(cellIDs);
	}
	else if (m_type==GMDS_QUAD){
		(*m_faces_container->m_Q2E)[m_type_id].allValues(cellIDs);
	}
	else if (m_type==GMDS_POLYGON){
		(*m_faces_container->m_P2E)[m_type_id].allValues(cellIDs);
	}
	else
		throw GMDSException("Not yet implemented");

	ACells.resize(cellIDs.size());
	for(unsigned int i=0;i<cellIDs.size();i++){
		ACells[i] =m_owner->m_edges_container->buildEdge(cellIDs[i]);
	}
}
/*----------------------------------------------------------------------------*/
void Face::delegateGetAll(std::vector<Face>&   ACells) const
{
	if(!m_owner->m_model.has(F2F))
		throw GMDSException("F2F adjacency is not supported by the mesh model");

	std::vector<TCellID> cellIDs;
	if(m_type==GMDS_TRIANGLE){
		(*m_faces_container->m_T2F)[m_type_id].allValues(cellIDs);
	}
	else if (m_type==GMDS_QUAD){
		(*m_faces_container->m_Q2F)[m_type_id].allValues(cellIDs);
	}
	else if (m_type==GMDS_POLYGON){
		(*m_faces_container->m_P2F)[m_type_id].allValues(cellIDs);
	}
	else
		throw GMDSException("Not yet implemented");

	ACells.resize(cellIDs.size());
	for(unsigned int i=0;i<cellIDs.size();i++){
		ACells[i] =m_faces_container->buildFace(cellIDs[i]);
	}
}
/*----------------------------------------------------------------------------*/
void Face::delegateGetAll(std::vector<Region>& ACells) const
{
	if(!m_owner->m_model.has(F2R))
		throw GMDSException("F2R adjacency is not supported by the mesh model");

	std::vector<TCellID> cellIDs;
	if(m_type==GMDS_TRIANGLE){
		(*m_faces_container->m_T2R)[m_type_id].allValues(cellIDs);
	}
	else if (m_type==GMDS_QUAD){
		(*m_faces_container->m_Q2R)[m_type_id].allValues(cellIDs);
	}
	else if (m_type==GMDS_POLYGON){
		(*m_faces_container->m_P2R)[m_type_id].allValues(cellIDs);
	}
	else
		throw GMDSException("Not yet implemented");

	ACells.resize(cellIDs.size());
	for(unsigned int i=0;i<cellIDs.size();i++){
		ACells[i] =m_owner->m_regions_container->buildRegion(cellIDs[i]);
	}
}
/*----------------------------------------------------------------------------*/
void Face::delegateGetAllNodeIDs(std::vector<TCellID>&   ACells) const
{
        if(!m_owner->m_model.has(F2N))
                throw GMDSException("F2N adjacency is not supported by the mesh model");

        if(m_type==GMDS_TRIANGLE){
                (*m_faces_container->m_T2N)[m_type_id].allValues(ACells);
        }
        else if (m_type==GMDS_QUAD){
                (*m_faces_container->m_Q2N)[m_type_id].allValues(ACells);
        }
        else if (m_type==GMDS_POLYGON){
                (*m_faces_container->m_P2N)[m_type_id].allValues(ACells);
        }
        else
                throw GMDSException("Not yet implemented");
}
/*----------------------------------------------------------------------------*/
void Face::delegateGetAllEdgeIDs(std::vector<TCellID>&   ACells) const
{
        if(!m_owner->m_model.has(F2E))
                throw GMDSException("F2E adjacency is not supported by the mesh model");
                
        if(m_type==GMDS_TRIANGLE){
                (*m_faces_container->m_T2E)[m_type_id].allValues(ACells);
        }       
        else if (m_type==GMDS_QUAD){
                (*m_faces_container->m_Q2E)[m_type_id].allValues(ACells);
        }
        else if (m_type==GMDS_POLYGON){
                (*m_faces_container->m_P2E)[m_type_id].allValues(ACells);
        }
        else
                throw GMDSException("Not yet implemented");
}
/*----------------------------------------------------------------------------*/
void Face::delegateGetAllFaceIDs(std::vector<TCellID>&   ACells) const
{
        if(!m_owner->m_model.has(F2F))
                throw GMDSException("F2F adjacency is not supported by the mesh model");
                
        if(m_type==GMDS_TRIANGLE){
                (*m_faces_container->m_T2F)[m_type_id].allValues(ACells);
        }       
        else if (m_type==GMDS_QUAD){
                (*m_faces_container->m_Q2F)[m_type_id].allValues(ACells);
        }
        else if (m_type==GMDS_POLYGON){
                (*m_faces_container->m_P2F)[m_type_id].allValues(ACells);
        }
        else
                throw GMDSException("Not yet implemented");
}
/*----------------------------------------------------------------------------*/
void Face::delegateGetAllRegionIDs(std::vector<TCellID>&   ACells) const
{
        if(!m_owner->m_model.has(F2R))
                throw GMDSException("F2R adjacency is not supported by the mesh model");
                
        if(m_type==GMDS_TRIANGLE){
                (*m_faces_container->m_T2R)[m_type_id].allValues(ACells);
        }       
        else if (m_type==GMDS_QUAD){
                (*m_faces_container->m_Q2R)[m_type_id].allValues(ACells);
        }
        else if (m_type==GMDS_POLYGON){
                (*m_faces_container->m_P2R)[m_type_id].allValues(ACells);
        }
        else
                throw GMDSException("Not yet implemented");
}
/*----------------------------------------------------------------------------*/
void Face::delegateSetNodeIDs(const std::vector<TCellID>&   ACells)
{
	if(!m_owner->m_model.has(F2N))
		throw GMDSException("F2N adjacency is not supported by the mesh model");

	TInt nb_cells = getNbNodes();
	if(m_type==GMDS_TRIANGLE){
		if(nb_cells!=ACells.size())
			throw GMDSException("Invalid number of adj. entities");

		(*m_faces_container->m_T2N)[m_type_id]=ACells;
	}
	else if (m_type==GMDS_QUAD){
		if(nb_cells!=ACells.size())
			throw GMDSException("Invalid number of adj. entities");

		(*m_faces_container->m_Q2N)[m_type_id]=ACells;
	}
	else if (m_type==GMDS_POLYGON){
		(*m_faces_container->m_P2N)[m_type_id]=ACells;
	}
}
/*----------------------------------------------------------------------------*/
void Face::delegateSetEdgeIDs(const std::vector<TCellID>&   ACells)
{
	if(!m_owner->m_model.has(F2E))
		throw GMDSException("F2E adjacency is not supported by the mesh model");

	TInt nb_cells = getNbEdges();
	if(m_type==GMDS_TRIANGLE){
		if(nb_cells!=ACells.size())
			throw GMDSException("Invalid number of adj. entities");

		(*m_faces_container->m_T2E)[m_type_id]=ACells;
	}
	else if (m_type==GMDS_QUAD){
		if(nb_cells!=ACells.size())
			throw GMDSException("Invalid number of adj. entities");

		(*m_faces_container->m_Q2E)[m_type_id]= ACells;
	}
	else if (m_type==GMDS_POLYGON){
		(*m_faces_container->m_P2E)[m_type_id]=ACells;
	}
}
/*----------------------------------------------------------------------------*/
void Face::delegateSetFaceIDs(const std::vector<TCellID>&   ACells)
{
	if(!m_owner->m_model.has(F2F))
		throw GMDSException("F2F adjacency is not supported by the mesh model");

	TInt nb_cells = getNbFaces();
	if(m_type==GMDS_TRIANGLE){
		if(nb_cells!=ACells.size())
			throw GMDSException("Invalid number of adj. entities");

		(*m_faces_container->m_T2F)[m_type_id]=ACells;
	}
	else if (m_type==GMDS_QUAD){
		if(nb_cells!=ACells.size())
			throw GMDSException("Invalid number of adj. entities");

		(*m_faces_container->m_Q2F)[m_type_id]=ACells;
	}
	else if (m_type==GMDS_POLYGON){

		(*m_faces_container->m_P2F)[m_type_id]=ACells;
	}
}
/*----------------------------------------------------------------------------*/
void Face::delegateSetRegionIDs(const std::vector<TCellID>& ACells)
{
	if(!m_owner->m_model.has(F2R))
		throw GMDSException("F2R adjacency is not supported by the mesh model");

	TInt nb_cells = getNbRegions();
	if(nb_cells!=ACells.size())
		throw GMDSException("Invalid number of adj. entities");

	if(m_type==GMDS_TRIANGLE){
		(*m_faces_container->m_T2R)[m_type_id]=ACells;
	}
	else if (m_type==GMDS_QUAD){
		(*m_faces_container->m_Q2R)[m_type_id]=ACells;
	}
	else if (m_type==GMDS_POLYGON){
		(*m_faces_container->m_P2R)[m_type_id]=ACells;
	}
}
/*----------------------------------------------------------------------------*/
void Face::delegateNodeAdd(TCellID AID)
{
	if(!m_owner->m_model.has(F2N))
		throw GMDSException("F2N adjacency is not supported by the mesh model");

	if(m_type==GMDS_TRIANGLE){
		(*m_faces_container->m_T2N)[m_type_id].add(AID);
	}
	else if (m_type==GMDS_QUAD){
		(*m_faces_container->m_Q2N)[m_type_id].add(AID);
	}
	else if (m_type==GMDS_POLYGON){
		(*m_faces_container->m_P2N)[m_type_id].add(AID);
	}
}
/*----------------------------------------------------------------------------*/
void Face::delegateEdgeAdd(TCellID AID)
{
	if(!m_owner->m_model.has(F2E))
		throw GMDSException("F2E adjacency is not supported by the mesh model");

	if(m_type==GMDS_TRIANGLE){
		(*m_faces_container->m_T2E)[m_type_id].add(AID);
	}
	else if (m_type==GMDS_QUAD){
		(*m_faces_container->m_Q2E)[m_type_id].add(AID);
	}
	else if (m_type==GMDS_POLYGON){
		(*m_faces_container->m_P2E)[m_type_id].add(AID);
	}
}
/*----------------------------------------------------------------------------*/
void Face::delegateFaceAdd(TCellID AID)
{
	if(!m_owner->m_model.has(F2F))
		throw GMDSException("F2F adjacency is not supported by the mesh model");

	if(m_type==GMDS_TRIANGLE){
		(*m_faces_container->m_T2F)[m_type_id].add(AID);
	}
	else if (m_type==GMDS_QUAD){
		(*m_faces_container->m_Q2F)[m_type_id].add(AID);
	}
	else if (m_type==GMDS_POLYGON){
		(*m_faces_container->m_P2F)[m_type_id].add(AID);
	}
}
/*----------------------------------------------------------------------------*/
void Face::delegateRegionAdd(TCellID AID)
{
	if(!m_owner->m_model.has(F2R))
		throw GMDSException("F2R adjacency is not supported by the mesh model");

	if(m_type==GMDS_TRIANGLE){
		(*m_faces_container->m_T2R)[m_type_id].add(AID);
	}
	else if (m_type==GMDS_QUAD){
		(*m_faces_container->m_Q2R)[m_type_id].add(AID);
	}
	else if (m_type==GMDS_POLYGON){
		(*m_faces_container->m_P2R)[m_type_id].add(AID);
	}
}
/*----------------------------------------------------------------------------*/
void Face::delegateNodeRemove(TCellID AID)
{
	if(!m_owner->m_model.has(F2N))
		throw GMDSException("F2N adjacency is not supported by the mesh model");


	if(m_type==GMDS_TRIANGLE){
		(*m_faces_container->m_T2N)[m_type_id].del(AID);
	}
	else if (m_type==GMDS_QUAD){
		(*m_faces_container->m_Q2N)[m_type_id].del(AID);
	}
	else if (m_type==GMDS_POLYGON){
		(*m_faces_container->m_P2N)[m_type_id].del(AID);
	}
}
/*----------------------------------------------------------------------------*/
void Face::delegateEdgeRemove(TCellID AID)
{
	if(!m_owner->m_model.has(F2E))
		throw GMDSException("F2E adjacency is not supported by the mesh model");


	if(m_type==GMDS_TRIANGLE){
		(*m_faces_container->m_T2E)[m_type_id].del(AID);
	}
	else if (m_type==GMDS_QUAD){
		(*m_faces_container->m_Q2E)[m_type_id].del(AID);
	}
	else if (m_type==GMDS_POLYGON){
		(*m_faces_container->m_P2E)[m_type_id].del(AID);
	}
}
/*----------------------------------------------------------------------------*/
void Face::delegateFaceRemove(TCellID AID)
{
	if(!m_owner->m_model.has(F2F))
		throw GMDSException("F2F adjacency is not supported by the mesh model");

	if(m_type==GMDS_TRIANGLE){
		(*m_faces_container->m_T2F)[m_type_id].del(AID);
	}
	else if (m_type==GMDS_QUAD){
		(*m_faces_container->m_Q2F)[m_type_id].del(AID);
	}
	else if (m_type==GMDS_POLYGON){
		(*m_faces_container->m_P2F)[m_type_id].del(AID);
	}
}
/*----------------------------------------------------------------------------*/
void Face::delegateRegionRemove(TCellID AID)
{
	if(!m_owner->m_model.has(F2R))
		throw GMDSException("F2R adjacency is not supported by the mesh model");

	if(m_type==GMDS_TRIANGLE){
		(*m_faces_container->m_T2R)[m_type_id].del(AID);
	}
	else if (m_type==GMDS_QUAD){
		(*m_faces_container->m_Q2R)[m_type_id].del(AID);
	}
	else if (m_type==GMDS_POLYGON){
		(*m_faces_container->m_P2R)[m_type_id].del(AID);
	}
}
/*----------------------------------------------------------------------------*/
void Face::delegateNodeReplace(TCellID AID1, TCellID AID2)
{
	if(!m_owner->m_model.has(F2N))
		throw GMDSException("F2N adjacency is not supported by the mesh model");

	if(m_type==GMDS_TRIANGLE){
		(*m_faces_container->m_T2N)[m_type_id].replace(AID1, AID2);
	}
	else if (m_type==GMDS_QUAD){
		(*m_faces_container->m_Q2N)[m_type_id].replace(AID1, AID2);
	}
	else if (m_type==GMDS_POLYGON){
		(*m_faces_container->m_P2N)[m_type_id].replace(AID1, AID2);
	}
}
/*----------------------------------------------------------------------------*/
void Face::delegateEdgeReplace(TCellID AID1, TCellID AID2)
{
	if(!m_owner->m_model.has(F2E))
		throw GMDSException("F2E adjacency is not supported by the mesh model");

	if(m_type==GMDS_TRIANGLE){
		(*m_faces_container->m_T2E)[m_type_id].replace(AID1, AID2);
	}
	else if (m_type==GMDS_QUAD){
		(*m_faces_container->m_Q2E)[m_type_id].replace(AID1, AID2);
	}
	else if (m_type==GMDS_POLYGON){
		(*m_faces_container->m_P2E)[m_type_id].replace(AID1, AID2);
	}
}
/*----------------------------------------------------------------------------*/
void Face::delegateFaceReplace(TCellID AID1, TCellID AID2)
{
	if(!m_owner->m_model.has(F2F))
		throw GMDSException("F2R adjacency is not supported by the mesh model");

	if(m_type==GMDS_TRIANGLE){
		(*m_faces_container->m_T2F)[m_type_id].replace(AID1, AID2);
	}
	else if (m_type==GMDS_QUAD){
		(*m_faces_container->m_Q2F)[m_type_id].replace(AID1, AID2);
	}
	else if (m_type==GMDS_POLYGON){
		(*m_faces_container->m_P2F)[m_type_id].replace(AID1, AID2);
	}
}
/*----------------------------------------------------------------------------*/
void Face::delegateRegionReplace(TCellID AID1, TCellID AID2)
{
	if(!m_owner->m_model.has(F2R))
		throw GMDSException("F2R adjacency is not supported by the mesh model");

	if(m_type==GMDS_TRIANGLE){
		(*m_faces_container->m_T2R)[m_type_id].replace(AID1, AID2);
	}
	else if (m_type==GMDS_QUAD){
		(*m_faces_container->m_Q2R)[m_type_id].replace(AID1, AID2);
	}
	else if (m_type==GMDS_POLYGON){
		(*m_faces_container->m_P2R)[m_type_id].replace(AID1, AID2);
	}
}
/*----------------------------------------------------------------------------*/
std::ostream & operator << (std::ostream & AStream, const Face & AF)
{
	AStream<<"Face "<<AF.getID()<<" - ";
	if(AF.getType()==GMDS_QUAD)
		AStream<<"Quad";
	else if(AF.getType()==GMDS_TRIANGLE)
		AStream<<"Triangle";
	else if(AF.getType()==GMDS_POLYGON)
		AStream<<"Polygon";
	else
		AStream<<"Invalid type";

	AStream<<std::endl<<"Nb Nodes "<<AF.getNbNodes();
	AStream<<std::endl<<"Nb Faces "<<AF.getNbFaces();
	return AStream;
}
/*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/
