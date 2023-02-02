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
/** \file    MesquiteMeshImplAdapterSurf.cpp
 *  \author  legoff
 *  \date    10/29/2015
 */
/*----------------------------------------------------------------------------*/
// stl File Headers
/*----------------------------------------------------------------------------*/
#include <set>
/*----------------------------------------------------------------------------*/
// External File Headers
/*----------------------------------------------------------------------------*/
// Mesquite
#include <MeshImplData.hpp>
#include <MsqError.hpp>
/*----------------------------------------------------------------------------*/
// CaGe File Headers
/*----------------------------------------------------------------------------*/
#include <GMDS/MesquiteMeshImplAdapterSurf.h>
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
MesquiteMeshImplAdapterSurf::MesquiteMeshImplAdapterSurf(
						 gmds::IGMesh& AMesh,
						 gmds::geom::GeomManager& AGeomManager,
						 int AMarkFixedNodes)
: m_mesh(AMesh), m_geomManager(AGeomManager)
{
	Mesquite::MsqError err;

	// We only keep faces and nodes on surfaces (and curves and vertices)
	std::set<gmds::Face> registeredFaces;
	std::set<gmds::Node> registeredNodes;

	if(!AMesh.doesGeometricClassificationExist(0) || !AMesh.doesGeometricClassificationExist(1) || !AMesh.doesGeometricClassificationExist(2)) {
		throw GMDSException("MesquiteMeshImplAdapterSurf::MesquiteMeshImplAdapterSurf call not permitted when there is no geom classification.");
	}

	gmds::IGMesh::face_iterator itf  = AMesh.faces_begin();
        for(;!itf.isDone();itf.next()) {
                gmds::Face current_face = itf.value();

		gmds::geom::GeomEntity* geomEntity = AMesh.getGeometricClassification(current_face);
		if(geomEntity != NULL) {
			registeredFaces.insert(current_face);		
		}
	}

std::cout<<"registeredFaces.size() "<<registeredFaces.size()<<std::endl;

	for(std::set<gmds::Face>::iterator it = registeredFaces.begin(); it != registeredFaces.end(); it++) {
		std::vector<gmds::Node> nodes = it->get<gmds::Node>();

		for(int iNode=0; iNode<nodes.size(); iNode++){
			registeredNodes.insert(nodes[iNode]);
		}
	}
	unsigned int vertexCount = registeredNodes.size();

	// We orient the faces so that their normals are outward
	std::vector<bool> isFaceOutward(registeredFaces.size(), true); 
	unsigned int facesIndex = 0;

	for(std::set<gmds::Face>::iterator it = registeredFaces.begin(); it != registeredFaces.end(); it++) {
		std::vector<gmds::Region> regions = it->get<gmds::Region>();	
		
		std::vector<gmds::Node> nodes = it->get<gmds::Node>();

		if(regions.size() == 1 || AMesh.getGeometricClassification(regions[0]) != NULL) {
			isFaceOutward[facesIndex] = regions[0].isFaceOrientedOutward(nodes);
		} else {
			isFaceOutward[facesIndex] = regions[1].isFaceOrientedOutward(nodes);
		}
		facesIndex++;
        }

	// Remplissage des structures pour Mesquite
	myMesh->allocate_vertices (vertexCount, err);
	MSQ_CHKERR (err);

	// table de correspondance entre noeuds Gmds et indices pour Mesquite
	std::map<gmds::TCellID, uint> num_insurf;

	unsigned int iVertexCount = 0;
	//gmds::IGMesh::node_iterator itn  = m_mesh.nodes_begin();
	//for(;!itn.isDone();itn.next()) {
	//	gmds::Node current_node = itn.value();
	for(std::set<gmds::Node>::iterator itn = registeredNodes.begin(); itn != registeredNodes.end(); itn++) {
		
		gmds::Node current_node = *itn;

		num_insurf[current_node.getID()] = iVertexCount;

		if(AMarkFixedNodes != -1 && m_mesh.isMarked(current_node,AMarkFixedNodes)) {
			myMesh->reset_vertex (
				iVertexCount,
				Mesquite::Vector3D (
						current_node.X(),
						current_node.Y(),
						current_node.Z()),
						true,
						err);
		} else {
			myMesh->reset_vertex (
                                iVertexCount,
                                Mesquite::Vector3D (
                                                current_node.X(),
                                                current_node.Y(),
                                                current_node.Z()),
                                                false,
                                                err);
		}
		MSQ_CHKERR (err);

		iVertexCount++;
	}

	{
		std::vector<size_t> verts;
		myMesh->all_vertices(verts,err);
		MSQ_CHKERR (err);

		std::vector<Mesquite::Mesh::VertexHandle> vertices;

		get_all_vertices(vertices,err);
		MSQ_CHKERR (err);

		unsigned int iVertexCount = 0;
		//gmds::IGMesh::node_iterator itn  = m_mesh.nodes_begin();

		//for(;!itn.isDone();itn.next()) {
			//gmds::Node current_node = itn.value();
		for(std::set<gmds::Node>::iterator itn = registeredNodes.begin(); itn != registeredNodes.end(); itn++) {

			gmds::Node current_node = *itn;

			m_mesquite2GMDSNodes[vertices[iVertexCount]] = current_node;
			iVertexCount++;
		}
	}

	//unsigned int elementCount = m_mesh.getNbRegions();
	unsigned int elementCount = registeredFaces.size();
	myMesh->allocate_elements (elementCount, err);
	MSQ_CHKERR (err);

	std::vector < size_t > vertices;
	Mesquite::EntityTopology elem_type;

	unsigned int iElementCount = 0;

	for(std::set<gmds::Face>::iterator itf = registeredFaces.begin(); itf != registeredFaces.end(); itf++) {

		gmds::Face current_face = *itf;
	
		std::vector<gmds::TCellID> nodes = current_face.getIDs<gmds::Node>();

		switch (current_face.getType()) {
		case gmds::GMDS_TRIANGLE:
			elem_type = Mesquite::TRIANGLE;
			vertices.resize (nodes.size());
			for (size_t iNode = 0; iNode < nodes.size(); iNode++) {
				if(isFaceOutward[iElementCount]) {
					vertices[iNode] = num_insurf[nodes[iNode]];
				} else {
					vertices[iNode] = num_insurf[nodes[nodes.size()-iNode-1]];
				}
			}
			break;
		case gmds::GMDS_QUAD:
			elem_type = Mesquite::QUADRILATERAL;
			vertices.resize (nodes.size());
			for (size_t iNode = 0; iNode < nodes.size(); iNode++) {
				if(isFaceOutward[iElementCount]) {
                                        vertices[iNode] = num_insurf[nodes[iNode]];
                                } else {
                                        vertices[iNode] = num_insurf[nodes[nodes.size()-iNode-1]];
                                }
			}
			break;
		default:
			throw GMDSException ("MesquiteMeshImplAdapter::MesquiteMeshImplAdapterSurf element type not recognized.");
			break;
		}

		myMesh->reset_element (iElementCount, vertices, elem_type, err);
		MSQ_CHKERR (err);
		iElementCount++;
	}

	{       
                std::vector<size_t> elmnts;
                myMesh->all_elements(elmnts,err);
                MSQ_CHKERR (err);
                
                std::vector<Mesquite::Mesh::ElementHandle> elements;
                
                get_all_elements(elements,err);
                MSQ_CHKERR (err);
                
                unsigned int iElementCount = 0;
                //gmds::IGMesh::node_iterator itn  = m_mesh.nodes_begin();
                
                //for(;!itn.isDone();itn.next()) {
                        //gmds::Node current_node = itn.value();
                for(std::set<gmds::Face>::iterator itf = registeredFaces.begin(); itf != registeredFaces.end(); itf++) {
                        
                        gmds::Face current_face = *itf;
                        
                        m_mesquite2GMDSCells2D[elements[iElementCount]] = current_face;
                        iElementCount++;
                }
        }

	MSQ_CHKERR (err);
#ifdef _DEBUG2
	std::cout << "iElementCount = "<<iElementCount<<std::endl;
#endif


}
/*----------------------------------------------------------------------------*/
MesquiteMeshImplAdapterSurf::~MesquiteMeshImplAdapterSurf()
{

}
/*----------------------------------------------------------------------------*/
std::map<Mesquite::Mesh::VertexHandle,gmds::Node>*
MesquiteMeshImplAdapterSurf::getMesquite2GMDSNodes()
{
	return &m_mesquite2GMDSNodes;
}
/*----------------------------------------------------------------------------*/
std::map<Mesquite::Mesh::ElementHandle,gmds::Face>*
MesquiteMeshImplAdapterSurf::getMesquite2GMDSCells2D()
{
    return &m_mesquite2GMDSCells2D;
}
/*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/
