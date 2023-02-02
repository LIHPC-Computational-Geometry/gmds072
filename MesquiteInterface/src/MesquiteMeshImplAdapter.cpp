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
/** \file    MesquiteMeshImplAdapter.cpp
 *  \author  legoff
 *  \date    27/05/2015
 */
/*----------------------------------------------------------------------------*/
// stl File Headers
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
// External File Headers
/*----------------------------------------------------------------------------*/
// Mesquite
#include <MeshImplData.hpp>
#include <MsqError.hpp>
/*----------------------------------------------------------------------------*/
// CaGe File Headers
/*----------------------------------------------------------------------------*/
#include <GMDS/MesquiteMeshImplAdapter.h>
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
MesquiteMeshImplAdapter::MesquiteMeshImplAdapter(
						 gmds::IGMesh& AMesh,
						 gmds::geom::GeomManager& AGeomManager,
						 int AMarkFixedNodes)
: m_mesh(AMesh), m_geomManager(AGeomManager)
{
	Mesquite::MsqError err;

	// Tous les noeuds sont pris
	unsigned int vertexCount = m_mesh.getNbNodes();

	// Remplissage des structures pour Mesquite
	myMesh->allocate_vertices (vertexCount, err);
	MSQ_CHKERR (err);

	// table de correspondance entre noeuds Gmds et indices pour Mesquite
	std::map<gmds::TCellID, uint> num_insurf;

	gmds::IGMesh::node_iterator itn  = m_mesh.nodes_begin();

	unsigned int iVertexCount = 0;
	for(;!itn.isDone();itn.next()) {
		gmds::Node current_node = itn.value();

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
		gmds::IGMesh::node_iterator itn  = m_mesh.nodes_begin();

		for(;!itn.isDone();itn.next()) {
			gmds::Node current_node = itn.value();

			m_mesquite2GMDSNodes[vertices[iVertexCount]] = current_node;
			iVertexCount++;
		}
	}

	unsigned int elementCount = m_mesh.getNbRegions();
	myMesh->allocate_elements (elementCount, err);
	MSQ_CHKERR (err);

	std::vector < size_t > vertices;
	Mesquite::EntityTopology elem_type;

	gmds::IGMesh::region_iterator itr  = m_mesh.regions_begin();

	unsigned int iElementCount = 0;
	for(;!itr.isDone();itr.next()) {
		gmds::Region current_region = itr.value();

		std::vector<gmds::TCellID> nodes = current_region.getIDs<gmds::Node>();

		switch (current_region.getType()) {
		case gmds::GMDS_HEX:
			elem_type = Mesquite::HEXAHEDRON;
			vertices.resize (nodes.size());
			for (size_t iNode = 0; iNode < nodes.size(); iNode++)
				vertices[iNode] = num_insurf[nodes[iNode]];
			break;
		case gmds::GMDS_TETRA:
			elem_type = Mesquite::TETRAHEDRON;
			vertices.resize (nodes.size());
			for (size_t iNode = 0; iNode < nodes.size(); iNode++)
				vertices[iNode] = num_insurf[nodes[iNode]];
			break;
		case gmds::GMDS_PYRAMID:
			elem_type = Mesquite::PYRAMID;
			vertices.resize (nodes.size());
			for (size_t iNode = 0; iNode < nodes.size(); iNode++)
				vertices[iNode] = num_insurf[nodes[iNode]];
			break;
		case gmds::GMDS_PRISM3:
			elem_type = Mesquite::PRISM;
			vertices.resize (nodes.size());
			for (size_t iNode = 0; iNode < nodes.size(); iNode++)
				vertices[iNode] = num_insurf[nodes[iNode]];
			break;
		default:
			throw GMDSException ("MesquiteMeshImplAdapter::MesquiteMeshImplAdapter element type not recognized.");
			break;
		}

		myMesh->reset_element (iElementCount, vertices, elem_type, err);
		MSQ_CHKERR (err);
		iElementCount++;
	}

	MSQ_CHKERR (err);
#ifdef _DEBUG2
	std::cout << "iElementCount = "<<iElementCount<<std::endl;
#endif


}
/*----------------------------------------------------------------------------*/
MesquiteMeshImplAdapter::~MesquiteMeshImplAdapter()
{

}
/*----------------------------------------------------------------------------*/
std::map<Mesquite::Mesh::VertexHandle,gmds::Node>*
MesquiteMeshImplAdapter::getMesquite2GMDSNodes()
{
	return &m_mesquite2GMDSNodes;
}
/*----------------------------------------------------------------------------*/
std::map<Mesquite::Mesh::ElementHandle,gmds::Face>*
MesquiteMeshImplAdapter::getMesquite2GMDSCells2D()
{
    return &m_mesquite2GMDSCells2D;
}
/*----------------------------------------------------------------------------*/
std::map<Mesquite::Mesh::ElementHandle,gmds::Region>*
MesquiteMeshImplAdapter::getMesquite2GMDSCells3D()
{
    return &m_mesquite2GMDSCells3D;
}
/*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/
