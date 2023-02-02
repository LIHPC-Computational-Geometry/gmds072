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
/*
 * \file   MesquiteDomainImplAdapter.h
 * \author legoff
   \date   02/06/2015
 */
/*----------------------------------------------------------------------------*/
#include <GMDS/MesquiteDomainImplAdapter.h>

#include <GMDS/CAD/GeomSurface.h>
/*----------------------------------------------------------------------------*/
namespace gmds {
/*----------------------------------------------------------------------------*/
MesquiteDomainImplAdapter::
MesquiteDomainImplAdapter(
		gmds::IGMesh& AMesh,
		gmds::geom::GeomManager& AGeomManager,
		Mesquite::MeshImpl& AMeshAdapter,
		std::map<Mesquite::Mesh::VertexHandle, gmds::Node>& AMesquite2GMDSNodes,
		std::map<Mesquite::Mesh::ElementHandle,gmds::Face>& AMesquite2GMDSCells2D,
		std::map<Mesquite::Mesh::ElementHandle,gmds::Region>& AMesquite2GMDSCells3D,
		bool AIsNormalOutward)
: m_mesh(AMesh), m_geomManager(AGeomManager), m_meshAdapter(AMeshAdapter),
  m_mesquite2GMDSNodes(AMesquite2GMDSNodes),
  m_mesquite2GMDSCells2D(AMesquite2GMDSCells2D),
  m_mesquite2GMDSCells3D(AMesquite2GMDSCells3D),
  m_isNormalOutward(AIsNormalOutward)
{

	gmds::Variable<gmds::geom::GeomEntity* >* nodesClassification = m_mesh.getGeometricClassification(0);

	gmds::IGMesh::node_iterator itn  = m_mesh.nodes_begin();

        for(;!itn.isDone();itn.next()) {
                gmds::Node current_node = itn.value();

		m_nodes2GeomEntity[current_node.getID()] = (*nodesClassification)[current_node.getID()];

		if((*nodesClassification)[current_node.getID()] != NULL) {

			switch((*nodesClassification)[current_node.getID()]->getDim()) {
			case 0:
				m_nodesIsOnVertex[current_node.getID()] = true;
				break;
			case 1:
				m_nodesIsOnCurve[current_node.getID()] = true;
                	        break;
			case 2:
				m_nodesIsOnSurface[current_node.getID()] = true;
                	        break;
			default:
				break;
			}
		}
	}

	gmds::Variable<gmds::geom::GeomEntity* >* facesClassification = m_mesh.getGeometricClassification(2);

	std::map<Mesquite::Mesh::ElementHandle,gmds::Face>::iterator it = AMesquite2GMDSCells2D.begin();
	for(;it!=AMesquite2GMDSCells2D.end(); it++) {
		gmds::Face current_face = it->second;
		m_faces2GeomEntity[current_face.getID()] = (*facesClassification)[current_face.getID()];
	}

}
/*----------------------------------------------------------------------------*/
MesquiteDomainImplAdapter::~MesquiteDomainImplAdapter()
{

}
/*----------------------------------------------------------------------------*/
void MesquiteDomainImplAdapter::snap_to (
		Mesquite::Mesh::EntityHandle AEntityHandle,
		Mesquite::Vector3D & coordinate) const
{
	gmds::TCellID node = m_mesquite2GMDSNodes[AEntityHandle].getID();

	if((m_nodesIsOnVertex.find(node)  != m_nodesIsOnVertex.end())
	|| (m_nodesIsOnCurve.find(node)   != m_nodesIsOnCurve.end())
	|| (m_nodesIsOnSurface.find(node) != m_nodesIsOnSurface.end()))
	{
		gmds::math::Point pt (coordinate.x(),coordinate.y(),coordinate.z());
		gmds::geom::GeomEntity* geomEntity = m_nodes2GeomEntity.find(node)->second;
		geomEntity->project(pt);
		coordinate.set(pt.X(), pt.Y(), pt.Z());
		return;

	} else {
		// in this case the node is free so coordinate is not modified
		return;
	}

}
/*----------------------------------------------------------------------------*/
void MesquiteDomainImplAdapter::vertex_normal_at(
		Mesquite::Mesh::EntityHandle AEntityHandle,
		Mesquite::Vector3D &coordinate) const
{
	gmds::TCellID node = m_mesquite2GMDSNodes[AEntityHandle].getID();

	if(m_nodesIsOnSurface.find(node) == m_nodesIsOnSurface.end()) {
		// in this case the node is not on a surface so the normal cannot be computed
		coordinate.set(0.,0.,0.);
	} else {
		gmds::math::Point pt (coordinate.x(),coordinate.y(),coordinate.z());
		gmds::geom::GeomEntity* geomEntity = m_nodes2GeomEntity.find(node)->second;
		gmds::geom::GeomSurface* surface = dynamic_cast<gmds::geom::GeomSurface*> (geomEntity);
		if(surface == NULL) {
			throw GMDSException("MesquiteDomainImplAdapter::vertex_normal_at dynamic_cast failed.");
		}

		gmds::math::Vector vec;
		surface->computeNormal(pt,vec);

		if(m_isNormalOutward) {
			coordinate.set(vec.X(), vec.Y(), vec.Z());
		} else {
			coordinate.set(-vec.X(), -vec.Y(), -vec.Z());
		}
	}	

}
/*----------------------------------------------------------------------------*/
void MesquiteDomainImplAdapter::element_normal_at(
		Mesquite::Mesh::ElementHandle AEntityHandle,
		Mesquite::Vector3D &coordinate) const
{
	gmds::Face face = m_mesquite2GMDSCells2D[AEntityHandle];
	
	std::vector<Mesquite::Mesh::VertexHandle> verticesHandles;
	std::vector<size_t> offsets;
	Mesquite::MsqError err;

	m_meshAdapter.elements_get_attached_vertices(
					&AEntityHandle,
					1,
					verticesHandles,
					offsets,
					err);
	MSQ_CHKERR (err);

	Mesquite::Mesh::VertexHandle verticesArray[verticesHandles.size()];
	for (unsigned int iVertex=0; iVertex<verticesHandles.size(); iVertex++) {
                verticesArray[iVertex] = verticesHandles[iVertex];
        }
	Mesquite::MsqVertex* coords = new Mesquite::MsqVertex[verticesHandles.size()];

	m_meshAdapter.vertices_get_coordinates(
					verticesArray,
					coords,
					verticesHandles.size(),
					err	
	);
	MSQ_CHKERR (err);

	gmds::math::Point pt(coords[0].x(),coords[0].y(),coords[0].z());

	delete[] coords;

	//gmds::math::Point pt = face.center();
	gmds::math::Vector vec;
	gmds::geom::GeomEntity* geomEntity = m_faces2GeomEntity.find(face.getID())->second;
        gmds::geom::GeomSurface* surface = dynamic_cast<gmds::geom::GeomSurface*> (geomEntity);		

	if(surface == NULL) {
                throw GMDSException("MesquiteDomainImplAdapter::element_normal_at dynamic_cast failed.");
        }

	surface->computeNormal(pt,vec);
	if(m_isNormalOutward) {
                coordinate.set(vec.X(), vec.Y(), vec.Z());
        } else {
                coordinate.set(-vec.X(), -vec.Y(), -vec.Z());
        }

}
/*----------------------------------------------------------------------------*/
void MesquiteDomainImplAdapter::vertex_normal_at(
		const Mesquite::Mesh::EntityHandle * m,
		Mesquite::Vector3D coords[],
		unsigned count,
		Mesquite::MsqError& ) const
{
  for (unsigned i = 0; i < count; ++i) {
	  vertex_normal_at(m[i], coords[i]);
  }

}
/*----------------------------------------------------------------------------*/
void MesquiteDomainImplAdapter::closest_point( Mesquite::Mesh::EntityHandle mesh_entity,
				       const Mesquite::Vector3D& position,
				       Mesquite::Vector3D& closest,
				       Mesquite::Vector3D& normal,
				       Mesquite::MsqError& ) const
{
  closest = position;
  normal = position;

  snap_to(mesh_entity, closest);
  vertex_normal_at(mesh_entity, normal);
}
/*----------------------------------------------------------------------------*/
void MesquiteDomainImplAdapter::domain_DoF( const Mesquite::Mesh::EntityHandle* handle_array,
		unsigned short* dof_array,
		size_t num_vertices,
		Mesquite::MsqError& err ) const
{
	for (unsigned i = 0; i < num_vertices; ++i) {
		gmds::TCellID node = m_mesquite2GMDSNodes[handle_array[i]].getID();

		if(m_nodesIsOnVertex.find(node)  != m_nodesIsOnVertex.end()) {
			dof_array[i] = 0;
			continue;
		} else if(m_nodesIsOnCurve.find(node)  != m_nodesIsOnCurve.end()) {
			dof_array[i] = 1;
			continue;
		} else if(m_nodesIsOnSurface.find(node)  != m_nodesIsOnSurface.end()) {
			dof_array[i] = 2;
			continue;
		}

		dof_array[i] = 3;
	}

}
/*----------------------------------------------------------------------------*/
} // end namespace gmds 
/*----------------------------------------------------------------------------*/
