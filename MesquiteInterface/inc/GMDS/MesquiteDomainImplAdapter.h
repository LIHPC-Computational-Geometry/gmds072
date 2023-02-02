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
#ifndef GMDS_MESQUITEDOMAINIMPLADAPTER_H_
#define GMDS_MESQUITEDOMAINIMPLADAPTER_H_
/*----------------------------------------------------------------------------*/
// stl File Headers
/*----------------------------------------------------------------------------*/
#include <map>
/*----------------------------------------------------------------------------*/
// External File Headers
/*----------------------------------------------------------------------------*/
// Mesquite
#include <MeshInterface.hpp>
#include <MeshImpl.hpp>
/*----------------------------------------------------------------------------*/
// Project File Headers
/*----------------------------------------------------------------------------*/
#include <GMDS/IG/IGMesh.h>
#include <GMDS/CAD/GeomManager.h>
/*----------------------------------------------------------------------------*/
namespace gmds {
/*----------------------------------------------------------------------------*/
/** @brief Classe pont entre une surface Magix3D et Mesquite
*/
class MesquiteDomainImplAdapter : public Mesquite::MeshDomain
{
public:
	MesquiteDomainImplAdapter(
			gmds::IGMesh& AMesh,
			gmds::geom::GeomManager& AGeomManager,
			Mesquite::MeshImpl& AMeshAdapter,
			std::map<Mesquite::Mesh::VertexHandle, gmds::Node>& AMesquite2GMDSNodes,
			std::map<Mesquite::Mesh::ElementHandle,gmds::Face>& AMesquite2GMDSCells2D,
			std::map<Mesquite::Mesh::ElementHandle,gmds::Region>& AMesquite2GMDSCells3D,
			bool AIsNormalOutward);
  ~MesquiteDomainImplAdapter();

  /** Modifies "coordinate" so that it lies on the
   domain to which "entity_handle" is constrained.
   The handle determines the domain.  The coordinate
   is the proposed new position on that domain.
   */
  virtual void snap_to(Mesquite::Mesh::EntityHandle,
		  Mesquite::Vector3D &coordinate) const ;

  /** Returns the normal of the domain to which
   "entity_handle" is constrained.  For non-planar surfaces,
   the normal is calculated at the point on the domain that
   is closest to the passed in value of "coordinate".  If the
   domain does not have a normal, or the normal cannot
   be determined, "coordinate" is set to (0,0,0).  Otherwise,
   "coordinate" is set to the domain's normal at the
   appropriate point.
   In summary, the handle determines the domain.  The coordinate
   determines the point of interest on that domain.

   User should see also PatchData::get_domain_normal_at_vertex and
   PatchData::get_domain_normal_at_element .
   */
  virtual void vertex_normal_at(Mesquite::Mesh::EntityHandle,
		  Mesquite::Vector3D &coordinate) const ;

  virtual void element_normal_at(Mesquite::Mesh::ElementHandle entity_handle,
		  Mesquite::Vector3D &coordinate) const;

  /**\brief evaluate surface normals
   *
   * Returns normals for a domain.
   *
   *\param handles       The domain evaluated is the one in which
   *                     this mesh entity is constrained.
   *\param coordinates   As input, a list of positions at which to
   *                     evaluate the domain.  As output, the resulting
   *                     domain normals.
   *\param count         The length of the coordinates array.
   *\param err           Erreur en retour
   */
  virtual void vertex_normal_at( const Mesquite::Mesh::EntityHandle  * handles,
		  Mesquite::Vector3D coordinates[],
			  unsigned int count,
			  Mesquite::MsqError& err ) const;

  /**\brief evaluate closest point and normal
   *
   * Given a position in space, return the closest
   * position in the domain and the domain normal
   * at that point.
   *
   *\param handle        Evaluate the subset of the domain contianing
   *                     this entity
   *\param position      Input position for which to evaluate
   *\param closest       Closest position in the domain.
   *\param normal        Domain normal at the location of 'closest'
   *\param err           Erreur en retour
   */
  virtual void closest_point( Mesquite::Mesh::EntityHandle handle,
			      const Mesquite::Vector3D& position,
			      Mesquite::Vector3D& closest,
			      Mesquite::Vector3D& normal,
			      Mesquite::MsqError& err ) const;

  /**\brief Get degrees of freedom in vertex movement.
   *
   * Given a vertex, return how the domain constrains the
   * location of that vertex as the number of degrees of
   * freedom in the motion of the vertex.  If the domain
   * is a geometric domain, the degrees of freedom for a
   * vertex is the dimension of the geometric entity the
   * vertex is constrained to lie on (e.g. point = 0, curve = 1,
   * surface = 2, volume = 3.)
   */
  virtual void domain_DoF( const Mesquite::Mesh::EntityHandle* handle_array,
			   unsigned short* dof_array,
			   size_t num_handles,
			   Mesquite::MsqError& err ) const;

private:

	gmds::IGMesh& m_mesh;
	gmds::geom::GeomManager& m_geomManager;

	Mesquite::MeshImpl& m_meshAdapter;


	std::map<Mesquite::Mesh::VertexHandle, gmds::Node>& m_mesquite2GMDSNodes;
	std::map<Mesquite::Mesh::ElementHandle,gmds::Face>& m_mesquite2GMDSCells2D;
	std::map<Mesquite::Mesh::ElementHandle,gmds::Region>& m_mesquite2GMDSCells3D;

	std::map<gmds::TCellID, gmds::geom::GeomEntity*> m_nodes2GeomEntity;
	std::map<gmds::TCellID, gmds::geom::GeomEntity*> m_faces2GeomEntity;
	std::map<gmds::TCellID, bool> m_nodesIsOnSurface;
	std::map<gmds::TCellID, bool> m_nodesIsOnCurve;
	std::map<gmds::TCellID, bool> m_nodesIsOnVertex;

	/* indicates whether the normal to a point on a surface is outward oriented */
	bool m_isNormalOutward;
};
/*----------------------------------------------------------------------------*/
} // end namespace gmds
/*----------------------------------------------------------------------------*/
#endif /* GMDS_MESQUITEDOMAINIMPLADAPTER_H_ */
