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
/** \file    MesquiteCaller.cpp
 *  \author  legoff
 *  \date    02/06/2015
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

#include "MeshInterface.hpp"
#include "QualityAssessor.hpp"
#include "ParallelMeshImpl.hpp"

#include "ShapeImprovementWrapper.hpp"
#include "LaplaceWrapper.hpp"
#include "ShapeImprover.hpp"

#include "PaverMinEdgeLengthWrapper.hpp"
#include "SizeAdaptShapeWrapper.hpp"

#include "UntangleWrapper.hpp"
/*----------------------------------------------------------------------------*/
// Project File Headers
/*----------------------------------------------------------------------------*/
#include <GMDS/MesquiteCaller.h>

#include <GMDS/MesquiteMeshImplAdapter.h>
#include <GMDS/MesquiteMeshImplAdapterSurf.h>
#include <GMDS/MesquiteDomainImplAdapter.h>
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
MesquiteCaller::MesquiteCaller(
						 gmds::IGMesh& AMesh,
						 gmds::geom::GeomManager& AGeomManager)
: m_mesh(AMesh), m_geomManager(AGeomManager)
{

}
/*----------------------------------------------------------------------------*/
MesquiteCaller::~MesquiteCaller()
{

}
/*----------------------------------------------------------------------------*/
void 
MesquiteCaller::exec(bool AIsNormalOutward, int AMarkFixedNodes)
{
	Mesquite::MsqError err;

	MesquiteMeshImplAdapter* meshAdapter = new MesquiteMeshImplAdapter(
			m_mesh,
			m_geomManager,
			AMarkFixedNodes);

	std::map<Mesquite::Mesh::VertexHandle,gmds::Node>* mesquite2GMDSNodes = meshAdapter->getMesquite2GMDSNodes();
	std::map<Mesquite::Mesh::VertexHandle,gmds::Face>* mesquite2GMDSCells2D = meshAdapter->getMesquite2GMDSCells2D();
	std::map<Mesquite::Mesh::ElementHandle,gmds::Region>* mesquite2GMDSCells3D = meshAdapter->getMesquite2GMDSCells3D();

	MesquiteDomainImplAdapter* domainAdapter = new MesquiteDomainImplAdapter(
			m_mesh,
			m_geomManager,
			*meshAdapter,
			*mesquite2GMDSNodes,
			*mesquite2GMDSCells2D,
			*mesquite2GMDSCells3D,
			AIsNormalOutward);

	Mesquite::MeshDomainAssoc* myAssoc = new Mesquite::MeshDomainAssoc(meshAdapter, domainAdapter, false, false, true);

	Mesquite::UntangleWrapper* untangleWrapper = new Mesquite::UntangleWrapper();
	untangleWrapper->run_instructions(myAssoc,err);
        MSQ_CHKERR (err);

//	Mesquite::ShapeImprovementWrapper* wrapper = new Mesquite::ShapeImprovementWrapper();
//	Mesquite::LaplaceWrapper* wrapper = new Mesquite::LaplaceWrapper();
	Mesquite::ShapeImprover* wrapper = new Mesquite::ShapeImprover();
//	Mesquite::PaverMinEdgeLengthWrapper* wrapper = new Mesquite::PaverMinEdgeLengthWrapper(10e-10);
//	Mesquite::SizeAdaptShapeWrapper* wrapper = new Mesquite::SizeAdaptShapeWrapper(10e-10);
//	Mesquite::Settings* dummySettings = new Mesquite::Settings();
//	Mesquite::QualityAssessor* qa = new Mesquite::QualityAssessor();
//	Mesquite::ParallelMeshImpl* poyop = new Mesquite::ParallelMeshImpl(meshAdapter);

	wrapper->run_instructions(myAssoc,err);
	MSQ_CHKERR (err);

//	wrapper->run_wrapper(
//			myAssoc,
//			poyop,
//			dummySettings,
//			qa,
//			err
//	);

	std::vector<Mesquite::Mesh::VertexHandle> vertices;
	meshAdapter->get_all_vertices(vertices,err);
	if(err.error()) {
		std::cout<<err.error_message()<<std::endl;
	}

	Mesquite::Mesh::VertexHandle verticesArray[vertices.size()];
	for (unsigned int iVertex=0; iVertex<vertices.size(); iVertex++) {
		verticesArray[iVertex] = vertices[iVertex];
	}

	Mesquite::MsqVertex* coords = new Mesquite::MsqVertex[vertices.size()];
	meshAdapter->vertices_get_coordinates(verticesArray,coords,vertices.size(),err);
	MSQ_CHKERR (err);

	gmds::IGMesh::node_iterator itn  = m_mesh.nodes_begin();

	for (unsigned int iVertex=0; iVertex<vertices.size(); iVertex++) {

//		if(mesquite2GMDSNodes->find(vertices[iVertex]) == mesquite2GMDSNodes->end()) {
//			std::cout<<"MeshImplementation::smooth pb with mesquite2GMDSNodes"<<std::endl;
//			throw TkUtil::Exception("MeshImplementation::smooth pb with mesquite2GMDSNodes");
//		}

//		gmds::Node* node = (*mesquite2GMDSNodes)[vertices[iVertex]];
		gmds::Node current_node = itn.value();
		current_node.setXYZ(coords[iVertex].x(),coords[iVertex].y(),coords[iVertex].z());

		itn.next();
	}

	delete[] coords;
}
/*----------------------------------------------------------------------------*/
void
MesquiteCaller::execSurf(int AMarkFixedNodes)
{
	// first build submesh containing only the faces associated to surfaces
	Mesquite::MsqError err;

        MesquiteMeshImplAdapterSurf* meshAdapter = new MesquiteMeshImplAdapterSurf(
                        m_mesh,
                        m_geomManager,
                        AMarkFixedNodes);

	std::map<Mesquite::Mesh::VertexHandle,gmds::Node>* mesquite2GMDSNodes = meshAdapter->getMesquite2GMDSNodes();
        std::map<Mesquite::Mesh::VertexHandle,gmds::Face>* mesquite2GMDSCells2D = meshAdapter->getMesquite2GMDSCells2D();
	std::map<Mesquite::Mesh::VertexHandle,gmds::Region>* mesquite2GMDSCells3D;

        MesquiteDomainImplAdapter* domainAdapter = new MesquiteDomainImplAdapter(
                        m_mesh,
                        m_geomManager,
			*meshAdapter,
                        *mesquite2GMDSNodes,
                        *mesquite2GMDSCells2D,
			*mesquite2GMDSCells3D,
			true
                        );

        Mesquite::MeshDomainAssoc* myAssoc = new Mesquite::MeshDomainAssoc(meshAdapter, domainAdapter, false, false, true);

        Mesquite::UntangleWrapper* untangleWrapper = new Mesquite::UntangleWrapper();
        untangleWrapper->run_instructions(myAssoc,err);
        MSQ_CHKERR (err);

	Mesquite::ShapeImprover* wrapper = new Mesquite::ShapeImprover();
	wrapper->run_instructions(myAssoc,err);

	std::vector<Mesquite::Mesh::VertexHandle> vertices;
        meshAdapter->get_all_vertices(vertices,err);
        if(err.error()) {
                std::cout<<err.error_message()<<std::endl;
        }

        Mesquite::Mesh::VertexHandle verticesArray[vertices.size()];
        for (unsigned int iVertex=0; iVertex<vertices.size(); iVertex++) {
                verticesArray[iVertex] = vertices[iVertex];
        }

        Mesquite::MsqVertex* coords = new Mesquite::MsqVertex[vertices.size()];
        meshAdapter->vertices_get_coordinates(verticesArray,coords,vertices.size(),err);
        MSQ_CHKERR (err);

        //gmds::IGMesh::node_iterator itn  = m_mesh.nodes_begin();

        for (unsigned int iVertex=0; iVertex<vertices.size(); iVertex++) {

              gmds::Node current_node = (*mesquite2GMDSNodes)[vertices[iVertex]];
              //gmds::Node current_node = itn.value();
              current_node.setXYZ(coords[iVertex].x(),coords[iVertex].y(),coords[iVertex].z());

              //itn.next();
        }	

	delete[] coords;
}
/*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/
