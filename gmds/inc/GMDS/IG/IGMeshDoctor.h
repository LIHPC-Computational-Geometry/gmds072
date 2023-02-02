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
 * IGMeshDoctor.h
 *
 *  Created on: 22 mai 2014
 *      Author: ledouxf
 */
/*----------------------------------------------------------------------------*/
#ifndef GMDS_IGMESHDOCTOR_H_
#define GMDS_IGMESHDOCTOR_H_
/*----------------------------------------------------------------------------*/
#include <GMDS/IG/IG.h>
#include <map>
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
/** \class IGMeshDoctor
 *
 *  \brief this class provides algorithm to clean, check and modify a mesh.
 *
 */
	class EXPORT_GMDS IGMeshDoctor
{
public:

	/*------------------------------------------------------------------------*/
	/** \brief Constructor
	 * \param AModel a mesh model defining available cells and  connectivities
	 */
	IGMeshDoctor(IGMesh* AMesh);

	/*------------------------------------------------------------------------*/
	/** \brief Destructor
	 */
	virtual ~IGMeshDoctor();

	/*------------------------------------------------------------------------*/
	/** \brief orient faces in the 2D case 
	 */
	int  orient2DFaces();
	/*------------------------------------------------------------------------*/
	/** \brief orient a 2D face
	 */
	bool  orient2DFace(Face& AF);

	/*------------------------------------------------------------------------*/
	/** \brief create faces and R2F adjacency
	 */
	void  buildFacesAndR2F() const;

	/*------------------------------------------------------------------------*/
	/** \brief create faces and R2F adjacency
	 */
	void  buildEdgesAndX2E() const;

	void  updateUpwardConnectivity() const;
    /*------------------------------------------------------------------------*/
    /** \brief create faces and R2F
     */
    void  buildFAndR2F() const;
    /*------------------------------------------------------------------------*/
    /** \brief create faces
     */
    void  buildF() const;

	/*------------------------------------------------------------------------*/
	/** \brief create faces
	 */
	void  buildE() const;

	/*------------------------------------------------------------------------*/
	/** \brief Fill up the X2Y connectivity for m_mesh
	 *  \param ARefModel model that indicates the adjacency relationship that
	 *  	   can be used.
	 */
	void  buildN2N(const MeshModel &ARefModel) const;
	void  buildN2E(const MeshModel &ARefModel) const;
	void  buildN2F(const MeshModel &ARefModel) const;
	void  buildN2R(const MeshModel &ARefModel) const;

	void  buildE2N(const MeshModel &ARefModel) const;
	void  buildE2E(const MeshModel &ARefModel) const;
	void  buildE2F(const MeshModel &ARefModel) const;
	void  buildE2R(const MeshModel &ARefModel) const;

	void  buildF2N(const MeshModel &ARefModel) const;
	void  buildF2E(const MeshModel &ARefModel) const;
	void  buildF2F(const MeshModel &ARefModel) const;
	void  buildF2R(const MeshModel &ARefModel) const;

	void  buildR2N(const MeshModel &ARefModel) const;
	void  buildR2E(const MeshModel &ARefModel) const;
	void  buildR2F(const MeshModel &ARefModel) const;
	void  buildR2R(const MeshModel &ARefModel) const;

	/*------------------------------------------------------------------------*/
	/** \brief Change the mesh algorithms will be applied on
	 */
	void setMesh(IGMesh* AMesh);

 protected:
	TCoord isLeft(Node& AN1, Node& AN2, Node& AN3);
	/*------------------------------------------------------------------------*/
        /** \brief Utilitary method to add a face if it exists
         *
         * \returns the id of the new or old face
         */
        TCellID  addFace(Face& AFace,
                                std::map<FakeFace::FaceID, TCellID>& AFakeFaceMap) const;

	/*------------------------------------------------------------------------*/
	/** \brief Utilitary method to add a face if it does not exist
     
     * \returns the id of the new or old face
     */
    TCellID  addFace(std::vector<TCellID>& ANodeIDs,
                     std::map<FakeFace::FaceID, TCellID>& AFakeFaceMap) const;
    /*------------------------------------------------------------------------*/
	/** \brief Utilitary method to add an edge if it does not exist
	 */
	void  addEdge(TCellID AN1, TCellID AN2,
				std::map<FakeEdge::EdgeID, TCellID>& AFakeEdgeMap) const;
private:
	IGMesh* m_mesh;


};
/*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/
#endif /* GMDS_IGMESHDOCTOR_H_ */
/*----------------------------------------------------------------------------*/

