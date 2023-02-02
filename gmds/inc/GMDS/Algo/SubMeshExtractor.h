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
/** \file    SubMeshExtractor.h
 *  \author  F. LEDOUX
 *  \date    08/31/2010
 */
/*----------------------------------------------------------------------------*/
#ifndef GMDS_EXTRACTOR_H_
#define GMDS_EXTRACTOR_H_
/*----------------------------------------------------------------------------*/
#include <GMDS/Utils/CommonTypes.h>
#include <GMDS/IG/IGMesh.h>
/*----------------------------------------------------------------------------*/
#include <vector>
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
/** \class SubMeshExtractor
 *  \brief Performs different operations in sequential or parallel in order
 *  	   to extract cells from a mesh.
 */
/*----------------------------------------------------------------------------*/

class SubMeshExtractor{

//public:
//
//    /*------------------------------------------------------------------------*/
//    /** \brief Constructor.
//     *
//     *  \param AMesh the mesh where we extract cells from
//     */
//	SubMeshExtractor(IGMesh& AMesh);
//
//    /*------------------------------------------------------------------------*/
//    /** \brief  Destructor.	*/
//	virtual ~SubMeshExtractor();
//
//    /*------------------------------------------------------------------------*/
//    /** \brief this operation builds submeshes according to the partition defined
//     * 		   in variable AName.
//     * 		   Downward connectivities are mandatory.
//     * 		   Same-level/upward connectivities are allowed, but not taken
//     * 		   into account.
//     *
//     *	\param ATypePartAlong the type of the partitioned cells (GMDS_NODE,
//     *					 GMDS_EDGE, GMDS_FACE, GMDS_REGION).
//     *					 Must be the top-level entities at the moment
//     *  \param AName     variable<TInt> name indicating the destination part of cells
//     *  \param SubMeshes container of the extracted submeshes, one for each part
//     */
//	void extractNoGhostWithNativeSharedInfoTopLevelRepartition(
//			const ECellType ATypePartAlong,
//			const std::string& AName,
//			std::vector<IGMesh >& SubMeshes);
//
//	void extractNoGhostWithNativeSharedInfoTopLevelRepartitionExp(
//			const ECellType ATypePartAlong,
//			const std::string& AName,
//			std::vector<IGMesh >& SubMeshes);
//
//	void extractNoGhostWithNativeSharedInfoTopLevelRepartitionExp2(
//			const ECellType ATypePartAlong,
//			const std::string& AName,
//			std::vector<IGMesh >& SubMeshes);
//
//	void extractWithGhostWithNativeSharedInfoTopLevelRepartitionExp2(
//			const ECellType ATypePartAlong,
//			const std::string& AName,
//			std::vector<IGMesh >& SubMeshes);
//
//	/*------------------------------------------------------------------------*/
//	/** \brief this operation builds submeshes according to the partition defined
//	 * 		   in variable AName.
//	 * 		   Same-level/upward connectivities will prompt the generation
//	 *         of a ghost layer; this ghost layer will be formed of the
//	 *         top-level entities and its descendants.
//	 *
//	 *	\param ATypePartAlong the type of the partitioned cells (GMDS_NODE,
//	 *					 GMDS_EDGE, GMDS_FACE, GMDS_REGION).
//	 *					 Must be the top-level entities at the moment
//	 *  \param AName     variable<TInt> name indicating the destination part of cells
//	 *  \param SubMeshes container of the extracted submeshes, one for each part
//	 */
//	void extractWithGhostWithNativeSharedInfoTopLevelRepartition(
//			const ECellType ATypePartAlong,
//			const std::string& AName,
//			std::vector<IGMesh >& SubMeshes);
//
//	// WARNING : do not use, here for testing purpose
//	void partitioning(const ECellType ATypePartAlong, const TInt nParts);
//
//protected:
//
//    /*------------------------------------------------------------------------*/
//    /** \brief  this operation extracts all the cells marked with AMark and
//     * 			fills with these cells the mesh ADestMesh. At the beginning,
//     * 			ADestMesh is cleaned to be empty.
//     *
//     * 			In order to get the cells to extract, this operation has to
//     * 			traverse all the entities of mesh_.
//     *
//     *  \param ADestMesh       the mesh where cells will be created
//     *  \param AMark 	       the boolean mark used to mark the cells to extract
//     *  \param PartID          the part we want to extract
//     *  \param AVarNodeLID     the nodes new local id on their original copy parts
//     *  \param AVarEdgeLID     the edges new local id on their original copy parts
//     *  \param AVarFaceLID     the faces new local id on their original copy parts
//     *  \param AVarRegionLID   the regions new local id on their original copy parts
//     *  \param AVarNodeOwner   the nodes new original copy parts
//     *  \param AVarEdgeOwner   the edges new original copy parts
//     *  \param AVarFaceOwner   the faces new original copy parts
//     *  \param AVarRegionOwner the regions new original copy parts
//     *  \param AMarkShared     the boolean mark used to know when to set master/slave
//     */
//	void getSubMeshRepartition(
//			IGMesh& ADestMesh,
//			const int AMark, const TInt PartID,
//			Variable<TCellID>* AVarNodeLID,
//			Variable<TCellID>* AVarEdgeLID,
//			Variable<TCellID>* AVarFaceLID,
//			Variable<TCellID>* AVarRegionLID,
//			const Variable<TInt>* AVarNodeOwner,
//			const Variable<TInt>* AVarEdgeOwner,
//			const Variable<TInt>* AVarFaceOwner,
//			const Variable<TInt>* AVarRegionOwner,
//			const int AMarkShared);
//
//	void getSubMeshRepartitionExp(
//			std::vector<IGMesh >& SubMeshes,
//			std::vector<TInt>& partsOrdering,
//			std::map<TInt,TInt>& partsOrder,
//			Variable<std::vector<TCellID> >* AVarNodeLIDs,
//			Variable<std::vector<TCellID> >* AVarEdgeLIDs,
//			Variable<std::vector<TCellID> >* AVarFaceLIDs,
//			Variable<std::vector<TCellID> >* AVarRegionLIDs,
//			Variable<std::vector<TInt> >* AVarNodeOwners,
//			Variable<std::vector<TInt> >* AVarEdgeOwners,
//			Variable<std::vector<TInt> >* AVarFaceOwners,
//			Variable<std::vector<TInt> >* AVarRegionOwners,
//			const int AMarkShared);
//
//	void getSubMeshRepartitionExp2(
//			std::vector<IGMesh >& SubMeshes,
//			std::vector<TInt>& partsOrdering,
//			std::map<TInt,TInt>& partsOrder,
//			bool* isSharedNode,
//			bool* isSharedEdge,
//			bool* isSharedFace,
//			bool* isSharedRegion,
//			id* nodes_LID,
//			id* edges_LID,
//			id* faces_LID,
//			id* regions_LID,
//			TInt* nodes_owner,
//			TInt* edges_owner,
//			TInt* faces_owner,
//			TInt* regions_owner,
//			std::vector<std::map<TInt,TCellID> >& nodes_owners2LIDs,
//			std::vector<std::map<TInt,TCellID> >& edges_owners2LIDs,
//			std::vector<std::map<TInt,TCellID> >& faces_owners2LIDs,
//			std::vector<std::map<TInt,TCellID> >& regions_owners2LIDs);
//
//    /*------------------------------------------------------------------------*/
//    /** \brief  this operation removes all the cells not marked with AMark.
//     *
//     * 			In order to get the cells to extract, this operation has to
//     * 			traverse all the entities of mesh_.
//     *
//     *  \param AMark 	       the boolean mark used to mark the cells to extract
//     *  \param PartID          the part we want to keep
//     *  \param AVarNodeLID     the nodes new local id on their original copy parts
//     *  \param AVarEdgeLID     the edges new local id on their original copy parts
//     *  \param AVarFaceLID     the faces new local id on their original copy parts
//     *  \param AVarRegionLID   the regions new local id on their original copy parts
//     *  \param AVarNodeOwner   the nodes new original copy parts
//     *  \param AVarEdgeOwner   the edges new original copy parts
//     *  \param AVarFaceOwner   the faces new original copy parts
//     *  \param AVarRegionOwner the regions new original copy parts
//     *  \param AMarkShared     the boolean mark used to know when to set master/slave
//     */
//	void getSubMeshRepartitionRemove(
//			const int AMark,
//			const TInt PartID,
//			Variable<TCellID>* AVarNodeLID,
//			Variable<TCellID>* AVarEdgeLID,
//			Variable<TCellID>* AVarFaceLID,
//			Variable<TCellID>* AVarRegionLID,
//			const Variable<TInt>* AVarNodeOwner,
//			const Variable<TInt>* AVarEdgeOwner,
//			const Variable<TInt>* AVarFaceOwner,
//			const Variable<TInt>* AVarRegionOwner,
//			const int AMarkShared);
//
//    /*------------------------------------------------------------------------*/
//    /** \brief  this operation marks ACell and its descendants with AMark, sets
//     * 			PartID as the original copy part if it is the first time ACell
//     *			is marked; marks with AMarkShared if it was already.
//	 *
//	 *			Depending on setLID, an original copy local ID is set.
//	 *			Depending on setPartID, an original copy part is set.
//     *
//     *  \param ACell	       the cell we will start the marking on
//     *  \param AMark 	       the boolean mark used to mark the cells
//     *  \param PartID          the part that will be set as owner when relevant
//     *  \param AVarNodeLID     the nodes new local id on their original copy parts
//     *  \param AVarEdgeLID     the edges new local id on their original copy parts
//     *  \param AVarFaceLID     the faces new local id on their original copy parts
//     *  \param AVarRegionLID   the regions new local id on their original copy parts
//     *  \param AVarNodeOwner   the nodes new original copy parts
//     *  \param AVarEdgeOwner   the edges new original copy parts
//     *  \param AVarFaceOwner   the faces new original copy parts
//     *  \param AVarRegionOwner the regions new original copy parts
//     *  \param AMarkShared     the boolean mark used to know whether ACell will be shared
//     *  \param setPartID       the boolean that determines if an owner ID must be set
//     *  \param setLID          the boolean that determines if a local ID must be set
//     */
//	void recursiveMarkWithSharedRepartition(
//			AbstractCell* ACell,
//			const int AMark,
//			TInt PartID,
//			Variable<TCellID>* AVarNodeLID,
//			Variable<TCellID>* AVarEdgeLID,
//			Variable<TCellID>* AVarFaceLID,
//			Variable<TCellID>* AVarRegionLID,
//			Variable<TInt>* AVarNodeOwner,
//			Variable<TInt>* AVarEdgeOwner,
//			Variable<TInt>* AVarFaceOwner,
//			Variable<TInt>* AVarRegionOwner,
//			const int AMarkShared,
//			const bool setPartID,
//			const bool setLID);
//
//	void recursiveMarkWithSharedRepartitionExp(
//			AbstractCell* ACell,
//			TInt PartID,
//			Variable<std::vector<TInt> >* AVarNodeOwners,
//			Variable<std::vector<TInt> >* AVarEdgeOwners,
//			Variable<std::vector<TInt> >* AVarFaceOwners,
//			Variable<std::vector<TInt> >* AVarRegionOwners,
//			const int AMarkShared);
//
//	void recursiveMarkWithSharedRepartitionExp2(
//			AbstractCell* ACell,
//			TInt PartID,
//			bool* isSharedNode,
//			bool* isSharedEdge,
//			bool* isSharedFace,
//			bool* isSharedRegion,
//			id* nodes_LID,
//			id* edges_LID,
//			id* faces_LID,
//			id* regions_LID,
//			TInt* nodes_owner,
//			TInt* edges_owner,
//			TInt* faces_owner,
//			TInt* regions_owner,
//			std::vector<std::map<TInt,TCellID> >& nodes_owners2LIDs,
//			std::vector<std::map<TInt,TCellID> >& edges_owners2LIDs,
//			std::vector<std::map<TInt,TCellID> >& faces_owners2LIDs,
//			std::vector<std::map<TInt,TCellID> >& regions_owners2LIDs);
//
//    /*------------------------------------------------------------------------*/
//    /** \brief  this operation returns true if ACell or at least one of its
//     *          descendants is marked with AMark.
//     *
//     *  \param ACell	the cell we will start the check on
//     *  \param AMark	the boolean mark looked for
//     */
//	bool recursiveIsMarked(AbstractCell* ACell, const int AMark);
//
//	void recursiveGhostMarkingExp2(
//			AbstractCell* ACell,
//			std::vector<TInt>& parts2add,
//			bool* isSharedNode,
//			bool* isSharedEdge,
//			bool* isSharedFace,
//			bool* isSharedRegion,
//			id* nodes_LID,
//			id* edges_LID,
//			id* faces_LID,
//			id* regions_LID,
//			TInt* nodes_owner,
//			TInt* edges_owner,
//			TInt* faces_owner,
//			TInt* regions_owner,
//			std::vector<std::map<TInt,TCellID> >& nodes_owners2LIDs,
//			std::vector<std::map<TInt,TCellID> >& edges_owners2LIDs,
//			std::vector<std::map<TInt,TCellID> >& faces_owners2LIDs,
//			std::vector<std::map<TInt,TCellID> >& regions_owners2LIDs);
//
//	/* a mesh */
//	IGMesh& m_mesh;
};
/*----------------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------------*/
#endif /* GMDS_EXTRACTOR_H_ */
/*----------------------------------------------------------------------------*/

