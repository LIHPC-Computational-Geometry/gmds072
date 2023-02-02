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
/** \file    MigrationTool.h
 *  \author  F. LEDOUX
 *  \date    06/16/2009
 */
/*----------------------------------------------------------------------------*/
#ifndef GMDS_MIGRATIONTOOL_H_
#define GMDS_MIGRATIONTOOL_H_
/*----------------------------------------------------------------------------*/
// STL File Headers
#include <string>
/*----------------------------------------------------------------------------*/
#include <GMDS/Utils/CommonTypes.h>
#include <GMDS/IG/IGMesh.h>
/*----------------------------------------------------------------------------*/
#ifdef GMDS_PARALLEL
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
/** \class MigrationTool
 *  \brief ????
 *
 *  \param TMask mesh mask
 */
/*----------------------------------------------------------------------------*/
class MigrationTool{
public:

	/*------------------------------------------------------------------------*/
    /** \brief  Default constructor.
     */
	MigrationTool();

	/*------------------------------------------------------------------------*/
    /** \brief  Destructor.
     */
	virtual ~MigrationTool();

	void exchangeCells(IGMesh& AMesh, const std::vector<TInt>& AFilter) const;

	void sequentialIDSplit(Mesh& AMesh) const;

	void newSequentialIDSplit(Mesh& AMesh) const;

	void getSlaveInfo(Mesh& AMesh,
			std::vector<TCellID>& ALocalNodeSlaves,
			std::vector<std::vector<std::pair<int, TCellID> > >& ADistNodeSlaves,
			std::vector<TCellID>& ALocalEdgeSlaves,
			std::vector<std::vector<std::pair<int, TCellID> > >& ADistEdgeSlaves,
			std::vector<TCellID>& ALocalFaceSlaves,
			std::vector<std::vector<std::pair<int, TCellID> > >& ADistFaceSlaves,
			std::vector<TCellID>& ALocalRegionSlaves,
			std::vector<std::vector<std::pair<int, TCellID> > >& ADistRegionSlaves) const;


	void getDestinations(const TCellID& ACelLID,
						 const std::map<TInt, std::vector<TInt> >& AFilters,
						 std::vector<TInt>& ADests) const;
	/*------------------------------------------------------------------------*/
	/** \brief Extracts from AFromMesh all the mesh faces defined by
	 * 		   AToExtract. All these faces are inserted in AToMesh and the
	 * 		   interfaces of AFromMesh and AToMesh are updated. The interface
	 * 		   between AFromMesh and AToMesh is the index of the resTCellIDent
	 * 		   processor.
	 *
	 * 		   This operation is local to a partition.
	 *
	 * 		   We suppose that every cell knows its adjacent nodes !!
	 *
	 *  \param AFromMesh  The mesh which cells are extracted from
	 *  \param AToMesh 	  The mesh where the cells are added
	 *  \param AToExtract The cells that must be moved.
	 *  \param AToPartitionID TCellID where the cells will be migrated.
	 * */
	void extractSubMesh(Mesh& AFromMesh, Mesh &AToMesh,
						std::vector<Face*>& ToExtract,
						const int& AToPartitionID) const;


	void extractSubMesh(Mesh& AFromMesh, Mesh &AToMesh,
						std::vector<Region*>& ToExtract,
						const int& AToPartitionID) const;

	/*------------------------------------------------------------------------*/
	/** \brief Insert mesh ASubMesh into AToMesh. Interfaces are recomputed and
	 * 		   ASubmesh is empty at the end of the process.
	 *
	 * 		   This operation is local to a partition.
	 *
	 *  \param AMesh 	The mesh where the cells are added
	 *  \param ASubMesh The mesh to insert into AToMesh
	 *  \param AFromPartitionID TCellID of the partition ASubmesh comes from
	 * */
	void insertSubMesh(Mesh& AMesh, Mesh &ASubMesh,
					   const int& AFromPartitionID) const;

	/*------------------------------------------------------------------------*/
	/** \brief Migrate
	 *
	 *  \param AToMigrated The cells that must be migrated.
	 * */
	void migrate( Mesh& AMesh,
			  	  std::map<int, std::vector<Face*> > AToMigrate,
				  std::vector<Mesh >& ASubMeshes) const;


private:

	void sequentialIDSplit2D(Mesh& AMesh) const;
	void newSequentialIDSplit2D(Mesh& AMesh) const;
	void newSequentialIDSplit3D(Mesh& AMesh) const;

	inline void insertNodeOnEdges(Node* AN1, Node* AN2,
								  std::map<FakeEdge::ID, TCellID>& AE2N,
								  std::vector<TCellID> ACentroTCellIDList);

	void sequentialIDSplit3D(Mesh& AMesh) const;

	void initFilters(Mesh& AMesh,
						 const std::vector<TInt>& AFilter,
						 const std::set<TInt>& ATargets,
						 std::vector<std::vector<TInt> >& ANodeFilters,
						 std::vector<std::vector<TInt> >& AEdgeFilters,
						 std::vector<Face*>& AExchangeFaces,
						 const int& AMoveMark) const;

	void initFilters(Mesh& AMesh,
						 const std::vector<TInt>& AFilter,
						 const std::set<TInt>& ATargets,
						 std::vector<std::vector<TInt> >& ANodeFilters,
						 std::vector<std::vector<TInt> >& AFaceFilters,
						 std::vector<Region*>& AExchangeRegions,
						 const int& AMoveMark) const;

	void updateInterfacesBeforeCellExchange(Mesh& AMesh,
			std::vector<std::vector<TInt> >& ANodeFilters,
			 std::map<TInt, std::set<Node*> >& ANodeInterfaces) const;

	void create3DInternalInterfaces(Mesh& AMesh,
			std::vector<std::vector<TInt> >& ANodeFilters,
			std::vector<std::vector<TInt> >& AFaceFilters) const;

	void createInternalInterfaces(Mesh& AMesh,
			std::vector<std::vector<TInt> >& ANodeFilters,
			std::vector<std::vector<TInt> >& AEdgeFilters) const;

	void prepareCellMigration(Mesh& 		AMesh,
		 const std::vector<TInt>& 				AFilter,
		 std::vector<Face*>& 					AExchangeFaces,
		 std::vector<std::vector<TInt> >& 	ANodeFilters,
		 std::vector<std::vector<TInt> >& 	AEdgeFilters,
		 std::map<TInt, std::set<Node*> >& ANodeInterfaces,
		 std::map<TInt, std::set<Edge*> >& AEdgeInterfaces,
		 DistributedMemoryManager::MigrationStruct* AMigInfo,
		 const TInt& ANbMigInfo,
		 const int& AMoveMark) const;

	void prepareCellMigration(Mesh& 		AMesh,
		 const std::vector<TInt>& 				AFilter,
		 std::vector<Region*>& 					AExchangeRegions,
		 std::vector<std::vector<TInt> >& 	ANodeFilters,
		 std::vector<std::vector<TInt> >& 	AFaceFilters,
		 std::map<TInt, std::set<Node*> >& ANodeInterfaces,
		 std::map<TInt, std::set<Face*> >& AFaceInterfaces,
		 DistributedMemoryManager::MigrationStruct* AMigInfo,
		 const TInt& ANbMigInfo,
		 const int& AMoveMark) const;

//	void performCellMigration( const TInt& ARank,
//		std::map<TInt, DistributedMemoryManager::MigrationStruct>& ASendInfo,
//		std::map<TInt, DistributedMemoryManager::MigrationStruct>& ARecvInfo) const;

	void performCellMigrationForIDSplit_Send( const TInt& ARank,TInt *ADestPartitions,
			DistributedMemoryManager::MigrationStruct* ASendInfo, const TInt& ANbSend) const;

	void performCellMigrationForIDSplit_Recv( const TInt& ARank,
		DistributedMemoryManager::MigrationStruct& ARecvInfo) const;

	void updateAfterCellMigration(Mesh& Amesh,
		DistributedMemoryManager::MigrationStruct& ARecvInfo) const;


	void updateAfterCellMigration3D(Mesh& Amesh,
		DistributedMemoryManager::MigrationStruct& ARecvInfo) const;

	void prepareInterfacesMigration(Mesh& AMesh,
			std::vector<std::vector<TInt> >& ANodeFilters,
		 std::map<TInt, std::set<Node*> >& ANodeInterfaces,
		 std::map<TInt, std::pair<std::vector<TCellID>, std::vector<TInt> > >& ADataToSend)
		const;
	void performInterfacesMigration(
		std::map<TInt, std::pair<std::vector<TCellID>, std::vector<TInt> > >& ASend,
		std::map<TInt, std::pair<std::vector<TCellID>, std::vector<TInt> > >& ARecv)
		const;
	void updateAfterInterfacesMigration(Mesh& AMesh,
		std::map<TInt, std::pair<std::vector<TCellID>, std::vector<TInt> > >& ARecv)
		const;

	void initProjectionMap(Mesh& AMesh, Mesh &ASubMesh,
						   std::map<TCellID,TCellID>& ANodeProj,
						   std::map<TCellID,TCellID>& AEdgeProj,
						   std::map<TCellID,TCellID>& AFaceProj,
						   std::map<TCellID,TCellID>& ARegionProj,
						   const int& AFromPartitionID,
						   const int& AInterfaceMark) const;


};
/*----------------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------------*/
#include "MigrationTool.t.h"
/*----------------------------------------------------------------------------*/
#endif // GMDS_PARALLEL
/*----------------------------------------------------------------------------*/
#endif /* GMDS_MIGRATIONTOOL_H_ */
/*----------------------------------------------------------------------------*/
