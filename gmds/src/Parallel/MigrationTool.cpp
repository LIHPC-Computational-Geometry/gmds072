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
/** \file    MigrationTool.t.h
 *  \author  F. LEDOUX
 *  \date    06/16/2009
 */
/*----------------------------------------------------------------------------*/
// GMDS File Headers
#include <GMDS/Parallel/MigrationTool.h>
#include <GMDS/Parallel/DistributedMemoryManager.h>
#include <GMDS/Utils/PlatformUtils.h>
#include <GMDS/Utils/Log.h>
#include <GMDS/Utils/Timer.h>

#include <GMDS/IG/IGMeshDoctor.h>
#include <GMDS/Algo/SubMeshExtractor.h>
//#include <GMDS/Algo/MeshFuse.h>
/*----------------------------------------------------------------------------*/
#include <fstream>
#include <iostream>
#include <iomanip>
///*----------------------------------------------------------------------------*/
//using namespace gmds;
///*----------------------------------------------------------------------------*/
//MigrationTool::MigrationTool()
//{}
///*----------------------------------------------------------------------------*/
//MigrationTool::~MigrationTool()
//{}
///*----------------------------------------------------------------------------*/
//void MigrationTool::initFilters(
//		IGMesh& AMesh,
//		const std::vector<TInt>& AFilter,
//		const std::set<TInt>& ATargets,
//		std::vector<std::vector<TInt> >& ANodeFilters,
//		std::vector<std::vector<TInt> >& AEdgeFilters,
//		std::vector<Face*>& AExchangeFaces,
//		const int& AMoveMark) const
//{
//	int partition_id = AMesh.getPartID();
//
//	TInt max_node_lid =AMesh.getMaxLocalID(0);
//	ANodeFilters.clear();
//	ANodeFilters.resize(max_node_lid+1); // item 0 is  used !!
//
//	TInt max_edge_lid =AMesh.getMaxLocalID(1);
//	AEdgeFilters.clear();
//	AEdgeFilters.resize(max_edge_lid+1); // item 0 is used !!
//
//	AExchangeFaces.clear();
//	/* We fill in the node and edge filters. For that we go through all the faces of the
//	 * current partition even if only a little number of cells are migrated. We must
//	 * do that to be generic. If a face must go to partition P, its nodes must go too.
//	 * Then the corresponding value in the filter will be 1.*/
//
//	/* we also create vector of faces, for the faces that must be migrated. This
//	 * operation is done to avoid to go accross all the faces afterward*/
//
//
//	typename Mesh::faces_iterator face_it 	 = AMesh.faces_begin();
//	typename Mesh::faces_iterator face_it_end = AMesh.faces_end();
//
//	for(;face_it!=face_it_end;face_it++){
//		TInt dest_id = AFilter[(*face_it)->getID()];
//
//		//this face must be moved
//		if(dest_id!=partition_id)
//			AExchangeFaces.push_back(*face_it);
//
//		/* TODO to extend in genericity
//		 * we suppose we have connection F2N */
//		std::vector<id> node_ids = (*face_it)->getNodeIDs();
//
////		//TODO template version
////		/* WARNING: we suppose we have connection F2E */
////		std::vector<id> edge_ids = (*face_it)->getEdgeIDs();
//
//		/* we indicate that each node of the current face must be in the
//		 * dest-id partition too. */
//		const int& nb_nodes = node_ids.size();
//		for(unsigned int i=0;i<nb_nodes;i++)
//		{
//			std::vector<TInt> node_filter = ANodeFilters[node_ids[i]];
//			if(std::find(node_filter.begin(),node_filter.end(),dest_id)==node_filter.end())
//			{
//				ANodeFilters[node_ids[i]].push_back(dest_id);
//				AMesh.mark(AMesh.getLNode(node_ids[i]),AMoveMark);
//			}
//		}
////		/* we indicate that each edge of the current face must be in the
////		 * dest-id partition too. */
////		const int& nb_edges = edge_ids.size();
////		for(unsigned int i=0;i<nb_edges;i++)
////		{
////			std::vector<TInt> edge_filter = AEdgeFilters[edge_ids[i]];
////			if(std::find(edge_filter.begin(),edge_filter.end(),dest_id)==edge_filter.end())
////			{
////				AEdgeFilters[edge_ids[i]].push_back(dest_id);
////				AMesh.mark(AMesh.getLEdge(edge_ids[i]),AMoveMark);
////			}
////		}
//	}
//}
///*----------------------------------------------------------------------------*/
//
//void MigrationTool::
//initFilters(Mesh& AMesh,
//				const std::vector<TInt>& AFilter,
//				const std::set<TInt>& ATargets,
//				std::vector<std::vector<TInt> >& ANodeFilters,
//				std::vector<std::vector<TInt> >& AFaceFilters,
//				std::vector<Region*>& AExchangeRegions,
//				const int& AMoveMark) const
//{
//	int partition_id = AMesh.getPartID();
//
//	TInt max_node_lid =AMesh.getMaxLocalID(0);
//	ANodeFilters.clear();
//	ANodeFilters.resize(max_node_lid+1); // item 0 is not used !!
//
//	TInt max_face_lid =AMesh.getMaxLocalID(2);
//	AFaceFilters.clear();
//	AFaceFilters.resize(max_face_lid+1); // item 0 is not used !!
//
//	/* We fill in the node and face filters. For that we go through all the regions of the
//	 * current partition even if only a little number of cells are migrated. We must
//	 * do that to be generic. If a region must go to partition P, its nodes must go too.
//	 * Then the corresponding value in the filter will be 1.*/
//
//	/* we also create vector of regions, for the regions that must be migrated. This
//	 * operation is done to avoid to go accross all the rgions afterward*/
//
//
//	typename Mesh::regions_iterator region_it 	 = AMesh.regions_begin();
//	typename Mesh::regions_iterator region_it_end = AMesh.regions_end();
//
//	for(;region_it!=region_it_end;region_it++){
//		TInt dest_id = AFilter[(*region_it)->getID()];
//
//		if(dest_id!=partition_id)
//			AExchangeRegions.push_back(*region_it);
//
//		/* we suppose we have connection R2N */
//		std::vector<id> node_ids = (*region_it)->getNodeIDs();
//		/* WARNING: we suppose we have connection R2F */
//		std::vector<id> face_ids = (*region_it)->getFaceIDs();
//
//		/* we indicate that each node of the current region must be in the
//		 * dest-id partition too. */
//		const int& nb_nodes = node_ids.size();
//		for(unsigned int i=0;i<nb_nodes;i++)
//		{
//			std::vector<TInt> node_filter = ANodeFilters[node_ids[i]];
//			if(std::find(node_filter.begin(),node_filter.end(),dest_id)==node_filter.end())
//			{
//				ANodeFilters[node_ids[i]].push_back(dest_id);
//				AMesh.mark(AMesh.getLNode(node_ids[i]),AMoveMark);
//			}
//		}
//		/* we indicate that each face of the current region must be in the
//		 * dest-id partition too. */
//		const int& nb_faces = face_ids.size();
//		for(unsigned int i=0;i<nb_faces;i++)
//		{
//			std::vector<TInt> face_filter = AFaceFilters[face_ids[i]];
//			if(std::find(face_filter.begin(),face_filter.end(),dest_id)==face_filter.end())
//			{
//				AFaceFilters[face_ids[i]].push_back(dest_id);
//				AMesh.mark(AMesh.getLFace(face_ids[i]),AMoveMark);
//			}
//		}
//	}
//}
///*----------------------------------------------------------------------------*/
//
//void MigrationTool::
//prepareInterfacesMigration(Mesh& AMesh,
//		std::vector<std::vector<TInt> >& ANodeFilters,
//		std::map<TInt, SmartVector<std::pair<id,id> > >& ANodeInterfaces,
//		std::map<TInt, std::vector<id> >& ADataToSend)
//		const
//{
//	/* For each neighbour, we know the slaved nodes */
//	std::map<TInt, SmartVector<std::pair<id,id> > >::iterator it = ANodeInterfaces.begin();
//	for(; it!=ANodeInterfaces.end();it++)
//	{
//		TInt neighbour_id = it->first;
//		SmartVector<std::pair<id,id> > shared_nodes = it->second;
//
//		std::vector<id> 	local_data1;
//
//		for(SmartVector<std::pair<id,id> >::iterator it_nodes=shared_nodes.begin();
//			it_nodes!=shared_nodes.end();it_nodes++)
//		{
//			id node_id = (*it_nodes).first;
//			std::vector<TInt> destination = ANodeFilters[node_id];
//			for(unsigned int i=0; i<destination.size();i++)
//			{
//				local_data1.push_back(node_id);
//			}
//		}
//		ADataToSend[neighbour_id] =local_data1;
//	}
//
//}
///*----------------------------------------------------------------------------*/
//
//void MigrationTool::
//performInterfacesMigration(
//	std::map<TInt, std::pair<std::vector<id>, std::vector<TInt> > >& ASend,
//	std::map<TInt, std::pair<std::vector<id>, std::vector<TInt> > >& ARecv)
//	const
//{
//	DistributedMemoryManager* mpi_manager = DistributedMemoryManager::instance();
//
//	TInt current_rank = mpi_manager->rank();
//	std::map<TInt, std::pair<std::vector<id>, std::vector<TInt> > >::iterator it;
//	for(it=ASend.begin();it!=ASend.end();it++)
//	{
//		TInt dest_rank = it->first;
//		std::vector<id>& send_cell_ids  = (it->second).first;
//		std::vector<TInt>&      send_cell_dest = (it->second).second;
//		std::vector<id> recv_cell_ids;
//		std::vector<TInt>      recv_cell_dest;
//
//		if(current_rank>dest_rank)
//		{
//			//send then recv
//			mpi_manager->sendPacked(dest_rank, send_cell_ids, send_cell_dest);
//			mpi_manager->receivePacked(dest_rank, recv_cell_ids, recv_cell_dest);
//		}
//		else {
//			//recv then send
//			mpi_manager->receivePacked(dest_rank, recv_cell_ids, recv_cell_dest);
//			mpi_manager->sendPacked(dest_rank, send_cell_ids, send_cell_dest);
//		}
//
//		ARecv[dest_rank] =
//		std::pair<std::vector<id>, std::vector<TInt> >(recv_cell_ids, recv_cell_dest);
//	}
//}
///*----------------------------------------------------------------------------*/
//
//void MigrationTool::
//updateAfterInterfacesMigration(Mesh& AMesh,
//	std::map<TInt, std::pair<std::vector<id>, std::vector<TInt> > >& ARecv)
//{
//	TInt current_rank = AMesh.getPartID();
//	std::map<TInt, std::pair<std::vector<id>, std::vector<TInt> > >::iterator it;
//	for(it=ARecv.begin();it!=ARecv.end();it++)
//	{
//		TInt from_rank = it->first;
//		SmartVector<std::pair<id,id> > local_interface;
//		SmartVector<std::pair<id,id> >::iterator it_interface;
//		AMesh.getInterface(local_interface,from_rank,0);
//
//		std::vector<id>& cell_ids  = (it->second).first;
//		std::vector<TInt>&      cell_dest = (it->second).second;
//		std::vector<id> to_keep;
//		for(unsigned int i=0; i<cell_ids.size();i++)
//		{
//			id cell_id = cell_ids[i];
//			/* if cell X goes from i to i or to this proc, nothing to do
//			 * we will just be sure to avoid to erase it by adding it again at
//			 * the end
//			 */
//			if( (from_rank==cell_dest[i]))// || (current_rank==cell_dest[i]) )
//			{
//				to_keep.push_back(cell_id);
//				continue;
//			}
//
//			Node* cell =0;
//			bool found_cell = false;
//			for(it_interface=local_interface.begin();
//				it_interface!=local_interface.end() && !found_cell ;
//				it_interface++)
//				if((*it_interface).first == cell_id)
//				{
//					found_cell = true;
//					cell = AMesh.getLNode(cell_id);
//				}
//
//			if(!found_cell)
//			{
//				std::cout<<current_rank<<") problem to find node "<<cell_id<<" with "<<cell_dest[i]<<std::endl;
//				throw GMDSException("Problem in exchangeCells");
//			}
//			AMesh.unsetSlave(cell,from_rank);
//			AMesh.setSlave(cell,NullID,cell_dest[i]);
//		}
//		for(unsigned int i=0;i<to_keep.size();i++)
//		{
//			bool found_cell = false;
//			Node* cell;
//			for(it_interface=local_interface.begin();
//				it_interface!=local_interface.end() && !found_cell ;
//				it_interface++)
//				if((*it_interface).first==to_keep[i])
//				{
//					found_cell = true;
//					cell = AMesh.getLNode(to_keep[i]);
//				}
//
//			if(!found_cell)
//			{
//				std::cout<<current_rank<<") problem to find node "<<to_keep[i]<<std::endl;
//				throw GMDSException("Problem in exchangeCells");
//			}
//			AMesh.setSlave(cell,NullID,from_rank);
//		}
//	}
//
//}
///*----------------------------------------------------------------------------*/
//
//void MigrationTool::
//updateInterfacesBeforeCellExchange(Mesh& AMesh,
//		std::vector<std::vector<TInt> >& ANodeFilters)
//{
//	Log::mng()<<"============ SLAVE NODE TO MIGRATE =============\n";
//
//
//	for(int k=0;k<ANodeFilters.size();k++)
//	{
//		std::vector<TInt> to_go = ANodeFilters[k];
//		if(to_go.size()!=0){
//			Log::mng()<<"Node "<<k<<": ";
//			for(int i=0;i<to_go.size();i++)
//				Log::mng()<<to_go[i]<<" ";
//			Log::mng()<<"\n";
//		}
//	}
//	/* now we have the necessary to update interfaces between partitions and
//	 * exchange cells. We do not create submeshes, it is not necessary. */
//	std::map<TInt, SmartVector<std::pair<id,id> > > initial_node_interfaces;
//
//	AMesh.getInterface(initial_node_interfaces,0);
//
//	/* Some cells of the AMesh interface are going to be moved or duplicated.
//	 * Thus the neighbours must be warned before the effective cell migration.
//	 * Otherwise, the partition interfaces could become invalid.
//	 */
//	std::map<TInt, std::vector<id> > data_to_send;
//	prepareInterfacesMigration(AMesh,ANodeFilters,initial_node_interfaces,data_to_send);
//
//
//	/* we have the full information to send to the neighbourhood in order
//	 * to update partition interfaces. The send/receive process is non blocking because
//	 * the communication are between all the neighbours in both direction. */
//
//	std::map<TInt, std::pair<std::vector<id>, std::vector<TInt> > > data_to_recv;
//	//performInterfacesMigration(data_to_send,data_to_recv);
//
////#ifdef _DEBUG
////	{
////		Log::mng()<<"test output\n";
////
////		Log::mng()<<"============ INIT =============\n";
////		std::map<TInt, std::set<Node*  > > output_interfaces;
////		AMesh.getInterfaces(output_interfaces);
////		std::map<TInt, std::set<Node*  > >::iterator it_out = output_interfaces.begin();
////		for(;it_out != output_interfaces.end();it_out++)
////		{
////
////			Log::mng()<<"INTERFACE WITH "<<it_out->first<<"\n";
////			std::set<Node*> nodes = it_out->second;
////			std::set<Node*>::iterator it_nodes = nodes.begin();
////			for(;it_nodes!=nodes.end();it_nodes++){
////				Log::mng()<<(*it_nodes)->getID()<<" ";
////			}
////			Log::mng()<<"\n";
////		}
////	}
////#endif
//
//
////#ifdef _DEBUG
//	Log::mng()<<"============ SLAVE NODE TO PROPAGATE =============\n";
//	std::map<TInt, std::vector<id> >::iterator it_data_to_send;
//	for(it_data_to_send=data_to_send.begin();it_data_to_send!=data_to_send.end();
//		it_data_to_send++)
//	{
//		TInt dest_rank = it_data_to_send->first;
//		std::vector<id>& cell_ids  = it_data_to_send->second;
//		Log::mng()<<"To "<<dest_rank<<"\n";
//		for(unsigned int i=0;i<cell_ids.size();i++)
//			Log::mng()<<cell_ids[i]<<" ";
//		Log::mng()<<"\n";
//
//	}
////#endif
//
//
//	/* At this time, the current mesh has received the complete information to
//	 * update its connectivity */
////	updateAfterInterfacesMigration(AMesh, data_to_recv);
//
////#ifdef _DEBUG
////	{
////		Log::mng()<<"============ AFTER PROPAG =============\n";
////		std::map<TInt, std::set<Node*  > > output_interfaces;
////		AMesh.getInterfaces(output_interfaces);
////		std::map<TInt, std::set<Node*  > >::iterator it_out = output_interfaces.begin();
////		for(;it_out != output_interfaces.end();it_out++)
////		{
////			Log::mng()<<"INTERFACE WITH "<<it_out->first<<"\n";
////			std::set<Node*> nodes = it_out->second;
////			std::set<Node*>::iterator it_nodes = nodes.begin();
////			for(;it_nodes!=nodes.end();it_nodes++){
////				Log::mng()<<(*it_nodes)->getID()<<" ";
////			}
////			Log::mng()<<"\n";
////		}
////	}
////#endif //_DEBUG
//
//}
///*----------------------------------------------------------------------------*/
//
//void MigrationTool::
//createInternalInterfaces(Mesh& AMesh,
//		std::vector<std::vector<TInt> >& ANodeFilters,
//		std::vector<std::vector<TInt> >& AEdgeFilters) const
//{
//	typename Mesh::nodes_iterator node_it 	 = AMesh.nodes_begin();
//	typename Mesh::nodes_iterator node_it_end = AMesh.nodes_end();
//	for(;node_it!=node_it_end;node_it++)
//	{
//		id cur_lid = (*node_it)->getID();
//		std::vector<TInt> sharing = ANodeFilters[cur_lid];
//		if(sharing.size()>1)
//			AMesh.shareWith(*node_it,sharing);
//	}
//
//	typename Mesh::edges_iterator edge_it 	 = AMesh.edges_begin();
//	typename Mesh::edges_iterator edge_it_end = AMesh.edges_end();
//	for(;edge_it!=edge_it_end;edge_it++)
//	{
//		id cur_lid = (*edge_it)->getID();
//		std::vector<TInt> sharing = AEdgeFilters[cur_lid];
//		if(sharing.size()>1)
//			AMesh.shareWith(*edge_it,sharing);
//	}
//#ifdef _DEBUG
//	{
////		Log::mng()<<"============ INTERNAL =============\n";
////
////		std::map<TInt, std::set<Node*> > output_interfaces;
////		AMesh.getInterfaces(output_interfaces);
////		std::map<TInt, std::set<Node*> >::iterator it_out = output_interfaces.begin();
////		for(;it_out != output_interfaces.end();it_out++)
////		{
////			Log::mng()<<"NODE INTERFACE WITH "<<it_out->first<<"\n";
////			std::set<Node*> nodes = it_out->second;
////			std::set<Node*>::iterator it_nodes = nodes.begin();
////			for(;it_nodes!=nodes.end();it_nodes++){
////				Log::mng()<<(*it_nodes)->getID()<<" ";
////			}
////			Log::mng()<<"\n";
////		}
////
////		std::map<TInt, std::set<Edge*> > edge_output_interfaces;
////		AMesh.getInterfaces(edge_output_interfaces);
////		std::map<TInt, std::set<Edge*> >::iterator it_edge_out = edge_output_interfaces.begin();
////		for(;it_edge_out != edge_output_interfaces.end();it_edge_out++)
////		{
////			Log::mng()<<"EDGE INTERFACE WITH "<<it_edge_out->first<<"\n";
////			std::set<Edge*> edges= it_edge_out->second;
////			std::set<Edge*>::iterator it_edges = edges.begin();
////			for(;it_edges!=edges.end();it_edges++){
////				Log::mng()<<(*it_edges)->getID()<<" ";
////			}
////			Log::mng()<<"\n";
////		}
//	}
//#endif
//
//}
///*----------------------------------------------------------------------------*/
//
//void MigrationTool::
//create3DInternalInterfaces(Mesh& AMesh,
//		std::vector<std::vector<TInt> >& ANodeFilters,
//		std::vector<std::vector<TInt> >& AFaceFilters) const
//{
//	typename Mesh::nodes_iterator node_it 	 = AMesh.nodes_begin();
//	typename Mesh::nodes_iterator node_it_end = AMesh.nodes_end();
//	for(;node_it!=node_it_end;node_it++)
//	{
//		id cur_lid = (*node_it)->getID();
//		std::vector<TInt> sharing = ANodeFilters[cur_lid];
//		if(sharing.size()>1)
//			AMesh.shareWith(*node_it,sharing);
//	}
//
//	typename Mesh::faces_iterator face_it 	 = AMesh.faces_begin();
//	typename Mesh::faces_iterator face_it_end = AMesh.faces_end();
//	for(;face_it!=face_it_end;face_it++)
//	{
//		id cur_lid = (*face_it)->getID();
//		std::vector<TInt> sharing = AFaceFilters[cur_lid];
//		if(sharing.size()>1)
//			AMesh.shareWith(*face_it,sharing);
//	}
//}
///*----------------------------------------------------------------------------*/
//
//void MigrationTool::
//prepareCellMigration(Mesh& AMesh,
//     const std::vector<TInt>& AFilter,
//	 std::vector<Face*>& AExchangeFaces,
//	 std::vector<std::vector<TInt> >& ANodeFilters,
//	 std::vector<std::vector<TInt> >& AEdgeFilters,
//	 std::map<TInt, std::set<Node*> >& ANodeInterfaces,
//	 std::map<TInt, std::set<Edge*> >& AEdgeInterfaces,
//	 DistributedMemoryManager::MigrationStruct* AMigInfo,
//	 const TInt& ANbMigInfo,
//	 const int& AMoveMark) const
//{
//	TInt partition_id = AMesh.getPartID();
//
//	/* we begin by initializing the migration info. Note that each mesh send and
//	 * receive from all its neighbours in order to avoid blocking cases.
//	 */
//	// INUTILE POUR LE DECOUPAGE
////	for(std::map<TInt, std::set<Node*> >::iterator it = ANodeInterfaces.begin();
////		it!=ANodeInterfaces.end();it++){
////		DistributedMemoryManager::MigrationStruct local_struct;
////		local_struct.init();
////		AMigInfo[it->first] = local_struct;
////
////	}
//
//	/* we have an id filter on the nodes to get a local reference to the nodes for every
//	 * partition we export cells. Expensive because depending on the number of partitions*/
//	std::vector<TInt> id_filter[ANbMigInfo];
//	TInt id_increment[ANbMigInfo];
//	for(int i=0;i<ANbMigInfo;i++){
//		AMesh.createFilter(0,id_filter[i]);
//		id_increment[i]=0;
//	}
//
//	/* we create a filter on the nodes giving for each node which partitions it belongs to*/
//	std::vector<std::vector<TInt> > nodes_to_sharing_partitions;
//	TInt max_node_lid = AMesh.getMaxLocalID(0);
//	nodes_to_sharing_partitions.resize(max_node_lid+1);
//	std::map<TInt, std::set<Node*> > node_partitions;
//	AMesh.getInterfaces(node_partitions);
//	std::map<TInt, std::set<Node*> >::iterator node_partitions_it     = node_partitions.begin();
//	std::map<TInt, std::set<Node*> >::iterator node_partitions_it_end = node_partitions.end  ();
//	for(;node_partitions_it!=node_partitions_it_end;node_partitions_it++)
//	{
//		std::set<Node*> nodes = (node_partitions_it)->second;
//		TInt current_pid      = (node_partitions_it)->first;
////		Node* n = *(nodes.begin());
////		std::vector<Node*> vec;
////		vec.resize(nodes.size());
////		const int vec_size = vec.size();
////		for(unsigned int i=0;i<vec_size;i++)
////			vec[i] = n;
//		std::set<Node*>::iterator node_it     = nodes.begin();
//		std::set<Node*>::iterator node_it_end = nodes.end();
//
//		for(;node_it!=node_it_end;node_it++)
//		{
//			nodes_to_sharing_partitions[(*node_it)->getID()].push_back(current_pid);
//		}
//	}
//
//	/* we create a filter on the edges giving for each edge which partitions it belongs to*/
//	std::vector<std::vector<TInt> > edges_to_sharing_partitions;
//	TInt max_edge_lid = AMesh.getMaxLocalID(1);
//	edges_to_sharing_partitions.resize(max_edge_lid+1);
//	std::map<TInt, std::set<Edge*> > edge_partitions;
//	AMesh.getInterfaces(edge_partitions);
//	std::map<TInt, std::set<Edge*> >::iterator edge_partitions_it     = edge_partitions.begin();
//	std::map<TInt, std::set<Edge*> >::iterator edge_partitions_it_end = edge_partitions.end  ();
//	for(;edge_partitions_it!=edge_partitions_it_end;edge_partitions_it++)
//	{
//		std::set<Edge*> edges = (edge_partitions_it)->second;
//		TInt current_pid      = (edge_partitions_it)->first;
////		Edge* n = *(edges.begin());
////		std::vector<Edge*> vec;
////		vec.resize(edges.size());
////		const int vec_size = vec.size();
////		for(unsigned int i=0;i<vec_size;i++)
////			vec[i] = n;
//
//		std::set<Edge*>::iterator edge_it     = edges.begin();
//		std::set<Edge*>::iterator edge_it_end = edges.end();
//
//		for(;edge_it!=edge_it_end;edge_it++)
//		{
//			edges_to_sharing_partitions[(*edge_it)->getID()].push_back(current_pid);
//		}
//	}
//	/* We go through the nodes that can be moved. For that we go through all the nodes and
//	 * we only consider nodes marked with AMoveMark */
//
//
//	typename Mesh::nodes_iterator node_it 	 = AMesh.nodes_begin();
//	typename Mesh::nodes_iterator node_it_end = AMesh.nodes_end();
//	for(;node_it!=node_it_end;node_it++)
//	{
//		/* it is not a node to be moving */
//		if(!AMesh.isMarked(*node_it,AMoveMark))
//			continue;
//
//		Node* current_node = *node_it;
//
//
//		/* on commence par regarder avec qui le noeud est partage et s'il doit etre
//		 * conserve localement */
//		std::vector<TInt> sharing_p = ANodeFilters[(*node_it)->getID()];
//
//
//		const int sharing_p_size = sharing_p.size();
//		for(unsigned int i_dest = 0; i_dest<sharing_p_size;i_dest++)
//		{
//			if(sharing_p[i_dest]==partition_id)
//				continue;
//
//			TInt destination = sharing_p[i_dest];
//			(id_filter[destination])[current_node->getID()] = (id_increment[destination])++;
//			AMigInfo[destination].nb_nodes++;
//			AMigInfo[destination].node_ids.push_back(current_node->getID());
//
//			AMigInfo[destination].node_coords.push_back(current_node->getX());
//			AMigInfo[destination].node_coords.push_back(current_node->getY());
//			AMigInfo[destination].node_coords.push_back(current_node->getZ());
//
//
//			/* on transmet les informations concernant les noeuds partagés. Un noeud
//			 * est partagé s'il appartient au moins à un filtre, ou qu'il appartient
//			 * à un filtre et une interface du maillage*/
////			if(sharing_p.size()>1)
////			{
//				for(unsigned int j_dest = 0; j_dest<sharing_p.size();j_dest++)
//					AMigInfo[destination].node_sharing[sharing_p[j_dest]].push_back(current_node->getID());
////			}
//			if(AMesh.isShared(current_node))
//			{
//				std::vector<TInt> shared_partitions = nodes_to_sharing_partitions[current_node->getID()];
////				AMesh.getPartitionsSharing(current_node,shared_partitions);
//				//migration_info[destination].node_sharing[sharing_p[0]].push_back(current_node->getID());
//				for(unsigned int j_dest = 0; j_dest<shared_partitions.size();j_dest++)
//					AMigInfo[destination].node_sharing[shared_partitions[j_dest]].push_back(current_node->getID());
//			}
//		}
//
//
//	}
//	typename Mesh::edges_iterator edge_it 	 = AMesh.edges_begin();
//	typename Mesh::edges_iterator edge_it_end = AMesh.edges_end();
//	for(;edge_it!=edge_it_end;edge_it++)
//	{
//		/* it is not an edge to be moving */
//		if(!AMesh.isMarked(*edge_it,AMoveMark))
//			continue;
//
//		Edge* current_edge= *edge_it;
//		bool keep_local =false;
//		/* on commence par regarder avec qui l'arete est partagee et si elle doit etre
//		 * conserve localement */
//		std::vector<TInt> sharing_p = AEdgeFilters[current_edge->getID()];
//
//		for(unsigned int i_dest = 0; i_dest<sharing_p.size();i_dest++)
//		{
//			if(sharing_p[i_dest]==partition_id)
//			{
//				keep_local = true;
//				continue;
//			}
//
//			TInt destination = sharing_p[i_dest];
//			AMigInfo[destination].nb_edges++;
//			AMigInfo[destination].edge_ids.push_back(current_edge->getID());
//
//			std::vector<id> edge_node_ids = current_edge->getNodeIDs();
//			AMigInfo[destination].edge_2N.push_back((id_filter[destination])[edge_node_ids[0]]);
//			AMigInfo[destination].edge_2N.push_back((id_filter[destination])[edge_node_ids[1]]);
//
//
//			/* on transmet les informations concernant les aretes partagés. Une arete
//			 * est partagé si elle appartient au moins à un filtre, ou qu'il appartient
//			 * à un filtre et une interface du maillage*/
////			if(sharing_p.size()>1)
////			{
//				for(unsigned int j_dest = 0; j_dest<sharing_p.size();j_dest++)
//					AMigInfo[destination].edge_sharing[sharing_p[j_dest]].push_back(current_edge->getID());
////			}
//			if(AMesh.isShared(current_edge))
//			{
//				std::vector<TInt> shared_partitions = edges_to_sharing_partitions[current_edge->getID()];
//
//				//migration_info[destination].node_sharing[sharing_p[0]].push_back(current_node->getID());
//				for(unsigned int j_dest = 0; j_dest<shared_partitions.size();j_dest++)
//					AMigInfo[destination].edge_sharing[shared_partitions[j_dest]].push_back(current_edge->getID());
//			}
//		}
//
//		/* the edge can be removed from the mesh now */
//		if(!keep_local)
//			AMesh.deleteEdge(current_edge);
//
//	}
//
//
//	for(unsigned i_f=0;i_f<AExchangeFaces.size();i_f++)
//	{
//		Face* current_face = AExchangeFaces[i_f];
//		TInt destination = AFilter[current_face->getID()];
//		AMigInfo[destination].nb_faces++;
//		AMigInfo[destination].face_ids.push_back(current_face->getID());
//		AMigInfo[destination].face_type.push_back(current_face->getType());
//		std::vector<id> face_node_ids = current_face->getNodeIDs();
//		std::vector<TInt> lids;
//		lids.reserve(face_node_ids.size());
//		for(unsigned int i_n=0;i_n<face_node_ids.size();i_n++)
//			lids.push_back((id_filter[destination])[face_node_ids[i_n] ]);
//
//		AMigInfo[destination].face_2N.push_back(lids);
//
//		/* the face can be removed from the mesh now */
//		AMesh.deleteFace(current_face);
//	}
//
//	for(node_it= AMesh.nodes_begin();node_it!=node_it_end;node_it++)
//	{
//		/* it is not a node to be moving */
//		if(!AMesh.isMarked(*node_it,AMoveMark))
//			continue;
//
//		bool keep_local = false;
//		Node* current_node = *node_it;
//		/* on regarde si le noeud doit etre conserve localement */
//		std::vector<TInt> dest_partitions = ANodeFilters[current_node->getID()];
//		if(std::find(dest_partitions.begin(),dest_partitions.end(),partition_id)==dest_partitions.end())
//			AMesh.deleteNode(current_node);
//	}
//
//	/* before sending we remove the artificial interface with itself */
//	std::set<Node*> artificial_node_interface;
//	AMesh.getInterface(artificial_node_interface,partition_id);
//	for(std::set<Node*>::iterator it = artificial_node_interface.begin();
//		it!=artificial_node_interface.end();it++)
//	{
//		AMesh.removeFromInterface(*it,partition_id);
//	}
//
//	/* before sending we remove the artificial interface with itself */
//	std::set<Edge*> artificial_edge_interface;
//	AMesh.getInterface(artificial_edge_interface,partition_id);
//	for(std::set<Edge*>::iterator it = artificial_edge_interface.begin();
//		it!=artificial_edge_interface.end();it++)
//	{
//		AMesh.removeFromInterface(*it,partition_id);
//	}
//
//#ifdef _DEBUG
//	{
////		Log::mng()<<"============ BEFORE SENDING =============\n";
////		std::map<TInt, std::set<Node*> > output_interfaces;
////		AMesh.getInterfaces(output_interfaces);
////		std::map<TInt, std::set<Node*> >::iterator it_out = output_interfaces.begin();
////		for(;it_out != output_interfaces.end();it_out++)
////		{
////			Log::mng()<<"INTERFACE WITH "<<it_out->first<<"\n";
////			std::set<Node*> nodes = it_out->second;
////			std::set<Node*>::iterator it_nodes = nodes.begin();
////			for(;it_nodes!=nodes.end();it_nodes++){
////				Log::mng()<<(*it_nodes)->getID()<<" ";
////			}
////			Log::mng()<<"\n";
////		}
//	}
//#endif //_DEBUG
//
//}
///*----------------------------------------------------------------------------*/
//
//void MigrationTool::
//prepareCellMigration(Mesh& AMesh,
//     const std::vector<TInt>& AFilter,
//	 std::vector<Region*>& AExchangeRegions,
//	 std::vector<std::vector<TInt> >& ANodeFilters,
//	 std::vector<std::vector<TInt> >& AFaceFilters,
//	 std::map<TInt, std::set<Node*> >& ANodeInterfaces,
//	 std::map<TInt, std::set<Face*> >& AFaceInterfaces,
//	 DistributedMemoryManager::MigrationStruct* AMigInfo,
//	 const TInt& ANbMigInfo,
//	 const int& AMoveMark) const
//{
//	TInt partition_id = AMesh.getPartID();
//
//	/* we have an id filter on the nodes to get a local reference to the nodes for every
//	 * partition we export cells. Expensive because depending on the number of partitions*/
//	std::vector<TInt> id_filter[ANbMigInfo];
//	TInt id_increment[ANbMigInfo];
//	for(int i=0;i<ANbMigInfo;i++){
//		AMesh.createFilter(0,id_filter[i]);
//		id_increment[i]=0;
//	}
//
//	/* we create a filter on the nodes giving for each node which partitions it belongs to*/
//	std::vector<std::vector<TInt> > nodes_to_sharing_partitions;
//	TInt max_node_lid = AMesh.getMaxLocalID(0);
//	nodes_to_sharing_partitions.resize(max_node_lid+1);
//	std::map<TInt, std::set<Node*> > node_partitions;
//	AMesh.getInterfaces(node_partitions);
//	std::map<TInt, std::set<Node*> >::iterator node_partitions_it     = node_partitions.begin();
//	std::map<TInt, std::set<Node*> >::iterator node_partitions_it_end = node_partitions.end  ();
//	for(;node_partitions_it!=node_partitions_it_end;node_partitions_it++)
//	{
//		std::set<Node*> nodes = (node_partitions_it)->second;
//		TInt current_pid      = (node_partitions_it)->first;
//		std::set<Node*>::iterator node_it     = nodes.begin();
//		std::set<Node*>::iterator node_it_end = nodes.end();
//
//		for(;node_it!=node_it_end;node_it++)
//			nodes_to_sharing_partitions[(*node_it)->getID()].push_back(current_pid);
//	}
//
//	/* we create a filter on the faces giving for each face which partitions it belongs to*/
//	std::vector<std::vector<TInt> > faces_to_sharing_partitions;
//	TInt max_face_lid = AMesh.getMaxLocalID(2);
//	faces_to_sharing_partitions.resize(max_face_lid+1);
//	std::map<TInt, std::set<Face*> > face_partitions;
//	AMesh.getInterfaces(face_partitions);
//	std::map<TInt, std::set<Face*> >::iterator face_partitions_it     = face_partitions.begin();
//	std::map<TInt, std::set<Face*> >::iterator face_partitions_it_end = face_partitions.end  ();
//	for(;face_partitions_it!=face_partitions_it_end;face_partitions_it++)
//	{
//		std::set<Face*> faces = (face_partitions_it)->second;
//		TInt current_pid      = (face_partitions_it)->first;
//
//		std::set<Face*>::iterator face_it     = faces.begin();
//		std::set<Face*>::iterator face_it_end = faces.end();
//
//		for(;face_it!=face_it_end;face_it++)
//		{
//			faces_to_sharing_partitions[(*face_it)->getID()].push_back(current_pid);
//		}
//	}
//	/* We go through the nodes that can be moved. For that we go through all the nodes and
//	 * we only consider nodes marked with AMoveMark */
//
//
//	typename Mesh::nodes_iterator node_it 	 = AMesh.nodes_begin();
//	typename Mesh::nodes_iterator node_it_end = AMesh.nodes_end();
//	for(;node_it!=node_it_end;node_it++)
//	{
//		/* it is not a node to be moving */
//		if(!AMesh.isMarked(*node_it,AMoveMark))
//			continue;
//
//		Node* current_node = *node_it;
//
//
//		/* on commence par regarder avec qui le noeud est partage et s'il doit etre
//		 * conserve localement */
//		std::vector<TInt> sharing_p = ANodeFilters[(*node_it)->getID()];
//
//
//		const int sharing_p_size = sharing_p.size();
//		for(unsigned int i_dest = 0; i_dest<sharing_p_size;i_dest++)
//		{
//			if(sharing_p[i_dest]==partition_id)
//				continue;
//
//			TInt destination = sharing_p[i_dest];
//			(id_filter[destination])[current_node->getID()] = (id_increment[destination])++;
//			AMigInfo[destination].nb_nodes++;
//			AMigInfo[destination].node_ids.push_back(current_node->getID());
//
//			AMigInfo[destination].node_coords.push_back(current_node->getX());
//			AMigInfo[destination].node_coords.push_back(current_node->getY());
//			AMigInfo[destination].node_coords.push_back(current_node->getZ());
//
//
//			/* on transmet les informations concernant les noeuds partagés. Un noeud
//			 * est partagé s'il appartient au moins à un filtre, ou qu'il appartient
//			 * à un filtre et une interface du maillage*/
//
//			for(unsigned int j_dest = 0; j_dest<sharing_p.size();j_dest++)
//				AMigInfo[destination].node_sharing[sharing_p[j_dest]].push_back(current_node->getID());
//
//			if(AMesh.isShared(current_node))
//			{
//				std::vector<TInt> shared_partitions = nodes_to_sharing_partitions[current_node->getID()];
//				for(unsigned int j_dest = 0; j_dest<shared_partitions.size();j_dest++)
//					AMigInfo[destination].node_sharing[shared_partitions[j_dest]].push_back(current_node->getID());
//			}
//
//			/* si le noeud est dans un groupe on indique aussi le groupe auquel il appartient. */
////			if(AMesh.isInAGroup(current_node))
//			{
//				std::vector<std::string> group_names;
//				AMesh.getGroups(current_node,group_names);
//				for(unsigned int i_name=0;i_name<group_names.size();i_name++)
//				{
//					AMigInfo[destination].node_groups[group_names[i_name]].push_back(current_node->getID());
//					std::cout<<"Send to: "<<group_names[i_name]<<" <- "<<(current_node->getID())<<std::endl;
//				}
//			}
//
//		}
//
//
//	}
//	typename Mesh::faces_iterator face_it 	 = AMesh.faces_begin();
//	typename Mesh::faces_iterator face_it_end = AMesh.faces_end();
//	for(;face_it!=face_it_end;face_it++)
//	{
//		/* it is not an edge to be moving */
//		if(!AMesh.isMarked(*face_it,AMoveMark))
//			continue;
//
//		Face* current_face= *face_it;
//		bool keep_local =false;
//		/* on commence par regarder avec qui la face est partagee et si elle doit etre
//		 * conserve localement */
//		std::vector<TInt> sharing_p = AFaceFilters[current_face->getID()];
//
//		for(unsigned int i_dest = 0; i_dest<sharing_p.size();i_dest++)
//		{
//			if(sharing_p[i_dest]==partition_id)
//			{
//				keep_local = true;
//				continue;
//			}
//
//			TInt destination = sharing_p[i_dest];
//			AMigInfo[destination].nb_faces++;
//			AMigInfo[destination].face_ids.push_back(current_face->getID());
//			AMigInfo[destination].face_type.push_back(current_face->getType());
//			std::vector<id> face_node_ids = current_face->getNodeIDs();
//			std::vector<TInt>     F2N_local_numbering;
//			F2N_local_numbering.reserve(face_node_ids.size());
//			for(unsigned int i_node=0;i_node<face_node_ids.size();i_node++)
//			{
//				F2N_local_numbering.push_back((id_filter[destination])[face_node_ids[i_node]]);
//			}
//			AMigInfo[destination].face_2N.push_back(F2N_local_numbering);
//
//			/* on transmet les informations concernant les faces partagées. Une face
//			 * est partagé si elle appartient au moins à un filtre, ou qu'il appartient
//			 * à un filtre et une interface du maillage*/
//			for(unsigned int j_dest = 0; j_dest<sharing_p.size();j_dest++)
//				AMigInfo[destination].face_sharing[sharing_p[j_dest]].push_back(current_face->getID());
//			if(AMesh.isShared(current_face))
//			{
//				std::vector<TInt> shared_partitions = faces_to_sharing_partitions[current_face->getID()];
//
//				for(unsigned int j_dest = 0; j_dest<shared_partitions.size();j_dest++)
//					AMigInfo[destination].face_sharing[shared_partitions[j_dest]].push_back(current_face->getID());
//			}
//			/* si la face est dans un groupe on indique aussi le groupe auquel il appartient. */
//			if(AMesh.isInAGroup(current_face))
//			{
//				std::vector<std::string> group_names;
//				AMesh.getGroups(current_face,group_names);
//				for(unsigned int i_name=0;i_name<group_names.size();i_name++)
//					AMigInfo[destination].face_groups[group_names[i_name]].push_back(current_face->getID());
//			}
//
//
//		}
//
//		/* the face can be removed from the mesh now */
//		if(!keep_local)
//			AMesh.deleteFace(current_face);
//
//	}
//
//
//	for(unsigned i_r=0;i_r<AExchangeRegions.size();i_r++)
//	{
//		Region* current_region = AExchangeRegions[i_r];
//		TInt destination = AFilter[current_region->getID()];
//		AMigInfo[destination].nb_regions++;
//		AMigInfo[destination].region_ids.push_back(current_region->getID());
//		AMigInfo[destination].region_type.push_back(current_region->getType());
//		std::vector<id> region_node_ids = current_region->getNodeIDs();
//		std::vector<TInt> lids;
//		lids.reserve(region_node_ids.size());
//		for(unsigned int i_n=0;i_n<region_node_ids.size();i_n++)
//			lids.push_back((id_filter[destination])[region_node_ids[i_n] ]);
//
//		AMigInfo[destination].region_2N.push_back(lids);
//
//		/* the region can be removed from the mesh now */
//		AMesh.deleteRegion(current_region);
//	}
//
//	for(node_it= AMesh.nodes_begin();node_it!=node_it_end;node_it++)
//	{
//		/* it is not a node to be moving */
//		if(!AMesh.isMarked(*node_it,AMoveMark))
//			continue;
//
//		bool keep_local = false;
//		Node* current_node = *node_it;
//		/* on regarde si le noeud doit etre conserve localement */
//		std::vector<TInt> dest_partitions = ANodeFilters[current_node->getID()];
//		if(std::find(dest_partitions.begin(),dest_partitions.end(),partition_id)==dest_partitions.end())
//			AMesh.deleteNode(current_node);
//	}
//
//	/* before sending we remove the artificial interface with itself */
//	std::set<Node*> artificial_node_interface;
//	AMesh.getInterface(artificial_node_interface,partition_id);
//	for(std::set<Node*>::iterator it = artificial_node_interface.begin();
//		it!=artificial_node_interface.end();it++)
//	{
//		AMesh.removeFromInterface(*it,partition_id);
//	}
//
//	/* before sending we remove the artificial interface with itself */
//	std::set<Face*> artificial_face_interface;
//	AMesh.getInterface(artificial_face_interface,partition_id);
//	for(std::set<Face*>::iterator it = artificial_face_interface.begin();
//		it!=artificial_face_interface.end();it++)
//	{
//		AMesh.removeFromInterface(*it,partition_id);
//	}
//
//#ifdef _DEBUG
//	{
////		Log::mng()<<"============ BEFORE SENDING =============\n";
////		std::map<TInt, std::set<Node*> > output_interfaces;
////		AMesh.getInterfaces(output_interfaces);
////		std::map<TInt, std::set<Node*> >::iterator it_out = output_interfaces.begin();
////		for(;it_out != output_interfaces.end();it_out++)
////		{
////			Log::mng()<<"INTERFACE WITH "<<it_out->first<<"\n";
////			std::set<Node*> nodes = it_out->second;
////			std::set<Node*>::iterator it_nodes = nodes.begin();
////			for(;it_nodes!=nodes.end();it_nodes++){
////				Log::mng()<<(*it_nodes)->getID()<<" ";
////			}
////			Log::mng()<<"\n";
////		}
//	}
//#endif //_DEBUG
//
//}
///*----------------------------------------------------------------------------*/
///*
// * Ancienne version où l'on calcule toutes les migrations une seule fois et on
// * les stocke. Avantage un seul parcours, inconvénient stockage mémoire
// */
////
////void MigrationTool::
////prepareCellMigration(Mesh& AMesh,
////     const std::vector<TInt>& AFilter,
////	 std::vector<Face*>& AExchangeFaces,
////	 std::vector<std::vector<TInt> >& ANodeFilters,
////	 std::vector<std::vector<TInt> >& AEdgeFilters,
////	 std::map<TInt, std::set<Node*> >& ANodeInterfaces,
////	 std::map<TInt, std::set<Edge*> >& AEdgeInterfaces,
////	 DistributedMemoryManager::MigrationStruct* AMigInfo,
////	 const TInt& ANbMigInfo,
////	 const int& AMoveMark) const
////{
////	TInt partition_id = AMesh.getPartID();
////
////	/* we begin by initializing the migration info. Note that each mesh send and
////	 * receive from all its neighbours in order to avoid blocking cases.
////	 */
////	// INUTILE POUR LE DECOUPAGE
//////	for(std::map<TInt, std::set<Node*> >::iterator it = ANodeInterfaces.begin();
//////		it!=ANodeInterfaces.end();it++){
//////		DistributedMemoryManager::MigrationStruct local_struct;
//////		local_struct.init();
//////		AMigInfo[it->first] = local_struct;
//////
//////	}
////
////	/* we have an id filter on the nodes to get a local reference to the nodes for every
////	 * partition we export cells. Expensive because depending on the number of partitions*/
////	std::vector<TInt> id_filter[ANbMigInfo];
////	TInt id_increment[ANbMigInfo];
////	for(int i=0;i<ANbMigInfo;i++){
////		AMesh.createFilter(0,id_filter[i]);
////		id_increment[i]=0;
////	}
////
////	/* We go through the nodes that can be moved. For that we go through all the nodes and
////	 * we only consider nodes marked with AMoveMark */
////	typename Mesh::nodes_iterator node_it 	 = AMesh.nodes_begin();
////	typename Mesh::nodes_iterator node_it_end = AMesh.nodes_end();
////	for(;node_it!=node_it_end;node_it++)
////	{
////		/* it is not a node to be moving */
////		if(!AMesh.isMarked(*node_it,AMoveMark))
////			continue;
////
////		Node* current_node = *node_it;
////
////
////		/* on commence par regarder avec qui le noeud est partage et s'il doit etre
////		 * conserve localement */
////		std::vector<TInt> sharing_p = ANodeFilters[(*node_it)->getID()];
////
////		const int sharing_p_size = sharing_p.size();
////		for(unsigned int i_dest = 0; i_dest<sharing_p_size;i_dest++)
////		{
////			if(sharing_p[i_dest]==partition_id)
////				continue;
////
////			TInt destination = sharing_p[i_dest];
////			(id_filter[destination])[current_node->getID()] = (id_increment[destination])++;
////			AMigInfo[destination].nb_nodes++;
////			AMigInfo[destination].node_ids.push_back(current_node->getID());
////
////			AMigInfo[destination].node_coords.push_back(current_node->getX());
////			AMigInfo[destination].node_coords.push_back(current_node->getY());
////			AMigInfo[destination].node_coords.push_back(current_node->getZ());
////
////
////			/* on transmet les informations concernant les noeuds partagés. Un noeud
////			 * est partagé s'il appartient au moins à un filtre, ou qu'il appartient
////			 * à un filtre et une interface du maillage*/
//////			if(sharing_p.size()>1)
//////			{
////				for(unsigned int j_dest = 0; j_dest<sharing_p.size();j_dest++)
////					AMigInfo[destination].node_sharing[sharing_p[j_dest]].push_back(current_node->getID());
//////			}
////			if(AMesh.isShared(current_node))
////			{
////				std::vector<TInt> shared_partitions;
////				AMesh.getPartitionsSharing(current_node,shared_partitions);
////				//migration_info[destination].node_sharing[sharing_p[0]].push_back(current_node->getID());
////				for(unsigned int j_dest = 0; j_dest<shared_partitions.size();j_dest++)
////					AMigInfo[destination].node_sharing[shared_partitions[j_dest]].push_back(current_node->getID());
////			}
////		}
////
////
////	}
////	typename Mesh::edges_iterator edge_it 	 = AMesh.edges_begin();
////	typename Mesh::edges_iterator edge_it_end = AMesh.edges_end();
////	for(;edge_it!=edge_it_end;edge_it++)
////	{
////		/* it is not an edge to be moving */
////		if(!AMesh.isMarked(*edge_it,AMoveMark))
////			continue;
////
////		Edge* current_edge= *edge_it;
////		bool keep_local =false;
////		/* on commence par regarder avec qui l'arete est partagee et si elle doit etre
////		 * conserve localement */
////		std::vector<TInt> sharing_p = AEdgeFilters[current_edge->getID()];
////
////		for(unsigned int i_dest = 0; i_dest<sharing_p.size();i_dest++)
////		{
////			if(sharing_p[i_dest]==partition_id)
////			{
////				keep_local = true;
////				continue;
////			}
////
////			TInt destination = sharing_p[i_dest];
////			AMigInfo[destination].nb_edges++;
////			AMigInfo[destination].edge_ids.push_back(current_edge->getID());
////
////			std::vector<id> edge_node_ids = current_edge->getNodeIDs();
////			AMigInfo[destination].edge_2N.push_back((id_filter[destination])[edge_node_ids[0]]);
////			AMigInfo[destination].edge_2N.push_back((id_filter[destination])[edge_node_ids[1]]);
////
////
////			/* on transmet les informations concernant les aretes partagés. Une arete
////			 * est partagé si elle appartient au moins à un filtre, ou qu'il appartient
////			 * à un filtre et une interface du maillage*/
//////			if(sharing_p.size()>1)
//////			{
////				for(unsigned int j_dest = 0; j_dest<sharing_p.size();j_dest++)
////					AMigInfo[destination].edge_sharing[sharing_p[j_dest]].push_back(current_edge->getID());
//////			}
////			if(AMesh.isShared(current_edge))
////			{
////				std::vector<TInt> shared_partitions;
////				AMesh.getPartitionsSharing(current_edge,shared_partitions);
////				//migration_info[destination].node_sharing[sharing_p[0]].push_back(current_node->getID());
////				for(unsigned int j_dest = 0; j_dest<shared_partitions.size();j_dest++)
////					AMigInfo[destination].edge_sharing[shared_partitions[j_dest]].push_back(current_edge->getID());
////			}
////		}
////
////		/* the edge can be removed from the mesh now */
////		if(!keep_local)
////			AMesh.deleteEdge(current_edge);
////
////	}
////
////
////	for(unsigned i_f=0;i_f<AExchangeFaces.size();i_f++)
////	{
////		Face* current_face = AExchangeFaces[i_f];
////		TInt destination = AFilter[current_face->getID()];
////		AMigInfo[destination].nb_faces++;
////		AMigInfo[destination].face_ids.push_back(current_face->getID());
////		AMigInfo[destination].face_type.push_back(current_face->getType());
////		std::vector<id> face_node_ids = current_face->getNodeIDs();
////		std::vector<TInt> lids;
////		lids.reserve(face_node_ids.size());
////		for(unsigned int i_n=0;i_n<face_node_ids.size();i_n++)
////			lids.push_back((id_filter[destination])[face_node_ids[i_n] ]);
////
////		AMigInfo[destination].face_2N.push_back(lids);
////
////		/* the face can be removed from the mesh now */
////		AMesh.deleteFace(current_face);
////	}
////
////	for(node_it= AMesh.nodes_begin();node_it!=node_it_end;node_it++)
////	{
////		/* it is not a node to be moving */
////		if(!AMesh.isMarked(*node_it,AMoveMark))
////			continue;
////
////		bool keep_local = false;
////		Node* current_node = *node_it;
////		/* on regarde si le noeud doit etre conserve localement */
////		std::vector<TInt> dest_partitions = ANodeFilters[current_node->getID()];
////		if(std::find(dest_partitions.begin(),dest_partitions.end(),partition_id)==dest_partitions.end())
////			AMesh.deleteNode(current_node);
////	}
////
////	/* before sending we remove the artificial interface with itself */
////	std::set<Node*> artificial_node_interface;
////	AMesh.getInterface(artificial_node_interface,partition_id);
////	for(std::set<Node*>::iterator it = artificial_node_interface.begin();
////		it!=artificial_node_interface.end();it++)
////	{
////		AMesh.removeFromInterface(*it,partition_id);
////	}
////
////	/* before sending we remove the artificial interface with itself */
////	std::set<Edge*> artificial_edge_interface;
////	AMesh.getInterface(artificial_edge_interface,partition_id);
////	for(std::set<Edge*>::iterator it = artificial_edge_interface.begin();
////		it!=artificial_edge_interface.end();it++)
////	{
////		AMesh.removeFromInterface(*it,partition_id);
////	}
////
////#ifdef _DEBUG
////	{
//////		Log::mng()<<"============ BEFORE SENDING =============\n";
//////		std::map<TInt, std::set<Node*> > output_interfaces;
//////		AMesh.getInterfaces(output_interfaces);
//////		std::map<TInt, std::set<Node*> >::iterator it_out = output_interfaces.begin();
//////		for(;it_out != output_interfaces.end();it_out++)
//////		{
//////			Log::mng()<<"INTERFACE WITH "<<it_out->first<<"\n";
//////			std::set<Node*> nodes = it_out->second;
//////			std::set<Node*>::iterator it_nodes = nodes.begin();
//////			for(;it_nodes!=nodes.end();it_nodes++){
//////				Log::mng()<<(*it_nodes)->getID()<<" ";
//////			}
//////			Log::mng()<<"\n";
//////		}
////	}
////#endif //_DEBUG
////
////}
///*----------------------------------------------------------------------------*/
////
////void MigrationTool::
////performCellMigration(const TInt& ARank,
////	std::map<TInt, DistributedMemoryManager::MigrationStruct>& ASendInfo,
////	std::map<TInt, DistributedMemoryManager::MigrationStruct>& ARecvInfo) const
////{
////	std::map<TInt, DistributedMemoryManager::MigrationStruct>::iterator it_mig = ASendInfo.begin();
////
////	for(;it_mig!=ASendInfo.end();it_mig++)
////	{
////		DistributedMemoryManager::MigrationStruct to_receive, to_send;
////		TInt dest_rank = it_mig->first;
////#ifdef _DEBUG
////		to_send = it_mig->second;
////		Log::mng()<<"================================== \n";
////		Log::mng()<<"Send to "<<dest_rank<<"\n";
////		Log::mng()<<"nb faces, nb nodes = "<<it_mig->second.nb_faces<<", "<<it_mig->second.nb_nodes<<"\n";
////		Log::mng()<<"nodes sent"<<"\n";
////		std::vector<id> send_node_ids = it_mig->second.node_ids;
////		for(unsigned int i=0; i<send_node_ids.size();i++)
////			Log::mng()<<" "<<send_node_ids[i]<<"\n";
////		Log::mng()<<"\n";
////		Log::mng()<<"sharing node info sent"<<"\n";
////		std::map<TInt, std::set<Node*> >::iterator it_new_node_interface;
////		for(it_new_node_interface  = it_mig->second.node_sharing.begin();
////			it_new_node_interface != it_mig->second.node_sharing.end();
////			it_new_node_interface++)
////		{
////			TInt share_partition = it_new_node_interface->first;
////			Log::mng()<<" "<<share_partition<<" -> "<<"\n";
////			std::vector<id> share_nodes = it_new_node_interface->second;
////			for(unsigned int i=0; i<share_nodes.size();i++){
////				Log::mng()<<share_nodes[i]<<" ";
////			}
////			Log::mng()<<"\n";
////		}
////#endif
////		DistributedMemoryManager* mpi_manager = DistributedMemoryManager::instance();
////		if(ARank>dest_rank)
////		{
////			//send then recv
////			mpi_manager->sendPacked(dest_rank, it_mig->second);
////			mpi_manager->receivePacked(dest_rank,to_receive);
////		}
////		else {
////			//recv then send
////			mpi_manager->receivePacked(dest_rank, to_receive);
////			mpi_manager->sendPacked(dest_rank, it_mig->second);
////		}
////		ARecvInfo[dest_rank] = to_receive;
////
////#ifdef _DEBUG
////
////		Log::mng()<<"          ----------              "<<"\n";
////		Log::mng()<<"Received from "<<dest_rank<<"\n";
////		Log::mng()<<"nb faces, nb nodes = "<<to_receive.nb_faces<<", "<<to_receive.nb_nodes<<"\n";
////		Log::mng()<<"\n";
////		Log::mng()<<"nodes received"<<"\n";
////		std::vector<id> recv_node_ids = to_receive.node_ids;
////		for(unsigned int i=0; i<recv_node_ids.size();i++)
////			Log::mng()<<" "<<recv_node_ids[i]<<"\n";
////		Log::mng()<<"\n";
////		Log::mng()<<"sharing node info received"<<"\n";
////		for(it_new_node_interface  = to_receive.node_sharing.begin();
////			it_new_node_interface != to_receive.node_sharing.end();
////			it_new_node_interface++)
////		{
////			TInt share_partition = it_new_node_interface->first;
////			Log::mng()<<" "<<share_partition<<" -> "<<"\n";
////			std::vector<id> share_nodes = it_new_node_interface->second;
////			for(unsigned int i=0; i<share_nodes.size();i++){
////				Log::mng()<<share_nodes[i]<<" ";
////			}
////			Log::mng()<<"\n";
////		}
////#endif // _DEBUG
////	}
////}
///*----------------------------------------------------------------------------*/
//
//void MigrationTool::
//performCellMigrationForIDSplit_Send(const TInt& ARank,
//		TInt *ADestPartitions,
//	DistributedMemoryManager::MigrationStruct* ASendInfo, const TInt& ANbSend) const
//{
//	DistributedMemoryManager* mpi_manager= DistributedMemoryManager::instance();
//	TInt rank 		   = mpi_manager->rank();
//	for(int i =0; i<ANbSend;i++)
//	{
//		DistributedMemoryManager::MigrationStruct  to_send=ASendInfo[i];
//		TInt dest_rank = ADestPartitions[i];
//		if(dest_rank==rank)
//			continue;
//#ifdef _DEBUG
////		to_send = it_mig->second;
////		Log::mng()<<"================================== \n";
////		Log::mng()<<"Send to "<<dest_rank<<"\n";
////		Log::mng()<<"nb faces, nb nodes = "<<it_mig->second.nb_faces<<", "<<it_mig->second.nb_nodes<<"\n";
////		Log::mng()<<"faces sent"<<"\n";
////		std::vector<id> send_face_ids = it_mig->second.face_ids;
////		for(unsigned int i=0; i<send_face_ids.size();i++)
////			Log::mng()<<" "<<send_face_ids[i]<<"\n";
////		Log::mng()<<"\n";
////		Log::mng()<<"nodes sent"<<"\n";
////		std::vector<id> send_node_ids = it_mig->second.node_ids;
////		for(unsigned int i=0; i<send_node_ids.size();i++)
////			Log::mng()<<" "<<send_node_ids[i]<<"\n";
////		Log::mng()<<"\n";
////		Log::mng()<<"sharing node info sent"<<"\n";
////		std::map<TInt, std::vector<id> >::iterator it_new_node_interface;
////		for(it_new_node_interface  = it_mig->second.node_sharing.begin();
////			it_new_node_interface != it_mig->second.node_sharing.end();
////			it_new_node_interface++)
////		{
////			TInt share_partition = it_new_node_interface->first;
////			Log::mng()<<" "<<share_partition<<" -> "<<"\n";
////			std::vector<id> share_nodes = it_new_node_interface->second;
////			for(unsigned int i=0; i<share_nodes.size();i++){
////				Log::mng()<<share_nodes[i]<<" ";
////			}
////			Log::mng()<<"\n";
////		}
//#endif
//		DistributedMemoryManager* mpi_manager = DistributedMemoryManager::instance();
//		mpi_manager->sendPacked(dest_rank, to_send);
//	}
//}
///*----------------------------------------------------------------------------*/
//
//void MigrationTool::
//performCellMigrationForIDSplit_Recv(const TInt& ARank,
//	DistributedMemoryManager::MigrationStruct& ARecvInfo) const
//{
//	DistributedMemoryManager* mpi_manager = DistributedMemoryManager::instance();
//	std::cout<<mpi_manager->rank()<<" - start recv"<<std::endl;
//	mpi_manager->receivePacked(ARank,ARecvInfo);
//	std::cout<<mpi_manager->rank()<<" - end recv"<<std::endl;
//
//#ifdef _DEBUG
////
////		Log::mng()<<"          ----------              "<<"\n";
////		Log::mng()<<"Received from "<<ARank<<"\n";
////		Log::mng()<<"nb faces, nb nodes = "<<ARecvInfo.nb_faces<<", "<<ARecvInfo.nb_nodes<<"\n";
////		Log::mng()<<"\n";
////		Log::mng()<<"faces received"<<"\n";
////		std::vector<id> send_face_ids = ARecvInfo.face_ids;
////		for(unsigned int i=0; i<send_face_ids.size();i++)
////			Log::mng()<<" "<<send_face_ids[i]<<"\n";
////		Log::mng()<<"\n";
////
////		Log::mng()<<"nodes received"<<"\n";
////		std::vector<id> recv_node_ids = ARecvInfo.node_ids;
////		for(unsigned int i=0; i<recv_node_ids.size();i++)
////			Log::mng()<<" "<<recv_node_ids[i]<<"\n";
////		Log::mng()<<"\n";
////		Log::mng()<<"sharing node info received"<<"\n";
////		std::map<TInt, std::vector<id> >::iterator it_new_node_interface;
////		for(it_new_node_interface  = ARecvInfo.node_sharing.begin();
////			it_new_node_interface != ARecvInfo.node_sharing.end();
////			it_new_node_interface++)
////		{
////			TInt share_partition = it_new_node_interface->first;
////			Log::mng()<<" "<<share_partition<<" -> "<<"\n";
////			std::vector<id> share_nodes = it_new_node_interface->second;
////			for(unsigned int i=0; i<share_nodes.size();i++){
////				Log::mng()<<share_nodes[i]<<" ";
////			}
////			Log::mng()<<"\n";
////		}
//#endif // _DEBUG
//
//}
//
///*----------------------------------------------------------------------------*/
//
//void MigrationTool::
//updateAfterCellMigration(Mesh& AMesh,
//	DistributedMemoryManager::MigrationStruct& ARecvInfo) const
//{
//
//	TInt partition_id = AMesh.getPartID();
//	std::cout<<partition_id<<" - start update"<<std::endl;
//	std::map<id, std::vector<TInt> > new_nodes_sharing_info;
//
//	std::map<TInt, std::vector<id> >::iterator it_new_node_interface;
//	for(it_new_node_interface  = ARecvInfo.node_sharing.begin();
//		it_new_node_interface != ARecvInfo.node_sharing.end();
//		it_new_node_interface++)
//	{
//		TInt share_partition = it_new_node_interface->first;
//		std::vector<id> share_nodes = it_new_node_interface->second;
//		for(unsigned int i=0; i<share_nodes.size();i++)
//			new_nodes_sharing_info[share_nodes[i]].push_back(share_partition);
//	}
//
//
//	/* Creation of the imported nodes and faces if they do not still exist on this proc */
//	std::set<Node*> node_interfaces;
//	std::set<Edge*> edge_interfaces;
//	// WARNING REPLACE 0 BY SENDER RANK
//	AMesh.getInterface(node_interfaces, 0);
//	AMesh.getInterface(edge_interfaces, 1);
//
//	std::vector<id> link_to_new_nodes;
//	link_to_new_nodes.resize(ARecvInfo.nb_nodes);
//
//	for(unsigned int i=0;i<ARecvInfo.nb_nodes;i++)
//	{
//		id gid = ARecvInfo.node_ids[i];
//		/* We only add the nodes that are not yet in this domain. Warning, it should
//		 * be possible to check only the interface between this domain and the coming
//		 * from domain. */
//		bool found = false;
//		for(std::set<Node*>::iterator it_loc=node_interfaces.begin();
//			it_loc != node_interfaces.end() && !found; it_loc++)
//			if((*it_loc)->getID()==gid)
//			{
//				found = true;
//				link_to_new_nodes[i] = (*it_loc)->getID();
//				std::vector<TInt> shared_partitions = new_nodes_sharing_info[gid];
//				for(unsigned int i=0;i<shared_partitions.size();i++)
//				{
//					if(shared_partitions[i]!=partition_id)
////								AMesh.removeFromInterface(*it_loc,recv_from);
////							else
//						AMesh.shareWith(*it_loc,shared_partitions[i]);
//				}
//			}
//
//		if(!found)
//		{
//			TCoord x = ARecvInfo.node_coords[3*i];
//			TCoord y = ARecvInfo.node_coords[3*i+1];
//			TCoord z = ARecvInfo.node_coords[3*i+2];
//			Node* n = AMesh.newNode(x,y,z);
//	//		n->setID(gid);
//			link_to_new_nodes[i] = n->getID();
//			std::vector<TInt> shared_partitions = new_nodes_sharing_info[gid];
//			for(unsigned int i=0;i<shared_partitions.size();i++)
//			{
//				if(shared_partitions[i]!=partition_id)
////						AMesh.shareWith(n,recv_from);
////					else
//					AMesh.shareWith(n,shared_partitions[i]);
//			}
//		}
//
//	}
//
//
//	std::map<id, std::vector<TInt> > new_edges_sharing_info;
//
//	std::map<TInt, std::vector<id> >::iterator it_new_edge_interface;
//	for(it_new_edge_interface  = ARecvInfo.edge_sharing.begin();
//		it_new_edge_interface != ARecvInfo.edge_sharing.end();
//		it_new_edge_interface++)
//	{
//		TInt share_partition = it_new_edge_interface->first;
//		std::vector<id> share_edges = it_new_edge_interface->second;
//		for(unsigned int i=0; i<share_edges.size();i++)
//			new_edges_sharing_info[share_edges[i]].push_back(share_partition);
//	}
//
//	for(unsigned int i = 0; i < ARecvInfo.nb_edges;i++)
//	{
//		id gid = ARecvInfo.edge_ids[i];
//		/* We only add the nodes that are not yet in this domain. Warning, it should
//		 * be possible to check only the interface between this domain and the coming
//		 * from domain. */
//		bool found = false;
//		Edge* e =0;
//		for(std::set<Edge*>::iterator it_loc=edge_interfaces.begin();
//			it_loc != edge_interfaces.end() && !found; it_loc++)
//			if((*it_loc)->getID()==gid)
//			{
//				found = true;
//				e = *it_loc;
//			}
//
//		if(!found)
//		{
//			e = AMesh.newEdge(link_to_new_nodes[ARecvInfo.edge_2N[2*i  ] ],
//							  link_to_new_nodes[ARecvInfo.edge_2N[2*i+1] ]);
////			e->setID(gid);
//		}
//
//		std::vector<TInt> shared_partitions = new_edges_sharing_info[gid];
//		for(unsigned int i=0;i<shared_partitions.size();i++)
//		{
//			if(shared_partitions[i]!=partition_id)
//				AMesh.shareWith(e,shared_partitions[i]);
//		}
//
//	}
//
//
//	for(unsigned int i = 0; i < ARecvInfo.nb_faces;i++)
//	{
//		id gid = ARecvInfo.face_ids[i];
//		ECellType typ = ARecvInfo.face_type[i];
//		std::vector<TInt> loc_node_ids = ARecvInfo.face_2N[i];
//		std::vector<Node*> nodes_of_face;
//		nodes_of_face.clear();
//		for(unsigned int i_n = 0; i_n <loc_node_ids.size(); i_n++)
//			nodes_of_face.push_back(AMesh.getLNode(link_to_new_nodes[loc_node_ids[i_n] ]));
//
//		Face* f = AMesh.newFace(nodes_of_face);
//	//	f->setID(gid);
//	}
//
//
//
//	MeshDoctor doctor(AMesh);
//	doctor.buildEdgesAndX2E();
//
//#ifdef _DEBUG
//	{
////		Log::mng()<<"==== AFTER UPDATE ====\n";
////		std::map<TInt, std::set<Node*> > output_interfaces;
////		AMesh.getInterfaces(output_interfaces);
////		std::map<TInt, std::set<Node*> >::iterator it_out = output_interfaces.begin();
////		for(;it_out != output_interfaces.end();it_out++)
////		{
////			Log::mng()<<"NODE INTERFACE WITH "<<it_out->first<<"\n";
////			std::set<Node*> nodes = it_out->second;
////			std::set<Node*>::iterator it_nodes = nodes.begin();
////			for(;it_nodes!=nodes.end();it_nodes++){
////				Log::mng()<<(*it_nodes)->getID()<<" ";
////			}
////			Log::mng()<<"\n";
////		}
////		std::map<TInt, std::set<Edge*> > edge_output_interfaces;
////		AMesh.getInterfaces(edge_output_interfaces);
////		std::map<TInt, std::set<Edge*> >::iterator it_edge_out = edge_output_interfaces.begin();
////		for(;it_edge_out != edge_output_interfaces.end();it_edge_out++)
////		{
////			Log::mng()<<"EDGE INTERFACE WITH "<<it_edge_out->first<<"\n";
////			std::set<Edge*> edges= it_edge_out->second;
////			std::set<Edge*>::iterator it_edges = edges.begin();
////			for(;it_edges!=edges.end();it_edges++){
////				Log::mng()<<(*it_edges)->getID()<<" ";
////			}
////			Log::mng()<<"\n";
////		}
////
//	}
//
//#endif
//	std::cout<<partition_id<<" - end update"<<std::endl;
//}
///*----------------------------------------------------------------------------*/
//
//void MigrationTool::
//updateAfterCellMigration3D(Mesh& AMesh,
//	DistributedMemoryManager::MigrationStruct& ARecvInfo) const
//{
//
//	TInt partition_id = AMesh.getPartID();
//	std::cout<<partition_id<<" - start update"<<std::endl;
//	std::map<id, std::vector<TInt> > new_nodes_sharing_info;
//
//	std::map<TInt, std::vector<id> >::iterator it_new_node_interface;
//	for(it_new_node_interface  = ARecvInfo.node_sharing.begin();
//		it_new_node_interface != ARecvInfo.node_sharing.end();
//		it_new_node_interface++)
//	{
//		TInt share_partition = it_new_node_interface->first;
//		std::vector<id> share_nodes = it_new_node_interface->second;
//		for(unsigned int i=0; i<share_nodes.size();i++)
//			new_nodes_sharing_info[share_nodes[i]].push_back(share_partition);
//	}
//
//
//	/* Creation of the imported nodes and faces if they do not still exist on this proc */
//	std::set<Node*> node_interfaces;
//	std::set<Face*> face_interfaces;
//	// WARNING REPLACE 0 BY SENDER RANK
//	AMesh.getInterface(node_interfaces, 0);
//	AMesh.getInterface(face_interfaces, 2);
//
//	std::vector<id> link_to_new_nodes;
//	link_to_new_nodes.resize(ARecvInfo.nb_nodes);
//
//	for(unsigned int i=0;i<ARecvInfo.nb_nodes;i++)
//	{
//		id gid = ARecvInfo.node_ids[i];
//		/* We only add the nodes that are not yet in this domain. Warning, it should
//		 * be possible to check only the interface between this domain and the coming
//		 * from domain. */
//		bool found = false;
//		for(std::set<Node*>::iterator it_loc=node_interfaces.begin();
//			it_loc != node_interfaces.end() && !found; it_loc++)
//			if((*it_loc)->getID()==gid)
//			{
//				found = true;
//				link_to_new_nodes[i] = (*it_loc)->getID();
//				std::vector<TInt> shared_partitions = new_nodes_sharing_info[gid];
//				for(unsigned int i=0;i<shared_partitions.size();i++)
//				{
//					if(shared_partitions[i]!=partition_id)
//						AMesh.shareWith(*it_loc,shared_partitions[i]);
//				}
//			}
//
//		if(!found)
//		{
//			TCoord x = ARecvInfo.node_coords[3*i];
//			TCoord y = ARecvInfo.node_coords[3*i+1];
//			TCoord z = ARecvInfo.node_coords[3*i+2];
//			Node* n = AMesh.newNode(x,y,z);
////			n->setID(gid);
//			link_to_new_nodes[i] = n->getID();
//			std::vector<TInt> shared_partitions = new_nodes_sharing_info[gid];
//			for(unsigned int i=0;i<shared_partitions.size();i++)
//			{
//				if(shared_partitions[i]!=partition_id)
//					AMesh.shareWith(n,shared_partitions[i]);
//			}
//		}
//
//	}
//
//
//	std::map<id, std::vector<TInt> > new_faces_sharing_info;
//
//	std::map<TInt, std::vector<id> >::iterator it_new_face_interface;
//	for(it_new_face_interface  = ARecvInfo.face_sharing.begin();
//		it_new_face_interface != ARecvInfo.face_sharing.end();
//		it_new_face_interface++)
//	{
//		TInt share_partition = it_new_face_interface->first;
//		std::vector<id> share_faces = it_new_face_interface->second;
//		for(unsigned int i=0; i<share_faces.size();i++)
//			new_faces_sharing_info[share_faces[i]].push_back(share_partition);
//	}
//
//	for(unsigned int i = 0; i < ARecvInfo.nb_faces;i++)
//	{
//		id gid = ARecvInfo.face_ids[i];
//		/* We only add the faces that are not yet in this domain. Warning, it should
//		 * be possible to check only the interface between this domain and the coming
//		 * from domain. */
//		bool found = false;
//		Face* f =0;
//		for(std::set<Face*>::iterator it_loc=face_interfaces.begin();
//			it_loc != face_interfaces.end() && !found; it_loc++)
//			if((*it_loc)->getID()==gid)
//			{
//				found = true;
//				f = *it_loc;
//			}
//
//		if(!found)
//		{
//			/* we create a new face */
//			ECellType typ = ARecvInfo.face_type[i];
//			std::vector<TInt> loc_node_ids = ARecvInfo.face_2N[i];
//			std::vector<Node*> nodes_of_face;
//			nodes_of_face.clear();
//			for(unsigned int i_n = 0; i_n <loc_node_ids.size(); i_n++)
//				nodes_of_face.push_back(AMesh.getLNode(link_to_new_nodes[loc_node_ids[i_n] ]));
//
//			if(nodes_of_face.empty())
//				std::cout<<"EMPTY "<<loc_node_ids.size()<<std::endl;
//			f = AMesh.newFace(nodes_of_face);
////			f->setID(gid);
//
//		}
//
//		std::vector<TInt> shared_partitions = new_faces_sharing_info[gid];
//		for(unsigned int i=0;i<shared_partitions.size();i++)
//		{
//			if(shared_partitions[i]!=partition_id)
//				AMesh.shareWith(f,shared_partitions[i]);
//		}
//	}
//
//
//	for(unsigned int i = 0; i < ARecvInfo.nb_regions;i++)
//	{
//		id gid = ARecvInfo.region_ids[i];
//		ECellType typ = ARecvInfo.region_type[i];
//		std::vector<TInt> loc_node_ids = ARecvInfo.region_2N[i];
//		std::vector<Node*> nodes_of_face;
//		nodes_of_face.clear();
//		for(unsigned int i_n = 0; i_n <loc_node_ids.size(); i_n++)
//			nodes_of_face.push_back(AMesh.getLNode(link_to_new_nodes[loc_node_ids[i_n] ]));
//
//		Region* r=0;
//		switch(typ){
//		case GMDS_HEX:
//			r = AMesh.newHex(nodes_of_face[0],nodes_of_face[1],nodes_of_face[2],nodes_of_face[3],
//							 nodes_of_face[4],nodes_of_face[5],nodes_of_face[6],nodes_of_face[7]);
//			break;
//		default:
//			throw GMDSException("Cell type not handled!");
//		}
////		r->setID(gid);
//	}
//
//	/* mise a jour des groupes de cellules */
//	std::map<std::string, std::vector<id> >::iterator it_group;
//
//	/* groupes de noeuds */
//	for(it_group  = ARecvInfo.node_groups.begin();it_group!= ARecvInfo.node_groups.end();
//		it_group++)
//	{
//		try{
//		typename Mesh::cloud& gp= AMesh.getCloud(it_group->first);
//		std::vector<id> group_nodes = it_group->second;
//		const unsigned int group_nodes_size = group_nodes.size();
//		for(unsigned int i=0; i<group_nodes_size;i++)
//			gp.add(AMesh.getLNode(AMesh.geid(group_nodes[i],0)));
//		}
//		catch (GMDSException& e){
//			typename Mesh::cloud&  gp = AMesh.newCloud(it_group->first);
//			std::vector<id> group_nodes = it_group->second;
//			std::cout<<gp.name()<<std::endl;
//			const unsigned int group_nodes_size = group_nodes.size();
//			for(unsigned int i=0; i<group_nodes_size;i++)
//				std::cout<<group_nodes[i]<<" ";
//			std::cout<<std::endl;
//
//			for(unsigned int i=0; i<group_nodes_size;i++)
//			{
//				id index = group_nodes[i];
//				id id = AMesh.geid(index,0);
//				Node *n = AMesh.getLNode(id);
//				gp.add(n);
//			}
//		}
//	}
//
//	/* groupes de faces */
//	for(it_group  = ARecvInfo.face_groups.begin();it_group!= ARecvInfo.face_groups.end();
//		it_group++)
//	{
//		try{
//		typename Mesh::surface& gp= AMesh.getSurface(it_group->first);
//		std::vector<id> group_faces = it_group->second;
//		const unsigned int group_faces_size = group_faces.size();
//		for(unsigned int i=0; i<group_faces_size;i++)
//			gp.add(AMesh.getLFace(AMesh.geid(group_faces[i],2)));
//		}
//		catch (GMDSException& e){
//			typename Mesh::surface&  gp = AMesh.newSurface(it_group->first);
//			std::vector<id> group_faces = it_group->second;
//			const unsigned int group_faces_size = group_faces.size();
//			for(unsigned int i=0; i<group_faces_size;i++)
//				gp.add(AMesh.getLFace(AMesh.geid(group_faces[i],2)));
//		}
//	}
//
//	MeshDoctor doctor(AMesh);
//	doctor.buildFacesAndR2F();
//
//#ifdef _DEBUG
//	{
////		Log::mng()<<"==== AFTER UPDATE ====\n";
////		std::map<TInt, std::set<Node*> > output_interfaces;
////		AMesh.getInterfaces(output_interfaces);
////		std::map<TInt, std::set<Node*> >::iterator it_out = output_interfaces.begin();
////		for(;it_out != output_interfaces.end();it_out++)
////		{
////			Log::mng()<<"NODE INTERFACE WITH "<<it_out->first<<"\n";
////			std::set<Node*> nodes = it_out->second;
////			std::set<Node*>::iterator it_nodes = nodes.begin();
////			for(;it_nodes!=nodes.end();it_nodes++){
////				Log::mng()<<(*it_nodes)->getID()<<" ";
////			}
////			Log::mng()<<"\n";
////		}
////		std::map<TInt, std::set<Edge*> > edge_output_interfaces;
////		AMesh.getInterfaces(edge_output_interfaces);
////		std::map<TInt, std::set<Edge*> >::iterator it_edge_out = edge_output_interfaces.begin();
////		for(;it_edge_out != edge_output_interfaces.end();it_edge_out++)
////		{
////			Log::mng()<<"EDGE INTERFACE WITH "<<it_edge_out->first<<"\n";
////			std::set<Edge*> edges= it_edge_out->second;
////			std::set<Edge*>::iterator it_edges = edges.begin();
////			for(;it_edges!=edges.end();it_edges++){
////				Log::mng()<<(*it_edges)->getID()<<" ";
////			}
////			Log::mng()<<"\n";
////		}
////
//	}
//
//#endif
//	std::cout<<partition_id<<" - end update"<<std::endl;
//}
///*----------------------------------------------------------------------------*/
//
//void MigrationTool::
//exchangeCells(Mesh& AMesh, const std::vector<TInt>& AFilter)
//{
//	std::set<TInt> 		destination_partitions;
//	std::vector<Face*> 	faces_to_be_moved;
//
//
//	for(unsigned int i =0; i<AFilter.size();i++)
//		destination_partitions.insert(AFilter[i]);
//
//	/* we get the partition number of the current mesh */
//	int partition_id = AMesh.getPartID();
//
//	/* for each target partition, we create a filter indicating
//	 * if this cell must be migrated to the target partition.
//	 * 	- the mak key is the partition number
//	 * 	- in the associated vector is store 0 in item i if it must not go there,
//	 * 	  1 otherwise
//	 */
//
//	int mv = AMesh.getNewMark();
//	std::vector<std::vector<TInt> > node_filters;
//	std::vector<std::vector<TInt> > edge_filters;
//	initFilters(AMesh, AFilter, destination_partitions,
//				node_filters, edge_filters, faces_to_be_moved, mv);
//
//	/* Partitions interface update */
//	updateInterfacesBeforeCellExchange(AMesh,node_filters);
////
////
////	/* creation of the internal interfaces */
////	/* we traverse the mesh nodes and  we share the involved node */
////	createInternalInterfaces(AMesh,node_filters);
////
////
////	/* interfaces are now correct, we can send the cells */
////	std::map<TInt, DistributedMemoryManager::MigrationStruct> migration_info;
////
////	prepareCellMigration(AMesh,AFilter, faces_to_be_moved, nodes_to_be_moved,
////						 node_filters, initial_node_interfaces, migration_info);
////	/* Now we can send our structure. Be careful to send nothing to all the neighbours even if
////	 * the message is empty. Otherwise, you will be stuck !! */
////	std::map<TInt, DistributedMemoryManager::MigrationStruct> recv_info;
////
////	performCellMigration(partition_id,migration_info,recv_info);
////
////	/* we get the data to modify mesh interfaces */
////	updateAfterCellMigration(AMesh,recv_info);
//
//}
///*----------------------------------------------------------------------------*/
//
//void MigrationTool::
//sequentialIDSplit(Mesh& AMesh) const
//{
//	if(MeshDescriptor::dimension==2)
//		sequentialIDSplit2D(AMesh);
//	else if(MeshDescriptor::dimension==3)
//		sequentialIDSplit3D(AMesh);
//	else
//		throw GMDSException("ERROR: invalid mesh dimension");
//
//}
///*----------------------------------------------------------------------------*/
//
//void MigrationTool::
//newSequentialIDSplit(Mesh& AMesh) const
//{
//	if(MeshDescriptor::dimension==2)
//		newSequentialIDSplit2D(AMesh);
//	else if(MeshDescriptor::dimension==3)
//		newSequentialIDSplit3D(AMesh);
//	else
//		throw GMDSException("ERROR: invalid mesh dimension");
//
//}
///*----------------------------------------------------------------------------*/
//
//void MigrationTool::
//newSequentialIDSplit2D(Mesh& AMesh) const
//{
//	Timer time_start;
//	DistributedMemoryManager* mpi_manager= DistributedMemoryManager::instance();
//
//	TInt nb_partitions = mpi_manager->nbPartitions();
//	TInt rank 		   = mpi_manager->rank();
//	TInt partition_id = rank;
//
//#ifdef _DEBUG
//	Log::mng()<<"NB partitions= "<<nb_partitions<<"\n";
//#endif
//	Log::mng()<<"Sequential cut on "<<rank<<"\n";
//	if (!platform::getEnvironmentVariable("INFO_CONSO").empty())
//		Log::mng()<<rank<<" - Consommation memoire en octets: "<<platform::getMemoryUsed()<<"\n";
//
//	if(rank==0){
//		Timer t1;
//		// id partitioning variable
//		Variable<TInt> *var_dest = AMesh.newVariable<TInt>(GMDS_FACE,"part");
//		int nb_part = mpi_manager->nbPartitions();
//		typename Mesh::cells_iterator it_f = AMesh.cells_begin();
//		typename Mesh::cells_iterator it_fe = AMesh.cells_end();
//		int cell_index=0;
//
//		int nb_cells_per_part = AMesh.getNbCells()/nb_part;
//		int part_index=0;
//		for(;part_index<nb_part;part_index++)
//			for(int i=0;i<nb_cells_per_part;i++){
//				(*var_dest)[(*it_f)->getID()] = part_index;
//				it_f++;
//			}
//
//		part_index--;
//		for(;it_f!=it_fe;it_f++){
//			(*var_dest)[(*it_f)->getID()] = part_index;
//		}
//
//		/* MESH EXTRACTION ON PART 0	 */
//		if (!platform::getEnvironmentVariable("INFO_CONSO").empty())
//			Log::mng()<<rank<<" - Consommation memoire en octets: "<<platform::getMemoryUsed()<<"\n";
//
//		Extractor extract(AMesh);
//		std::vector<Mesh > subMeshes;
//		Timer t2;
//		Log::mng()<<rank<<"/ Time to initialize the partitioning varible: "<<t2-t1<<"\n";
//		extract.extractNoGhostWithNativeSharedInfoTopLevelRepartitionExp2(GMDS_FACE,"part",subMeshes);
//		Timer t3;
//		Log::mng()<<rank<<"/ Time to extract submeshes: "<<t3-t2<<"\n";
//		if (!platform::getEnvironmentVariable("INFO_CONSO").empty())
//			Log::mng()<<rank<<" - Consommation memoire en octets: "<<platform::getMemoryUsed()<<"\n";
//
////		for(int i=0;i<subMeshes.size();i++){
////			LimaWriter writer(subMeshes[i]);
////			//	VTKWriter writer(subMeshes[i]);
////			std::ostringstream file_name;
////			file_name<<"homePOYOP/mesh/new_parallel_test2D_sub"<<i<<".mli";
////			writer.write(file_name.str(),N|F);
////		}
//
//		/* SUBMESHES SERIALIZATION */
//		for(int i=0;i<subMeshes.size();i++){
//			Mesh& subM = subMeshes[i];
//			std::ostringstream sub_stream;
//			subM.serialize(sub_stream);
//
//			std::string str = sub_stream.str();
//
//			int sub_size   = str.size()*sizeof(char);
//			const char* sub_char = str.c_str();
//
//			mpi_manager->send(i+1,0, (void*)sub_char,sub_size);
//
//		}
//		Timer t4;
//		Log::mng()<<rank<<"/ Time to serialize all submeshes: "<<t4-t3<<"\n";
//		if (!platform::getEnvironmentVariable("INFO_CONSO").empty())
//			Log::mng()<<rank<<" - Consommation memoire en octets: "<<platform::getMemoryUsed()<<"\n";
//
//
//	}
//	else{
//		if (!platform::getEnvironmentVariable("INFO_CONSO").empty())
//			Log::mng()<<rank<<" - Consommation memoire en octets: "<<platform::getMemoryUsed()<<"\n";
//
//		void *s;
//		int t=0;
//		s = mpi_manager->receive(0,0,t);
//		std::cout<<rank<<" - receive "<<t<<std::endl;
//
//		std::istrstream sub_stream2((char*)s,t);
//
//		AMesh.unserialize(sub_stream2);
//	}
//	if (!platform::getEnvironmentVariable("INFO_CONSO").empty())
//		Log::mng()<<rank<<" - Consommation memoire en octets: "<<platform::getMemoryUsed()<<"\n";
//
//}
///*----------------------------------------------------------------------------*/
//
//void MigrationTool::
//newSequentialIDSplit3D(Mesh& AMesh) const
//{
//	Timer time_start;
//	DistributedMemoryManager* mpi_manager= DistributedMemoryManager::instance();
//
//	TInt nb_partitions = mpi_manager->nbPartitions();
//	TInt rank 		   = mpi_manager->rank();
//	TInt partition_id = rank;
//
//#ifdef _DEBUG
//	Log::mng()<<"NB partitions= "<<nb_partitions<<"\n";
//#endif
//	Log::mng()<<"Sequential cut on "<<rank<<"\n";
//	if (!platform::getEnvironmentVariable("INFO_CONSO").empty())
//		Log::mng()<<rank<<" - Consommation memoire en octets: "<<platform::getMemoryUsed()<<"\n";
//
//	if(rank==0){
//		Timer t1;
//		// id partitioning variable
//		Variable<TInt> *var_dest = AMesh.template newVariable<TInt>(GMDS_REGION,"part");
//		int nb_part = mpi_manager->nbPartitions();
//		typename Mesh::cells_iterator it_f = AMesh.cells_begin();
//		typename Mesh::cells_iterator it_fe = AMesh.cells_end();
//		int cell_index=0;
//
//		int nb_cells_per_part = AMesh.getNbCells()/nb_part;
//		int part_index=0;
//		for(;part_index<nb_part;part_index++)
//			for(int i=0;i<nb_cells_per_part;i++){
//				(*var_dest)[(*it_f)->getID()] = part_index;
//				it_f++;
//			}
//
//		part_index--;
//		for(;it_f!=it_fe;it_f++){
//			(*var_dest)[(*it_f)->getID()] = part_index;
//		}
//
//		/* MESH EXTRACTION ON PART 0	 */
//		if (!platform::getEnvironmentVariable("INFO_CONSO").empty())
//			Log::mng()<<rank<<" - Consommation memoire en octets: "<<platform::getMemoryUsed()<<"\n";
//
//		Extractor extract(AMesh);
//		std::vector<Mesh > subMeshes;
//		Timer t2;
//		Log::mng()<<rank<<"/ Time to initialize the partitioning varible: "<<t2-t1<<"\n";
//		extract.extractWithGhostWithNativeSharedInfoTopLevelRepartitionExp2(GMDS_REGION,"part",subMeshes);
//		Timer t3;
//		Log::mng()<<rank<<"/ Time to extract submeshes: "<<t3-t2<<"\n";
//		if (!platform::getEnvironmentVariable("INFO_CONSO").empty())
//			Log::mng()<<rank<<" - Consommation memoire en octets: "<<platform::getMemoryUsed()<<"\n";
//
////		for(int i=0;i<subMeshes.size();i++){
////			LimaWriter writer(subMeshes[i]);
////			//	VTKWriter writer(subMeshes[i]);
////			std::ostringstream file_name;
////			file_name<<"homePOYOP/mesh/new_parallel_test2D_sub"<<i<<".mli";
////			writer.write(file_name.str(),N|F);
////		}
//
//		/* SUBMESHES SERIALIZATION */
//		for(int i=0;i<subMeshes.size();i++){
//			Mesh& subM = subMeshes[i];
//			std::ostringstream sub_stream;
//
//			subM.serialize(sub_stream);
//
//			std::string str = sub_stream.str();
//
//			int sub_size   = str.size()*sizeof(char);
//			const char* sub_char = str.c_str();
//
//			mpi_manager->send(i+1,0, (void*)sub_char,sub_size);
//
//		}
//		Timer t4;
//		Log::mng()<<rank<<"/ Time to serialize all submeshes: "<<t4-t3<<"\n";
//		if (!platform::getEnvironmentVariable("INFO_CONSO").empty())
//			Log::mng()<<rank<<" - Consommation memoire en octets: "<<platform::getMemoryUsed()<<"\n";
//
//
//	}
//	else{
//		if (!platform::getEnvironmentVariable("INFO_CONSO").empty())
//			Log::mng()<<rank<<" - Consommation memoire en octets: "<<platform::getMemoryUsed()<<"\n";
//
//		void *s;
//		int t=0;
//		s = mpi_manager->receive(0,0,t);
//		std::cout<<rank<<" - receive "<<t<<std::endl;
//
//		std::istrstream sub_stream2((char*)s,t);
//
//		AMesh.unserialize(sub_stream2);
//	}
//	if (!platform::getEnvironmentVariable("INFO_CONSO").empty())
//		Log::mng()<<rank<<" - Consommation memoire en octets: "<<platform::getMemoryUsed()<<"\n";
//
//}
///*----------------------------------------------------------------------------*/
//
//void MigrationTool::
//getSlaveInfo(Mesh& AMesh,
//		std::vector<id>& ALocalNodeSlaves,
//		std::vector<std::vector<std::pair<int, id> > >& ADistNodeSlaves,
//		std::vector<id>& ALocalEdgeSlaves,
//		std::vector<std::vector<std::pair<int, id> > >& ADistEdgeSlaves,
//		std::vector<id>& ALocalFaceSlaves,
//		std::vector<std::vector<std::pair<int, id> > >& ADistFaceSlaves,
//		std::vector<id>& ALocalRegionSlaves,
//		std::vector<std::vector<std::pair<int, id> > >& ADistRegionSlaves) const
//{
//	Timer time_start;
//	DistributedMemoryManager* mpi_manager= DistributedMemoryManager::instance();
//
//	TInt nb_partitions = mpi_manager->nbPartitions();
//	TInt rank 		   = mpi_manager->rank();
//	TInt partition_id = rank;
//
//#ifdef _DEBUG
//	Log::mng()<<"NB partitions= "<<nb_partitions<<"\n";
//#endif
//	Log::mng()<<"Get slavery info on "<<rank<<"\n";
//	if (!platform::getEnvironmentVariable("INFO_CONSO").empty())
//		Log::mng()<<rank<<" - Consommation memoire en octets: "<<platform::getMemoryUsed()<<"\n";
//
//	/* build the list of queries to the neighbourhood */
//	std::set<int> nodeNeighbourhood = AMesh.getNodeNeighbourhood();
//	for(std::set<int>::iterator it = nodeNeighbourhood.begin();it!=nodeNeighbourhood.end();it++)
//	{
//
//		int dest_part = *it;
////		if(dest_part<rank){ //SEND then RECV
////		mpi_manager->send(i+1,0, (void*)sub_char,sub_size);
////		}
////		else { //RECV then SEND
////
////		}
//
//	}
//
//}
///*----------------------------------------------------------------------------*/
//
//void MigrationTool::
//sequentialIDSplit2D(Mesh& AMesh) const
//{
//	Timer time_start;
//	/* the aim of this method is to split mesh_ into n mesh partitionned onto
//	 * n processors. The initial mesh is considered as being on processor 0.*/
//	DistributedMemoryManager* mpi_manager= DistributedMemoryManager::instance();
//
//	TInt nb_partitions = mpi_manager->nbPartitions();
//	TInt rank 		   = mpi_manager->rank();
//	TInt partition_id = rank;
//#ifdef _DEBUG
//	Log::mng()<<"NB partitions= "<<nb_partitions<<"\n";
//#endif
//	Log::mng()<<"Sequential cut on "<<rank<<"\n";
//	if (!platform::getEnvironmentVariable("INFO_CONSO").empty())
//		Log::mng()<<rank<<" - Consommation memoire en octets: "<<platform::getMemoryUsed()<<"\n";
//
//
//	/* interfaces are now correct, we can send the cells */
//	DistributedMemoryManager::MigrationStruct *migration_info, recv_info;
//	recv_info.init();
//	if(rank==0)
//	{
//		Timer time1;
//		std::vector<TInt> cellFilter;
//		int filter_size = AMesh.getNbFaces();
//		cellFilter.resize(filter_size+1);
//		//AMesh.createFilter<TInt>(2,cellFilter);
//		TInt maxNbCells = cellFilter.size()-1;
//#ifdef _DEBUG
//		Log::mng()<<"NB cells= "<<maxNbCells<<"\n";
//#endif
//		/* we compute the number of cells to keep on proc 0 and the number of cells
//		 * to send to each other partition */
//		TInt nbCellToSend = maxNbCells/nb_partitions;
//		TInt nbCellToKeep = nbCellToSend;
//		if (nbCellToSend*nb_partitions!=maxNbCells)
//			nbCellToKeep++;
//#ifdef _DEBUG
//		Log::mng()<<"to keep, to send = "<<nbCellToKeep<<", "<<nbCellToSend<<"\n";
//#endif
//		/* initialisation of the filter */
//		TInt filter_item = 1;
//		for(int i=0; i <nbCellToKeep;i++)
//			cellFilter[filter_item++]= 0;
//		for(int p = 1;p<nb_partitions;p++)
//			for(int i=0; i <nbCellToSend;i++)
//					cellFilter[filter_item++]= p;
//
//
//		/* WARNING there is no cell in the item 0 of a filter */
//		std::set<TInt> destination_partitions;
//		std::vector<Face*> faces_to_be_moved;
//
//		for(unsigned int i =0; i<nb_partitions;i++)
//			destination_partitions.insert(i);
//
//		/*we get the number of destination partitions */
//		TInt nb_migration_info = destination_partitions.size();
//		migration_info = new DistributedMemoryManager::MigrationStruct[nb_migration_info];
//		TInt *dest_tab = new TInt[nb_migration_info];
//		int dest_tab_index=0;
//		for(std::set<TInt>::iterator it= destination_partitions.begin();
//			it!= destination_partitions.end();it++)
//			dest_tab[dest_tab_index++] = *it;
//
//		/* for each target partition, we create a filter indicating
//		 * if this cell must be migrated to the target partition. */
//		std::vector<std::vector<TInt> > node_filters;
//		std::vector<std::vector<TInt> > edge_filters;
//		int move_mark = AMesh.getNewMark();
//
//		Timer time2;
//		Log::mng()<<"Init after: "<<time2-time1<<"\n";
//		if (!platform::getEnvironmentVariable("INFO_CONSO").empty())
//			Log::mng()<<rank<<" - Consommation memoire en octets: "<<platform::getMemoryUsed()<<"\n";
//
//		initFilters(AMesh,cellFilter,destination_partitions,node_filters,
//					edge_filters,faces_to_be_moved, move_mark);
//
//		Timer time3;
//		Log::mng()<<"After initFilters: "<<time3-time2<<"\n";
//		if (!platform::getEnvironmentVariable("INFO_CONSO").empty())
//			Log::mng()<<rank<<" - Consommation memoire en octets: "<<platform::getMemoryUsed()<<"\n";
//
//		/* creation of the internal interfaces */
//		/* we traverse the mesh nodes and  we share the involved node */
//		createInternalInterfaces(AMesh,node_filters, edge_filters);
//		Timer time4;
//		Log::mng()<<"After creation of internal interfaces: "<<time4-time3<<"\n";
//		if (!platform::getEnvironmentVariable("INFO_CONSO").empty())
//			Log::mng()<<rank<<" - Consommation memoire en octets: "<<platform::getMemoryUsed()<<"\n";
//
//
//
//		/* the initial interface is empty in this case */
//		std::map<TInt, std::set<Node*> > initial_node_interfaces;
//		std::map<TInt, std::set<Edge*> > initial_edge_interfaces;
//		//AMesh.getInterfaces(initial_node_interfaces);
//
//		prepareCellMigration(AMesh,cellFilter, faces_to_be_moved,
//							 node_filters, edge_filters,
//							 initial_node_interfaces, initial_edge_interfaces,
//							 migration_info, nb_migration_info, move_mark);
//		Timer time5;
//		Log::mng()<<"After preparing cell migration: "<<time5-time4<<"\n";
//		if (!platform::getEnvironmentVariable("INFO_CONSO").empty())
//			Log::mng()<<rank<<" - Consommation memoire en octets: "<<platform::getMemoryUsed()<<"\n";
//
//		/* Now we can send our structure. Be careful to send nothing to all the neighbours even if
//		 * the message is empty. Otherwise, you will be stuck !! */
//		performCellMigrationForIDSplit_Send(partition_id,dest_tab,migration_info,nb_migration_info);
//		Timer time6;
//		Log::mng()<<"After send: "<<time6-time5<<"\n";
//		if (!platform::getEnvironmentVariable("INFO_CONSO").empty())
//			Log::mng()<<rank<<" - Consommation memoire en octets: "<<platform::getMemoryUsed()<<"\n";
//
//		delete[] dest_tab;
//		AMesh.unmarkNodes(move_mark);
//		AMesh.unmarkEdges(move_mark);
//		AMesh.freeMark(move_mark);
//
//	}
//	else
//	{
//		Timer time5;
//		performCellMigrationForIDSplit_Recv(0,recv_info);
//		Timer time6;
//		Log::mng()<<"After recv: "<<time6-time5<<"\n";
//		if (!platform::getEnvironmentVariable("INFO_CONSO").empty())
//			Log::mng()<<rank<<" - Consommation memoire en octets: "<<platform::getMemoryUsed()<<"\n";
//
//	}
//
//	Timer time7;
//	/* we get the data to modify mesh interfaces */
//	updateAfterCellMigration(AMesh,recv_info);
//	Timer time8;
//	Log::mng()<<"After updating: "<<time8-time7<<"\n";
//	if (!platform::getEnvironmentVariable("INFO_CONSO").empty())
//		Log::mng()<<rank<<" - Consommation memoire en octets: "<<platform::getMemoryUsed()<<"\n";
//
////	if(migration_info!=0)
////		delete[] migration_info;
//
//
//}
///*----------------------------------------------------------------------------*/
//
//void MigrationTool::
//sequentialIDSplit3D(Mesh& AMesh) const
//{
//	Timer time_start;
//	/* the aim of this method is to split mesh_ into n mesh partitionned onto
//	 * n processors. The initial mesh is considered as being on processor 0.*/
//	DistributedMemoryManager* mpi_manager= DistributedMemoryManager::instance();
//
//	TInt nb_partitions = mpi_manager->nbPartitions();
//	TInt rank 		   = mpi_manager->rank();
//	TInt partition_id = rank;
//#ifdef _DEBUG
//	Log::mng()<<"NB partitions= "<<nb_partitions<<"\n";
//#endif
//	Log::mng()<<"Sequential cut on "<<rank<<"\n";
//	if (!platform::getEnvironmentVariable("INFO_CONSO").empty())
//		Log::mng()<<rank<<" - Consommation memoire en octets: "<<platform::getMemoryUsed()<<"\n";
//
//
//	/* interfaces are now correct, we can send the cells */
//	DistributedMemoryManager::MigrationStruct *migration_info, recv_info;
//	recv_info.init();
//	if(rank==0)
//	{
//		Timer time1;
//		std::vector<TInt> cellFilter;
//		int filter_size = AMesh.getNbRegions();
//		cellFilter.resize(filter_size+1);
//		//AMesh.createFilter<TInt>(2,cellFilter);
//		TInt maxNbCells = cellFilter.size()-1;
//#ifdef _DEBUG
//		Log::mng()<<"NB cells= "<<maxNbCells<<"\n";
//#endif
//		/* we compute the number of cells to keep on proc 0 and the number of cells
//		 * to send to each other partition */
//		TInt nbCellToSend = maxNbCells/nb_partitions;
//		TInt nbCellToKeep = nbCellToSend;
//		if (nbCellToSend*nb_partitions!=maxNbCells)
//			nbCellToKeep++;
//#ifdef _DEBUG
//		Log::mng()<<"to keep, to send = "<<nbCellToKeep<<", "<<nbCellToSend<<"\n";
//#endif
//		/* initialisation of the filter */
//		TInt filter_item = 1;
//		for(int i=0; i <nbCellToKeep;i++)
//			cellFilter[filter_item++]= 0;
//		for(int p = 1;p<nb_partitions;p++)
//			for(int i=0; i <nbCellToSend;i++)
//					cellFilter[filter_item++]= p;
//
//
//		/* WARNING there is no cell in the item 0 of a filter */
//		std::set<TInt> destination_partitions;
//		std::vector<Region*> regions_to_be_moved;
//
//		for(unsigned int i =0; i<nb_partitions;i++)
//			destination_partitions.insert(i);
//
//		/*we get the number of destination partitions */
//		TInt nb_migration_info = destination_partitions.size();
//		migration_info = new DistributedMemoryManager::MigrationStruct[nb_migration_info];
//		TInt *dest_tab = new TInt[nb_migration_info];
//		int dest_tab_index=0;
//		for(std::set<TInt>::iterator it= destination_partitions.begin();
//			it!= destination_partitions.end();it++)
//			dest_tab[dest_tab_index++] = *it;
//
//		/* for each target partition, we create a filter indicating
//		 * if this cell must be migrated to the target partition. */
//		std::vector<std::vector<TInt> > node_filters;
//		std::vector<std::vector<TInt> > face_filters;
//		int move_mark = AMesh.getNewMark();
//
//		Timer time2;
//		Log::mng()<<"Init after: "<<time2-time1<<"\n";
//		if (!platform::getEnvironmentVariable("INFO_CONSO").empty())
//			Log::mng()<<rank<<" - Consommation memoire en octets: "<<platform::getMemoryUsed()<<"\n";
//
//		initFilters(AMesh,cellFilter,destination_partitions,node_filters,
//					face_filters,regions_to_be_moved, move_mark);
//
//		Timer time3;
//		Log::mng()<<"After initFilters: "<<time3-time2<<"\n";
//		if (!platform::getEnvironmentVariable("INFO_CONSO").empty())
//			Log::mng()<<rank<<" - Consommation memoire en octets: "<<platform::getMemoryUsed()<<"\n";
//
//		/* creation of the internal interfaces */
//		/* we traverse the mesh nodes and  we share the involved node */
//		create3DInternalInterfaces(AMesh,node_filters, face_filters);
//		Timer time4;
//		Log::mng()<<"After creation of internal interfaces: "<<time4-time3<<"\n";
//		if (!platform::getEnvironmentVariable("INFO_CONSO").empty())
//			Log::mng()<<rank<<" - Consommation memoire en octets: "<<platform::getMemoryUsed()<<"\n";
//
//
//
//		/* the initial interface is empty in this case */
//		std::map<TInt, std::set<Node*> > initial_node_interfaces;
//		std::map<TInt, std::set<Face*> > initial_face_interfaces;
//		//AMesh.getInterfaces(initial_node_interfaces);
//
//		prepareCellMigration(AMesh,cellFilter, regions_to_be_moved,
//							 node_filters, face_filters,
//							 initial_node_interfaces, initial_face_interfaces,
//							 migration_info, nb_migration_info, move_mark);
//		Timer time5;
//		Log::mng()<<"After preparing cell migration: "<<time5-time4<<"\n";
//		if (!platform::getEnvironmentVariable("INFO_CONSO").empty())
//			Log::mng()<<rank<<" - Consommation memoire en octets: "<<platform::getMemoryUsed()<<"\n";
//
//		/* Now we can send our structure. Be careful to send nothing to all the neighbours even if
//		 * the message is empty. Otherwise, you will be stuck !! */
//		performCellMigrationForIDSplit_Send(partition_id,dest_tab,migration_info,nb_migration_info);
//		Timer time6;
//		Log::mng()<<"After send: "<<time6-time5<<"\n";
//		if (!platform::getEnvironmentVariable("INFO_CONSO").empty())
//			Log::mng()<<rank<<" - Consommation memoire en octets: "<<platform::getMemoryUsed()<<"\n";
//
//		delete[] dest_tab;
//		AMesh.unmarkNodes(move_mark);
//		AMesh.unmarkEdges(move_mark);
//		AMesh.freeMark(move_mark);
//
//	}
//	else
//	{
//		Timer time5;
//		performCellMigrationForIDSplit_Recv(0,recv_info);
//		Timer time6;
//		Log::mng()<<"After recv: "<<time6-time5<<"\n";
//		if (!platform::getEnvironmentVariable("INFO_CONSO").empty())
//			Log::mng()<<rank<<" - Consommation memoire en octets: "<<platform::getMemoryUsed()<<"\n";
//
//	}
//
//	Timer time7;
//	/* we get the data to modify mesh interfaces */
//	updateAfterCellMigration3D(AMesh,recv_info);
//	Timer time8;
//	Log::mng()<<"After updating: "<<time8-time7<<"\n";
//	if (!platform::getEnvironmentVariable("INFO_CONSO").empty())
//		Log::mng()<<rank<<" - Consommation memoire en octets: "<<platform::getMemoryUsed()<<"\n";
//
////	if(migration_info!=0)
////		delete[] migration_info;
//
//
//}
///*----------------------------------------------------------------------------*/
//
//void MigrationTool::
//getDestinations(const id& ACellID,
//				const std::map<TInt, std::vector<TInt> >& AFilters,
//				std::vector<TInt>& ADests) const
//{
//	ADests.clear();
//
//	std::map<TInt, std::vector<TInt> >::const_iterator it = AFilters.begin();
//	for(;it!=AFilters.end();it++)
//		if((it->second)[ACellID]==1)
//			ADests.push_back(it->first);
//}
///*----------------------------------------------------------------------------*/
//
//
//void MigrationTool::
//extractSubMesh(Mesh& AFromMesh, Mesh &AToMesh,
//			   std::vector<Face*>& AToExtract, const int& AToPartitionID) const
//{
//	MeshDescriptor meshDescriptor;
//
//	int migrationMark = AFromMesh.getNewMark();
//
//	std::set<id> nodesToCopy;
//	std::set<id> edgesToCopy;
//
//	if(meshDescriptor.hasNodes)
//	{
//		for(unsigned int i=0;i<AToExtract.size();i++)
//		{
//			/* we suppose that the relation to nodes is always available*/
//			AFromMesh.mark(AToExtract[i],migrationMark);
//			std::vector<id> local_nodes = AToExtract[i]->getNodeIDs();
//			for(unsigned int node_index = 0; node_index<local_nodes.size();
//				node_index++)
//			{
//				nodesToCopy.insert(local_nodes[node_index]);
//				AFromMesh.mark(AFromMesh.getLNode(local_nodes[node_index]),migrationMark);
//			}
//		}
//	}
//
//	if(meshDescriptor.hasEdges)
//	{
//		/*it is the most safe way to get the right result even if the connection are
//		 * not well initialized */
//		typename Mesh::edges_iterator edge_it 	 = AFromMesh.edges_begin();
//		typename Mesh::edges_iterator edge_it_end = AFromMesh.edges_end();
//		for(;edge_it!=edge_it_end;edge_it++)
//		{
//			std::vector<Node*> local_nodes= (*edge_it)->getNodes();
//			if (AFromMesh.isMarked(local_nodes[0],migrationMark) &&
//				AFromMesh.isMarked(local_nodes[1],migrationMark))
//				edgesToCopy.insert((*edge_it)->getID());
//		}
//
//	}
//
//	/* Now, we know the cells that must be migrated*/
//
//	/* NODES' MIGRATION */
//	std::set<id>::iterator node_it ;
//	std::map<id,id> map_nodes;
//
//	for(node_it= nodesToCopy.begin(); node_it!=nodesToCopy.end();node_it++)
//	{
//		Node* from_node = AFromMesh.getLNode(*node_it);
//		Node* to_node = AToMesh.newNode(from_node->getX(),from_node->getY(),
//										from_node->getZ());
//	//	to_node->setID(from_node->getID());
//		map_nodes[*node_it]=to_node->getID();
//
//		/*migration of the sharing behaviour*/
//		if(AFromMesh.isShared(from_node))
//		{
//			std::vector<TInt> sharing_partition;
//			AFromMesh.getPartitionsSharing(from_node,sharing_partition);
//			AToMesh.shareWith(to_node,sharing_partition);
//		}
//	}
//
//	/* EDGES' MIGRATION */
//	std::set<id>::iterator edge_it ;
//	std::map<id,id> map_edges;
//
//	for(edge_it= edgesToCopy.begin(); edge_it!=edgesToCopy.end();edge_it++)
//	{
//		Edge* from_edge = AFromMesh.getLEdge(*edge_it);
//		std::vector<id> local_nodes = from_edge->getNodeIDs();
//
//		Edge* to_edge   = AToMesh.newEdge(map_nodes[local_nodes[0]],map_nodes[local_nodes[1]]);
//
//		//to_edge->setID(from_edge->getID());
//		map_edges[*edge_it]=to_edge->getID();
//
//		/*migration of the sharing behaviour*/
//		if(AFromMesh.isShared(from_edge))
//		{
//			std::vector<TInt> sharing_partition;
//			AFromMesh.getPartitionsSharing(from_edge,sharing_partition);
//			AToMesh.shareWith(to_edge,sharing_partition);
//
//		}
//	}
//	/* FACES' MIGRATION */
//	for(unsigned int i=0;i<AToExtract.size();i++)
//	{
//		std::vector<id> fromNodeIds = AToExtract[i]->getNodeIDs();
//		std::vector<id> toNodeIds;
//		toNodeIds.resize(fromNodeIds.size());
//		for(unsigned int j=0;j<fromNodeIds.size();j++)
//			toNodeIds[j]=map_nodes[fromNodeIds[j]];
//		Face* to_face = AToMesh.newFace(toNodeIds);
//	//	to_face->setID(AToExtract[i]->getID());
//	}
//
//	/*--------- CELLS' EXTRACTION --------*/
//
//	/* we remove the extracted cells from AFromMesh*/
//	/* Some of the nodes to copy must remain in AFromMesh too (interface nodes).
//	 * The others are nodes of nodesToMove ant they are migrated without local copy */
//	std::set<id> nodesToMove;
//	std::vector<id> nodeAdjacency;
//	nodeAdjacency.resize(AFromMesh.getMaxLocalID(0)+1);
//	typename Mesh::faces_iterator face_it 	= AFromMesh.faces_begin();
//	typename Mesh::faces_iterator face_it_end = AFromMesh.faces_end();
//	for(;face_it!=face_it_end;face_it++)
//		if(!AFromMesh.isMarked(*face_it,migrationMark))
//		{
//			std::vector<id> node_ids = (*face_it)->getNodeIDs();
//			for(unsigned int id=0;id<node_ids.size();id++)
//				nodeAdjacency[node_ids[id]]++;
//		}
//
//
//	/* before removing the nodes, we mark them to remove some edges and
//	 * determinate interface edges */
//	int nodesToMoveMark    = AFromMesh.getNewMark();
//	int nodesOnInterfaceMark = AFromMesh.getNewMark();
//
//
//
//	/* we must remove some nodes from the mesh*/
//	for(unsigned int id=1;id<nodeAdjacency.size();id++)
//		if(nodeAdjacency[id]==0 && AFromMesh.isAValidNodeLocalID(id))
//		{
//			nodesToMove.insert(id);
//			AFromMesh.mark(AFromMesh.getLNode(id),nodesToMoveMark);
//		}
//
//	/* Now we define the nodes on the interface*/
//	std::set<id>::iterator it_nodesToCopy = nodesToCopy.begin();
//	std::vector<id> nodesOnInterface;
//	for(;it_nodesToCopy!=nodesToCopy.end();it_nodesToCopy++)
//	{
//		id lid = *it_nodesToCopy;
//		/* A REMPLACER PAR UN TEST SUR LA MARQUE BOOLEENE nodesToMoveMark*/
//		std::set<id>::iterator it_nodesToMove = nodesToMove.begin();
//		bool found = false;
//		for(;it_nodesToMove!=nodesToMove.end() && !found; it_nodesToMove++)
//			if(lid==*it_nodesToMove)
//				found = true;
//
//		if(!found)
//		{
//			AFromMesh.shareWith(AFromMesh.getLNode(lid),AToPartitionID);
//			AToMesh.shareWith(AToMesh.getLNode(map_nodes[lid]),
//							  AFromMesh.getPartID());
//			AFromMesh.mark(AFromMesh.getLNode(lid),nodesOnInterfaceMark);
//			nodesOnInterface.push_back(lid);
//		}
//	}
//
//	/* Now we define the edges on the interface and the edges to remove*/
//	/* Edges of edgesToMove are migrated without local copy */
//	if(meshDescriptor.hasEdges)
//	{
//		typename Mesh::edges_iterator edge_mesh_it 	 = AFromMesh.edges_begin();
//		typename Mesh::edges_iterator edge_mesh_it_end = AFromMesh.edges_end();
//		for(;edge_mesh_it!=edge_mesh_it_end;edge_mesh_it++)
//		{
//			std::vector<Node*> nodes = (*edge_mesh_it)->getNodes();
//			if(AFromMesh.isMarked(nodes[0],nodesToMoveMark) || AFromMesh.isMarked(nodes[1],nodesToMoveMark))
//				AFromMesh.deleteEdge(*edge_mesh_it);
//			else if (AFromMesh.isMarked(nodes[0],nodesOnInterfaceMark) &&
//					 AFromMesh.isMarked(nodes[1],nodesOnInterfaceMark) )
//			{
//				AFromMesh.shareWith(*edge_mesh_it,AToPartitionID);
//				AToMesh.shareWith(AToMesh.getLEdge(map_edges[(*edge_mesh_it)->getID()]),
//								  AFromMesh.getPartID());
//			}
//		}
//	}
//
//	/* now we definitively remove the nodes */
//	std::set<id>::iterator it_nodesToMove = nodesToMove.begin();
//	for(;it_nodesToMove!=nodesToMove.end();it_nodesToMove++)
//		AFromMesh.deleteNode(*it_nodesToMove);
//
//	/*FACES' REMOVAL*/
//	for(unsigned int i=0;i<AToExtract.size();i++)
//		AFromMesh.deleteFace(AToExtract[i]);
//
//	/* We unmark the nodes on the boundary before making the corresponding mark
//	 * free
//	 */
//	for(unsigned int i=0;i<nodesOnInterface.size();i++)
//		AFromMesh.unmark(AFromMesh.getLNode(nodesOnInterface[i]),nodesOnInterfaceMark);
//
//
//	AFromMesh.freeMark(migrationMark);
//	AFromMesh.freeMark(nodesToMoveMark);
//	AFromMesh.freeMark(nodesOnInterfaceMark);
//}
///*----------------------------------------------------------------------------*/
//
//void MigrationTool::
//extractSubMesh(Mesh& AFromMesh, Mesh &AToMesh,
//			   std::vector<Region*>& AToExtract, const int& AToPartitionID) const
//{
//	MeshDescriptor meshDescriptor;
//
//	int migrationMark = AFromMesh.getNewMark();
//
//	std::set<id> nodesToCopy;
//	std::set<id> edgesToCopy;
//	std::set<id> facesToCopy;
//
//
//	if(meshDescriptor.hasNodes)
//	{
//		for(unsigned int i=0;i<AToExtract.size();i++)
//		{
//			/* we suppose that the relation to nodes is always available*/
//			AFromMesh.mark(AToExtract[i],migrationMark);
//			std::vector<id> local_nodes = AToExtract[i]->getNodeIDs();
//			for(unsigned int node_index = 0; node_index<local_nodes.size();
//				node_index++)
//			{
//				nodesToCopy.insert(local_nodes[node_index]);
//				AFromMesh.mark(AFromMesh.getLNode(local_nodes[node_index]),migrationMark);
//			}
//		}
//	}
//	if(meshDescriptor.hasEdges)
//	{
//		/*it is the most safe way to get the right result even if the connection are
//		 * not well initialized */
//		typename Mesh::edges_iterator edge_it 	 = AFromMesh.edges_begin();
//		typename Mesh::edges_iterator edge_it_end = AFromMesh.edges_end();
//		for(;edge_it!=edge_it_end;edge_it++)
//		{
//			std::vector<Node*> local_nodes= (*edge_it)->getNodes();
//			if (AFromMesh.isMarked(local_nodes[0],migrationMark) &&
//				AFromMesh.isMarked(local_nodes[1],migrationMark))
//				edgesToCopy.insert((*edge_it)->getID());
//		}
//
//	}
//	if(meshDescriptor.hasFaces)
//	{
//		/*it is the most safe way to get the right result even if the connection are
//		 * not well initialized */
//		typename Mesh::faces_iterator face_it 	 = AFromMesh.faces_begin();
//		typename Mesh::faces_iterator face_it_end = AFromMesh.faces_end();
//		for(;face_it!=face_it_end;face_it++)
//		{
//			std::vector<Node*> local_nodes= (*face_it)->getNodes();
//			bool marked = true;
//			for(unsigned int i=0; i<local_nodes.size() && marked;i++)
//				if (!AFromMesh.isMarked(local_nodes[i],migrationMark))
//					marked=false;
//
//			if(marked)
//				facesToCopy.insert((*face_it)->getID());
//		}
//
//	}
//	/* Now, we know the cells that must be migrated*/
//
//	/* NODES' MIGRATION */
//	std::set<id>::iterator node_it ;
//	std::map<id,id> map_nodes;
//
//	for(node_it= nodesToCopy.begin(); node_it!=nodesToCopy.end();node_it++)
//	{
//		Node* from_node = AFromMesh.getLNode(*node_it);
//		Node* to_node = AToMesh.newNode(from_node->getX(),from_node->getY(),
//										from_node->getZ());
//		//to_node->setID(from_node->getID());
//		map_nodes[*node_it]=to_node->getID();
//
//		/*migration of the sharing behaviour*/
//		if(AFromMesh.isShared(from_node))
//		{
//			std::vector<TInt> sharing_partition;
//			AToMesh.shareWith(to_node,AFromMesh.getPartitionsSharing(from_node,sharing_partition));
//		}
//	}
//
//	/* EDGES' MIGRATION */
//	std::set<id>::iterator edge_it ;
//	std::map<id,id> map_edges;
//
//	for(edge_it= edgesToCopy.begin(); edge_it!=edgesToCopy.end();edge_it++)
//	{
//		Edge* from_edge = AFromMesh.getLEdge(*edge_it);
//		std::vector<id> local_nodes = from_edge->getNodeIDs();
//
//		Edge* to_edge   = AToMesh.newEdge(map_nodes[local_nodes[0]],map_nodes[local_nodes[1]]);
//
//	//	to_edge->setID(from_edge->getID());
//		map_edges[*edge_it]=to_edge->getID();
//
//		/*migration of the sharing behaviour*/
//		if(AFromMesh.isShared(from_edge))
//		{
//			std::vector<TInt> sharing_partition;
//			AFromMesh.getPartitionsSharing(from_edge,sharing_partition);
//			AToMesh.shareWith(to_edge,sharing_partition);
//		}
//	}
//
//	/* FACES' MIGRATION */
//	std::set<id>::iterator face_it ;
//	std::map<id,id> map_faces;
//	for(face_it= facesToCopy.begin(); face_it!=facesToCopy.end();face_it++)
//	{
//		Face* from_face = AFromMesh.getLFace(*face_it);
//		std::vector<id> fromNodeIds = from_face->getNodeIDs();
//		std::vector<id> toNodeIds;
//		toNodeIds.resize(fromNodeIds.size());
//		for(unsigned int j=0;j<fromNodeIds.size();j++)
//			toNodeIds[j]=map_nodes[fromNodeIds[j]];
//
//		Face* to_face = AToMesh.newFace(toNodeIds);
//
//	//	to_face->setID(from_face->getID());
//		map_faces[*face_it]=to_face->getID();
//
//		/*migration of the sharing behaviour*/
//		if(AFromMesh.isShared(from_face))
//		{
//			std::vector<TInt> sharing_partition;
//			AFromMesh.getPartitionsSharing(from_face,sharing_partition);
//			AToMesh.shareWith(to_face,sharing_partition);
//		}
//	}
//
//	/* REGION' MIGRATION */
//	for(unsigned int i=0;i<AToExtract.size();i++)
//	{
//		std::vector<id> fromNodeIds = AToExtract[i]->getNodeIDs();
//		std::vector<id> toNodeIds;
//		toNodeIds.resize(fromNodeIds.size());
//		for(unsigned int j=0;j<fromNodeIds.size();j++)
//			toNodeIds[j]=map_nodes[fromNodeIds[j]];
//
//		Region* region_new;
//		switch(AToExtract[i]->getType())
//		{
//		case GMDS_HEX:
//			region_new = AToMesh.newHex(toNodeIds[0],toNodeIds[1],toNodeIds[2],toNodeIds[3],
//									  toNodeIds[4],toNodeIds[5],toNodeIds[6],toNodeIds[7]);
//			break;
//		case GMDS_TETRA:
//			region_new = AToMesh.newTet(toNodeIds[0],toNodeIds[1],toNodeIds[2],toNodeIds[3]);
//			break;
//		case GMDS_PRISM3:
//			region_new = AToMesh.newPrism3(toNodeIds[0],toNodeIds[1],toNodeIds[2],toNodeIds[3],
//										 toNodeIds[4],toNodeIds[5]);
//			break;
//		case GMDS_PRISM5:
//			region_new = AToMesh.newPrism5(toNodeIds[0],toNodeIds[1],toNodeIds[2],toNodeIds[3],
//										 toNodeIds[4],toNodeIds[5],toNodeIds[6],toNodeIds[7],
//										 toNodeIds[8],toNodeIds[9]);
//			break;
//		case GMDS_PRISM6:
//			region_new = AToMesh.newPrism6(toNodeIds[0],toNodeIds[1],toNodeIds[2],toNodeIds[3],
//										 toNodeIds[4],toNodeIds[5],toNodeIds[6],toNodeIds[7],
//										 toNodeIds[8],toNodeIds[9],toNodeIds[10],toNodeIds[11]);
//			break;
//		case GMDS_PYRAMID:
//			region_new = AToMesh.newPyramid(toNodeIds[0],toNodeIds[1],toNodeIds[2],toNodeIds[3],
//										  toNodeIds[4]);
//			break;
//		case GMDS_POLYHEDRA:
//			std::cout<<"Copy of polyhedra is not implemented yet"<<std::endl;
//			continue;
//			break;
//
//		}
//		//region_new->setID(AToExtract[i]->getID());
//
//	}
//
//	/*--------- CELLS' EXTRACTION --------*/
//
//	/* we remove the extracted cells from AFromMesh*/
//	/* Some of the nodes to copy must remain in AFromMesh too (interface nodes).
//	 * The others are nodes of nodesToMove ant they are migrated without local copy */
//	std::set<id> nodesToMove;
//	std::vector<id> nodeAdjacency;
//	std::vector<id> nodesOnInterface;
//	nodeAdjacency.resize(AFromMesh.getMaxLocalID(0)+1);
//
//	typename Mesh::regions_iterator region_it 	 = AFromMesh.regions_begin();
//	typename Mesh::regions_iterator region_it_end = AFromMesh.regions_end();
//	for(;region_it!=region_it_end;region_it++)
//		if(!AFromMesh.isMarked(*region_it,migrationMark))
//		{
//			std::vector<id> node_ids = (*region_it)->getNodeIDs();
//			for(unsigned int id=0;id<node_ids.size();id++)
//				nodeAdjacency[node_ids[id]]++;
//		}
//
//
//	/* before removing the nodes, we mark them to remove some edges and faces and
//	 * determinate interface edges and faces */
//	int nodesToMoveMark    = AFromMesh.getNewMark();
//	int nodesOnInterfaceMark = AFromMesh.getNewMark();
//
//	/* we remove some nodes from the mesh*/
//	for(unsigned int id=1;id<nodeAdjacency.size();id++)
//		if(nodeAdjacency[id]==0 && AFromMesh.isAValidNodeLocalID(id))
//		{
//			nodesToMove.insert(id);
//			AFromMesh.mark(AFromMesh.getLNode(id),nodesToMoveMark);
//		}
//
//	/* Now we define the nodes on the interface*/
//	std::set<id>::iterator it_nodesToCopy = nodesToCopy.begin();
//	for(;it_nodesToCopy!=nodesToCopy.end();it_nodesToCopy++)
//	{
//		id lid = *it_nodesToCopy;
//		std::set<id>::iterator it_nodesToMove = nodesToMove.begin();
//		bool found = false;
//		for(;it_nodesToMove!=nodesToMove.end() && !found; it_nodesToMove++)
//			if(lid==*it_nodesToMove)
//				found = true;
//
//		if(!found)
//		{
//			AFromMesh.shareWith(AFromMesh.getLNode(lid),AToPartitionID);
//			AToMesh.shareWith(AToMesh.getLNode(map_nodes[lid]),
//							  AFromMesh.getPartID());
//			AFromMesh.mark(AFromMesh.getLNode(lid),nodesOnInterfaceMark);
//			nodesOnInterface.push_back(lid);
//		}
//	}
//
//	/* Now we define the edges on the interface and the edges to remove*/
//	if(meshDescriptor.hasEdges)
//	{
//		typename Mesh::edges_iterator edge_mesh_it 	 = AFromMesh.edges_begin();
//		typename Mesh::edges_iterator edge_mesh_it_end = AFromMesh.edges_end();
//		for(;edge_mesh_it!=edge_mesh_it_end;edge_mesh_it++)
//		{
//			std::vector<Node*> nodes = (*edge_mesh_it)->getNodes();
//			if(AFromMesh.isMarked(nodes[0],nodesToMoveMark) || AFromMesh.isMarked(nodes[1],nodesToMoveMark))
//				AFromMesh.deleteEdge(*edge_mesh_it);
//			else if (AFromMesh.isMarked(nodes[0],nodesOnInterfaceMark) &&
//					 AFromMesh.isMarked(nodes[1],nodesOnInterfaceMark) )
//			{
//				AFromMesh.shareWith(*edge_mesh_it,AToPartitionID);
//				AToMesh.shareWith(AToMesh.getLEdge(map_edges[(*edge_mesh_it)->getID()]),
//								  AFromMesh.getPartID());
//			}
//		}
//	}
//
//	/* Now we define the faces on the interface and the faces to remove*/
//	if(meshDescriptor.hasFaces)
//	{
//		typename Mesh::faces_iterator face_mesh_it 	  = AFromMesh.faces_begin();
//		typename Mesh::faces_iterator face_mesh_it_end = AFromMesh.faces_end();
//		for(;face_mesh_it!=face_mesh_it_end;face_mesh_it++)
//		{
//			std::vector<Node*> nodes = (*face_mesh_it)->getNodes();
//			bool faceOnInterface = true;
//			bool faceToMove = false;
//			for(unsigned int i=0;i<nodes.size() && !faceToMove;i++)
//			{
//				if(AFromMesh.isMarked(nodes[i],nodesToMoveMark) )
//					faceToMove=true;
//				if (!AFromMesh.isMarked(nodes[i],nodesOnInterfaceMark))
//					faceOnInterface=false;
//			}
//			if (faceToMove)
//				AFromMesh.deleteFace(*face_mesh_it);
//			else if (faceOnInterface)
//			{
//				AFromMesh.shareWith(*face_mesh_it,AToPartitionID);
//				AToMesh.shareWith(AToMesh.getLFace(map_faces[(*face_mesh_it)->getID()]),
//								  AFromMesh.getPartID());
//			}
//		}
//	}
//
//	/* now we definitively remove the nodes */
//	std::set<id>::iterator it_nodesToMove = nodesToMove.begin();
//	for(;it_nodesToMove!=nodesToMove.end();it_nodesToMove++)
//		AFromMesh.deleteNode(*it_nodesToMove);
//
//
//	for(unsigned int i=0;i<AToExtract.size();i++)
//		AFromMesh.deleteRegion(AToExtract[i]);
//
//	/* We unmark the nodes on the boundary before making the corresponding mark
//	 * free
//	 */
//	for(unsigned int i=0;i<nodesOnInterface.size();i++)
//		AFromMesh.unmark(AFromMesh.getLNode(nodesOnInterface[i]),nodesOnInterfaceMark);
//
//	AFromMesh.freeMark(migrationMark);
//	AFromMesh.freeMark(nodesToMoveMark);
//	AFromMesh.freeMark(nodesOnInterfaceMark);
//}
///*----------------------------------------------------------------------------*/
//
//void MigrationTool::
//insertSubMesh(Mesh& AMesh, Mesh &ASubMesh,
//			  const int& AFromPartitionID) const
//{
//	MeshDescriptor meshDescriptor;
//	std::set<TInt> partitionIds;
//	AMesh.getPartitions(partitionIds);
//
//	/** all the cells which are copied in AMesh are marked in ASubmesh. As a
//	 * consequence, all the cells in ASubmesh are marked by the end of this
//	 * algorithm.
//	 */
//	TInt mark = ASubMesh.getNewMark();
//
//	std::map<id,id> nodes_projection;
//	std::map<id,id> edges_projection;
//	std::map<id,id> faces_projection;
//	std::map<id,id> regions_projection;
//	/* we look for the cells which are shared between AMesh and ASubMesh*/
//	if(partitionIds.find(AFromPartitionID)!=partitionIds.end())
//		initProjectionMap(AMesh, ASubMesh,
//						  nodes_projection,edges_projection,
//						  faces_projection,regions_projection,
//						  AFromPartitionID, mark);
//
//
//	/*================== CELL COPY ==================*/
//
//	/* copy of the necessary nodes from ASubMesh to AMesh*/
//	typename Mesh::nodes_iterator it_n_submesh= ASubMesh.nodes_begin();
//	for(;it_n_submesh!=ASubMesh.nodes_end();it_n_submesh++)
//		if(!ASubMesh.isMarked(*it_n_submesh,mark))
//		{
//			Node* node_sub = *it_n_submesh;
//			ASubMesh.mark(node_sub,mark);
//
//			Node* node_new = AMesh.newNode(node_sub->getX(),node_sub->getY(), node_sub->getZ());
//		//	node_new->setID(node_sub->getID());
//
//			nodes_projection[(*it_n_submesh)->getID()]= node_new->getID();
//
//			/* partitionning information*/
//			std::vector<TInt> sharing_partition;
//			ASubMesh.getPartitionsSharing(node_sub,sharing_partition);
//			AMesh.shareWith(node_new,sharing_partition);
//		}
//
//	/* copy of the necessary edges from ASubMesh to AMesh*/
//	typename Mesh::edges_iterator it_e_submesh= ASubMesh.edges_begin();
//	for(;it_e_submesh!=ASubMesh.edges_end();it_e_submesh++)
//		if(!ASubMesh.isMarked(*it_e_submesh,mark))
//		{
//			Edge* edge_sub = *it_e_submesh;
//			ASubMesh.mark(edge_sub,mark);
//			std::vector<Node*> nodes_of_edge_sub = edge_sub->getNodes();
//			Edge* edge_new = AMesh.newEdge(nodes_projection[nodes_of_edge_sub[0]->getID()],
//										   nodes_projection[nodes_of_edge_sub[1]->getID()]);
//		//	edge_new->setID(edge_sub->getID());
//
//			edges_projection[(*it_e_submesh)->getID()]= edge_new->getID();
//			/* partitionning information*/
//			std::vector<TInt> sharing_partition;
//			ASubMesh.getPartitionsSharing(edge_sub,sharing_partition);
//			AMesh.shareWith(edge_new,sharing_partition);
//		}
//
//	/* copy of the necessary faces from ASubMesh to AMesh*/
//	typename Mesh::faces_iterator it_f_submesh= ASubMesh.faces_begin();
//	for(;it_f_submesh!=ASubMesh.faces_end();it_f_submesh++)
//		if(!ASubMesh.isMarked(*it_f_submesh,mark))
//		{
//			Face* face_sub = *it_f_submesh;
//			ASubMesh.mark(face_sub,mark);
//
//			std::vector<id> meshNodeIds, subMeshNodeIds = face_sub->getNodeIDs();
//			meshNodeIds.reserve(subMeshNodeIds.size());
//
//			for(unsigned int i =0; i< subMeshNodeIds.size();i++)
//				meshNodeIds.push_back(nodes_projection[subMeshNodeIds[i] ]);
//
//			Face* face_new = AMesh.newFace(meshNodeIds);
//		//	face_new->setID(face_sub->getID());
//
//			faces_projection[(*it_f_submesh)->getID()]= face_new->getID();
//			/* partitionning information*/
//			std::vector<TInt> sharing_partition;
//			ASubMesh.getPartitionsSharing(face_sub,sharing_partition);
//			AMesh.shareWith(face_new,sharing_partition);
//
//		}
//
//	/* copy of the necessary regions from ASubMesh to AMesh*/
//	typename Mesh::regions_iterator it_r_submesh= ASubMesh.regions_begin();
//	for(;it_r_submesh!=ASubMesh.regions_end();it_r_submesh++)
//		if(!ASubMesh.isMarked(*it_r_submesh,mark))
//		{
//			Region* region_sub = *it_r_submesh;
//			ASubMesh.mark(region_sub,mark);
//			std::vector<Node*> nodes_of_region_sub = region_sub->getNodes();
//
//			Region* region_new;
//			switch(region_sub->getType())
//			{
//			case GMDS_HEX:
//				region_new = AMesh.newHex(nodes_projection[nodes_of_region_sub[0]->getID()],
//										  nodes_projection[nodes_of_region_sub[1]->getID()],
//										  nodes_projection[nodes_of_region_sub[2]->getID()],
//										  nodes_projection[nodes_of_region_sub[3]->getID()],
//										  nodes_projection[nodes_of_region_sub[4]->getID()],
//										  nodes_projection[nodes_of_region_sub[5]->getID()],
//										  nodes_projection[nodes_of_region_sub[6]->getID()],
//										  nodes_projection[nodes_of_region_sub[7]->getID()]);
//				break;
//			case GMDS_TETRA:
//				region_new = AMesh.newTet(nodes_projection[nodes_of_region_sub[0]->getID()],
//										  nodes_projection[nodes_of_region_sub[1]->getID()],
//										  nodes_projection[nodes_of_region_sub[2]->getID()],
//										  nodes_projection[nodes_of_region_sub[3]->getID()]);
//				break;
//			case GMDS_PRISM3:
//				region_new = AMesh.newPrism3(nodes_projection[nodes_of_region_sub[0]->getID()],
//											 nodes_projection[nodes_of_region_sub[1]->getID()],
//											 nodes_projection[nodes_of_region_sub[2]->getID()],
//											 nodes_projection[nodes_of_region_sub[3]->getID()],
//											 nodes_projection[nodes_of_region_sub[4]->getID()],
//											 nodes_projection[nodes_of_region_sub[5]->getID()]);
//				break;
//			case GMDS_PRISM5:
//				region_new = AMesh.newPrism5(nodes_projection[nodes_of_region_sub[0]->getID()],
//											 nodes_projection[nodes_of_region_sub[1]->getID()],
//											 nodes_projection[nodes_of_region_sub[2]->getID()],
//											 nodes_projection[nodes_of_region_sub[3]->getID()],
//											 nodes_projection[nodes_of_region_sub[4]->getID()],
//											 nodes_projection[nodes_of_region_sub[5]->getID()],
//											 nodes_projection[nodes_of_region_sub[6]->getID()],
//											 nodes_projection[nodes_of_region_sub[7]->getID()],
//											 nodes_projection[nodes_of_region_sub[8]->getID()],
//											 nodes_projection[nodes_of_region_sub[9]->getID()]);
//				break;
//			case GMDS_PRISM6:
//				region_new = AMesh.newPrism6(nodes_projection[nodes_of_region_sub[0]->getID()],
//											 nodes_projection[nodes_of_region_sub[1]->getID()],
//											 nodes_projection[nodes_of_region_sub[2]->getID()],
//											 nodes_projection[nodes_of_region_sub[3]->getID()],
//											 nodes_projection[nodes_of_region_sub[4]->getID()],
//											 nodes_projection[nodes_of_region_sub[5]->getID()],
//											 nodes_projection[nodes_of_region_sub[6]->getID()],
//											 nodes_projection[nodes_of_region_sub[7]->getID()],
//											 nodes_projection[nodes_of_region_sub[8]->getID()],
//											 nodes_projection[nodes_of_region_sub[9]->getID()],
//											 nodes_projection[nodes_of_region_sub[10]->getID()],
//											 nodes_projection[nodes_of_region_sub[11]->getID()]);
//				break;
//			case GMDS_PYRAMID:
//				region_new = AMesh.newPyramid(nodes_projection[nodes_of_region_sub[0]->getID()],
//											  nodes_projection[nodes_of_region_sub[1]->getID()],
//											  nodes_projection[nodes_of_region_sub[2]->getID()],
//											  nodes_projection[nodes_of_region_sub[3]->getID()],
//											  nodes_projection[nodes_of_region_sub[4]->getID()]);
//				break;
//			case GMDS_POLYHEDRA:
//				std::cout<<"Copy of polyhedra is not implemented yet"<<std::endl;
//				continue;
//				break;
//
//			}
//		//	region_new->setID(region_sub->getID());
//
//			regions_projection[(*it_r_submesh)->getID()]= region_new->getID();
//			/* partitionning information*/
//			std::vector<TInt> sharing_partition;
//			ASubMesh.getPartitionsSharing(region_sub,sharing_partition);
//			AMesh.shareWith(region_new,sharing_partition);
//
//		}
//
//	/*we clear the submesh now*/
//	ASubMesh.clear();
//}
///*----------------------------------------------------------------------------*/
//
//void MigrationTool::
//migrate( Mesh& AMesh, std::map<int, std::vector<Face*> > AToMigrate,
//		std::vector<Mesh >& ASubMeshes) const
//{
//	int current_partition_ID = AMesh.getPartID();
//	std::map<int, std::vector<Face*> >::iterator it_dest = AToMigrate.begin();
//
////	std::vector<Mesh > submeshes;
////	std::vector<Tint>  dest_partition;
////
////	submeshes.resize(AToMigrate.size());
////	dest_partition.resize(AToMigrate.size());
////
////	/* STEP 1
////	 * All the cells to migrate are extracted from the original mesh
////	 * and gathered in destination meshes
////	 */
////	for(;it_dest != AToMigrate.end();it_dest++)
////	{
////		extractSubMesh(AMesh, /*sub*/ASubmeshes[i], it_dest->second,it_dest->first);
////		dest_partition[i] = it_dest->first;
////	}
//
////	/* STEP 2
////	 * For all submeshes to export, we get the ids of the cells on the interface
////	 * and where they are going to be migrated
////	 */
////	std::map<id, std::vector<TInt> > node_to_partitions;
////	std::map<id, std::vector<TInt> > edge_to_partitions;
////	std::map<id, std::vector<TInt> > face_to_partitions;
////	std::map<id, std::vector<TInt> > region_to_partitions;
////
////	std::map<TInt, std::vector<id> > partition_to_nodes;
////	std::map<TInt, std::vector<id> > partition_to_edges;
////	std::map<TInt, std::vector<id> > partition_to_faces;
////	std::map<TInt, std::vector<id> > partition_to_regions;
////
////	for(unsigned int i=0; i<submeshes.size();i++)
////	{
////		Mesh current_submesh = submeshes[i];
////
////		std::map<TInt, std::set<Node*  > >& interface_nodes;
////		std::map<TInt, std::set<Edge*  > >& interface_edges;
////		std::map<TInt, std::set<Face*  > >& interface_faces;
////		std::map<TInt, std::set<Region*> >& interface_regions;
////
////		/* we get the interfaces of the current submesh*/
////		current_submesh.getInterfaces(interface_nodes);
////		current_submesh.getInterfaces(interface_edges);
////		current_submesh.getInterfaces(interface_faces);
////		current_submesh.getInterfaces(interface_regions);
////
////		/* we traverse the nodes participating to this interface */
////		std::map<TInt, std::set<Node*  > >::iterator it_nodes = interface_nodes.begin();
////		for(;it_nodes!=interface_nodes.end();it_nodes++)
////			if(it_nodes->first!=current_partition_ID){
////				std::set<Node*> nodes = it_nodes->second;
////				std::set<Node*>::iterator it = nodes.begin();
////				for(;it!=nodes.end();it++)
////				{
////					std::vector<TInt> sharing_partitions;
////					current_submesh.getPartitionsSharing(*it,sharing_partitions);
////					node_to_partitions[it->getID()]= sharing_partitions;
////					for(unsigned int i = 0; i<sharing_partitions.size(); i++)
////						partition_to_node[sharing_partitions[i]].push_back(*it->getID());
////				}
////			}
////
////		/* todo EDGES, FACES, REGIONS */
////	}
////
////	/* STEP 3
////	 * We know all the cells to move (submeshes) and the information to get to neighbourhood.
////	 * Now we tranfer this information.
////	 */
////	std::set<TInt> neighbours;
////	AMesh.getPartitions(neighbours);
////	std::set<TInt>::iterator it_neighbours= neighbours.begin();
////	for(;it_neighbours!=neighbours.end(); it_neighbours++){
////		updateBoundary(it_neighbours, neighbours,
////					   sharing_nodes, sharing_edges,
////					   sharing_faces, sharing_regions);
////	}
////
//
//
////	ASubMeshes.clear();
////	ASubMeshes.resize(submeshes.size());
////	for(unsigned int i=0; i<submeshes.size();i++)
////		ASubMeshes[i]=submeshes[i];
//}
///*----------------------------------------------------------------------------*/
//
//void MigrationTool::
//initProjectionMap(Mesh& AMesh, Mesh &ASubMesh,
//			  std::map<id,id>& ANodeProj,
//			  std::map<id,id>& AEdgeProj,
//			  std::map<id,id>& AFaceProj,
//			  std::map<id,id>& ARegionProj,
//			  const int& AFromPartitionID,
//			  const int& AInterfaceMark) const
//
//{
//	MeshDescriptor meshDescriptor;
//	std::set<Node*>   nodes_mesh  , nodes_submesh;
//	std::set<Edge*>   edges_mesh  , edges_submesh;
//	std::set<Face*>   faces_mesh  , faces_submesh;
//	std::set<Region*> regions_mesh, regions_submesh;
//
//
//	AMesh.getInterface(nodes_mesh,AFromPartitionID);
//
//	ASubMesh.getInterface(nodes_submesh,AMesh.getPartID());
//	/* we begin to build the projection fonction from local id in ASubmesh to
//	 * local id in AMesh*/
//	std::set<Node*>::const_iterator it_n_submesh= nodes_submesh.begin();
//	for(;it_n_submesh!=nodes_submesh.end();it_n_submesh++)
//	{
//		std::set<Node*>::const_iterator it_n_mesh= nodes_mesh.begin();
//		bool found = false;
//		for(;it_n_mesh!=nodes_mesh.end() && !found;it_n_mesh++)
//		{
//			if((*it_n_mesh)->getID()==(*it_n_submesh)->getID())
//			{
//				found=true;
//				ANodeProj[(*it_n_submesh)->getID()]= (*it_n_mesh)->getID();
//				ASubMesh.mark(*it_n_submesh,AInterfaceMark);
//				if(!ASubMesh.isSharedWith(*it_n_submesh,AFromPartitionID))
//				{
//					AMesh.removeFromInterface(*it_n_mesh,AFromPartitionID);
//				}
//			}
//		}
//	}
//
//	if(meshDescriptor.hasEdges)
//	{
//		AMesh.getInterface(edges_mesh,AFromPartitionID);
//		ASubMesh.getInterface(edges_submesh,AMesh.getPartID());
//
//		/* we begin to build the projection fonction from local id in ASubmesh to
//		 * local id in AMesh*/
//		std::set<Edge*>::const_iterator it_e_submesh= edges_submesh.begin();
//		for(;it_e_submesh!=edges_submesh.end();it_e_submesh++)
//		{
//			std::set<Edge*>::const_iterator it_e_mesh= edges_mesh.begin();
//			bool found = false;
//			for(;it_e_mesh!=edges_mesh.end() && !found;it_e_mesh++)
//			{
//				if((*it_e_mesh)->getID()==(*it_e_submesh)->getID())
//				{
//					found=true;
//					AEdgeProj[(*it_e_submesh)->getID()]= (*it_e_mesh)->getID();
//					ASubMesh.mark(*it_e_submesh,AInterfaceMark);
//					if(!ASubMesh.isSharedWith(*it_e_submesh,AFromPartitionID))
//					{
//						AMesh.removeFromInterface(*it_e_mesh,AFromPartitionID);
//					}
//
//				}
//			}
//		}
//	}
//
//	if(meshDescriptor.hasFaces)
//	{
//		AMesh.getInterface(faces_mesh,AFromPartitionID);
//		ASubMesh.getInterface(faces_submesh,AMesh.getPartID());
//
//		/* we begin to build the projection fonction from local id in ASubmesh to
//		 * local id in AMesh*/
//		std::set<Face*>::const_iterator it_f_submesh= faces_submesh.begin();
//		for(;it_f_submesh!=faces_submesh.end();it_f_submesh++)
//		{
//			std::set<Face*>::const_iterator it_f_mesh= faces_mesh.begin();
//			bool found = false;
//			for(;it_f_mesh!=faces_mesh.end() && !found;it_f_mesh++)
//			{
//				if((*it_f_mesh)->getID()==(*it_f_submesh)->getID())
//				{
//					found=true;
//					AFaceProj[(*it_f_submesh)->getID()]= (*it_f_mesh)->getID();
//					ASubMesh.mark(*it_f_submesh,AInterfaceMark);
//					if(!ASubMesh.isSharedWith(*it_f_submesh,AFromPartitionID))
//					{
//						AMesh.removeFromInterface(*it_f_mesh,AFromPartitionID);
//					}
//
//				}
//			}
//		}
//	}
//
//	if(meshDescriptor.hasRegions && meshDescriptor.RKnowR)
//	{
//		AMesh.getInterface(regions_mesh,AFromPartitionID);
//		ASubMesh.getInterface(regions_submesh,AMesh.getPartID());
//
//		/* we begin to build the projection fonction from local id in ASubmesh to
//		 * local id in AMesh*/
//		std::set<Region*>::const_iterator it_r_submesh= regions_submesh.begin();
//		for(;it_r_submesh!=regions_submesh.end();it_r_submesh++)
//		{
//			std::set<Region*>::const_iterator it_r_mesh= regions_mesh.begin();
//			bool found = false;
//			for(;it_r_mesh!=regions_mesh.end() && !found;it_r_mesh++)
//			{
//				if((*it_r_mesh)->getID()==(*it_r_submesh)->getID())
//				{
//					found=true;
//					ARegionProj[(*it_r_submesh)->getID()]= (*it_r_mesh)->getID();
//					ASubMesh.mark(*it_r_submesh,AInterfaceMark);
//					if(!ASubMesh.isSharedWith(*it_r_submesh,AFromPartitionID))
//					{
//						AMesh.removeFromInterface(*it_r_mesh,AFromPartitionID);
//					}
//				}
//			}
//		}
//	}
//
//}
/*----------------------------------------------------------------------------*/
