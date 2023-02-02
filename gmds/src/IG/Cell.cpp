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
 * Cell.cpp
 *
 *  Created on: 5 févr. 2014
 *      Author: ledouxf
 */
/*----------------------------------------------------------------------------*/
#include <GMDS/IG/Cell.h>
#include <GMDS/IG/Node.h>
#include <GMDS/IG/Edge.h>
#include <GMDS/IG/Face.h>
#include <GMDS/IG/Region.h>
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
Cell::Cell(IGMesh* AMesh,const ECellType& AType, const TCellID& AID)
:m_owner(AMesh), m_type(AType), m_id(AID)
{;}
/*----------------------------------------------------------------------------*/
TCellID Cell::getID() const
{
	return m_id;
}
/*----------------------------------------------------------------------------*/
ECellType  Cell::getType() const
{
	return m_type;
}
/*----------------------------------------------------------------------------*/
void
Cell::computeBoundingBox(TCoord minXYZ[3], TCoord maxXYZ[3]) const
{
	if(GMDS_NODE == getType()) {
		throw GMDSException("Cell::computeBoundingBox not available for nodes.");
	}

        minXYZ[0] =  HUGE_VALF;
        minXYZ[1] =  HUGE_VALF;
        minXYZ[2] =  HUGE_VALF;
        maxXYZ[0] = -HUGE_VALF;
        maxXYZ[1] = -HUGE_VALF;
        maxXYZ[2] = -HUGE_VALF;

	std::vector<Node> nodes = get<Node>();

        for(int iNode=0; iNode<nodes.size(); iNode++) {
                if(nodes[iNode].X() < minXYZ[0]) minXYZ[0] = nodes[iNode].X();
                if(nodes[iNode].Y() < minXYZ[1]) minXYZ[1] = nodes[iNode].Y();
                if(nodes[iNode].Z() < minXYZ[2]) minXYZ[2] = nodes[iNode].Z();
                if(nodes[iNode].X() > maxXYZ[0]) maxXYZ[0] = nodes[iNode].X();
                if(nodes[iNode].Y() > maxXYZ[1]) maxXYZ[1] = nodes[iNode].Y();
                if(nodes[iNode].Z() > maxXYZ[2]) maxXYZ[2] = nodes[iNode].Z();
        }
}
/*------------------------------------------------------------------------*/
void Cell::serializeCellData(std::ostream& AStr) const{
	//the mesh owner pointer is not serialize
	AStr.write((char*)&m_type, sizeof(ECellType));
	AStr.write((char*)&m_id	 , sizeof(TCellID)	);
}
/*------------------------------------------------------------------------*/
void Cell::unserializeCellData(std::istream& AStr) {
	AStr.read((char*)&m_type, sizeof(ECellType)	);
	AStr.read((char*)&m_id	, sizeof(TCellID)	);

//	/* we have to free the indirect connectivity before creating new ones. To
//	 * do that, we remove all the connections (direct and indirect).
//	 */
//	removeAllConnections();
//
//	AStr.read((char*)&adjDirect_  [0],descriptor::sizeDirect  *sizeof(id));
//
//	/* now we create new indirection connections:
//	 * we do not care about which type of connections it is (to nodes, edges,
//	 * faces, or regions. We just allocate new memory for that.
//	 */
//	LIDVectorAllocator<GChunkSize>* allocator =
//					this->mesh_->getUndefinedSizeConnectivityAllocator();
//	for(int i=0;i<descriptor::sizeIndirect;i++){
//		int nb_ids=0;
//		// get the number of ids stored in the ith connection
//		AStr.read((char*)&nb_ids,sizeof(id));
//		if(nb_ids!=0){
//			id* p = allocator->allocate(nb_ids);
//			/* p[0] stores nb_ids ids */
//
//			AStr.read((char*)&p[1],nb_ids*sizeof(id));
//			adjIndirect_[i]=p;
//		}
//	}
}

/*----------------------------------------------------------------------------*/
template<> void Cell::getIDs<Node>(std::vector<TCellID>& ACells) const {
	delegateGetNodeIDs(ACells);
}
/*----------------------------------------------------------------------------*/
template<> void Cell::getIDs<Edge>(std::vector<TCellID>& ACells) const {
	delegateGetEdgeIDs(ACells);
}
/*----------------------------------------------------------------------------*/
template<> void Cell::getIDs<Face>(std::vector<TCellID>& ACells) const {
	delegateGetFaceIDs(ACells);
}
/*----------------------------------------------------------------------------*/
template<> void Cell::getIDs<Region>(std::vector<TCellID>& ACells) const {
	delegateGetRegionIDs(ACells);
}
/*----------------------------------------------------------------------------*/
template<> void Cell::getAllIDs<Node>(std::vector<TCellID>& ACells) const {
        delegateGetAllNodeIDs(ACells);
}
/*----------------------------------------------------------------------------*/
template<> void Cell::getAllIDs<Edge>(std::vector<TCellID>& ACells) const {
        delegateGetAllEdgeIDs(ACells);
}
/*----------------------------------------------------------------------------*/
template<> void Cell::getAllIDs<Face>(std::vector<TCellID>& ACells) const {
        delegateGetAllFaceIDs(ACells);
}
/*----------------------------------------------------------------------------*/
template<> void Cell::getAllIDs<Region>(std::vector<TCellID>& ACells) const {
        delegateGetAllRegionIDs(ACells);
}
/*----------------------------------------------------------------------------*/
template<> void Cell::set<Node>(const std::vector<TCellID>& ACells)  {
	delegateSetNodeIDs(ACells);
}
/*----------------------------------------------------------------------------*/
template<> void Cell::set<Edge>(const std::vector<TCellID>& ACells)  {
	delegateSetEdgeIDs(ACells);
}
/*----------------------------------------------------------------------------*/
template<> void Cell::set<Face>(const std::vector<TCellID>& ACells)  {
	delegateSetFaceIDs(ACells);
}
/*----------------------------------------------------------------------------*/
template<> void Cell::set<Region>(const std::vector<TCellID>& ACells)  {
	delegateSetRegionIDs(ACells);
}
/*----------------------------------------------------------------------------*/
template<> void Cell::set<Node>(const std::vector<Node>& ACells)  {
	delegateSetNodeIDs(convertCellToID<Node>(ACells));
}
/*----------------------------------------------------------------------------*/
template<> void Cell::set<Edge>(const std::vector<Edge>& ACells)  {
	delegateSetEdgeIDs(convertCellToID<Edge>(ACells));
}
/*----------------------------------------------------------------------------*/
template<> void Cell::set<Face>(const std::vector<Face>& ACells)  {
	delegateSetFaceIDs(convertCellToID<Face>(ACells));
}
/*----------------------------------------------------------------------------*/
template<> void Cell::set<Region>(const std::vector<Region>& ACells)  {
	delegateSetRegionIDs(convertCellToID<Region>(ACells));
}

/*----------------------------------------------------------------------------*/
template<> void Cell::add<Node>(TCellID AElt) {
	delegateNodeAdd(AElt);
}
/*----------------------------------------------------------------------------*/
template<> void Cell::add<Edge>(TCellID AElt) {
	delegateEdgeAdd(AElt);
}
/*----------------------------------------------------------------------------*/
template<> void Cell::add<Face>(TCellID AElt) {
	delegateFaceAdd(AElt);
}
/*----------------------------------------------------------------------------*/
template<> void Cell::add<Region>(TCellID AElt) {
	delegateRegionAdd(AElt);
}
/*----------------------------------------------------------------------------*/
template<> void Cell::remove<Node>(TCellID AElt) {
	delegateNodeRemove(AElt);
}
/*----------------------------------------------------------------------------*/
template<> void Cell::remove<Edge>(TCellID AElt) {
	delegateEdgeRemove(AElt);
}
/*----------------------------------------------------------------------------*/
template<> void Cell::remove<Face>(TCellID AElt) {
	delegateFaceRemove(AElt);
}
/*----------------------------------------------------------------------------*/
template<> void Cell::remove<Region>(TCellID AElt) {
	delegateRegionRemove(AElt);
}
/*----------------------------------------------------------------------------*/
template<> void Cell::replace<Node>(TCellID AID1, TCellID AID2) {
	delegateNodeReplace(AID1,AID2);
}
/*----------------------------------------------------------------------------*/
template<> void Cell::replace<Edge>(TCellID AID1, TCellID AID2) {
	delegateEdgeReplace(AID1,AID2);
}
/*----------------------------------------------------------------------------*/
template<> void Cell::replace<Face>(TCellID AID1, TCellID AID2) {
	delegateFaceReplace(AID1,AID2);
}
/*----------------------------------------------------------------------------*/
template<> void Cell::replace<Region>(TCellID AID1, TCellID AID2) {
	delegateRegionReplace(AID1,AID2);
}
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/
