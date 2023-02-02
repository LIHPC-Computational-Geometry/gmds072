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
/** \file    MeshSplitter3Refinement.h
 *  \author  legoff
 *  \date    23/01/2015
 */
/*----------------------------------------------------------------------------*/
#ifndef GMDS_MESHSPLITTER3REFINEMENT_H_
#define GMDS_MESHSPLITTER3REFINEMENT_H_
/*----------------------------------------------------------------------------*/
#include <set>
/*----------------------------------------------------------------------------*/
#include "GMDS/IG/IGMesh.h"
#include "GMDS/CAD/GeomManager.h"
#include "GeomMeshIntersectionService.h"
/*----------------------------------------------------------------------------*/
#include "MeshSplitter.h"
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
const int MESH_SPLITTER_3_REFINEMENT_HEX_ORDERED_FACES[6][4] = {
		{0,3,2,1}, // top
		{4,7,6,5}, // bottom
		{4,0,1,5}, // left
		{7,3,2,6}, // right
		{4,7,3,0}, // front
		{5,6,2,1}  // back
};
/*----------------------------------------------------------------------------*/
struct RegionIDs{

	gmds::TCellID ids[8];
};
struct FaceIDs{

	gmds::TCellID ids[4];
};
struct NodeHexIJK{

	gmds::TCellID ids[3];
};
struct HexIJK{

	NodeHexIJK ijk[8];
};


/*----------------------------------------------------------------------------*/
template<int TSize>
class tableMarkedNodes {
public :
	tableMarkedNodes() {
		for(unsigned int iMark = 0; iMark<TSize; iMark++) {
			mark_[iMark] = 0;
		}
	};

	tableMarkedNodes(const bool& AMark1, const bool& AMark2, const bool& AMark3, const bool& AMark4) {

		mark_[0] = AMark1;
		mark_[1] = AMark2;
		mark_[2] = AMark3;
		mark_[3] = AMark4;
	};

	tableMarkedNodes(const bool& AMark1, const bool& AMark2, const bool& AMark3, const bool& AMark4,
					 const bool& AMark5, const bool& AMark6, const bool& AMark7, const bool& AMark8) {

		mark_[0] = AMark1;
		mark_[1] = AMark2;
		mark_[2] = AMark3;
		mark_[3] = AMark4;
		mark_[4] = AMark5;
		mark_[5] = AMark6;
		mark_[6] = AMark7;
		mark_[7] = AMark8;
	};

	tableMarkedNodes(const bool AMarks[TSize]) {
		for(unsigned int iMark = 0; iMark<TSize; iMark++) {
			mark_[iMark] = AMarks[iMark];
		}
	};

	tableMarkedNodes(const tableMarkedNodes<TSize>& ATableMarkedNodes) {
		for(unsigned int iMark=0; iMark<TSize; iMark++) {
			mark_[iMark] = ATableMarkedNodes.mark_[iMark];
		}
	};

	~tableMarkedNodes() {

	};

	tableMarkedNodes<TSize>& operator=(const tableMarkedNodes<TSize>& ATableMarkedNodes) {
		for(unsigned int iMark=0; iMark<TSize; iMark++) {
			mark_[iMark] = ATableMarkedNodes.mark_[iMark];
		}

		return *this;
	};

	bool operator<(const tableMarkedNodes<TSize>& ATableMarkedNodes) const {
		for(unsigned int iMark=0; iMark<TSize; iMark++) {
			if(!mark_[iMark] && ATableMarkedNodes.mark_[iMark]) {
				return true;
			}
			if(mark_[iMark] && !ATableMarkedNodes.mark_[iMark]) {
				return false;
			}
		}

		return false;
	};

	bool operator>(const tableMarkedNodes<TSize>& ATableMarkedNodes) const {
		for(unsigned int iMark=0; iMark<TSize; iMark++) {
			if(mark_[iMark] && !ATableMarkedNodes.mark_[iMark]) {
				return true;
			}
			if(!mark_[iMark] && ATableMarkedNodes.mark_[iMark]) {
				return false;
			}
		}

		return false;
	};

	bool operator!=(const tableMarkedNodes<TSize>& ATableMarkedNodes) const {
		for(unsigned int iMark=0; iMark<TSize; iMark++) {
			if(mark_[iMark] != ATableMarkedNodes.mark_[iMark]) {
				return true;
			}
		}

		return false;
	};

	bool operator==(const tableMarkedNodes<TSize>& ATableMarkedNodes) const {
		for(unsigned int iMark=0; iMark<TSize; iMark++) {
			if(mark_[iMark] != ATableMarkedNodes.mark_[iMark]) {
				return false;
			}
		}

		return true;
	};

	bool isMarked(const int& AIndex) const {
		return mark_[AIndex];
	};

	bool isEqual(const tableMarkedNodes<TSize>& ATableMarkedNodes) const {
		return (*this==ATableMarkedNodes);
	};

private :
	bool mark_[TSize];
};
/*----------------------------------------------------------------------------*/
template<int TSize>
class compare_tableMarkedNodes {
public :
	bool operator()(const tableMarkedNodes<TSize>& ATableMarkedNodes1, const tableMarkedNodes<TSize>& ATableMarkedNodes2) {
		return (ATableMarkedNodes1<ATableMarkedNodes2);
	}
};
/*----------------------------------------------------------------------------*/
const bool splittingConvert[256][8] =
{
		{0,0,0,0,0,0,0,0}, // 00000000 //   0
		{1,0,0,0,0,0,0,0}, // 10000000
		{0,1,0,0,0,0,0,0}, // 01000000
		{1,1,0,0,0,0,0,0}, // 11000000
		{0,0,1,0,0,0,0,0}, // 00100000
		{1,1,1,1,0,0,0,0}, // 10100000
		{0,1,1,0,0,0,0,0}, // 01100000
		{1,1,1,1,0,0,0,0}, // 11100000
		{0,0,0,1,0,0,0,0}, // 00010000
		{1,0,0,1,0,0,0,0}, // 10010000
		{1,1,1,1,0,0,0,0}, // 01010000 //  10
		{1,1,1,1,0,0,0,0}, // 11010000
		{0,0,1,1,0,0,0,0}, // 00110000
		{1,1,1,1,0,0,0,0}, // 10110000
		{1,1,1,1,0,0,0,0}, // 01110000
		{1,1,1,1,0,0,0,0}, // 11110000 //  15



		{0,0,0,0,1,0,0,0}, // 00001000
		{0,0,0,0,0,1,0,0}, // 00000100
		{0,0,0,0,0,0,1,0}, // 00000010
		{0,0,0,0,0,0,0,1}, // 00000001
};
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/**
 *  \class MeshSplitter3Refinement
 *
 *  \brief Interface
 */
/*----------------------------------------------------------------------------*/
class MeshSplitter3Refinement: public MeshSplitter {

public:

    /*------------------------------------------------------------------------*/
    /** \brief Constructor.
     *
     *  \param AMesh the mesh to split.
     */
	MeshSplitter3Refinement(
			gmds::IGMesh& AMesh,
			gmds::geom::GeomManager& AManager,
			gmds::geom::GeomMeshIntersectionService& AService);

    /*------------------------------------------------------------------------*/
    /** \brief  Destructor.
     *
     */
	virtual ~MeshSplitter3Refinement();

	virtual void exportVTK(const std::string& AFile);

	std::set<Region> getSplitRegions() const;
	std::set<Face> getSplitFaces() const;
	std::set<Edge> getSplitEdges() const;
	std::set<Region> getCreatedRegions() const;
	std::set<Face> getCreatedFaces() const;
	std::set<Edge> getCreatedEdges() const;
	std::set<Node> getCreatedNodes() const;

protected:

	/*------------------------------------------------------------------------*/
	/** \brief mark the regions that will be split.
	 *
	 */
	virtual void markRegionsToSplit(
			const int& AMarkRegionToSplit,
			const int& AMarkFaceToSplit,
			const int& AMarkEdgeToSplit,
			const int& AMarkNodeToSplit);

	/*------------------------------------------------------------------------*/
        /** \brief mark the entities that will be split in order to add edges around 
 	 *	a node.
         *
         */
        virtual void markEntitiesToSplitForAddEdges(
			gmds::Node ANode,
			gmds::Face AFace,
                        const int& AMarkRegionToSplit,
                        const int& AMarkFaceToSplit,
                        const int& AMarkEdgeToSplit,
                        const int& AMarkNodeToSplit);

	/*------------------------------------------------------------------------*/
	/** \brief mark the regions that will be split in order to avoid bad
	 * 		   configurations depending on the refinement scheme used.
	 *
	 */
	virtual void avoidBadConfigurations(
			const int& AMarkRegionToSplit,
			const int& AMarkFaceToSplit,
			const int& AMarkEdgeToSplit,
			const int& AMarkNodeToSplit);

	/*------------------------------------------------------------------------*/
	/** \brief split the mesh.
	 *
	 */
	virtual void split(
			const int& AMarkRegionToSplit,
			const int& AMarkFaceToSplit,
			const int& AMarkEdgeToSplit,
			const int& AMarkNodeToSplit);

	/*------------------------------------------------------------------------*/
	/** \brief split an edge.
	 *
	 */
	virtual bool splitEdge(const int& AMarkNodeToSplit, const int& AMarkEdgeToSplit, Edge& AEdge);

	/*------------------------------------------------------------------------*/
	/** \brief split a face.
	 *
	 */
	virtual bool splitFace(const int& AMarkNodeToSplit, const int& AMarkFaceToSplit, Face& AFace);

	/*------------------------------------------------------------------------*/
	/** \brief split a region.
	 *
	 */
	virtual void splitRegion(const int& AMarkNodeToSplit, Region& ARegion);

	/*------------------------------------------------------------------------*/
	/** \brief get node of edge by id.
	 *
	 *  \param AEdge an edge
	 *  \param AId an id
	 */
	virtual Node getNodeOfEdge(Edge& AEdge, const TCellID& AId);

	/*------------------------------------------------------------------------*/
	/** \brief get node of face by id.
	 *
	 *  \param AFace a face
	 *  \param AId an id
	 */
	virtual Node getNodeOfFace(Face& AFace, const TCellID& AId);

	/*------------------------------------------------------------------------*/
	/** \brief get node id of face from index I and J.
	 *
	 *  \param AOffset
	 *  \param AIndexI
	 *  \param AIndexJ
	 *
	 *  \return a TCellID
	 */
	virtual TCellID getNodeIdOfFace(const TCellID& AOffset, const bool& AIsDirect, const TCellID& AIndexI, const TCellID& AIndexJ);

	/*------------------------------------------------------------------------*/
	/** \brief check whether ARegion is intersected twice by the model.
	 *		   Intersected by two surfaces that are not neighbors?
	 *
	 */
	virtual bool regionIsIntersectedTwiceByModel(Region& ARegion);

	/*------------------------------------------------------------------------*/
	/** \brief check whether ARegion is intersected by the model.
	 *
	 */
	virtual bool regionIsIntersectedByModel(Region& ARegion);

	/*------------------------------------------------------------------------*/
	/** \brief Based on nodes being in or out between neighbor regions
	 *
	 */
	virtual void markNeighborRegionsInOutModelAABBTree(
                        const int& AMarkRegionToSplit);

	virtual void markNeighborRegionsInOutModel(
			const int& AMarkRegionToSplit);

	/*------------------------------------------------------------------------*/
	/** \brief Based on nodes being in or out between neighbor regions
	 *
	 */
	virtual void markRegionsIntersectedTwiceByModelAABBTree(
			const int& AMarkRegionToSplit);

	virtual void markRegionsIntersectedTwiceByModel(
				const int& AMarkRegionToSplit);

	/*------------------------------------------------------------------------*/
	/** \brief Based on nodes being in or out between neighbor regions
	 *
	 */
	virtual void markRegionsIntersectedByModel(
			const int& AMarkRegionToSplit);

	/*------------------------------------------------------------------------*/
	/** \brief build surfaces neighbors
	 *
	 */
	virtual void buildSurfacesNeighbors();


	virtual void initialization();

	/*------------------------------------------------------------------------*/
//	/* a mesh */
//	Mesh<TMask>& mesh_;
//
//	/* a geometric model */
//	gmds::geom::FacetedGeomManager<TBase>& manager_;
//
//	/* service associated to the geometric model */
//	gmds::geom::GeomMeshIntersectionService<TBase>& service_;

//	std::vector<std::set<gmds::geom::GeomSurface<TBase>* > > surfacesNeighbors_;
	std::map<gmds::geom::GeomSurface*, std::set<gmds::geom::GeomSurface*> > surfacesNeighbors_;

	//std::map<bool[8], bool[8]> lookupNodes_;
	//std::map<int, bool[8]> lookupNodes_;
	std::map<tableMarkedNodes<8>, tableMarkedNodes<8>, compare_tableMarkedNodes<8> > lookupNodes_;

	std::map<tableMarkedNodes<4>, unsigned int> nbNodesOnFaceToBuild_;
	std::map<tableMarkedNodes<4>, unsigned int> nbFacesOnFaceToBuild_;
	std::map<tableMarkedNodes<4>, std::vector<TCellID> > nodesOnFaceToBuild_;
	std::map<tableMarkedNodes<4>, std::vector<FaceIDs> > facesOnFaceToBuild_;

	std::map<tableMarkedNodes<8>, RegionIDs> nodePermutToUnitHex_;
	std::map<tableMarkedNodes<8>, unsigned int> nbNodesOnRegionToBuild_;
	std::map<tableMarkedNodes<8>, std::vector<NodeHexIJK> > nodesOnRegionToBuild_;
	std::map<tableMarkedNodes<8>, unsigned int> nbRegionsOnRegionToBuild_;
	std::map<tableMarkedNodes<8>, std::vector<HexIJK> > regionsOnRegionToBuild_;

//	std::map<tableMarkedNodes<8>, std::vector<id[8]> > hexToBuild_;

	std::map<Edge,std::vector<Node> > edgeSplitNodes_;
	std::map<Face,std::vector<Node> > faceSplitNodes_;

	std::set<Region> splitRegions_;
	std::set<Face  > splitFaces_;
	std::set<Edge  > splitEdges_;
	std::set<Region> splitRegionsCreated_;
	std::set<Face  > splitFacesCreated_;
	std::set<Edge  > splitEdgesCreated_;
	std::set<Node  > splitNodesCreated_;

	// Axis-Aligned Bounding Box tree for the regions
	// The bounded value is the CellID
	GNode* aabbRegionsTree_;

};
/*----------------------------------------------------------------------------*/
} // end namespace gmds
/*----------------------------------------------------------------------------*/
#endif /* GMDS_MESHSPLITTER3REFINEMENT_H_ */
/*----------------------------------------------------------------------------*/
