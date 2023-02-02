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
/** \file    MeshInsertDetailInOut.h
 *  \author  N. LE GOFF
 *  \date    27/08/2012
 */
/*----------------------------------------------------------------------------*/
#ifndef GMDS_MESHINSERTDETAILINOUT_H_
#define GMDS_MESHINSERTDETAILINOUT_H_
/*----------------------------------------------------------------------------*/
#include <set>
#include <gts.h>
/*----------------------------------------------------------------------------*/
#include "GMDS/IG/IGMesh.h"
#include "GMDS/CAD/GeomManager.h"
#include "GMDS/CAD/FacetedSurface.h"
#include "GMDS/Utils/Exception.h"
#include "GMDS/IG/IG.h"

//#include "ModelGraph.h"
//#include "LaplacianSmoothingGeomClassificationNew.h"
//#include "GMDSMeshTools/SheetOperator.h"
#include "MeshInsertDetail.h"
#include "GeomMeshIntersectionService.h"
//#include "MeshGraph.h"
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
// contains the nodes ids of the faces ordered with an outward normal
const int MESH_INSERT_DETAIL_INOUT_HEX_ORDERED_FACES[6][4] = {
		{0,3,2,1},
		{7,4,5,6},
		{5,4,0,1},
		{7,6,2,3},
		{4,7,3,0},
		{6,5,1,2}
};

// number of points used to discretize a loop curve
const unsigned int MESH_INSERT_DETAIL_INOUT_CURVE_NBPOINTS = 3;

// number of candidate nodes that will be studied for each point classification
//const unsigned int MESH_INSERT_DETAIL_INOUT_NB_NODES_KEPT = 10;
const unsigned int MESH_INSERT_DETAIL_INOUT_NB_NODES_KEPT = 1;
/*----------------------------------------------------------------------------*/
/**
 *  \class MeshInsertDetailInOut
 *
 *  \brief Interface
 */
/*----------------------------------------------------------------------------*/
class MeshInsertDetailInOut: public MeshInsertDetail {
public:

    /*------------------------------------------------------------------------*/
    /** \brief Constructor.
     *
     *  \param AMesh the mesh to build.
     */
	MeshInsertDetailInOut(
		IGMesh& AMesh,
		gmds::geom::GeomManager& AManager,
		gmds::geom::GeomMeshIntersectionService& AService);

    /*------------------------------------------------------------------------*/
    /** \brief  Destructor.	*/
	virtual ~MeshInsertDetailInOut();

//    /*------------------------------------------------------------------------*/
//	/** \brief  check the validity of the mesh model and the geometry model
//	 * 			for the desired algorithm.
//	 *
//	 * 			Currently empty; always returns true.
//	 *
//	 *  \return a boolean
//	 */
//	virtual bool checkValidity();
//
	/*------------------------------------------------------------------------*/
	/** \brief  insert the geometric detail
	 *
	 */
	virtual void autoInsert(gmds::geom::GeomVolume* AVol);

	/*------------------------------------------------------------------------*/
	/** \brief  project the nodes on the model depending on their classification
	 *
	 */
	virtual void project();
//
//	/*------------------------------------------------------------------------*/
//	/** \brief  project the nodes on the model depending on their classification
//	 */
//	void insertFunSheets(gmds::geom::GeomVolume* vol);

	virtual void exportCurvesEdgesVTK(const std::string& AFile);

protected:

//	void exportVTKBis(gmds::geom::GeomManager& model,
//			IGMesh& mesh,
//			std::string filename);
	/*------------------------------------------------------------------------*/
	/** \brief  associate nodes/edges/faces/regions to the geometric model
	 *
	 */
	virtual void associateGeometryClassification(gmds::geom::GeomVolume* AVol);

	/*------------------------------------------------------------------------*/
	/** \brief  associate nodes to the model points
	 *
	 */
	virtual void associateGeometryClassificationPoints(gmds::geom::GeomVolume* AVol);

	/*------------------------------------------------------------------------*/
	/** \brief  associate a node to a model point
	 *
	 */
	virtual void associateGeometryClassificationPoint(
		gmds::geom::GeomPoint* APoint,
		gmds::geom::GeomVolume* AVol);

	/*------------------------------------------------------------------------*/
	/** \brief  select the best node that fits a model point.
	 *
	 */
	virtual Node associateGeometryClassificationPointSelectBestNode(
		std::vector<Node>& ANodes,
		gmds::geom::GeomPoint* APoint,
		gmds::geom::GeomVolume* AVol);

	/*------------------------------------------------------------------------*/
	/** \brief  compute the cost (along which criteria?)
	 * 			if ANode were to be associated with APoint.
	 *
	 */
	virtual double computeNodePointClassificationCost(
			gmds::Node ANode,
			gmds::geom::GeomPoint* APoint,
			gmds::geom::GeomVolume* AVol,
			std::vector<Edge>& AEdges);

	/*------------------------------------------------------------------------*/
	/** \brief  return direct ordered edges adjacent to ANode
	 *
	 */
	virtual void getOrderedDirectEdges(
			std::vector<Edge>& AEdges,
			Node ANode,
			gmds::geom::GeomPoint* APoint,
			gmds::geom::GeomVolume* AVol);

	/*------------------------------------------------------------------------*/
	/** \brief  associate edges to the curves at the model points.
	 *
	 * 			Edges are ordered, same as the curves.
	 *
	 */
	virtual void associateGeometryClassificationCurvesAtPoints(gmds::geom::GeomVolume* AVol);

	/*------------------------------------------------------------------------*/
	/** \brief  associate edges to the model curves
	 *
	 */
	virtual void associateGeometryClassificationCurves(gmds::geom::GeomVolume* AVol);

	/*------------------------------------------------------------------------*/
	/** \brief  remove edge association except
	 *
	 */
	virtual void removeGeometryClassificationCurves(gmds::geom::GeomVolume* AVol);

	/*------------------------------------------------------------------------*/
	/** \brief  associate edges to a model curve
	 *
	 *  \param ACurve a curve
	 *  \param AMarkNodeConstraint a mark that prohibits selecting edges
	 *  	   owning one of those nodes
	 */
	virtual bool associateGeometryClassificationEdge(
			double& ADistance,
			gmds::geom::GeomCurve& ACurve,
			const int& AMarkNodeConstraint
			);

	/*------------------------------------------------------------------------*/
	/** \brief  associate regions to the model volume.
	 *
	 * 			Only one volume at the moment, and it keeps the inside and intersecting
	 * 			hexahedra.
	 *
	 */
	virtual void associateGeometryClassificationCells(gmds::geom::GeomVolume* AVol);

	/*------------------------------------------------------------------------*/
	/** \brief  associate faces to the model surfaces.
	 *
	 * 			Works using an face propagation stopped by the edges classification,
	 * 			and characterizes the surface by its curves.
	 * 			Should not work when there are no curves (sphere) are only one (two
	 * 			overlapping spheres).
	 *
	 */
	virtual void associateGeometryClassificationSurfaces(gmds::geom::GeomVolume* AVol);

	/*------------------------------------------------------------------------*/
	/** \brief  marks faces delimited by an edge mark among a pool of faces
	 *
	 *  \param AFace a starting face
	 *  \param AMark a mark that will be propagated
	 *  \param AMarkAcceptable a mark designating the faces considered for propagation.
	 *  	   Once mark with AMark, faces are removed from the AMarkAcceptable pool.
	 *  \param AMarkBarrier a mark (carried by edges) delimiting the propagation
	 *  \param AEdgeClassification the edges geometric classification.
	 *  	   It is used only for keeping track of the curves (associated to the
	 *  	   delimiting edges) encountered.
	 *  \param ACurvesSet a set of curves.
	 *  	   It is filled in order to later characterize the surface AFace is part of.
	 *
	 */
	virtual void markRecurse(
		Face AFace, const int& AMark, const int& AMarkAcceptable, const int& AMarkBarrier,
		Variable<geom::GeomEntity* >* AEdgeClassification,
		std::set<geom::GeomEntity* >& ACurvesSet);
//
////    /*------------------------------------------------------------------------*/
////	/** \brief Find and split hexahedra that contains a node given in parameter.
////	 *
////	 * 		   It will split around other nodes, nodes opposite the given node
////	 * 		   (opposite in the hexahedra).
////	 * 		   It is so that we will add edges incident to the given node.
////	 *
////	 *  \param AID a node we want to split cells around of.
////	 */
////	virtual void brickSplitFind(const id& AID);
////
////    /*------------------------------------------------------------------------*/
////	/** \brief Find and split hexahedra that contains a node given in parameter
////	 *
////	 *  \param AID a node we want to split cells around of.
////	 */
////	virtual void brickSplitDo(const id& AID);
////
////	/*------------------------------------------------------------------------*/
////    /** \brief Remove the mesh entities that were deleted from the mesh from the
////     * 		   connectivities.
////     *
////     */
////	void purgeMesh();
////
////	virtual void splitMeshAroundNodes();
//
    /*------------------------------------------------------------------------*/
	/** \brief Propagates a mark on the mesh, face-wise. Propagation is stopped
	 * 		   by a barrier mark and fails if it encounters a foe mark.
	 *
	 *  \param ACell the origin cell.
	 *  \param AMark the mark that will be propagated.
	 *  \param AMarkBarrier the mark that acts as a barrier. There should be a
	 *  		layer of such cells between a mark area and a foe mark area.
	 *  \param AMarkFoe the mark that must not be encountered.
	 *  */
	void propagateMark(Region ARegion, const int& AMark, const int& AMarkBarrier, const int& AMarkFoe);
//
//	/*------------------------------------------------------------------------*/
//	/** \brief refine mesh near geometric details in order to avoid robustness issues.
//	 *
//	 */
//	void refineMeshNearGeometry();
//
//	/*------------------------------------------------------------------------*/
//	/** \brief clean a few
//	 *
//	 */
//	void cleanAFewThings();
//
//	/*------------------------------------------------------------------------*/
//	/** \brief  reassociate regions to the model volume.
//	 *
//	 */
////	virtual void reassociateVolume(gmds::geom::GeomVolume<TBase>* AVol);
//	virtual void reassociateVolumeUsingNodeSurfaceCriteria(gmds::geom::GeomVolume* AVol);
//
	/*------------------------------------------------------------------------*/
	/** \brief  propagate across regions that are listed as "inside"
	 *
	 */
	virtual void propagateVisited(gmds::Region& ARegion, std::map<Region,bool>& AIsInRegionMap, std::map<Region,bool>& AWasVisited);

	/*------------------------------------------------------------------------*/
	/** \brief  propagate across faces that are listed as "inside"
	 *
	 */
	virtual void propagateVisited(gmds::Face& AFace, std::map<Face,bool>& AIsInFaceMap, std::map<Face,bool>& AWasVisited);

	/*------------------------------------------------------------------------*/
	/** \brief check whether the sub-mesh formed by the regions associated to AVol
	 * 	       is non-manifold. It is based on edges.
	 *
	 *  \param AVol a volume.
	 *
	 *  \return a boolean.
	 */
	virtual bool checkNonManifoldnessEdges(gmds::geom::GeomVolume* AVol);

//	/*------------------------------------------------------------------------*/
//	/** \brief check whether the sub-mesh formed by the regions associated to AVol
//	 * 	       is manifold. It is based on nodes.
//	 *
//	 *  \param AVol a volume.
//	 *
//	 *  \return a boolean.
//	 */
//	virtual bool checkNonManifoldnessNodes(gmds::geom::GeomVolume* AVol);
//
//	/*------------------------------------------------------------------------*/
//	/** \brief check whether the sub-mesh formed by the regions associated to AVol
//	 * 	       is manifold. It is based on edges, and returns the first edge encountered
//	 * 	       that highlights this manifoldness.
//	 *
//	 *  \param AVol a volume.
//	 *
//	 *  \return an edge.
//	 */
//	virtual Edge* checkNonManifoldnessWithReturnEdges(gmds::geom::GeomVolume* AVol);
//
//	/*------------------------------------------------------------------------*/
//	/** \brief check whether the sub-mesh formed by the regions associated to AVol
//	 * 	       is manifold. It is based on nodes, and returns the first node encountered
//	 * 	       that highlights this manifoldness.
//	 *
//	 *  \param AVol a volume.
//	 *
//	 *  \return an edge.
//	 */
//	virtual Node* checkNonManifoldnessWithReturnNodes(gmds::geom::GeomVolume* AVol);

	/*------------------------------------------------------------------------*/
	/** \brief build the boundary information. That means marking mesh elements
	 * 		   with the markIsInside_ and markIsOnBoundary_ marks attributes.
	 *
	 *  \param AVol a volume.
	 */
	virtual void buildBoundaryInfo(gmds::geom::GeomVolume* AVol);

	/*------------------------------------------------------------------------*/
	/** \brief build entities (edges and nodes) on the boundary.
	 * 		   It is solely based on markFaceIsOnBoundary.
	 *
	 *  \param AVol a volume.
	 */
	virtual void buildBoundaryEntities();

	/*------------------------------------------------------------------------*/
	/** \brief Associate nodes to curve, surface and volume of AVol.
	 *
	 *  \param AVol a volume.
	 */
	virtual void completeNodeAssociation(gmds::geom::GeomVolume* AVol);

//	/*------------------------------------------------------------------------*/
//	/** \brief correct manifoldness
//	 *
//	 *  \param AVol a volume.
//	 *
//	 *  \return a boolean true if the resulting mesh classified to AVol is non-manifold.
//	 *  		It will return false only if it fails to produce a non-manifold
//	 *  		AVol classification, most likely if the whole mesh itself is manifold.
//	 */
//	virtual bool correctManifoldness(gmds::geom::GeomVolume* AVol);
//
	/*------------------------------------------------------------------------*/
	/** \brief count manifold edges
	 *
	 *  \param AVol a volume.
	 *
	 *  \return the number of manifold edges
	 */
	virtual unsigned int countManifoldEdges(gmds::geom::GeomVolume* AVol);


	virtual void initialization();

	virtual unsigned int intersectedCellsDetectionAABBTree(std::map<gmds::TCellID, std::vector<gmds::math::Triangle> >& AInOutTriangles);
	virtual unsigned int intersectedCellsDetection(std::map<gmds::TCellID, std::vector<gmds::math::Triangle> >& AInOutTriangles);


//	/* a mesh */
//	Mesh<TMask>& mesh_;
//
//	/* a geometric model */
//	gmds::geom::FacetedGeomManager<TBase>& manager_;
//
//	/* service associated to the geometric model */
//	gmds::geom::GeomMeshIntersectionService<TBase>& service_;

	/* association between GeomNode and mesh node*/
	std::map<gmds::geom::GeomPoint*,Node> geom2MeshNode_;

	/* association between GeomCurve and mesh edge*/
	std::map<gmds::geom::GeomCurve*,std::vector<Edge> > geom2MeshEdge_;

	/* association between GeomSurface and mesh face*/
	std::map<gmds::geom::GeomSurface*,std::vector<Face> > geom2MeshFace_;

	std::map<gmds::geom::GeomPoint*,std::vector<Edge> > geomCurvesAtPoints2MeshEdge_;
//	std::vector<std::vector<std::vector<id> > > geomSurfacesAtPoints2MeshFace_;
//
//
	// mark carried by inside entities
	int markRegionIsInside_;
	int markFaceIsInside_;
	int markEdgeIsInside_;
	int markNodeIsInside_;


	// mark carried by entities on the boundary of the model
	int markRegionIsOnBoundary_;
	int markFaceIsOnBoundary_;
	int markEdgeIsOnBoundary_;
	int markNodeIsOnBoundary_;

	//
	Variable<geom::GeomEntity* >* previousVolumeClassification_;

	Variable<double>* insideRatio_;

	// Axis-Aligned Bounding Box tree for the surfaces triangles
	GNode* aabbSurfacesTrianglesTree_;

	// Axis-Aligned Bounding Box tree for the surfaces triangles;
	// one tree for each surface
	std::map<gmds::geom::GeomSurface*,GNode*> aabbSurfacesTrianglesTrees_;

	// Axis-Aligned Bounding Box tree for the regions
	// The bounded value is the CellID
	GNode* aabbRegionsTree_;
};
/*----------------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------------*/
#endif /* GMDS_MESHINSERTDETAILINOUT_H_ */
/*----------------------------------------------------------------------------*/
