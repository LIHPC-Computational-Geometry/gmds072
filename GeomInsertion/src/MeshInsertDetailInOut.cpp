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
/** \file    MeshInsertDetailInOut.cpp
 *  \author  legoff
 *  \date    09/07/2014
 */
/*----------------------------------------------------------------------------*/
// stl File Headers
/*----------------------------------------------------------------------------*/
#include <algorithm>
#include <set>
/*----------------------------------------------------------------------------*/
// GMDS File Headers
/*----------------------------------------------------------------------------*/
#include <GMDS/CAD/GeomSurface.h>
#include <GMDS/Math/Hexahedron.h>
#include <GMDS/IG/IGMeshDoctor.h>
#include <GMDS/IO/VTKWriter.h>
#include <GMDS/MeshInsertDetailInOut.h>
#include <GMDS/MeshGraph.h>
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
MeshInsertDetailInOut::MeshInsertDetailInOut(
		IGMesh& AMesh,
		gmds::geom::GeomManager& AManager,
		gmds::geom::GeomMeshIntersectionService& AService)
: MeshInsertDetail(AMesh,AManager,AService)
{

}
/*----------------------------------------------------------------------------*/
//MeshInsertDetailInOut::MeshInsertDetailInOut(const MeshInsertDetailInOut& AMeshInsertDetailInOut)
//{
//
//}
/*----------------------------------------------------------------------------*/
MeshInsertDetailInOut::~MeshInsertDetailInOut()
{

}
///*----------------------------------------------------------------------------*/
//MeshInsertDetailInOut&
//MeshInsertDetailInOut::operator=(const MeshInsertDetailInOut& AMeshInsertDetailInOut)
//{
//	return *this;
//}
/*----------------------------------------------------------------------------*/
void
MeshInsertDetailInOut::autoInsert(gmds::geom::GeomVolume* AVol)
{
//	this->exportModelVTK("/homePOYOP/travail/workspaces/gscc_workspace/Cube/src/CaGe2/debug_model.mli");
//	this->exportModelVTK("/homePOYOP/travail/workspaces/gscc_workspace/Cube/src/CaGe2/debug_model");

	std::cout<<"autoInsert"<<std::endl;
	this->associateGeometryClassification(AVol);
}
/*----------------------------------------------------------------------------*/
void MeshInsertDetailInOut::
associateGeometryClassification(gmds::geom::GeomVolume* AVol)
{
	std::cout<<"=======> associateGeometryClassification"<<std::endl;

	std::vector<gmds::geom::GeomSurface*> surfaces;
	AVol->get(surfaces);
	this->service_.initialization(surfaces);
	this->initialization();
	this->aabbSurfacesTrianglesTree_ = this->service_.buildAABBSurfacesTriangulationTree(surfaces,this->aabbSurfacesTrianglesTrees_);

	associateGeometryClassificationCells(AVol);
	std::cout<<"=======> cell classification DONE"<<std::endl;
//#ifdef GMDS_DEBUG_INSERTION_FILES
//	this->exportVTK("debug_cells.mli");
//#endif
//
	//bool isNonManifold = checkNonManifoldnessEdges(AVol);
	//if(!isNonManifold) {
	//	throw GMDSException("MeshInsertDetailInOut::associateGeometryClassification" 
	//		"Selected cells form a manifold volume 1.");
	//}
//	std::cout<<"=======> checkNonManifoldness first DONE"<<std::endl;
//
	associateGeometryClassificationPoints(AVol);
//	std::cout<<"=======> point classification DONE"<<std::endl;
//	// maybe split hexahedra around selected nodes?
////	splitMeshAroundNodes();
////	isNonManifold = checkNonManifoldnessEdges(AVol);
////	if(!isNonManifold) {
////		throw GMDSException("Selected cells form a manifold volume 2.");
////	}
////	std::cout<<"=======> checkNonManifoldness second DONE"<<std::endl;
//
	associateGeometryClassificationCurvesAtPoints(AVol);
//#ifdef GMDS_DEBUG_INSERTION_FILES
	this->exportMeshVTK("debug_curvesatpoints.mli",gmds::R|gmds::F|gmds::E|gmds::N);
	this->exportCurvesEdgesVTK("debug_curvesatpoints_curves");
//#endif
////	isNonManifold = checkNonManifoldnessEdges(AVol);
////	if(!isNonManifold) {
////		throw GMDSException("Selected cells form a manifold volume 3.");
////	}
	std::cout<<"=======> curves at point classification DONE"<<std::endl;
	associateGeometryClassificationCurves(AVol);
//#ifdef GMDS_DEBUG_INSERTION_FILES
	this->exportMeshVTK("debug_edges.mli",gmds::R|gmds::F|gmds::E|gmds::N);
	this->exportCurvesEdgesVTK("debug_edges_curves");
//#endif
	std::cout<<"=======> curves classification DONE"<<std::endl;
//
////	this->exportVTK("debug_beforereassociate.unf");
////	reassociateVolumeUsingNodeSurfaceCriteria(AVol);
////	this->exportVTK("debug_afterreassociate.unf");
////	std::cout<<"=======> volume cleaning DONE"<<std::endl;
//
////	isNonManifold = checkNonManifoldnessEdges(AVol);
////	if(!isNonManifold) {
////		throw GMDSException("Selected cells form a manifold volume 4.");
////	}
	associateGeometryClassificationSurfaces(AVol);
//#ifdef GMDS_DEBUG_INSERTION_FILES
	this->exportMeshVTK("debug_surfaces.mli",gmds::R|gmds::F|gmds::E|gmds::N);
//#endif
//	std::cout<<"=======> surfaces classification DONE"<<std::endl;
//
	completeNodeAssociation(AVol);
//	cleanAFewThings();
//	std::cout<<"=======> cleaning DONE"<<std::endl;
}
/*----------------------------------------------------------------------------*/
void MeshInsertDetailInOut::
associateGeometryClassificationCells(gmds::geom::GeomVolume* AVol)
{
	// associate Cells
	Variable<geom::GeomEntity* >* volumeClassification = this->mesh_.getGeometricClassification(3);

	// fill previous classification
	{
		previousVolumeClassification_ =
				this->mesh_.newVariable<gmds::geom::GeomEntity*>(GMDS_REGION,"previousVolumeClassification");
		IGMesh::region_iterator it  = this->mesh_.regions_begin();

		for(;!it.isDone();it.next()){
			(*previousVolumeClassification_)[it.value().getID()] = (*volumeClassification)[it.value().getID()];
		}
	}

	// several marks that cells will carry
	int isInner = this->mesh_.getNewMark<gmds::Region>();
	int isOuter = this->mesh_.getNewMark<gmds::Region>();
	int isOnIntersection = this->mesh_.getNewMark<gmds::Region>();

	unsigned int nbIsOnIntersection = 0;
	int nbInner = 0;
	int nbOuter = 0;

	std::vector<Region> intersectedRegions;
	std::map<gmds::TCellID, std::vector<gmds::math::Triangle> > in_out_triangles;

	std::vector<gmds::geom::GeomSurface* > surfaces;
	AVol->get(surfaces);

	this->service_.initialization(surfaces);

	std::cout<<"associateGeometryClassificationCells nb surfaces "<<surfaces.size()<<std::endl;

	nbIsOnIntersection = intersectedCellsDetectionAABBTree(in_out_triangles);
	{
		std::map<gmds::TCellID, std::vector<gmds::math::Triangle> >::iterator it = in_out_triangles.begin();
		for(; it!=in_out_triangles.end(); it++) {

			gmds::Region current_region = this->mesh_.get<gmds::Region>(it->first);

			this->mesh_.mark(current_region,isOnIntersection);
			intersectedRegions.push_back(current_region);
		}
	}
	std::cout<<"Intersected cells are resolved!!! "<<nbIsOnIntersection<<std::endl;

	if(0 == nbIsOnIntersection) {
		throw GMDSException("MeshInsertDetailInOut::associateGeometryClassificationCells "
				"no cell was intersected.");
	}

	// debug output
	{
		IGMesh::volume& intersection = this->mesh_.newVolume("intersection");

		for(unsigned int k=0; k<intersectedRegions.size();k++) {
			intersection.add(intersectedRegions[k]);
		}

		gmds::VTKWriter<gmds::IGMesh> meshWriter(this->mesh_);
		meshWriter.write("debug_intersected_cells.mli",gmds::R|gmds::N);

//#ifdef GMDS_DEBUG_INSERTION_FILES
//	this->exportVTK("debug_intersected_cells.mli");
//#endif
		this->mesh_.deleteVolume(intersection);
	}

	// now marking inner and outer cells
	{
	IGMesh::region_iterator it  = this->mesh_.regions_begin();
	it  = this->mesh_.regions_begin();

	for(;!it.isDone();it.next()){

		gmds::Region current_region = it.value();

		// this cell was already treated
		if(this->mesh_.isMarked(current_region,isInner) || this->mesh_.isMarked(current_region,isOuter) || this->mesh_.isMarked(current_region,isOnIntersection))
			continue;

		std::cout<<"inside treating cell "<<current_region.getID()<<std::endl;

		// check whether the hex is inside or outside the model
		std::vector<gmds::geom::GeomSurface* > surfaces;
		AVol->get(surfaces);

		if(this->service_.isInsideMC(surfaces,current_region))
		{
			this->mesh_.mark(current_region,isInner);
		}
		else
		{
			this->mesh_.mark(current_region,isOuter);
		}

		// propagate the inside/outside mark, propagation face_wise that stops at intersecting cells.
		// Intersecting cells should form a water-tight envelope
		if(this->mesh_.isMarked(current_region,isInner)) {
			propagateMark(current_region,isInner,isOnIntersection,isOuter);
		}
		if(this->mesh_.isMarked(current_region,isOuter)) {
			propagateMark(current_region,isOuter,isOnIntersection,isInner);
		}

	} // for(;it!=ite;it++){
	}

	// classify all regions marked by the isInner mark
	{
		IGMesh::region_iterator it  = this->mesh_.regions_begin();

		for(;!it.isDone();it.next()) {
			if ((this->mesh_.isMarked(it.value(),isInner))) {
				(*volumeClassification)[(it.value()).getID()] = AVol;
			}
		}
	}

	std::cout<<"INNER cells are resolved!!!"<<std::endl;
//#ifdef GMDS_DEBUG_INSERTION_FILES
	this->exportMeshVTK("debug_inner_cells.mli",gmds::N|gmds::R);
//#endif

	this->insideRatio_ = this->mesh_.newVariable<double>(GMDS_REGION,"MeshInsertDetailInOut_insideRatio");

	// SELECTION OF INTERSECTED CELLS TO BE IN OR OUT
	{
		int iPercent = 0;

		for(unsigned int ir=0;ir<intersectedRegions.size();ir++)
		{
			// display progression
			if(ir == ((intersectedRegions.size()/10)*iPercent)) {
				std::cout<<"intersected cell traversed: "<<ir<<" of "<<intersectedRegions.size()<<" "<<iPercent*10<<"%"<<std::endl;
				iPercent++;
			}

			Region current_region = intersectedRegions[ir];
			std::vector<gmds::math::Triangle> current_triangles = in_out_triangles[current_region.getID()];

			// check whether the region is mainly inside or outside the model
			double insideRatio;
			if(this->service_.isMainlyInsideMC(current_triangles,current_region,insideRatio))
			{
				(*volumeClassification)[current_region.getID()] = AVol;
			}
			(*insideRatio_)[current_region.getID()] = insideRatio;
		}
	}
	std::cout<<"Classification of intersected cells is done!!!"<<std::endl;

	this->exportMeshVTK("debug_isMainlyInside_cells.mli",gmds::N|gmds::R);

	// count the number of selected cells
	{
		unsigned int nbSelectedCells = 0;

		IGMesh::region_iterator it  = this->mesh_.regions_begin();

		for(;!it.isDone();it.next()) {
			if((*volumeClassification)[(it.value()).getID()] == AVol) {
				nbSelectedCells++;
			}
		}

		if(0 == nbSelectedCells) {
			throw GMDSException("MeshInsertDetailInOut::associateGeometryClassificationCells no cell was selected "
					"to be part of the volume");
		}
	}
//
//#ifdef GMDS_DEBUG_INSERTION_FILES
//	this->exportVTK("debug_manifold_cells_before.mli");
//#endif
//
	unsigned int nbManifoldEdges = countManifoldEdges(AVol);
	std::cout<<"nbManifoldEdges "<<nbManifoldEdges<<std::endl;
//
//	// correct manifoldness
//	if(!correctManifoldness(AVol)) {
//		throw GMDSException("MeshInsertDetailInOut::associateGeometryClassificationCells "
//				"failed to correct manifoldness.");
//	}
//
//#ifdef GMDS_DEBUG_INSERTION_FILES
//	this->exportVTK("debug_manifold_cells_after.mli");
//#endif
//
	this->mesh_.unmarkAll<Region>(isInner);
	this->mesh_.unmarkAll<Region>(isOuter);
	this->mesh_.unmarkAll<Region>(isOnIntersection);
	this->mesh_.freeMark<Region>(isInner);
	this->mesh_.freeMark<Region>(isOuter);
	this->mesh_.freeMark<Region>(isOnIntersection);
//
//	// mark which entities are on the boundary
	buildBoundaryInfo(AVol);
//
//#ifdef GMDS_DEBUG_INSERTION_FILES
	this->exportMeshVTK("debug_in_out_cells.mli",gmds::N|gmds::E|gmds::F|gmds::R);
//#endif
}
/*----------------------------------------------------------------------------*/
void
MeshInsertDetailInOut::propagateMark(
		Region ARegion,
		const int& AMark, const int& AMarkBarrier, const int& AMarkFoe)
{
	std::vector<Face> faces = ARegion.get<Face>();

	for(size_t iFace=0; iFace<faces.size(); iFace++) {
		std::vector<Region> regions = faces[iFace].get<Region>();

		for(size_t iRegion=0; iRegion<regions.size(); iRegion++) {
			if(this->mesh_.isMarked(regions[iRegion],AMark))
				continue;
			if(this->mesh_.isMarked(regions[iRegion],AMarkBarrier))
				continue;
			if(this->mesh_.isMarked(regions[iRegion],AMarkFoe)) {
				throw GMDSException("Error a Foe Mark was encountered during mark propagation!");
				return;
			}

			this->mesh_.mark(regions[iRegion],AMark);
			propagateMark(regions[iRegion],AMark,AMarkBarrier,AMarkFoe);

		}
	}
}
/*----------------------------------------------------------------------------*/
void
MeshInsertDetailInOut::associateGeometryClassificationPoints(gmds::geom::GeomVolume* AVol)
{
	// associate Points
	Variable<geom::GeomEntity* >* pointClassification = this->mesh_.getGeometricClassification(0);

	// find nearest node of the grid
	std::vector<gmds::geom::GeomPoint*> points;
	AVol->get(points);

	std::cout<<"number of sharp points = "<<points.size()<<std::endl;

	for(size_t iPoint=0; iPoint<points.size(); iPoint++)
	{
		associateGeometryClassificationPoint(points[iPoint],AVol);
	}

	this->exportMeshVTK("debug_points.mli",gmds::N|gmds::E|gmds::F|gmds::R);
}
/*----------------------------------------------------------------------------*/
void
MeshInsertDetailInOut::associateGeometryClassificationPoint(
		gmds::geom::GeomPoint* APoint,
		gmds::geom::GeomVolume* AVol)
{
	// associate Points
	Variable<geom::GeomEntity* >* pointClassification = this->mesh_.getGeometricClassification(0);

	// find nearest nodes of the grid
	// We keep only the NB_NODES_KEPT nearest nodes

	gmds::math::Point point = APoint->getPoint();

	// selecting possible candidates
	std::vector<double> distCandidates(MESH_INSERT_DETAIL_INOUT_NB_NODES_KEPT,HUGE_VALF);
	std::vector<double> distCandidates_tmp(MESH_INSERT_DETAIL_INOUT_NB_NODES_KEPT,HUGE_VALF);
	std::vector<gmds::Node> nodesCandidates(MESH_INSERT_DETAIL_INOUT_NB_NODES_KEPT,gmds::Node());
	std::vector<gmds::Node> nodesCandidates_tmp(MESH_INSERT_DETAIL_INOUT_NB_NODES_KEPT,gmds::Node());


	IGMesh::node_iterator it  = this->mesh_.nodes_begin();

	for(;!it.isDone();it.next()) {

		Node current_node = it.value();

		// node must be on the volume boundary
		if(this->mesh_.isMarked(current_node,this->markNodeIsOnBoundary_)) {

			// node must not be already associated to a point
			if((*pointClassification)[current_node.getID()] == NULL ||
				((*pointClassification)[current_node.getID()])->getDim()!=0	) {

				// if APoint is included in a model entity, a candidate
				// node must be part of the mesh entities associated to said entity.

///////////////POYOP need to do that
///////////////POYOP need to do that

				// we limit the search to nodes that have enough adjacent boundary edges
				std::vector<Edge> edges = current_node.get<Edge>();
				unsigned int nbBoundaryEdges = 0;
				for(unsigned int iEdge=0; iEdge<edges.size(); iEdge++) {
					if(this->mesh_.isMarked(edges[iEdge],this->markEdgeIsOnBoundary_)) {
						nbBoundaryEdges++;
					}
				}
				std::vector<gmds::geom::GeomCurve* > curves;
				APoint->get(curves);
				if(curves.size() > nbBoundaryEdges) {
					continue;
				}

				// check if current_node is among the NB_NODES_KEPT nearest nodes.
				// If it is, we insert it in nodes_array.
				double dist_tmp = point.distance(current_node.getPoint());

				if(distCandidates[MESH_INSERT_DETAIL_INOUT_NB_NODES_KEPT-1]>dist_tmp) {

					// look for the index where current_node will be inserted
					distCandidates_tmp = distCandidates;
					nodesCandidates_tmp = nodesCandidates;

					unsigned int iNode = 0;
					for(;distCandidates[iNode]<dist_tmp; iNode++) {
					}
					distCandidates[iNode] = dist_tmp;
					nodesCandidates[iNode] = current_node;
					iNode++;
					for(; iNode<MESH_INSERT_DETAIL_INOUT_NB_NODES_KEPT; iNode++) {
						distCandidates[iNode] = distCandidates_tmp[iNode-1];
						nodesCandidates[iNode] = nodesCandidates_tmp[iNode-1];
					}
				}

			} else {
				// node can be associated to two points if one of these points
				// is included in the other.
				// In such a case we throw an exception because we do not know
				// what to do at the moment.
				if(APoint->isIncludedIn()) {
					if((*pointClassification)[current_node.getID()] == APoint->getIncludedGeomEntity()) {
						throw GMDSException(""
								"MeshInsertDetailInOut::associateGeometryClassificationPoint "
								"one model point is included in another");
					}
				}
			}
		}
	}

	// Mirror Mirror, who is the nearest of them all?
	// at the moment we select the nearest node among the nearest
	Node selectedNode = associateGeometryClassificationPointSelectBestNode(nodesCandidates,APoint,AVol);

	if(selectedNode.getID() == gmds::NullID)
		throw GMDSException("A GeomPoint did not find a corresponding grid node");

	(*pointClassification)[selectedNode.getID()] = APoint;
	geom2MeshNode_[APoint] = selectedNode;
}
/*----------------------------------------------------------------------------*/
Node MeshInsertDetailInOut::
associateGeometryClassificationPointSelectBestNode(
		std::vector<Node>& ANodes,
		gmds::geom::GeomPoint* APoint,
		gmds::geom::GeomVolume* AVol)
{
	// associate Points
	Variable<geom::GeomEntity* >* pointClassification = this->mesh_.getGeometricClassification(0);

	// compute cost for each candidate nodes
	std::vector<double> nodesCosts(ANodes.size(),HUGE_VALF);

	for(unsigned int iNode=0; iNode<ANodes.size(); iNode++) {

		gmds::Node current_node = ANodes[iNode];

		std::vector<Edge> edges_tmp;
		nodesCosts[iNode] = computeNodePointClassificationCost(current_node,APoint,AVol,edges_tmp);
	}

	int indexBestNode = -1;
	double maxNodeCost = -HUGE_VALF;
	for(unsigned int iNode=0; iNode<ANodes.size(); iNode++) {
		if(maxNodeCost<nodesCosts[iNode]) {
			maxNodeCost = nodesCosts[iNode];
			indexBestNode = iNode;
		}
	}

	if(-1 == indexBestNode) {
		throw GMDSException("MeshInsertDetailInOut::associateGeometryClassificationPointSelectBestNode "
				"could not find a valid node?");
	}

//	return ANodes[indexBestNode];
	return ANodes[0];
}
/*----------------------------------------------------------------------------*/
double
MeshInsertDetailInOut::computeNodePointClassificationCost
(
		gmds::Node ANode,
		gmds::geom::GeomPoint* APoint,
		gmds::geom::GeomVolume* AVol,
		std::vector<gmds::Edge>& AEdges)
{
	std::vector<gmds::geom::GeomCurve*> curves;
	APoint->getOrderedDirect(curves);

	std::vector<gmds::math::Vector> curvesVector;
	curvesVector.resize(curves.size());

	for(unsigned int iCurve=0; iCurve<curves.size(); iCurve++)
	{
		curves[iCurve]->computeVector(*APoint,curvesVector[iCurve]);
	}

	// we will get the regions adjacent to the node
	std::vector<gmds::Region> regions;
	std::vector<gmds::Face> faces;
	std::vector<gmds::Edge> edges;

	regions = ANode.get<Region>();
	faces = ANode.get<Face>();
	edges = ANode.get<Edge>();

	// we keep only edges marked as boundary
	std::vector<Edge> edge_tmp;
	{
		for(unsigned int iEdge=0; iEdge<edges.size(); iEdge++) {
			if(this->mesh_.isMarked(edges[iEdge],this->markEdgeIsOnBoundary_)) {
				edge_tmp.push_back(edges[iEdge]);
			}
		}
	}

	// if there are less edges kept than curves, change node candidate
	if(edge_tmp.size() < curves.size()) {
		return -HUGE_VALF;
	}

	// we will order the edges in the direct order, direct related to the inside regions
	getOrderedDirectEdges(edge_tmp,ANode,APoint,AVol);

	std::vector<gmds::math::Vector> edgesVector;
	edgesVector.resize(edge_tmp.size());
	int permutEdges[edge_tmp.size()];

	std::cout<<"egdesVector "<<std::endl;

	for(unsigned int iEdge=0; iEdge<edge_tmp.size(); iEdge++)
	{
		// initialize the array that will be permuted
		permutEdges[iEdge] = iEdge;

		// compute the edges vectors
		std::vector<gmds::Node> nodes;
		nodes = edge_tmp[iEdge].get<Node>();
		if(ANode == nodes[0])
		{
			edgesVector[iEdge] = gmds::math::Vector (nodes[0].getPoint(),nodes[1].getPoint());
		}
		else
		{
			edgesVector[iEdge] = gmds::math::Vector (nodes[1].getPoint(),nodes[0].getPoint());
		}

		edgesVector[iEdge].normalize();
//		std::cout<<"edge_tmp[iEdge]->getID() "<<edge_tmp[iEdge].getID()<<std::endl;
//		std::cout<<edgesVector[iEdge]<<std::endl;
	}

	std::sort(permutEdges,permutEdges+edge_tmp.size());

	double Uij_max = 0.;
	int permutEdges_max[curves.size()];
	bool permutFound = false;

	do{
		// we check that this permutation is ordered increasingly.
		bool oneRound = false;
		bool orderedIncrease = true;
		for(unsigned int iCurve=0; iCurve<curves.size(); iCurve++) {
			if(oneRound && (permutEdges[(iCurve+1)%curves.size()] > permutEdges[0])) {
				orderedIncrease = false;
				break;
			}
			if(permutEdges[(iCurve+1)%curves.size()] < permutEdges[(iCurve)%curves.size()]) {
				if(oneRound) {
					orderedIncrease = false;
					break;
				} else {
					if(permutEdges[(iCurve+1)%curves.size()] > permutEdges[0]) {
						orderedIncrease = false;
						break;
					} else {
						oneRound = true;
					}
				}
			}
		}
		if(!orderedIncrease) {
			continue;
		}

		double Uij = 0.;

		for(unsigned int iCurve=0; iCurve<curves.size(); iCurve++)
		{
			Uij += 1. + curvesVector[iCurve].dot(edgesVector[permutEdges[iCurve]]);
		}
		Uij /=2.;

//			std::cout<<"permutEdges ";
//			for(unsigned int iEdge=0; iEdge<edge_tmp.size(); iEdge++) {
//				std::cout<<permutEdges[iEdge]<<" ";
//			}
//			std::cout<<endl;
//			std::cout<<"Uij Uij_max "<<Uij<<" "<<Uij_max<<std::endl;

		if(Uij>=Uij_max)
		{
			Uij_max = Uij;
			for(unsigned int iCurve=0; iCurve<curves.size(); iCurve++)
			{
				permutEdges_max[iCurve] = permutEdges[iCurve];
			}
//			isInRegionMap_max = isInRegionMap;
//			isInRegion_max = isInRegion;
//			edges_max = edge_tmp;
			permutFound = true;
		}
	}
	while(std::next_permutation(permutEdges,permutEdges+edge_tmp.size())); // do{

	if(!permutFound) {
		throw GMDSException("TBase MeshInsertDetailInOut::computeNodePointClassificationCost "
				"no suitable permutation was found.");
	}

	AEdges.clear();
	AEdges.resize(curves.size(),Edge());
	for(unsigned int iEdge=0; iEdge<curves.size(); iEdge++) {
		AEdges[iEdge] = edge_tmp[permutEdges_max[iEdge]];
	}

	return Uij_max;
}
/*----------------------------------------------------------------------------*/
void
MeshInsertDetailInOut::getOrderedDirectEdges
(
		std::vector<Edge>& AEdges,
		Node ANode,
		gmds::geom::GeomPoint* APoint,
		gmds::geom::GeomVolume* AVol)
{
	std::vector<gmds::geom::GeomCurve*> curves;
	APoint->getOrderedDirect(curves);

	std::vector<gmds::math::Vector> curvesVector;
	curvesVector.resize(curves.size());

	for(unsigned int iCurve=0; iCurve<curves.size(); iCurve++)
	{
		curves[iCurve]->computeVector(*APoint,curvesVector[iCurve]);
	}

	// we will get the regions adjacent to the node
	std::vector<gmds::Region> regions;
	std::vector<gmds::Face> faces;
	std::vector<gmds::Edge> edges;
	regions = ANode.get<Region>();
	faces = ANode.get<Face>();
	edges = ANode.get<Edge>();

	// we keep only edges marked as boundary
	std::vector<Edge> edge_tmp;
	{
		for(unsigned int iEdge=0; iEdge<edges.size(); iEdge++) {
			if(this->mesh_.isMarked(edges[iEdge],this->markEdgeIsOnBoundary_)) {
				edge_tmp.push_back(edges[iEdge]);
			}
		}
	}

	// we will order the edges in the direct order, direct related to the inside regions
	{
		std::vector<Edge> edges_ordered;

		Edge current_edge = edge_tmp[0];

		// keep only the two faces on the boundary
		std::vector<Face> faces = current_edge.get<Face>();
		for(std::vector<Face>::iterator itf  = faces.begin(); itf!=faces.end();) {
			if(!this->mesh_.isMarked(*itf,this->markFaceIsOnBoundary_)) {
				itf = faces.erase(itf);
			} else {
				itf++;
			}
		}
		if(faces.size() != 2)
			throw GMDSException("there should be two faces here. first.");

		// select which face is the one in the direct order
		Face current_face = faces[0];

		std::vector<Node> nodes = current_face.get<Node>();
		std::vector<Region> regions = current_face.get<Region>();

		Region current_region;
		if(this->mesh_.isMarked(regions[0],this->markRegionIsInside_)) {
			current_region = regions[0];
		} else {
			current_region = regions[1];
		}

		// find which face of the region corresponds to current_face.
		bool faceDirect;

		std::vector<Node> nodes_region = current_region.get<Node>();
		std::vector<Node> nodes_face = current_face.get<Node>();

		std::vector<std::vector<Node> > orderedNodesFaces = current_region.getOrderedNodesFaces();

		for(unsigned int iFace=0; iFace<orderedNodesFaces.size(); iFace++) {

			// compare the size of the faces
			if(orderedNodesFaces[iFace].size() != nodes_face.size()) {
				continue;
			}

			unsigned int nbCommonNodes = 0;
			for(unsigned int iNode1=0; iNode1<nodes_face.size(); iNode1++) {
				for(unsigned int iNode2=0; iNode2<orderedNodesFaces[iFace].size(); iNode2++) {
					if(nodes_face[iNode1] == orderedNodesFaces[iFace][iNode2]) {
						nbCommonNodes++;
					}
				}
			}

			if(nbCommonNodes == nodes_face.size()) {

				for(unsigned int iNode=0; iNode<nodes_face.size(); iNode++)
				{
					if(nodes_face[iNode] == orderedNodesFaces[iFace][0]) {
						if(nodes_face[(iNode+1)%(nodes_face.size())] == orderedNodesFaces[iFace][1]) {
							faceDirect = true;
							break;
						} else {
							faceDirect = false;
							break;
						}
					}
				}
			}

		} // for(unsigned int iFace=0; iFace<orderedNodesFaces.size(); iFace++) {

		// direct order is face current_face or the other one?
		Node node_other;
		std::vector<Node> nodes_edge = current_edge.get<Node>();
		if(nodes_edge[0] == ANode) {
			node_other = nodes_edge[1];
		} else {
			node_other = nodes_edge[0];
		}

		Face face_next = Face();

		for(unsigned int iNode=0; iNode<nodes_face.size(); iNode++) {
			if(nodes_face[iNode] == ANode) {
				if(faceDirect) {
					if(nodes_face[(iNode+1)%nodes_face.size()] == node_other) {
						face_next = current_face;
						break;
					} else {
						face_next = faces[1];
						break;
					}
				} else {
					if(nodes_face[(iNode+nodes_face.size()-1)%nodes_face.size()] == node_other) {
						face_next = current_face;
						break;
					} else {
						face_next = faces[1];
						break;
					}
				}
			}
		}

		// now that we have a first edge and a first face in direct order,
		// propagate across the edges.
		Edge starting_edge = current_edge;
		Edge next_edge = current_edge;
		Face next_face = face_next;

		do {
			current_edge = next_edge;
			current_face = next_face;

			edges_ordered.push_back(current_edge);

			std::vector<Edge> edges_face = current_face.get<Edge>();

			for(unsigned int iEdge=0; iEdge<edges_face.size(); iEdge++) {
				if(edges_face[iEdge].getID() != current_edge.getID()) {
					std::vector<Node> nodes_egde = edges_face[iEdge].get<Node>();
					if((nodes_egde[0] == ANode) || (nodes_egde[1] == ANode)) {
						next_edge = edges_face[iEdge];
						break;
					}
				}
			}
			// keep only the two faces on the boundary
			std::vector<Face> faces = next_edge.get<Face>();
			for(std::vector<Face>::iterator itf  = faces.begin(); itf!=faces.end();) {
				if(!this->mesh_.isMarked(*itf,this->markFaceIsOnBoundary_)) {
					itf = faces.erase(itf);
				} else {
					itf++;
				}
			}
			if(faces.size() != 2)
				throw GMDSException("MeshInsertDetailInOut::getOrderedDirectEdges there should be two faces here. second.");

			if(faces[0].getID() == current_face.getID()) {
				next_face = faces[1];
			} else {
				next_face = faces[0];
			}

		}
		while(next_edge.getID() != starting_edge.getID());

		if(edge_tmp.size() != edges_ordered.size()) {
			throw GMDSException("MeshInsertDetailInOut::getOrderedDirectEdges we missed or added an edge, this is baaaad!");
		}

		AEdges.clear();
		for(unsigned int iEdge=0; iEdge<edges_ordered.size(); iEdge++) {
			AEdges.push_back(edges_ordered[iEdge]);
		}
	}

}
/*----------------------------------------------------------------------------*/
void
MeshInsertDetailInOut::buildBoundaryInfo(
		gmds::geom::GeomVolume* AVol)
{
	Variable<geom::GeomEntity* >* volumeClassification = this->mesh_.getGeometricClassification(3);

	this->markRegionIsInside_ = this->mesh_.getNewMark<Region>();
	this->markFaceIsOnBoundary_ = this->mesh_.getNewMark<Face>();
	this->markEdgeIsOnBoundary_ = this->mesh_.getNewMark<Edge>();
	this->markNodeIsOnBoundary_ = this->mesh_.getNewMark<Node>();
	this->mesh_.unmarkAll<Region>(this->markRegionIsInside_);
	this->mesh_.unmarkAll<Face>(this->markFaceIsOnBoundary_);
	this->mesh_.unmarkAll<Edge>(this->markEdgeIsOnBoundary_);
	this->mesh_.unmarkAll<Node>(this->markNodeIsOnBoundary_);

	// marking regions
	IGMesh::region_iterator itr  = this->mesh_.regions_begin();

	for(;!itr.isDone();itr.next()){
		Region current_region = itr.value();

		if((*volumeClassification)[current_region.getID()] == AVol) {
			this->mesh_.mark(current_region,this->markRegionIsInside_);
		}
	}

	// marking faces
	IGMesh::face_iterator itf  = this->mesh_.faces_begin();

	for(;!itf.isDone();itf.next()){
		Face current_face = itf.value();

		std::vector<Region> regions = current_face.get<Region>();

		if(regions.size() == 2) {

			if((this->mesh_.isMarked(regions[0],this->markRegionIsInside_) && !this->mesh_.isMarked(regions[1],this->markRegionIsInside_))
					|| (this->mesh_.isMarked(regions[1],this->markRegionIsInside_) && !this->mesh_.isMarked(regions[0],this->markRegionIsInside_))) {

				this->mesh_.mark(current_face,markFaceIsOnBoundary_);
			}
		} else {
			if(regions.size() != 1) {
				throw GMDSException("MeshInsertDetailInOut::buildBoundaryInfo F2R incoherent.");
			}

			// in this case the face is on the boundary of the domain.
			if(this->mesh_.isMarked(regions[0],this->markRegionIsInside_)) {

				this->mesh_.mark(current_face,markFaceIsOnBoundary_);
			}
		} // if(regions.size() == 2) {

	} // for(;itf!=itfe;itf++){

	// created edges and adjacencies
	buildBoundaryEntities();

	// marking edges and nodes
	itf  = this->mesh_.faces_begin();

	for(;!itf.isDone();itf.next()){
		Face current_face = itf.value();

		if(this->mesh_.isMarked(current_face,markFaceIsOnBoundary_)) {

			std::vector<Edge> edges = current_face.get<Edge>();
			for(unsigned int iEdge=0; iEdge<edges.size(); iEdge++) {
				this->mesh_.mark(edges[iEdge],markEdgeIsOnBoundary_);
			}

			std::vector<Node> nodes = current_face.get<Node>();
			for(unsigned int iNode=0; iNode<nodes.size(); iNode++) {
				this->mesh_.mark(nodes[iNode],markNodeIsOnBoundary_);
			}
		}
	}
}
/*----------------------------------------------------------------------------*/
void
MeshInsertDetailInOut::buildBoundaryEntities()
{
	// change model
    gmds::MeshModel mod = gmds::DIM3|gmds::N|gmds::F|gmds::R|gmds::F2N|gmds::R2N|gmds::R2F|gmds::F2R|gmds::N2F|gmds::N2R|
    		gmds::E|gmds::E2N|gmds::F2E|gmds::E2F|gmds::N2E;
    this->mesh_.changeModel(mod,false);

	std::map<FakeEdge,Edge> fakeEdgesMap;

	IGMesh::face_iterator itf  = this->mesh_.faces_begin();

	for(;!itf.isDone();itf.next()) {
		Face current_face = itf.value();

		if(this->mesh_.isMarked(current_face,markFaceIsOnBoundary_)) {

			std::vector<Node> nodes = current_face.get<Node>();

			for(int iNode=0; iNode<nodes.size(); iNode++) {
				this->mesh_.mark(nodes[iNode],markNodeIsOnBoundary_);
			}

			// create the edges or retrieve the previously created edges
			for(int iNode=0; iNode<nodes.size(); iNode++) {
				Node n1 = nodes[iNode];
				Node n2 = nodes[(iNode+1)%nodes.size()];

				FakeEdge fakeEdge(n1.getID(),n2.getID());
				Edge newEdge;
				if(fakeEdgesMap.find(fakeEdge) == fakeEdgesMap.end()) {
					newEdge = this->mesh_.newEdge(n1,n2);
					fakeEdgesMap[fakeEdge] = newEdge;
					this->mesh_.mark(newEdge,markEdgeIsOnBoundary_);
					// filling the N2E adjacency
					n1.add<Edge>(newEdge);
					n2.add<Edge>(newEdge);
				} else {
					newEdge = fakeEdgesMap[fakeEdge];
				}

				// filling the F2E E2F adjacency
				current_face.add<Edge>(newEdge);
				newEdge.add<Face>(current_face);
			}
		}
	}
}
/*----------------------------------------------------------------------------*/
void
MeshInsertDetailInOut::associateGeometryClassificationCurvesAtPoints(gmds::geom::GeomVolume* AVol)
{
	// loop on Geom Points
	std::vector<gmds::geom::GeomPoint*> points;
	AVol->get(points);

	for(unsigned int iPoint=0; iPoint<points.size(); iPoint++)
	{
		gmds::Node current_node = geom2MeshNode_[points[iPoint]];

		std::cout<<"matching edges curves at node "<<current_node.getID()<<std::endl;

		std::vector<gmds::geom::GeomCurve*> curves;
		points[iPoint]->getOrderedDirect(curves);

		geomCurvesAtPoints2MeshEdge_[points[iPoint]].resize(curves.size());

		std::vector<Edge> edges_tmp;
		computeNodePointClassificationCost(current_node,points[iPoint],AVol,edges_tmp);

		Variable<geom::GeomEntity* >* edgeClassification = this->mesh_.getGeometricClassification(1);

		for(unsigned int iCurve=0; iCurve<curves.size(); iCurve++)
		{
			geomCurvesAtPoints2MeshEdge_[points[iPoint]][iCurve] = edges_tmp[iCurve];
			(*edgeClassification)[edges_tmp[iCurve].getID()] = curves[iCurve];
		}

	} // for(unsigned int iPoint=0; iPoint<points.size(); iPoint++)

}
/*----------------------------------------------------------------------------*/
void
MeshInsertDetailInOut::propagateVisited(gmds::Region& ARegion, std::map<Region,bool>& AIsInRegionMap, std::map<Region,bool>& AWasVisited)
{
	if(!AWasVisited[ARegion]) {
		if(AIsInRegionMap[ARegion]) {
			AWasVisited[ARegion] = true;

			std::vector<Face> faces = ARegion.get<Face>();

			for(unsigned int iFace=0; iFace<faces.size(); iFace++) {
				std::vector<Region> regions = faces[iFace].get<Region>();

				if(regions.size() == 1) {
					continue;
				}

				Region next_region;

				if(regions[0] == ARegion) {
					next_region = regions[1];
				} else {
					next_region = regions[0];
				}

				if(AIsInRegionMap.find(next_region) != AIsInRegionMap.end()) {
					if(AIsInRegionMap[next_region]) {
						propagateVisited(next_region,AIsInRegionMap,AWasVisited);
					}
				}
			}
		}
	}

}
/*----------------------------------------------------------------------------*/
void
MeshInsertDetailInOut::propagateVisited(gmds::Face& AFace, std::map<Face,bool>& AIsInFaceMap, std::map<Face,bool>& AWasVisited)
{
	if(!AWasVisited[AFace]) {
		if(AIsInFaceMap[AFace]) {
			AWasVisited[AFace] = true;

			std::vector<Edge> edges = AFace.get<Edge>();

			for(unsigned int iEdge=0; iEdge<edges.size(); iEdge++) {
				std::vector<Face> faces = edges[iEdge].get<Face>();

				Face next_face;

				// here it works because such an edge has only two faces registered as "inside"
				for(unsigned int iFace=0; iFace<faces.size(); iFace++) {
					if(AIsInFaceMap.find(faces[iFace]) != AIsInFaceMap.end()) {
						if(faces[iFace] != AFace) {
							if(AIsInFaceMap[faces[iFace]]) {
								next_face = faces[iFace];
								break;
							}
						}
					}
				}

				propagateVisited(next_face,AIsInFaceMap,AWasVisited);
			}
		}
	}

}
/*----------------------------------------------------------------------------*/
void
MeshInsertDetailInOut::associateGeometryClassificationCurves(gmds::geom::GeomVolume* AVol)
{
	std::vector<gmds::geom::GeomCurve*> curves;
	AVol->get(curves);

	// this mark will be used to avoid non-permissible configurations.
	// For example, (self-)intersecting curves
	int markNodeConstraint = this->mesh_.getNewMark<Node>();
	this->mesh_.unmarkAll<Node>(markNodeConstraint);

	std::vector<int*> permutPool;
	int permutCurvesFirst[curves.size()];
	for(unsigned int iCurve=0; iCurve<curves.size(); iCurve++)
	{
		permutCurvesFirst[iCurve] = iCurve;
	}
	permutPool.push_back(permutCurvesFirst);

	std::map<gmds::geom::GeomCurve*, bool> alreadyShuffledCurves;

	bool fullAssociationSuccess = true;
	int nbTry = 0;

	do {
		// we mark the nodes associated to points
		std::vector<gmds::geom::GeomPoint* > points;
		AVol->get(points);
		for(unsigned int iPoint=0; iPoint<points.size(); iPoint++)
		{
			this->mesh_.mark(geom2MeshNode_[points[iPoint]],markNodeConstraint);

			// add in restricted nodes the nodes of the edges associated at points
			for(unsigned int iCurveAtPoint=0; iCurveAtPoint<geomCurvesAtPoints2MeshEdge_[points[iPoint]].size(); iCurveAtPoint++) {
				Edge current_edge = geomCurvesAtPoints2MeshEdge_[points[iPoint]][iCurveAtPoint];

				std::vector<Node> nodes = current_edge.get<Node>();
				this->mesh_.mark(nodes[0],markNodeConstraint);
				this->mesh_.mark(nodes[1],markNodeConstraint);

				this->mesh_.mark(current_edge,markNodeConstraint);
			}
		}

		fullAssociationSuccess = true;
		double distanceTot = 0.;

		std::cout<<"permutPool.size() "<<permutPool.size()<<std::endl;
		for(unsigned int iCurve=0; iCurve<curves.size(); iCurve++) {
			std::cout<<(permutPool.back())[iCurve]<<" ";
		}
		std::cout<<std::endl;

		for(size_t iCurve=0; iCurve<curves.size(); iCurve++)
		{
			gmds::geom::GeomCurve* current_curve =  curves[(permutPool.back())[iCurve]];

			double distance = 0.;
			bool associationSuccess = associateGeometryClassificationEdge(
					distance,
					*current_curve,
					markNodeConstraint);
#ifdef GMDS_DEBUG_INSERTION_FILES
//			this->exportVTK("debug_acceptable_edges.mli");
#endif
//			if(nbTry<4) {
//				associationSuccess = false;
//				nbTry++;
//			}

			if(!associationSuccess) {
				std::cout<<"Fail to associate curve "<<current_curve->getName()<<std::endl;
				this->exportMeshVTK("debug_edges_beforefail",gmds::N|gmds::E|gmds::F|gmds::R);
				this->exportCurvesEdgesVTK("debug_edges_beforefail_curves");
				exit(-1);
				removeGeometryClassificationCurves(AVol);
				this->mesh_.unmarkAll<Node>(markNodeConstraint);
				fullAssociationSuccess = false;

				permutPool.pop_back();

				// check that current_curve has not already been scheduled for reshuffling
				if(alreadyShuffledCurves.find(current_curve) == alreadyShuffledCurves.end()) {

					// we get the curves adjacent to the end-points of current_curve
					std::vector<gmds::geom::GeomPoint*> points;
					current_curve->get(points);

					std::vector<gmds::geom::GeomCurve*> curves_pb;
					std::vector<gmds::geom::GeomCurve*> curves2;

					points[0]->get(curves_pb);
					if(!current_curve->isALoop()) {
						points[1]->get(curves2);
					}

					for(unsigned int iCurve2=0; iCurve2<curves2.size(); iCurve2++) {
						// check that the a curve is not added twice in curves_pb
						bool found = false;
						for(unsigned int iCurve_pb=0; iCurve_pb<curves_pb.size(); iCurve_pb++) {
							if(curves_pb[iCurve_pb] == curves2[iCurve2]) {
								found = true;
								break;
							}
						}
						if(!found) {
							curves_pb.push_back(curves2[iCurve2]);
						}
					}

					// now for each possible permutation between those curves
					// we add one entry to permutPool
					unsigned int origCurveIndex[curves_pb.size()];

					for(unsigned int iCurve_pb=0; iCurve_pb<curves_pb.size(); iCurve_pb++) {
						for(unsigned int iCurveOrig=0; iCurveOrig<curves.size(); iCurveOrig++) {
							if(curves_pb[iCurve_pb] == curves[iCurveOrig]) {
								origCurveIndex[iCurve_pb] = iCurveOrig;
							}
						}
					}

					int permutCurves_pb[curves_pb.size()];

					for(unsigned int iCurve_pb=0; iCurve_pb<curves_pb.size(); iCurve_pb++) {
						permutCurves_pb[iCurve_pb] = origCurveIndex[iCurve_pb];
					}
					std::sort(permutCurves_pb,permutCurves_pb+curves_pb.size());

					int nbPermut = 0;
					do {
						nbPermut++;
						if(nbPermut==1) {
							continue;
						}

						int* permutCurvesNew = new int[curves.size()];
						for(unsigned int iCurve1=0; iCurve1<curves.size(); iCurve1++) {
							permutCurvesNew[iCurve1] = iCurve1;
						}

						for(unsigned int iCurve_pb=0; iCurve_pb<curves_pb.size(); iCurve_pb++) {
							permutCurvesNew[origCurveIndex[iCurve_pb]] = permutCurves_pb[iCurve_pb];
						}

						permutPool.push_back(permutCurvesNew);
					}
					while(std::next_permutation(permutCurves_pb,permutCurves_pb+curves_pb.size()));

					alreadyShuffledCurves[current_curve] = true;

				} // if(alreadyShuffledCurves.find(current_curve) != alreadyShuffledCurves.end()) {

				break;

			} // if(!associationSuccess) {

			distanceTot += distance;

		} // for(size_t iCurve=0; iCurve<curves.size(); iCurve++)

		if(fullAssociationSuccess) {
			break;
		}
	}
	while(!permutPool.empty()); // do {

	if(!fullAssociationSuccess) {
		throw GMDSException("MeshInsertDetailInOut::associateGeometryClassificationCurves "
				"no valid edges curves association found.");
	}

//	int permutCurves[curves.size()];
//	for(unsigned int iCurve=0; iCurve<curves.size(); iCurve++)
//	{
//		permutCurves[iCurve] = iCurve;
//	}
//
//	std::sort(permutCurves,permutCurves+curves.size());
//	bool fullAssociationSuccess = true;
//	double distanceTotMin = HUGE_VALF;
//	int permutCurvesWithMinDistance[curves.size()];
//	unsigned int tryNum = 0;
//
//	do {
//		// we mark the nodes associated to points
//		std::vector<gmds::geom::GeomPoint<TBase>* > points;
//		AVol->get(points);
//		for(unsigned int iPoint=0; iPoint<points.size(); iPoint++)
//		{
//			this->mesh_.mark(geom2MeshNode_[points[iPoint]],markNodeConstraint);
//
//			// add in restricted nodes the nodes of the edges associated at points
//			for(unsigned int iCurveAtPoint=0; iCurveAtPoint<geomCurvesAtPoints2MeshEdge_[points[iPoint]].size(); iCurveAtPoint++) {
//				Edge* current_edge = geomCurvesAtPoints2MeshEdge_[points[iPoint]][iCurveAtPoint];
//
//				std::vector<Node*> nodes = current_edge->getNodes();
//				this->mesh_.mark(nodes[0],markNodeConstraint);
//				this->mesh_.mark(nodes[1],markNodeConstraint);
//
//				this->mesh_.mark(current_edge,markNodeConstraint);
//			}
//		}
//
//		std::cout<<"permutCurves "<<tryNum<<std::endl;
//		for(unsigned int iCurve=0; iCurve<curves.size(); iCurve++) {
//			std::cout<<permutCurves[iCurve]<<" ";
//		}
//		std::cout<<std::endl;
//
//		fullAssociationSuccess = true;
//		double distanceTot = 0.;
//
//		for(size_t iCurve=0; iCurve<curves.size(); iCurve++)
//		{
//			double distance = 0.;
//			bool associationSuccess = associateGeometryClassificationEdge(
//					distance,
//					*(curves[permutCurves[iCurve]]),
//					markNodeConstraint);
//#ifdef GMDS_DEBUG_INSERTION_FILES
////			this->exportVTK("debug_acceptable_edges.mli");
//#endif
//			if(!associationSuccess) {
//				removeGeometryClassificationCurves(AVol);
//				this->mesh_.unmarkAll(markNodeConstraint);
//				fullAssociationSuccess = false;
//				break;
//			}
//
//			distanceTot += distance;
//		}
//
//		std::cout<<"fullAssociationSuccess "<<fullAssociationSuccess<<std::endl;
//		std::cout<<"distanceTot "<<distanceTot<<" distanceTotMin "<<distanceTotMin<<std::endl;
//
//		if(fullAssociationSuccess && (distanceTot < distanceTotMin)) {
//			distanceTotMin = distanceTot;
//			for(unsigned int iCurve=0; iCurve<curves.size(); iCurve++) {
//				permutCurvesWithMinDistance[iCurve] = permutCurves[iCurve];
//			}
//
//			this->mesh_.unmarkAll(markNodeConstraint);
//			break;
//		}
//
//		this->mesh_.unmarkAll(markNodeConstraint);
//		tryNum++;
//	}
//	while(std::next_permutation(permutCurves,permutCurves+curves.size()));// && !fullAssociationSuccess); // do {
//
//	if(!fullAssociationSuccess) {
//		throw GMDSException("MeshInsertDetailInOut::associateGeometryClassificationCurves "
//				"no valid edges curves association found.");
//	}
//
//	removeGeometryClassificationCurves(AVol);
//
//	// we mark the nodes associated to points
//	std::vector<gmds::geom::GeomPoint<TBase>* > points;
//	AVol->get(points);
//	for(unsigned int iPoint=0; iPoint<points.size(); iPoint++)
//	{
//		this->mesh_.mark(geom2MeshNode_[points[iPoint]],markNodeConstraint);
//
//		// add in restricted nodes the nodes of the edges associated at points
//		for(unsigned int iCurveAtPoint=0; iCurveAtPoint<geomCurvesAtPoints2MeshEdge_[points[iPoint]].size(); iCurveAtPoint++) {
//			Edge* current_edge = geomCurvesAtPoints2MeshEdge_[points[iPoint]][iCurveAtPoint];
//
//			std::vector<Node*> nodes = current_edge->getNodes();
//			this->mesh_.mark(nodes[0],markNodeConstraint);
//			this->mesh_.mark(nodes[1],markNodeConstraint);
//
//			this->mesh_.mark(current_edge,markNodeConstraint);
//		}
//	}
//
//	for(size_t iCurve=0; iCurve<curves.size(); iCurve++) {
//		double distance = 0.;
//		bool associationSuccess = associateGeometryClassificationEdge(
//				distance,
//				*(curves[permutCurvesWithMinDistance[iCurve]]),
//				markNodeConstraint);
//		if(!associationSuccess) {
//			throw GMDSException("MeshInsertDetailInOut::associateGeometryClassificationCurves "
//					"problem of non-valid path");
//		}
//	}

//	for(size_t iCurve=0; iCurve<curves.size(); iCurve++)
//	{
//		std::cout<<"matching curve "<<iCurve<<std::endl;
//		bool associationSuccess = associateGeometryClassificationEdge(*(curves[iCurve]),markNodeConstraint);
//		if(!associationSuccess) {
//			removeGeometryClassificationEdge();
//		}
//#ifdef GMDS_DEBUG_INSERTION_FILES
//		this->exportVTK("debug_acceptable_edges.mli");
//#endif
//	}

	this->mesh_.unmarkAll<Node>(markNodeConstraint);
	this->mesh_.freeMark<Node>(markNodeConstraint);
}
/*----------------------------------------------------------------------------*/
bool
MeshInsertDetailInOut::associateGeometryClassificationEdge(
		double& ADistance,
		gmds::geom::GeomCurve& ACurve,
		const int& AMarkNodeContraint)
{
	ADistance = 0.;

	std::vector<gmds::geom::GeomPoint*> points;

	ACurve.get(points);
	if(points.size() != 2 && !ACurve.isALoop())
		throw GMDSException("Cannot work with a curve that does not have two points"
				"or that is not a loop");

	Variable<geom::GeomEntity*>* edgeClassification = this->mesh_.getGeometricClassification(1);

	Node starting_node;
	Node ending_node;
	Node curve_nodes[MESH_INSERT_DETAIL_INOUT_CURVE_NBPOINTS];

	// case when the curve is a loop
	if(ACurve.isALoop())
	{
		// get several points on the curve
		gmds::math::Point curve_points[MESH_INSERT_DETAIL_INOUT_CURVE_NBPOINTS];
		ACurve.getMultiplePoints(MESH_INSERT_DETAIL_INOUT_CURVE_NBPOINTS,curve_points);

		Variable<geom::GeomEntity*>* pointClassification = this->mesh_.getGeometricClassification(0);

		// now associate a node to each of these points but only if the loop has one point
		unsigned int iPointStart = 0;
		if(0 != points.size()) {
			curve_nodes[0] = geom2MeshNode_[points[0]];
			this->mesh_.mark(curve_nodes[0],AMarkNodeContraint);
			iPointStart++;
		}

		for(unsigned int iPoint=iPointStart; iPoint<MESH_INSERT_DETAIL_INOUT_CURVE_NBPOINTS; iPoint++) {

			double distmin = HUGE_VALF;

			IGMesh::node_iterator it  = this->mesh_.nodes_begin();

			for(;!it.isDone();it.next()) {
				if(((*pointClassification)[(it.value()).getID()] == NULL  ||
					((*pointClassification)[(it.value()).getID()])->getDim()!=0) &&
					!this->mesh_.isMarked(it.value(),AMarkNodeContraint) &&
					this->mesh_.isMarked(it.value(),markNodeIsOnBoundary_))
				{

					gmds::math::Point node_point = it.value().getPoint();

					double dist_tmp = node_point.distance(curve_points[iPoint]);
					if(dist_tmp < distmin) {
						curve_nodes[iPoint] = it.value();
						distmin = dist_tmp;
					}
				}
			}

			this->mesh_.mark(curve_nodes[iPoint],AMarkNodeContraint);
			std::cout<<"curve node id "<<curve_nodes[iPoint].getID()<<std::endl;
		}

	}
	else
	{
		gmds::geom::GeomPoint* point_k = points[0];
		gmds::geom::GeomPoint* point_k1 = points[1];

		starting_node = geom2MeshNode_[point_k];
		ending_node = geom2MeshNode_[point_k1];
		std::cout<<"starting_node "<<point_k->getId()-1<<" "<<starting_node.getID()<<std::endl;
		std::cout<<"ending_node "<<point_k1->getId()-1<<" "<<ending_node.getID()<<std::endl;

		// we offset starting_node and ending_nodes because we already have a starting edge and ending edge
		// at both ends of the curve.
		std::vector<gmds::geom::GeomCurve* > starting_curves;
		std::vector<gmds::geom::GeomCurve* > ending_curves;

		point_k->getOrderedDirect(starting_curves);
		point_k1->getOrderedDirect(ending_curves);

		// offsetting starting node
		for(unsigned int iCurve=0; iCurve<starting_curves.size(); iCurve++)
		{
			if(starting_curves[iCurve] == &ACurve)
			{
				gmds::Edge starting_edge = geomCurvesAtPoints2MeshEdge_[point_k][iCurve];

				geom2MeshEdge_[&ACurve].push_back(starting_edge);

				std::vector<gmds::Node> nodes;
				nodes = starting_edge.get<Node>();

				if(nodes[0] == starting_node)
				{
					starting_node = nodes[1];
				}
				else
				{
					starting_node = nodes[0];
				}

				// special case where there is only one edge between
				// the curve's two endpoints.
				if(starting_node == ending_node) {
					return true;
				}
			}
		}

		// offsetting ending node
		for(unsigned int iCurve=0; iCurve<ending_curves.size(); iCurve++)
		{
			if(ending_curves[iCurve] == &ACurve)
			{
				gmds::Edge ending_edge = geomCurvesAtPoints2MeshEdge_[point_k1][iCurve];

				geom2MeshEdge_[&ACurve].push_back(ending_edge);

				std::vector<gmds::Node> nodes;
				nodes = ending_edge.get<Node>();

				if(nodes[0] == ending_node)
				{
					ending_node = nodes[1];
				}
				else
				{
					ending_node = nodes[0];
				}
			}
		}

	}

	if(!ACurve.isALoop())
	{
		this->mesh_.unmark(starting_node,AMarkNodeContraint);
		this->mesh_.unmark(ending_node,AMarkNodeContraint);
	} else {
		for(unsigned int iPoint=0; iPoint<MESH_INSERT_DETAIL_INOUT_CURVE_NBPOINTS; iPoint++) {
			this->mesh_.unmark(curve_nodes[iPoint],AMarkNodeContraint);

			// case of the first point if the curve is a loop and a first edge was associated
			// to this curve
			if((0 != points.size()) && (0 == iPoint)) {
				Node starting_node = curve_nodes[iPoint];

				// offsetting starting node but only in the case where the loop has one point
				std::vector<gmds::geom::GeomCurve* > starting_curves;
				points[0]->getOrderedDirect(starting_curves);

				for(unsigned int iCurve=0; iCurve<starting_curves.size(); iCurve++)
				{
					if(&ACurve == starting_curves[iCurve])
					{
						gmds::Edge starting_edge = geomCurvesAtPoints2MeshEdge_[points[0]][iCurve];

						std::vector<gmds::Node> nodes;
						nodes = starting_edge.get<Node>();

						if(nodes[0] == curve_nodes[0]) {
							starting_node = nodes[1];
						} else {
							starting_node = nodes[0];
						}
					}
				}

				this->mesh_.unmark(starting_node,AMarkNodeContraint);
			} // if((1 == points.size()) && (0 == iPoint))
		}
	}

	// building a set of acceptable edges.
	// Edges that are owned by a hexahedron intersected by the curve.
	int markEdgeAcceptable = this->mesh_.getNewMark<Edge>();
	this->mesh_.unmarkAll<Edge>(markEdgeAcceptable);

	IGMesh::edge_iterator it_edges  = this->mesh_.edges_begin();

	for(;!it_edges.isDone();it_edges.next()){
		gmds::Edge current_edge = it_edges.value();
		if(this->mesh_.isMarked(current_edge,this->markEdgeIsOnBoundary_)) {
			if(!this->mesh_.isMarked(current_edge,AMarkNodeContraint)) {

				std::vector<Node> nodesOfEdge = current_edge.get<Node>();
				if(!this->mesh_.isMarked(nodesOfEdge[0],AMarkNodeContraint) && !this->mesh_.isMarked(nodesOfEdge[1],AMarkNodeContraint)) {
					this->mesh_.mark(current_edge,markEdgeAcceptable);

//					(*edgeClassification)[current_edge->getID()] = &ACurve;
				}
			}
//			(*edgeClassification)[current_edge->getID()] = &ACurve;
		}
	}
//	return;

//	{
//		typename Mesh<TMask>::regions_iterator it  = this->mesh_.regions_begin();
//
//		for(;!it->isDone();it->next()){
//			gmds::Region* current_region = it->currentItem();
//
//			std::vector<gmds::Node*> nodes;
//			current_region->getNodes(nodes);
//
//			GEPETO::Hexahedron<TBase> hexa(
//					nodes[4]->getPoint(),
//					nodes[7]->getPoint(),
//					nodes[6]->getPoint(),
//					nodes[5]->getPoint(),
//					nodes[0]->getPoint(),
//					nodes[3]->getPoint(),
//					nodes[2]->getPoint(),
//					nodes[1]->getPoint());
//
//			if(this->service_.intersects(ACurve,hexa))
//			{
//				std::vector<gmds::Face*> faces;
//				current_region->getFaces(faces);
//				for(unsigned int iFace=0; iFace<faces.size(); iFace++)
//				{
//					std::vector<gmds::Edge*> edges;
//					faces[iFace]->getEdges(edges);
//
//					for(unsigned int iEdge=0; iEdge<edges.size(); iEdge++)
//					{
//						if(this->mesh_.isMarked(edges[iEdge],this->markIsOnBoundary_)) {
//							std::vector<Node*> nodes = edges[iEdge]->getNodes();
//							if(!this->mesh_.isMarked(nodes[0],AMarkNodeContraint) && !this->mesh_.isMarked(nodes[1],AMarkNodeContraint)) {
//								this->mesh_.mark(edges[iEdge],markEdgeAcceptable);
////								if(curve->getId() == 7) {
////								(*edgeClassification)[edges[iEdge]->getID()] = &ACurve;
////								geom2MeshEdge_[curve].push_back(edges[iEdge]);
////								}
//							}
//
////							// dihedral angle criteria
////							// this criteria will be used only if we are not near
////							// starting_node and ending_node.
////							{
////								std::vector<Node*> nodes = edges[iEdge]->getNodes();
////							if(starting_node != nodes[0] && starting_node != nodes[1] && ending_node != nodes[0] && ending_node != nodes[1]) {
////
////							TBase angle_surf = curve->computeDihedralAngle();
////							if(fabs(angle_surf.toDouble()) > 10.) {
////								std::vector<Face*> faces = edges[iEdge]->getFaces();
////								Face* face0 = NULL;
////								Face* face1 = NULL;
////								for(unsigned int iFace=0; iFace<faces.size(); iFace++) {
////									if(this->mesh_.isMarked(faces[iFace],this->markIsOnBoundary_)) {
////										if(face0 == NULL) {
////											face0 = faces[iFace];
////										} else {
////											face1 = faces[iFace];
////										}
////									}
////								}
////
////								GEPETO::Vector<3,TBase> vector0 = face0->getNormal();
////								vector0.normalize();
////								GEPETO::Vector<3,TBase> vector1 = face1->getNormal();
////								vector1.normalize();
////
////								TBase absDotFacesNormals = fabs((vector0.dot(vector1)).toDouble());
////								if(absDotFacesNormals > 0.9) {
////									this->mesh_.unmark(edges[iEdge],markEdgeAcceptable);
////								}
////							}
////							} // if(angle_surf.abs() > 10.) {
////							}
//					}
//					}
//				}
//
//				// we also add the edges of the regions adjacent to the intersecting region.
//				// adjacent by edges
//				std::vector<Edge*> edges = current_region->getEdges();
//				for(unsigned int iEdge=0; iEdge<edges.size(); iEdge++) {
//					std::vector<Region*> regions = edges[iEdge]->getRegions();
//
//					for(unsigned int iRegion=0; iRegion<regions.size(); iRegion++) {
//						std::vector<Edge*> edgestoMark = regions[iRegion]->getEdges();
//						for(unsigned int iEdgeToMark=0; iEdgeToMark<edgestoMark.size(); iEdgeToMark++) {
//							if(this->mesh_.isMarked(edgestoMark[iEdgeToMark],this->markIsOnBoundary_)) {
//								std::vector<Node*> nodes = edgestoMark[iEdgeToMark]->getNodes();
//								if(!this->mesh_.isMarked(nodes[0],AMarkNodeContraint) && !this->mesh_.isMarked(nodes[1],AMarkNodeContraint)) {
//									this->mesh_.mark(edgestoMark[iEdgeToMark],markEdgeAcceptable);
////									if(curve->getId() == 7) {
////									(*edgeClassification)[edgestoMark[iEdgeToMark]->getID()] = &ACurve;
////									geom2MeshEdge_[curve].push_back(edgestoMark[iEdgeToMark]);
////									}
//								}
////								// dihedral angle criteria
////								// this criteria will be used only if we are not near
////								// starting_node and ending_node.
////								{
////									std::vector<Node*> nodes = edgestoMark[iEdgeToMark]->getNodes();
////								if(starting_node != nodes[0] && starting_node != nodes[1] && ending_node != nodes[0] && ending_node != nodes[1]) {
////
////								TBase angle_surf = curve->computeDihedralAngle();
////								if(fabs(angle_surf.toDouble()) > 10.) {
////									std::vector<Face*> faces = edgestoMark[iEdgeToMark]->getFaces();
////									Face* face0 = NULL;
////									Face* face1 = NULL;
////									for(unsigned int iFace=0; iFace<faces.size(); iFace++) {
////										if(this->mesh_.isMarked(faces[iFace],this->markIsOnBoundary_)) {
////											if(face0 == NULL) {
////												face0 = faces[iFace];
////											} else {
////												face1 = faces[iFace];
////											}
////										}
////									}
////
////									GEPETO::Vector<3,TBase> vector0 = face0->getNormal();
////									vector0.normalize();
////									GEPETO::Vector<3,TBase> vector1 = face1->getNormal();
////									vector1.normalize();
////
////									TBase absDotFacesNormals = fabs((vector0.dot(vector1).toDouble()));
////									if(absDotFacesNormals > 0.9) {
////										this->mesh_.unmark(edgestoMark[iEdgeToMark],markEdgeAcceptable);
////									}
////								} // if(angle_surf.abs() > 10.) {
////								}
////								}
//							}
//						}
//					}
//				}
//
//			}
//		}
//	}
//	return;
//	if(curve->getId() == 5) {
//		this->exportVTK("poyop.mli");
//		return;
//	}


	if(!ACurve.isALoop())
	{
		std::cout<<"looking for curve "<<starting_node.getID()<<" "<<ending_node.getID()<<std::endl;
	} else {
		std::cout<<"looking for curve ";
		for(unsigned int iPoint=0; iPoint<MESH_INSERT_DETAIL_INOUT_CURVE_NBPOINTS; iPoint++) {
			std::cout<<curve_nodes[iPoint].getID()<<" ";
		}
		std::cout<<std::endl;
	}

	if((starting_node == ending_node) && (!ACurve.isALoop())) {
		// nothing to do; this is the case when the curve is associated to two edges.
	} else {

	if(starting_node != ending_node) {

		MeshGraph meshGraph(this->mesh_,markEdgeAcceptable);
		std::vector<Edge> edges_path;


//		meshGraph.findShortestPath(starting_node,ending_node,edges_path);
		double distance = 0.;
		bool pathFound = meshGraph.findShortestPathHaussdorf(distance,starting_node,ending_node,edges_path,&ACurve);
		if(!pathFound) {
			std::cout<<"path not found"<<std::endl;
			this->mesh_.unmarkAll<Edge>(markEdgeAcceptable);
			this->mesh_.freeMark<Edge>(markEdgeAcceptable);
			return false;
		}
		ADistance = distance;

		for(unsigned int iEdge=0; iEdge<edges_path.size(); iEdge++) {
			(*edgeClassification)[edges_path[iEdge].getID()] = &ACurve;
			geom2MeshEdge_[&ACurve].push_back(edges_path[iEdge]);

			std::vector<Node> nodes = edges_path[iEdge].get<Node>();
			this->mesh_.mark(nodes[0],AMarkNodeContraint);
			this->mesh_.mark(nodes[1],AMarkNodeContraint);
		}
	} else {
		for(unsigned int iPoint=0; iPoint<MESH_INSERT_DETAIL_INOUT_CURVE_NBPOINTS; iPoint++) {

			Node starting_node;
			Node ending_node;

			starting_node = curve_nodes[iPoint];

			// offsetting starting node but only in the case where the loop has one point
			if((1 == points.size()) && (0 == iPoint)) {
				std::vector<gmds::geom::GeomCurve* > starting_curves;
				points[0]->getOrderedDirect(starting_curves);

				for(unsigned int iCurve=0; iCurve<starting_curves.size(); iCurve++)
				{
					if(&ACurve == starting_curves[iCurve])
					{
						gmds::Edge starting_edge = geomCurvesAtPoints2MeshEdge_[points[0]][iCurve];

						geom2MeshEdge_[&ACurve].push_back(starting_edge);

						std::vector<gmds::Node> nodes;
						nodes = starting_edge.get<Node>();

						if(nodes[0] == curve_nodes[0]) {
							starting_node = nodes[1];
						} else {
							starting_node = nodes[0];
						}
					}
				}
			} // if(0 == iPoint) {

			ending_node = curve_nodes[(iPoint+1)%MESH_INSERT_DETAIL_INOUT_CURVE_NBPOINTS];

			MeshGraph meshGraph(this->mesh_,markEdgeAcceptable);
			std::vector<Edge> edges_path;
			std::cout<<"looking for subcurve "<<starting_node.getID()<<" "<<ending_node.getID()<<std::endl;

//			meshGraph.findShortestPathHaussdorfSubCurve(MESH_INSERT_DETAIL_INOUT_CURVE_NBPOINTS,iPoint,curve_nodes[iPoint],curve_nodes[(iPoint+1)%MESH_INSERT_DETAIL_INOUT_CURVE_NBPOINTS],edges_path,&ACurve);
			double distance = 0.;
			bool pathFound = meshGraph.findShortestPathHaussdorfSubCurve(distance,MESH_INSERT_DETAIL_INOUT_CURVE_NBPOINTS,iPoint,starting_node,ending_node,edges_path,&ACurve);
			if(!pathFound) {
				std::cout<<"path not found"<<std::endl;
				this->mesh_.unmarkAll<Edge>(markEdgeAcceptable);
				this->mesh_.freeMark<Edge>(markEdgeAcceptable);
				return false;
			}
			ADistance += distance;

			for(unsigned int iEdge=0; iEdge<edges_path.size(); iEdge++) {
				(*edgeClassification)[edges_path[iEdge].getID()] = &ACurve;
				geom2MeshEdge_[&ACurve].push_back(edges_path[iEdge]);

				std::vector<Node> nodes = edges_path[iEdge].get<Node>();
				this->mesh_.mark(nodes[0],AMarkNodeContraint);
				this->mesh_.mark(nodes[1],AMarkNodeContraint);

				if((nodes[0] != curve_nodes[iPoint]) && (nodes[0] != curve_nodes[(iPoint+1)%MESH_INSERT_DETAIL_INOUT_CURVE_NBPOINTS])) {
					std::vector<Edge> edges = nodes[0].get<Edge>();
					for(unsigned int iEdge=0; iEdge<edges.size(); iEdge++) {
						this->mesh_.unmark(edges[iEdge],markEdgeAcceptable);
					}
				}
				if((nodes[1] != curve_nodes[iPoint]) && (nodes[1] != curve_nodes[(iPoint+1)%MESH_INSERT_DETAIL_INOUT_CURVE_NBPOINTS])) {
					std::vector<Edge> edges = nodes[1].get<Edge>();
					for(unsigned int iEdge=0; iEdge<edges.size(); iEdge++) {
						this->mesh_.unmark(edges[iEdge],markEdgeAcceptable);
					}
				}

				this->mesh_.unmark(edges_path[iEdge],markEdgeAcceptable);
			}
		}
	} // if(starting_node != ending_node)

	} // if((starting_node == ending_node) && (!ACurve.isALoop()))

//	for(unsigned int iEdge=0; iEdge<edges_path.size(); iEdge++) {
//		(*edgeClassification)[edges_path[iEdge]->getID()] = &ACurve;
//		geom2MeshEdge_[curve->getId()-1].push_back(edges_path[iEdge]);
//
//		std::vector<Node*> nodes = edges_path[iEdge]->getNodes();
//		this->mesh_.mark(nodes[0],AMarkNodeContraint);
//		this->mesh_.mark(nodes[1],AMarkNodeContraint);
//	}

//	this->mesh_.mark(starting_node,AMarkNodeContraint);
//	this->mesh_.mark(ending_node,AMarkNodeContraint);

	std::cout<<std::endl;

	this->mesh_.unmarkAll<Edge>(markEdgeAcceptable);
	this->mesh_.freeMark<Edge>(markEdgeAcceptable);

	return true;
}
/*----------------------------------------------------------------------------*/
void
MeshInsertDetailInOut::removeGeometryClassificationCurves(gmds::geom::GeomVolume* AVol)
{
	std::vector<gmds::geom::GeomCurve*> curves;
	AVol->get(curves);

	Variable<geom::GeomEntity*>* edgeClassification = this->mesh_.getGeometricClassification(1);

	// remove all edge association
	IGMesh::edge_iterator it  = this->mesh_.edges_begin();
	for(;!it.isDone();it.next()) {
		(*edgeClassification)[(it.value()).getID()] = NULL;
	}

	// we re-associate the first edges to each point
	std::vector<gmds::geom::GeomPoint*> points;
	AVol->get(points);
	for(unsigned int iPoint=0; iPoint<points.size(); iPoint++)
	{
		std::vector<gmds::geom::GeomCurve*> curves;
		points[iPoint]->getOrderedDirect(curves);

		for(unsigned int iCurveAtPoint=0; iCurveAtPoint<geomCurvesAtPoints2MeshEdge_[points[iPoint]].size(); iCurveAtPoint++) {
			Edge current_edge = geomCurvesAtPoints2MeshEdge_[points[iPoint]][iCurveAtPoint];

			(*edgeClassification)[current_edge.getID()] = curves[iCurveAtPoint];
		}
	}
}
/*----------------------------------------------------------------------------*/
void
MeshInsertDetailInOut::associateGeometryClassificationSurfaces(gmds::geom::GeomVolume* AVol)
{
	Variable<geom::GeomEntity* >* edgeClassification = this->mesh_.getGeometricClassification(1);
	Variable<geom::GeomEntity* >* faceClassification = this->mesh_.getGeometricClassification(2);

	std::vector<gmds::geom::GeomSurface* > surfaces;

	AVol->get(surfaces);

	int markFaceAcceptable = this->mesh_.getNewMark<Face>();
	unsigned int nbFaceAcceptable = 0;

	{
		IGMesh::face_iterator it  = this->mesh_.faces_begin();

		for(;!it.isDone();it.next()){
			if(this->mesh_.isMarked(it.value(),this->markFaceIsOnBoundary_)) {
				this->mesh_.mark(it.value(),markFaceAcceptable);
				nbFaceAcceptable++;
//				(*faceClassification)[it.value().getID()] = surfaces[0];
			}
		}
	}

//	this->exportMeshVTK("/homePOYOP/travail/workspaces/gscc_workspace/Cube/src/CaGe2/debug_surf.mli",gmds::N|gmds::E|gmds::F|gmds::R);
//	exit(-1);

	std::cout<<"there are "<<nbFaceAcceptable<<" acceptable faces"<<std::endl;
//	return;

	int markEdgeBarrier = this->mesh_.getNewMark<Edge>();
	{
		IGMesh::edge_iterator it  = this->mesh_.edges_begin();

		for(;!it.isDone();it.next()){
			if((*edgeClassification)[it.value().getID()] != NULL) {
				this->mesh_.mark(it.value(),markEdgeBarrier);
			}

		}
	}

	// we will pick at random a face, propagate a mark till meeting associated edges
	// and associate a corresponding surface to this set of faces.
	while(nbFaceAcceptable>0) {
		IGMesh::face_iterator it  = this->mesh_.faces_begin();


		while(!this->mesh_.isMarked(it.value(),markFaceAcceptable)) {
			it.next();
		}

		if(it.isDone())
			throw GMDSException("there should still be at least one acceptable face.");

		// recursive
		int markFaceRecurse = this->mesh_.getNewMark<Face>();

		std::set<geom::GeomEntity* > curvesSet;
		markRecurse(it.value(),markFaceRecurse,markFaceAcceptable,markEdgeBarrier,edgeClassification,curvesSet);

		unsigned int nbFaceRecurse = 0;

		{
			IGMesh::face_iterator it  = this->mesh_.faces_begin();

			for(;!it.isDone();it.next()){
				if(this->mesh_.isMarked(it.value(),markFaceRecurse)) {
					nbFaceRecurse++;
					(*faceClassification)[it.value().getID()] = surfaces[1];
				}
			}
		}
//		this->exportMeshVTK("/homePOYOP/travail/workspaces/gscc_workspace/Cube/src/CaGe2/debug_surf.mli",gmds::N|gmds::E|gmds::F|gmds::R);
		std::cout<<"nb of faces recurse "<<nbFaceRecurse<<std::endl;
//		exit(-1);

		std::cout<<"set of faces with "<<curvesSet.size()<<" curves."<<std::endl;
		std::set<geom::GeomEntity* >::iterator its = curvesSet.begin();
		for(;its!=curvesSet.end();its++){
			geom::GeomEntity* e = *its;
			std::cout<<e<<" ";
		}
		std::cout<<std::endl;
		// finding which surface is concerned
		gmds::geom::GeomSurface* current_surface = NULL;

		for(unsigned int iSurface=0; iSurface<surfaces.size(); iSurface++) {
			std::vector<gmds::geom::GeomCurve* > curves;
			surfaces[iSurface]->get(curves);
			std::cout<<"Nb curves for surf "<<iSurface<<"= "<<curves.size()<<std::endl;
			bool found = true;
			if(curvesSet.size() != curves.size())
				found = false;
			for(unsigned int iCurve=0; iCurve<curves.size(); iCurve++) {
				std::cout<<curves[iCurve]<<" ";
			}
			std::cout<<std::endl;
			for(unsigned int iCurve=0; iCurve<curves.size(); iCurve++) {
				if(curvesSet.end() == curvesSet.find(curves[iCurve])) {
					found = false;
				}
			}
			if(found) {
				current_surface = surfaces[iSurface];
				break;
			}
		}

		if(current_surface == NULL)
			throw GMDSException("surface not found.");

		{
			IGMesh::face_iterator it  = this->mesh_.faces_begin();

			for(;!it.isDone(); it.next()) {
				if(this->mesh_.isMarked(it.value(),markFaceRecurse)) {
					this->mesh_.unmark(it.value(),markFaceAcceptable);
					nbFaceAcceptable--;
					(*faceClassification)[it.value().getID()] = current_surface;
					geom2MeshFace_[current_surface].push_back(it.value());
				}
			}
		}

		this->mesh_.unmarkAll<Face>(markFaceRecurse);
		this->mesh_.freeMark<Face>(markFaceRecurse);
		//break;
	}

	this->mesh_.unmarkAll<Edge>(markEdgeBarrier);
	this->mesh_.freeMark<Edge>(markEdgeBarrier);
	this->mesh_.unmarkAll<Face>(markFaceAcceptable);
	this->mesh_.freeMark<Face>(markFaceAcceptable);
}
/*----------------------------------------------------------------------------*/
void
MeshInsertDetailInOut::markRecurse(Face AFace, const int& AMark, const int& AMarkAcceptable, const int& AMarkBarrier,
		Variable<geom::GeomEntity* >* AEdgeClassification,
		std::set<geom::GeomEntity* >& ACurvesSet)
{
	// if the face is already marked, do nothing
	if(this->mesh_.isMarked(AFace,AMark)) {
		return;
	}

	this->mesh_.mark(AFace,AMark);

	std::vector<Edge> edges = AFace.get<Edge>();
	for(unsigned iEdge=0; iEdge<edges.size(); iEdge++) {
		if(!this->mesh_.isMarked(edges[iEdge],AMarkBarrier)) {
			std::vector<Face> faces = edges[iEdge].get<Face>();

			for(unsigned int iFace=0; iFace<faces.size(); iFace++) {
				if(!this->mesh_.isMarked(faces[iFace],AMark) && this->mesh_.isMarked(faces[iFace],AMarkAcceptable)) {
					markRecurse(faces[iFace],AMark,AMarkAcceptable,AMarkBarrier,AEdgeClassification,ACurvesSet);
				}
			}
		} else {
			if(ACurvesSet.end() == ACurvesSet.find((*AEdgeClassification)[edges[iEdge].getID()])) {
				ACurvesSet.insert((*AEdgeClassification)[edges[iEdge].getID()]);
			}
		}
	}

}
///*----------------------------------------------------------------------------*/
//bool
//MeshInsertDetailInOut::checkNonManifoldnessEdges(gmds::geom::GeomVolume* AVol)
//{
//	// marking regions
//	IGMesh::edge_iterator ite  = this->mesh_.edges_begin();
//
//	for(;!ite.isDone();ite.next()){
//		Edge current_edge = ite.value();
//
//		std::vector<Face> faces = current_edge.get<Face>();
//
//		// we only keep the boundary faces
//	}
//}
/*----------------------------------------------------------------------------*/
void
MeshInsertDetailInOut::project()
{
	Variable<geom::GeomEntity* >* nodeClassification   = this->mesh_.getGeometricClassification(0);
	Variable<geom::GeomEntity* >* edgeClassification   = this->mesh_.getGeometricClassification(1);
	Variable<geom::GeomEntity* >* faceClassification   = this->mesh_.getGeometricClassification(2);
	Variable<geom::GeomEntity* >* regionClassification = this->mesh_.getGeometricClassification(3);

	int markNodeDontMove = this->mesh_.getNewMark<gmds::Node>();

	{
		gmds::IGMesh::node_iterator it  = this->mesh_.nodes_begin();
		for(;!it.isDone(); it.next()) {
			gmds::Node current_node = it.value();

			gmds::geom::GeomEntity* ge =(*nodeClassification)[current_node.getID()];
			if(ge != NULL)
			{
				if(ge->getDim()==0)
				{
					//this->mesh_.mark(current_node,markNodeDontMove);
					gmds::geom::GeomPoint* point =  dynamic_cast<gmds::geom::GeomPoint* >(ge);
					current_node.setPoint(point->getPoint());
				}
				else if (ge->getDim()==1)
				{
					//this->mesh_.mark(current_node,markNodeDontMove);
					gmds::geom::GeomCurve* curve = dynamic_cast<gmds::geom::GeomCurve* >(ge);

					gmds::math::Point point = current_node.getPoint();
					curve->project(point);
					current_node.setPoint(point);
				}
				else if (ge->getDim()==2){
					//this->mesh_.mark(current_node,markNodeDontMove);
					gmds::geom::GeomSurface* surface = dynamic_cast<gmds::geom::GeomSurface* >(ge);

					gmds::math::Point point = current_node.getPoint();
					//surface->project(point);
					GNode* gnode = (GNode*) (this->aabbSurfacesTrianglesTrees_)[surface];
					this->service_.project(point,gnode);
					current_node.setPoint(point);
				}

				//	continue;
			}
			//	} else {


//			bool edgeFound = false;
//			std::vector<Edge> edges = current_node.get<Edge>();
//			for(unsigned int iEdge=0; iEdge<edges.size(); iEdge++) {
//				if((*edgeClassification)[edges[iEdge].getID()] != NULL) {
//					gmds::geom::GeomCurve* curve =
//							(dynamic_cast<gmds::geom::GeomCurve* >((*edgeClassification)[edges[iEdge].getID()]));
//
//					gmds::math::Point point = current_node.getPoint();
//					curve->project(point);
//					current_node.setPoint(point);
//					edgeFound = true;
//					break;
//				}
//			}
//			if(edgeFound)
//				continue;
//
//			bool faceFound = false;
//			std::vector<Face> faces = current_node.get<Face>();
//			for(unsigned int iFace=0; iFace<faces.size(); iFace++) {
//				if((*faceClassification)[faces[iFace].getID()] != NULL) {
//					gmds::geom::GeomSurface* surface =
//							(dynamic_cast<gmds::geom::GeomSurface* >((*faceClassification)[faces[iFace].getID()]));
//
//					gmds::math::Point point = current_node.getPoint();
//					surface->project(point);
//					current_node.setPoint(point);
//					faceFound = true;
//					break;
//				}
//			}
//			if(faceFound)
//				continue;

		}
	}

	this->mesh_.unmarkAll<Node>(markNodeDontMove);
	this->mesh_.freeMark<Node>(markNodeDontMove);
}
/*----------------------------------------------------------------------------*/
void
MeshInsertDetailInOut::completeNodeAssociation(gmds::geom::GeomVolume* AVol)
{
	std::vector<gmds::geom::GeomPoint*> points;
	std::vector<gmds::geom::GeomCurve*> curves;
	std::vector<gmds::geom::GeomSurface*> surfaces;

	AVol->get(points);
	AVol->get(curves);
	AVol->get(surfaces);

	Variable<geom::GeomEntity* >* nodeClassification   = this->mesh_.getGeometricClassification(0);
	Variable<geom::GeomEntity* >* edgeClassification   = this->mesh_.getGeometricClassification(1);
	Variable<geom::GeomEntity* >* faceClassification   = this->mesh_.getGeometricClassification(2);
	Variable<geom::GeomEntity* >* regionClassification = this->mesh_.getGeometricClassification(3);

	gmds::IGMesh::region_iterator itr  = this->mesh_.regions_begin();
	for(;!itr.isDone(); itr.next()) {
		gmds::Region current_region = itr.value();

		if((*regionClassification)[current_region.getID()]  == AVol) {

			std::vector<gmds::Node> nodes = current_region.get<Node>();

			for(int iNode=0; iNode<nodes.size(); iNode++) {
				if((*nodeClassification)[nodes[iNode].getID()] == NULL || (*nodeClassification)[nodes[iNode].getID()]->getDim()  != 0) {
					(*nodeClassification)[nodes[iNode].getID()] = AVol;
				}
			}
		}
	}

	for(int iSurf=0; iSurf<surfaces.size(); iSurf++) {

		gmds::geom::GeomSurface* surface = surfaces[iSurf];

		for(int iFace=0; iFace<this->geom2MeshFace_[surface].size(); iFace++) {
			std::vector<gmds::Node> nodes = geom2MeshFace_[surface][iFace].get<gmds::Node>();

			for(int iNode=0; iNode<nodes.size(); iNode++) {
				if((*nodeClassification)[nodes[iNode].getID()] == NULL || (*nodeClassification)[nodes[iNode].getID()]->getDim()  != 0) {
					(*nodeClassification)[nodes[iNode].getID()] = surface;
				}
			}
		}
	}

	for(int iCurve=0; iCurve<curves.size(); iCurve++) {

		gmds::geom::GeomCurve* curve = curves[iCurve];

		for(int iEdge=0; iEdge<this->geom2MeshEdge_[curve].size(); iEdge++) {
			std::vector<gmds::Node> nodes = geom2MeshEdge_[curve][iEdge].get<gmds::Node>();

			for(int iNode=0; iNode<nodes.size(); iNode++) {
				if((*nodeClassification)[nodes[iNode].getID()] == NULL || (*nodeClassification)[nodes[iNode].getID()]->getDim()  != 0) {
					(*nodeClassification)[nodes[iNode].getID()] = curve;
				}
			}
		}
	}

}
/*----------------------------------------------------------------------------*/
unsigned int
MeshInsertDetailInOut::intersectedCellsDetectionAABBTree(
		std::map<gmds::TCellID, std::vector<gmds::math::Triangle> >& AInOutTriangles)
{
	// first get the intersecting cells IDs bbox / triangles bbox pairs.
	std::map<gmds::TCellID, std::vector<gmds::math::Triangle*> > inOutTriangles_tmp;

	GNode* boxTree = this->aabbRegionsTree_;
	this->service_.intersects(boxTree,inOutTriangles_tmp);

	// then we keep only the intersecting cell/triangles pairs.
	unsigned int nbIsOnIntersection = 0;

	std::map<gmds::TCellID, std::vector<gmds::math::Triangle*> >::iterator itCellID = inOutTriangles_tmp.begin();

	for(; itCellID != inOutTriangles_tmp.end(); itCellID++) {

		gmds::Region current_region = this->mesh_.get<gmds::Region>(itCellID->first);
		bool isRegionIntersected = false;

		for(unsigned int iTriangle=0; iTriangle<itCellID->second.size(); iTriangle++) {

			if(this->service_.intersects(*(itCellID->second[iTriangle]),current_region)) {
				AInOutTriangles[current_region.getID()].push_back(*(itCellID->second[iTriangle]));
				isRegionIntersected = true;
			}
		}
		if(isRegionIntersected) {
			nbIsOnIntersection++;
		}
	}

	return nbIsOnIntersection;
}
/*----------------------------------------------------------------------------*/
unsigned int
MeshInsertDetailInOut::intersectedCellsDetection(
		std::map<gmds::TCellID, std::vector<gmds::math::Triangle> >& AInOutTriangles)
{
	unsigned int nbIsOnIntersection = 0;

	// used for progress display
	int iPercent = 0;
	unsigned int iCell = 0;

	// marking intersecting cells
	IGMesh::region_iterator it  = this->mesh_.regions_begin();
	for(;!it.isDone();it.next()){
		gmds::Region current_region = it.value();

		// display progression
		if(iCell == ((this->mesh_.getNbRegions()/10)*iPercent)) {
			std::cout<<"intersect treating cell "<<iCell<<" of "<<this->mesh_.getNbRegions()<<" "<<iPercent*10<<"%"<<std::endl;
			iPercent++;
		}
		iCell++;

	       // check whether the hex is inside, outside or intersects the model
		std::vector<gmds::math::Triangle> intersectedTriangles;

		if(this->service_.intersects(current_region, intersectedTriangles))
		{
			AInOutTriangles[current_region.getID()] = intersectedTriangles;
			nbIsOnIntersection++;
		}
	} // for(;it!=ite;it++){

	return nbIsOnIntersection;
}
/*----------------------------------------------------------------------------*/
void
MeshInsertDetailInOut::initialization()
{
    GSList* list = NULL;
	gmds::IGMesh::region_iterator itr  = this->mesh_.regions_begin();

	for(;!itr.isDone();itr.next()) {
		double minXYZ[3];
		double maxXYZ[3];

		gmds::Region current_region = itr.value();
		current_region.computeBoundingBox(minXYZ,maxXYZ);

		gpointer pointer = GINT_TO_POINTER(current_region.getID());
		GtsBBox* bbox = gts_bbox_new(
								gts_bbox_class (),
								pointer,
								minXYZ[0],minXYZ[1],minXYZ[2],
								maxXYZ[0],maxXYZ[1],maxXYZ[2]);

		list = g_slist_prepend(list,bbox);
	}
	this->aabbRegionsTree_ = gts_bb_tree_new(list);
}
/*----------------------------------------------------------------------------*/
void
MeshInsertDetailInOut::exportCurvesEdgesVTK(const std::string& AFile) 
{
	gmds::MeshModel mod = gmds::DIM3|gmds::N|gmds::F|gmds::F2N;
	gmds::IGMesh mesh(mod);

	gmds::Variable<gmds::geom::GeomEntity* >* curveClassification  = this->mesh_.getGeometricClassification(1);
	std::vector<gmds::geom::GeomCurve* > curves;
	this->manager_.getCurves(curves);

	gmds::Variable<int>* curvesIndex= mesh.newVariable<int>(GMDS_FACE,"exportCurvesEdgesVTK_curvesIndex");

	for(int iLine=0; iLine<this->manager_.getNbCurves(); iLine++)
        {
		gmds::geom::GeomCurve* current_curve = curves[iLine];
		gmds::IGMesh::surface& surf = mesh.newSurface(current_curve->getName());
		
		for(unsigned int iEdge=0; iEdge<geom2MeshEdge_[current_curve].size(); iEdge++) {
			std::vector<gmds::Node> nodes = geom2MeshEdge_[current_curve][iEdge].get<gmds::Node>();
			gmds::Node n0 = nodes[0];
			gmds::Node n1 = nodes[1];

			gmds::Node newN0 = mesh.newNode(n0.getPoint());
			gmds::Node newN1 = mesh.newNode(n1.getPoint());
				
			gmds::Face newF = mesh.newTriangle(newN0,newN0,newN1);
			surf.add(newF);
		
			(*curvesIndex)[newF.getID()] = iLine;			
		}
	
	}

	gmds::VTKWriter<gmds::IGMesh> w(mesh);
        w.write(AFile,gmds::N|gmds::F);
}
/*----------------------------------------------------------------------------*/
bool
MeshInsertDetailInOut::checkNonManifoldnessEdges(gmds::geom::GeomVolume* AVol)
{
  gmds::Variable<gmds::geom::GeomEntity* >* volumeClassification  = this->mesh_.getGeometricClassification(3);

  int markIsOnBoundaryFace = this->mesh_.getNewMark<gmds::Face>();

  // mark boundary faces
  gmds::IGMesh::face_iterator itf  = this->mesh_.faces_begin();
  for(;!itf.isDone();itf.next()) {
    gmds::Face current_face = itf.value();

    std::vector<gmds::Region> regions = current_face.get<gmds::Region>();

    if(1 == regions.size()) {
      if(AVol == (*volumeClassification)[regions[0].getID()]) {
	this->mesh_.mark(current_face,markIsOnBoundaryFace);
      }
    } else if(2 == regions.size()) {
      if((AVol == (*volumeClassification)[regions[0].getID()]) && 
	 (AVol != (*volumeClassification)[regions[1].getID()])) {
	this->mesh_.mark(current_face,markIsOnBoundaryFace);
      }
    } else {
      std::cout<<"regions.size() "<<regions.size()<<std::endl;
      this->mesh_.unmarkAll<gmds::Face>(markIsOnBoundaryFace);
      this->mesh_.freeMark<gmds::Face>(markIsOnBoundaryFace);
      
      throw GMDSException("MeshInsertDetailInOut::checkNonManifoldnessEdges bad F2R adjacency.");
      
    }
  }

  gmds::IGMesh::edge_iterator ite  = this->mesh_.edges_begin();
  for(;!ite.isDone();ite.next()) {
    gmds::Edge current_edge = ite.value();
    
    std::vector<gmds::Face> faces = current_edge.get<gmds::Face>();

    unsigned int nbBoundaryFaces = 0;

    for(int iFace=0; iFace<faces.size(); iFace++) {
      if(this->mesh_.isMarked(faces[iFace],markIsOnBoundaryFace)) {
	nbBoundaryFaces++;
      }
    }
    
    if(nbBoundaryFaces>2) {
      std::cout<<"found non-manifold edge "<<current_edge.getID()<<std::endl;
      this->mesh_.unmarkAll<gmds::Face>(markIsOnBoundaryFace);
      this->mesh_.freeMark<gmds::Face>(markIsOnBoundaryFace);

      return false;
    }
  }
   
  this->mesh_.unmarkAll<gmds::Face>(markIsOnBoundaryFace);
  this->mesh_.freeMark<gmds::Face>(markIsOnBoundaryFace);
    
  return true;
}
/*----------------------------------------------------------------------------*/
unsigned int 
MeshInsertDetailInOut::countManifoldEdges(gmds::geom::GeomVolume* AVol)
{
  gmds::Variable<gmds::geom::GeomEntity* >* volumeClassification  = this->mesh_.getGeometricClassification(3);

  int markIsOnBoundaryFace = this->mesh_.getNewMark<gmds::Face>();

  // mark boundary faces
  gmds::IGMesh::face_iterator itf  = this->mesh_.faces_begin();
  for(;!itf.isDone();itf.next()) {
    gmds::Face current_face = itf.value();

    std::vector<gmds::Region> regions = current_face.get<gmds::Region>();

    if(1 == regions.size()) {
      if(AVol == (*volumeClassification)[regions[0].getID()]) {
	this->mesh_.mark(current_face,markIsOnBoundaryFace);
      }
    } else if(2 == regions.size()) {
      if((AVol == (*volumeClassification)[regions[0].getID()]) && 
	 (AVol != (*volumeClassification)[regions[1].getID()])) {
	this->mesh_.mark(current_face,markIsOnBoundaryFace);
      }
    } else {
      std::cout<<"regions.size() "<<regions.size()<<std::endl;
      this->mesh_.unmarkAll<gmds::Face>(markIsOnBoundaryFace);
      this->mesh_.freeMark<gmds::Face>(markIsOnBoundaryFace);
      
      throw GMDSException("MeshInsertDetailInOut::checkNonManifoldnessEdges bad F2R adjacency.");   
    }
  }

  unsigned int nbManifoldEdges = 0;

  gmds::IGMesh::edge_iterator ite  = this->mesh_.edges_begin();
  for(;!ite.isDone();ite.next()) {
    gmds::Edge current_edge = ite.value();

    std::vector<gmds::Face> faces = current_edge.get<gmds::Face>();

    unsigned int nbBoundaryFaces = 0;

    for(int iFace=0; iFace<faces.size(); iFace++) {
      if(this->mesh_.isMarked(faces[iFace],markIsOnBoundaryFace)) {
	nbBoundaryFaces++;
      }
    }
    
    if(nbBoundaryFaces>2) {
      nbManifoldEdges++;
    }
  }
   
  this->mesh_.unmarkAll<gmds::Face>(markIsOnBoundaryFace);
  this->mesh_.freeMark<gmds::Face>(markIsOnBoundaryFace);
    
  return nbManifoldEdges;
}
/*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/
