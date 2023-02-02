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
/** \file    MeshInsertDetail.t.h
 *  \author  N. LE GOFF
 *  \date    27/06/2012
 */
/*----------------------------------------------------------------------------*/
// stl File Headers
/*----------------------------------------------------------------------------*/
#include <sstream>
#include <set>
#include <map>
/*----------------------------------------------------------------------------*/
// External File Headers
/*----------------------------------------------------------------------------*/
#include <GMDS/IG/IGMesh.h>
#include <GMDS/IO/VTKWriter.h>
#include <GMDS/CAD/GeomVolume.h>
#include <GMDS/CAD/GeomSurface.h>
#include <GMDS/CAD/FacetedCurve.h>
#include <GMDS/MeshInsertDetail.h>
#include <GMDS/Algo/SheetOperator.h>
#include <GMDS/Algo/Pillowing.h>
#include <GMDS/IG/IGMeshQualityEvaluation.h>
#include <GMDS/Utils/CommonTypes.h>
/*----------------------------------------------------------------------------*/
namespace gmds {
/*----------------------------------------------------------------------------*/
MeshInsertDetail::MeshInsertDetail(
		IGMesh& AMesh,
		gmds::geom::GeomManager& AManager,
		gmds::geom::GeomMeshIntersectionService& AService)
:mesh_(AMesh), manager_(AManager), service_(AService)
{}
/*----------------------------------------------------------------------------*/
MeshInsertDetail::~MeshInsertDetail()
{}
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void
MeshInsertDetail::exportMeshVTK(const std::string& AFile)
{
	gmds::VTKWriter<gmds::IGMesh> w(this->mesh_);
	w.write(AFile,gmds::R|gmds::N);
}
/*----------------------------------------------------------------------------*/
void
MeshInsertDetail::exportMeshVTK(const std::string& AFile, const gmds::MeshModel& AMeshModel)
{
	std::cout<<"begin MeshInsertDetail::exportMeshVTK "<<AFile<<std::endl;

	Variable<geom::GeomEntity* >* nodeClassification   = this->mesh_.getGeometricClassification(0);
	Variable<geom::GeomEntity* >* curveClassification  = this->mesh_.getGeometricClassification(1);
	Variable<geom::GeomEntity* >* surfaceClassification= this->mesh_.getGeometricClassification(2);
	Variable<geom::GeomEntity* >* volumeClassification = this->mesh_.getGeometricClassification(3);

	std::vector<gmds::geom::GeomVolume* > volumes;
	std::vector<gmds::geom::GeomSurface* > surfaces;
	std::vector<gmds::geom::GeomCurve* > curves;
	std::vector<gmds::geom::GeomPoint* > points;
	this->manager_.getVolumes(volumes);
	this->manager_.getSurfaces(surfaces);
	this->manager_.getCurves(curves);
	this->manager_.getPoints(points);


	for(int iVolume=0; iVolume<this->manager_.getNbVolumes(); iVolume++)
	{
		std::ostringstream  volume_name;
		volume_name<<"volume_"<<iVolume;

		gmds::IGMesh::volume& vol = this->mesh_.newVolume(volume_name.str());

		gmds::IGMesh::region_iterator it  = this->mesh_.regions_begin();

		for(;!it.isDone();it.next()){
			Region current_region = it.value();

			if((*volumeClassification)[current_region.getID()] != NULL)
			{
				if((*volumeClassification)[current_region.getID()] == volumes[iVolume])
				{
					vol.add(current_region);
				}
			}
		}
	}
	if(this->manager_.getNbVolumes() > 0) {
                gmds::Variable<int>* volVariable = this->mesh_.newVariable<int>(gmds::GMDS_REGION,"Volumes");

                for(int iVolume=0; iVolume<this->manager_.getNbVolumes(); iVolume++) {
                        gmds::IGMesh::region_iterator it  = this->mesh_.regions_begin();

                        for(;!it.isDone();it.next()){
                               Region current_region = it.value();

                                if((*volumeClassification)[current_region.getID()] != NULL)
                                {
                                        if((*volumeClassification)[current_region.getID()] == volumes[iVolume])
                                        {
                                                (*volVariable)[current_region.getID()] = iVolume+1; // +1 so as to discriminate
                                                                                                    // those that are not in a volume
                                        }
                                }
			}
		}
	}


	for(int iSurface=0; iSurface<this->manager_.getNbSurfaces(); iSurface++)
	{
		std::ostringstream surface_name;
		surface_name<<"surface_"<<iSurface;

		gmds::IGMesh::surface& surf = this->mesh_.newSurface(surface_name.str());

		gmds::IGMesh::face_iterator it  = this->mesh_.faces_begin();

		for(;!it.isDone();it.next()){
			Face current_face = it.value();

			if((*surfaceClassification)[current_face.getID()] != NULL)
			{
				if((*surfaceClassification)[current_face.getID()] == surfaces[iSurface])
				{
					surf.add(current_face);
				}
			}
		}
	}
	if(this->manager_.getNbSurfaces() > 0) {
		gmds::Variable<int>* surfVariable = this->mesh_.newVariable<int>(gmds::GMDS_FACE,"Surfaces");

		for(int iSurface=0; iSurface<this->manager_.getNbSurfaces(); iSurface++) {
			gmds::IGMesh::face_iterator it  = this->mesh_.faces_begin();

		        for(;!it.isDone();it.next()){
                 	       Face current_face = it.value();

                        	if((*surfaceClassification)[current_face.getID()] != NULL)
	                        {
        	                        if((*surfaceClassification)[current_face.getID()] == surfaces[iSurface])
                	                {
                        	                (*surfVariable)[current_face.getID()] = iSurface+1; // +1 so as to discriminate
												    // those that are not on a surface
                                	}
	                        }
        	        }
		}
	}
	

	for(int iLine=0; iLine<this->manager_.getNbCurves(); iLine++)
	{
		std::ostringstream line_name;
		line_name<<"line_"<<iLine;
		std::ostringstream cloud_name;
		cloud_name<<"cloud_line_"<<iLine;

		gmds::IGMesh::line& line_l = this->mesh_.newLine(line_name.str());
		gmds::IGMesh::cloud& cloud_l = this->mesh_.newCloud(cloud_name.str());

		gmds::IGMesh::edge_iterator it  = this->mesh_.edges_begin();
		int nbE=0, nbN=0;
		std::set<TCellID> line_nodes;
		for(;!it.isDone();it.next()){
			Edge current_edge = it.value();

			if((*curveClassification)[current_edge.getID()] != NULL)
			{
				if((*curveClassification)[current_edge.getID()] == curves[iLine])
				{
					line_l.add(current_edge);
					nbE++;
					std::vector<Node> current_nodes = current_edge.get<Node>();
					for(unsigned int in=0;in<current_nodes.size();in++){
						line_nodes.insert(current_nodes[in].getID());
					}
				}
			}
		}

		std::set<TCellID>::iterator it_nodes = line_nodes.begin();
		for(;it_nodes!=line_nodes.end();it_nodes++)
		{
			cloud_l.add(this->mesh_.get<Node>(*it_nodes));
			nbN++;
		}

		if(line_nodes.empty()){

			gmds::IGMesh::node_iterator it  = this->mesh_.nodes_begin();
			for(;!it.isDone();it.next()){
				Node current_node = it.value();

				if((*nodeClassification)[current_node.getID()] != NULL)
				{
					if( (*nodeClassification)[current_node.getID()]==	curves[iLine])
					{
						cloud_l.add(current_node);
						nbN++;
					}
				}
			}
		}


//		std::cout<<"line "<<iLine<<" - "<<nbE<<", "<<nbN<<std::endl;

	}
	if(this->manager_.getNbCurves() > 0) {
		
		gmds::Variable<int>* curvesVariable = this->mesh_.newVariable<int>(gmds::GMDS_NODE,"Curves");

		for(int iLine=0; iLine<this->manager_.getNbCurves(); iLine++) {
			gmds::IGMesh::edge_iterator it  = this->mesh_.edges_begin();

			for(;!it.isDone();it.next()){
                        Edge current_edge = it.value();

                        	if((*curveClassification)[current_edge.getID()] != NULL)
	                        {
        	                        if((*curveClassification)[current_edge.getID()] == curves[iLine])
                	                {
						std::vector<Node> current_nodes = current_edge.get<Node>();
                                	        for(unsigned int in=0;in<current_nodes.size();in++){
                                        	        (*curvesVariable)[current_nodes[in].getID()] = iLine+1; // +1 so as to discriminate
                                                                                                    		// those that are not on a line
                                        	}
                                	}
                        	}
	                }

		}
	}

	for(int iCloud=0; iCloud<this->manager_.getNbPoints(); iCloud++)
	{
		std::ostringstream cloud_name;
		cloud_name<<"point_"<<iCloud;

		gmds::IGMesh::cloud& cloud_l = this->mesh_.newCloud(cloud_name.str());
		gmds::IGMesh::node_iterator it  = this->mesh_.nodes_begin();
		int nbN=0;
		for(;!it.isDone();it.next()){
			Node current_node = it.value();

			if((*nodeClassification)[current_node.getID()] != NULL)
			{
				if((*nodeClassification)[current_node.getID()] == points[iCloud])
				{
					cloud_l.add(current_node);
					nbN++;
				}
			}
		}
//		std::cout<<"point "<<iCloud<<" - "<<nbN<<std::endl;

	}
        if(this->manager_.getNbPoints() > 0) {
                gmds::Variable<int>* verticesVariable = this->mesh_.newVariable<int>(gmds::GMDS_NODE,"Points");

                for(int iPoint=0; iPoint<this->manager_.getNbPoints(); iPoint++) {
                        gmds::IGMesh::node_iterator it  = this->mesh_.nodes_begin();

                        for(;!it.isDone();it.next()){
                               Node current_node = it.value();

                                if((*nodeClassification)[current_node.getID()] != NULL)
                                {
                                        if((*nodeClassification)[current_node.getID()] == points[iPoint])
                                        {
                                                (*verticesVariable)[current_node.getID()] = iPoint+1; // +1 so as to discriminate
                                                                                                      // those that are not on a surface
                                        }
				}
			}
		}
	}

	gmds::MeshModel mod = gmds::MeshModel::intersection(AMeshModel,gmds::N|gmds::E|gmds::F|gmds::R);
	int mask = 0;
	if(mod.has(gmds::N)) mask|=gmds::N;
	if(mod.has(gmds::E)) mask|=gmds::E;
	if(mod.has(gmds::F)) mask|=gmds::F;
	if(mod.has(gmds::R)) mask|=gmds::R;

	gmds::VTKWriter<gmds::IGMesh> w(mesh_);
	w.write(AFile,mask);

	// carefully remove every created groups on the mesh, else this function
	// can not be called several times.
	for(int iVolume=0; iVolume<this->manager_.getNbVolumes(); iVolume++)
	{
		std::ostringstream volume_name;
		volume_name<<"volume_"<<iVolume;

		gmds::IGMesh::volume& vol = this->mesh_.getVolume(volume_name.str());
		this->mesh_.deleteVolume(vol);
	}
	if(this->manager_.getNbVolumes() > 0) {
                this->mesh_.deleteVariable(gmds::GMDS_REGION,"Volumes");
        }

	for(int iSurface=0; iSurface<this->manager_.getNbSurfaces(); iSurface++)
	{
		std::ostringstream surface_name;
		surface_name<<"surface_"<<iSurface;

		gmds::IGMesh::surface& surf = this->mesh_.getSurface(surface_name.str());
		this->mesh_.deleteSurface(surf);
	}
	if(this->manager_.getNbSurfaces() > 0) {
		this->mesh_.deleteVariable(gmds::GMDS_FACE,"Surfaces");
	}

	for(int iLine=0; iLine<this->manager_.getNbCurves(); iLine++)
	{
		std::ostringstream line_name;
		line_name<<"line_"<<iLine;
		std::ostringstream cloud_name;
		cloud_name<<"cloud_line_"<<iLine;

		gmds::IGMesh::line& line_l = this->mesh_.getLine(line_name.str());
		this->mesh_.deleteLine(line_l);
		gmds::IGMesh::cloud& cloud_l = this->mesh_.getCloud(cloud_name.str());
		this->mesh_.deleteCloud(cloud_l);
	}
	if(this->manager_.getNbCurves() > 0) {
                this->mesh_.deleteVariable(gmds::GMDS_NODE,"Curves");
        }

	for(int iPnt=0; iPnt<this->manager_.getNbPoints(); iPnt++)
	{
		std::ostringstream cloud_name;
		cloud_name<<"point_"<<iPnt;

		gmds::IGMesh::cloud& cloud_l = this->mesh_.getCloud(cloud_name.str());
		this->mesh_.deleteCloud(cloud_l);
	}
	if(this->manager_.getNbPoints() > 0) {
                this->mesh_.deleteVariable(gmds::GMDS_NODE,"Points");
        }

	std::cout<<"end MeshInsertDetail::exporMeshtVTK "<<AFile<<std::endl;
}
/*----------------------------------------------------------------------------*/
void
MeshInsertDetail::exportModelVTK(const std::string& AFile)
{
    gmds::MeshModel mod = gmds::DIM3|gmds::N|gmds::E|gmds::F|gmds::F2N|gmds::E2N;
    gmds::IGMesh mesh(mod);

	std::vector<gmds::geom::GeomVolume* > volumes;
	std::vector<gmds::geom::GeomSurface* > surfaces;
	std::vector<gmds::geom::GeomCurve* > curves;
	std::vector<gmds::geom::GeomPoint* > points;
	this->manager_.getVolumes(volumes);
	this->manager_.getSurfaces(surfaces);
	this->manager_.getCurves(curves);
	this->manager_.getPoints(points);

	std::map<gmds::geom::GeomSurface*,std::vector<gmds::Face> > surfaces2Faces;
	std::map<gmds::geom::GeomCurve*,std::vector<gmds::Edge> > curves2Edges;

	for(int iSurface=0; iSurface<this->manager_.getNbSurfaces(); iSurface++)
	{
		std::ostringstream surface_name;
		surface_name<<"surface_"<<iSurface;

		gmds::IGMesh::surface& surf = mesh.newSurface(surface_name.str());

		std::vector<gmds::math::Triangle> triangles;
		surfaces[iSurface]->getTriangulation(triangles);


		for(int iTriangle=0; iTriangle<triangles.size(); iTriangle++) {
			gmds::Node n0 = mesh.newNode(triangles[iTriangle].getPoint(0));
			gmds::Node n1 = mesh.newNode(triangles[iTriangle].getPoint(1));
			gmds::Node n2 = mesh.newNode(triangles[iTriangle].getPoint(2));
			gmds::Face face = mesh.newTriangle(n0,n1,n2);
			surfaces2Faces[surfaces[iSurface]].push_back(face);

			surf.add(face);
		}

	}

	for(int iCurve=0; iCurve<this->manager_.getNbCurves(); iCurve++)
	{
		std::ostringstream line_name;
		line_name<<"line_"<<iCurve;

		gmds::IGMesh::cloud& line_l = mesh.newCloud(line_name.str());

		// arbitrary number of points that will represent a curve
		const int nbPnts = 10;
		gmds::math::Point pnts[nbPnts];
		curves[iCurve]->getMultiplePoints(nbPnts,pnts);

		for(int iPnt=0; iPnt<10; iPnt++) {
			gmds::Node node = mesh.newNode(pnts[iPnt]);
			line_l.add(node);
		}

//		std::vector<gmds::Edge> edges;
//		gmds::geom::FacetedCurve* facetedCurve = dynamic_cast<gmds::geom::FacetedCurve*>(curves[iCurve]);
//		if(facetedCurve == NULL) {
//			throw GMDSException("MeshInsertDetail::exportModelVTK GeomCurve not a FacetedCurve.");
//		}
//
//		facetedCurve->getMeshEdges(edges);
//
//		for(int iEdge=0; iEdge<edges.size(); iEdge++) {
//
//			std::vector<gmds::Node> nodes = edges[iEdge].get<Node>();
//			gmds::Node n0 = mesh.newNode(nodes[0].getPoint());
//			gmds::Node n1 = mesh.newNode(nodes[1].getPoint());
//			gmds::Edge edge = mesh.newEdge(n0,n1);
//			curves2Edges[curves[iCurve]].push_back(edge);
//
//			line_l.add(edge);
//		}

	}

	for(int iVolume=0; iVolume<this->manager_.getNbVolumes(); iVolume++)
	{
		std::ostringstream  volume_name;
		volume_name<<"volume_"<<iVolume;

		gmds::IGMesh::surface& vol = mesh.newSurface(volume_name.str());

		std::vector<gmds::geom::GeomSurface* > surfaces;
		volumes[iVolume]->get(surfaces);

		for(int iSurf=0; iSurf<surfaces.size(); iSurf++) {
			for(int iFace=0; iFace<surfaces2Faces[surfaces[iSurf]].size(); iFace++) {
				vol.add(surfaces2Faces[surfaces[iSurf]][iFace]);
			}
		}
	}

	for(int iPoint=0; iPoint<this->manager_.getNbPoints(); iPoint++)
	{
		std::ostringstream cloud_name;
		cloud_name<<"point_"<<iPoint;

		gmds::IGMesh::cloud& cloud_l = mesh.newCloud(cloud_name.str());

		gmds::Node node = mesh.newNode(points[iPoint]->getPoint());
		cloud_l.add(node);
	}

	gmds::VTKWriter<gmds::IGMesh> w(mesh);
	w.write(AFile,gmds::F|gmds::E|gmds::N);

	for(int iVolume=0; iVolume<this->manager_.getNbVolumes(); iVolume++)
	{
		std::ostringstream volume_name;
		volume_name<<"volume_"<<iVolume;

		gmds::IGMesh::surface& vol = mesh.getSurface(volume_name.str());
		this->mesh_.deleteSurface(vol);
	}

	for(int iSurface=0; iSurface<this->manager_.getNbSurfaces(); iSurface++)
	{
		std::ostringstream surface_name;
		surface_name<<"surface_"<<iSurface;

		gmds::IGMesh::surface& surf = mesh.getSurface(surface_name.str());
		this->mesh_.deleteSurface(surf);
	}

	for(int iLine=0; iLine<this->manager_.getNbCurves(); iLine++)
	{
		std::ostringstream line_name;
		line_name<<"line_"<<iLine;

		gmds::IGMesh::cloud& line_l = mesh.getCloud(line_name.str());
		this->mesh_.deleteCloud(line_l);

//		gmds::IGMesh::line& line_l = mesh.getLine(line_name.str());
//		this->mesh_.deleteLine(line_l);
	}

	for(int iPnt=0; iPnt<this->manager_.getNbPoints(); iPnt++)
	{
		std::ostringstream cloud_name;
		cloud_name<<"point_"<<iPnt;

		gmds::IGMesh::cloud& cloud_l = mesh.getCloud(cloud_name.str());
		this->mesh_.deleteCloud(cloud_l);
	}

}
/*----------------------------------------------------------------------------*/
void
MeshInsertDetail::pillow(gmds::geom::GeomVolume* AVol)
{
	std::vector<gmds::geom::GeomSurface* > surfaces;
	AVol->get(surfaces);

	// first get the list of faces we want to pillow
	// and any additional data required for the 
	// sheet operator
	Variable<geom::GeomEntity* >* surfaceClassification= this->mesh_.getGeometricClassification(2);
	Variable<geom::GeomEntity* >* volumeClassification= this->mesh_.getGeometricClassification(3);
	
	std::set<gmds::Region> regions2PillowSet;
	std::vector<gmds::Face> faces2Pillow;
	std::vector<gmds::Region> newRegions;

	for(int iSurface=0; iSurface<surfaces.size(); iSurface++) {
		
		gmds::IGMesh::face_iterator it  = this->mesh_.faces_begin();
                for(;!it.isDone();it.next()) {
                        gmds::Face current_face = it.value();
			
			if((*surfaceClassification)[current_face.getID()] == surfaces[iSurface]) {
				faces2Pillow.push_back(current_face);

				std::vector<gmds::Node> nodes = current_face.get<gmds::Node>();

				for(int iNode=0; iNode<nodes.size(); iNode++) {
					gmds::Node current_node = nodes[iNode];
					std::vector<gmds::Region> regions = current_node.get<gmds::Region>();
					for(int iRegion=0; iRegion<regions.size(); iRegion++) {
						if((*volumeClassification)[regions[iRegion].getID()] == AVol) {
							regions2PillowSet.insert(regions[iRegion]);
						}
					}
				}
			}
		}	
	}

	std::vector<gmds::Region> regions2Pillow;
	std::set<gmds::Region>::iterator itr = regions2PillowSet.begin();
	while(itr != regions2PillowSet.end()) {
		regions2Pillow.push_back(*itr);
		itr++;
	}

	gmds::SheetOperator sheetOp(this->mesh_);
	sheetOp.pillow(regions2Pillow,faces2Pillow,newRegions,true);

	// associate the newly created internal nodes after pillowing
	Variable<geom::GeomEntity* >* nodesClassification= this->mesh_.getGeometricClassification(0);
	std::set<gmds::Node> newCreatedNodes = sheetOp.getNewInternalNodes();
	for(std::set<gmds::Node>::iterator it = newCreatedNodes.begin(); it != newCreatedNodes.end(); it++) {
		(*nodesClassification)[(*it).getID()] = AVol;
	}

	// put the regions2Pillow and newRegions in variables in order to debug
	Variable<int>* regions2PillowVar = this->mesh_.newVariable<int>(GMDS_REGION,"regions2Pillow");
	Variable<int>* newRegionsPillowVar = this->mesh_.newVariable<int>(GMDS_REGION,"newRegionsPillow");
	for(int iRegion=0; iRegion<regions2Pillow.size(); iRegion++) {
		(*regions2PillowVar)[regions2Pillow[iRegion].getID()] = 1;
	}

	Variable<geom::GeomEntity* >* nodeClassification= this->mesh_.getGeometricClassification(0);
	for(int iRegion=0; iRegion<newRegions.size(); iRegion++) {
                (*newRegionsPillowVar)[newRegions[iRegion].getID()] = 1;
        }	
}
/*----------------------------------------------------------------------------*/
void
MeshInsertDetail::pillowExt(gmds::geom::GeomVolume* AVol)
{
	std::vector<gmds::geom::GeomSurface* > surfaces;
        AVol->get(surfaces);
	
	// first get the list of faces we want to pillow
	// and any additional data required for the
	// sheet operator
        Variable<geom::GeomEntity* >* surfaceClassification= this->mesh_.getGeometricClassification(2);
        Variable<geom::GeomEntity* >* volumeClassification= this->mesh_.getGeometricClassification(3);

	std::set<gmds::Region> regions2PillowSet;
        std::vector<gmds::Face> faces2Pillow;
        std::vector<gmds::Region> newRegions;

        for(int iSurface=0; iSurface<surfaces.size(); iSurface++) {

                gmds::IGMesh::face_iterator it  = this->mesh_.faces_begin();
                for(;!it.isDone();it.next()) {
                        gmds::Face current_face = it.value();

                        if((*surfaceClassification)[current_face.getID()] == surfaces[iSurface]) {
                                faces2Pillow.push_back(current_face);

				std::vector<gmds::Node> nodes = current_face.get<gmds::Node>();

                                for(int iNode=0; iNode<nodes.size(); iNode++) {
                                        gmds::Node current_node = nodes[iNode];
                                        std::vector<gmds::Region> regions = current_node.get<gmds::Region>();
                                        for(int iRegion=0; iRegion<regions.size(); iRegion++) {
                                                if((*volumeClassification)[regions[iRegion].getID()] != AVol) {
                                                        regions2PillowSet.insert(regions[iRegion]);
                                                }
                                        }
                                }
                        }
                }
        }

	std::vector<gmds::Region> regions2Pillow;
        std::set<gmds::Region>::iterator itr = regions2PillowSet.begin();
        while(itr != regions2PillowSet.end()) {
                regions2Pillow.push_back(*itr);
                itr++;
        }

        gmds::SheetOperator sheetOp(this->mesh_);
        sheetOp.pillow(regions2Pillow,faces2Pillow,newRegions,true);

	// associate the newly created external nodes after pillowing
        Variable<geom::GeomEntity* >* nodesClassification= this->mesh_.getGeometricClassification(0);
        std::set<gmds::Node> newCreatedNodes = sheetOp.getNewInternalNodes();
        for(std::set<gmds::Node>::iterator it = newCreatedNodes.begin(); it != newCreatedNodes.end(); it++) {
                (*nodesClassification)[(*it).getID()] = NULL;
        }

	// put the regions2Pillow and newRegions in variables in order to debug
        Variable<int>* regions2PillowVar = this->mesh_.newVariable<int>(GMDS_REGION,"regions2PillowEXT");
        Variable<int>* newRegionsPillowVar = this->mesh_.newVariable<int>(GMDS_REGION,"newRegionsPillowEXT");
        for(int iRegion=0; iRegion<regions2Pillow.size(); iRegion++) {
                (*regions2PillowVar)[regions2Pillow[iRegion].getID()] = 1;
        }

        Variable<geom::GeomEntity* >* nodeClassification= this->mesh_.getGeometricClassification(0);
        for(int iRegion=0; iRegion<newRegions.size(); iRegion++) {
                (*newRegionsPillowVar)[newRegions[iRegion].getID()] = 1;
        }
}
/*----------------------------------------------------------------------------*/
void
MeshInsertDetail::pillowIntWithoutVerticesCurves(gmds::geom::GeomVolume* AVol)
{
        std::vector<gmds::geom::GeomSurface* > surfaces;
        AVol->get(surfaces);

	Variable<geom::GeomEntity* >* nodesClassification= this->mesh_.getGeometricClassification(0);
        Variable<geom::GeomEntity* >* surfaceClassification= this->mesh_.getGeometricClassification(2);
        Variable<geom::GeomEntity* >* volumeClassification= this->mesh_.getGeometricClassification(3);

        std::vector<gmds::Region> regions2Pillow;
        std::vector<gmds::Face> faces2Pillow;
        std::vector<gmds::Region> newRegions;

	// first get the shrink set of regions for the pillowing
        gmds::IGMesh::region_iterator itr  = this->mesh_.regions_begin();
        for(;!itr.isDone();itr.next()) {
                gmds::Region current_region = itr.value();

                if((*volumeClassification)[current_region.getID()] == AVol) {

			// do not keep the regions adjacent to a vertex or a curve
			std::vector<gmds::TCellID> nodesIDs = current_region.getIDs<gmds::Node> ();
			bool toKeep = true;
			for(int iNode=0; iNode<nodesIDs.size(); iNode++) {
				if(((*nodesClassification)[nodesIDs[iNode]] != NULL) && ((*nodesClassification)[nodesIDs[iNode]])->getDim() < 2) {
					toKeep = false;
					break;
				}
			}

			// but keep it when region has two or more faces on the same surface.
			// Since they are associated to a surface, they actually exist
			if (!toKeep) {
				std::vector<gmds::TCellID> facesIDs = current_region.getIDs<gmds::Face> ();
				std::set<gmds::geom::GeomEntity*> surfacesFound;
				for(int iFace=0; iFace<facesIDs.size(); iFace++) {
					if((*surfaceClassification)[facesIDs[iFace]] != NULL) {
						if(surfacesFound.find((*surfaceClassification)[facesIDs[iFace]]) != surfacesFound.end()) {
							toKeep = true;
							break;
						} else {
							surfacesFound.insert((*surfaceClassification)[facesIDs[iFace]]);
						}
					}
				}
			}
			
			if(toKeep) {
                        	regions2Pillow.push_back(current_region);
			}
                }
        }

	// get the set of fakefaces (the skin) necessary to pass as parameter to the pillowing
	std::map<gmds::FakeFace::FaceID,int> facesOccurence;
	std::map<gmds::FakeFace::FaceID,gmds::Region> faces2Region;
	for(int iRegion=0; iRegion<regions2Pillow.size(); iRegion++) {
		gmds::Region current_region = regions2Pillow[iRegion];

		std::vector<std::vector<Node> > nodesFaces = current_region.getOrderedNodesFaces();
		
		for(int iFace=0; iFace<nodesFaces.size(); iFace++) {
			std::vector<gmds::TCellID> nodesIDs;

			for(int iNode=0; iNode<nodesFaces[iFace].size(); iNode++) {
				nodesIDs.push_back(nodesFaces[iFace][iNode].getID());
			}
			FakeFace f(nodesIDs);

			if(facesOccurence.find(f.getID()) == facesOccurence.end()) {
				facesOccurence[f.getID()] = 1;
			} else {
				facesOccurence[f.getID()]++;
			}

			// it is not strictly consistent but as we plan only on creating fakefaces that appear only once,
			// it is OK.
			faces2Region[f.getID()] = current_region;
		}
	}
	
	// now we will build the faces that do not already exist
	// R2F F2R adjacencies have to be consistent

	// first get all the existing faces
	std::map<gmds::FakeFace::FaceID, gmds::Face> fakeFace2Face;
	gmds::IGMesh::face_iterator itf  = this->mesh_.faces_begin();
        for(;!itf.isDone();itf.next()) {
                gmds::Face current_face = itf.value();
		std::vector<gmds::TCellID> nodesIDs = current_face.getIDs<gmds::Node> ();

		FakeFace f(nodesIDs);
		fakeFace2Face[f.getID()] = current_face;
	}

	// create the relevant faces
	std::map<gmds::FakeFace::FaceID,int>::iterator itff = facesOccurence.begin();
	for(; itff != facesOccurence.end(); itff++) {

		gmds::FakeFace::FaceID fID = itff->first;
		if(facesOccurence[fID] == 1) {
			if(fakeFace2Face.find(fID) == fakeFace2Face.end()) {
				gmds::Face newF = this->mesh_.newFace(fID.getNodeIDs());
				faces2Pillow.push_back(newF);

				faces2Region[fID].add(newF);
				newF.add(faces2Region[fID]);
			} else {
				faces2Pillow.push_back(fakeFace2Face[fID]);
			}
		}
	}

std::cout<<"POYOP "<<regions2Pillow.size()<<std::endl;
std::cout<<"POYOP "<<faces2Pillow.size()<<std::endl;

	// pillow
	gmds::SheetOperator sheetOp(this->mesh_);
        sheetOp.pillow(regions2Pillow,faces2Pillow,newRegions,true);

        // associate the newly created external nodes after pillowing
        std::set<gmds::Node> newCreatedNodes = sheetOp.getNewInternalNodes();
        for(std::set<gmds::Node>::iterator it = newCreatedNodes.begin(); it != newCreatedNodes.end(); it++) {
                (*nodesClassification)[(*it).getID()] = NULL;
        }

        // put the regions2Pillow and newRegions in variables in order to debug
        Variable<int>* regions2PillowVar = this->mesh_.newVariable<int>(GMDS_REGION,"regions2PillowINT");
        Variable<int>* newRegionsPillowVar = this->mesh_.newVariable<int>(GMDS_REGION,"newRegionsPillowINT");
        for(int iRegion=0; iRegion<regions2Pillow.size(); iRegion++) {
                (*regions2PillowVar)[regions2Pillow[iRegion].getID()] = 1;
        }

        for(int iRegion=0; iRegion<newRegions.size(); iRegion++) {
                (*newRegionsPillowVar)[newRegions[iRegion].getID()] = 1;
        }
}
/*----------------------------------------------------------------------------*/
void
MeshInsertDetail::pillowExtWithoutVerticesCurves(gmds::geom::GeomVolume* AVol)
{
        std::vector<gmds::geom::GeomSurface* > surfaces;
        AVol->get(surfaces);

	Variable<geom::GeomEntity* >* nodesClassification= this->mesh_.getGeometricClassification(0);
        Variable<geom::GeomEntity* >* surfaceClassification= this->mesh_.getGeometricClassification(2);
        Variable<geom::GeomEntity* >* volumeClassification= this->mesh_.getGeometricClassification(3);

        std::vector<gmds::Region> regions2Pillow;
        std::vector<gmds::Face> faces2Pillow;
        std::vector<gmds::Region> newRegions;

	// first get the shrink set of regions for the pillowing
        gmds::IGMesh::region_iterator itr  = this->mesh_.regions_begin();
        for(;!itr.isDone();itr.next()) {
                gmds::Region current_region = itr.value();

                if((*volumeClassification)[current_region.getID()] != AVol) {

			// do not keep the regions adjacent to a vertex or a curve
			std::vector<gmds::TCellID> nodesIDs = current_region.getIDs<gmds::Node> ();
			bool toKeep = true;
			for(int iNode=0; iNode<nodesIDs.size(); iNode++) {
				if(((*nodesClassification)[nodesIDs[iNode]] != NULL) && ((*nodesClassification)[nodesIDs[iNode]])->getDim() < 2) {
					toKeep = false;
					break;
				}
			}

			// but keep it when region has two or more faces on the same surface.
			// Since they are associated to a surface, they actually exist
			if (!toKeep) {
				std::vector<gmds::TCellID> facesIDs = current_region.getIDs<gmds::Face> ();
				std::set<gmds::geom::GeomEntity*> surfacesFound;
				for(int iFace=0; iFace<facesIDs.size(); iFace++) {
					if((*surfaceClassification)[facesIDs[iFace]] != NULL) {
						if(surfacesFound.find((*surfaceClassification)[facesIDs[iFace]]) != surfacesFound.end()) {
							toKeep = true;
							break;
						} else {
							surfacesFound.insert((*surfaceClassification)[facesIDs[iFace]]);
						}
					}
				}
			}
			
			if(toKeep) {
                        	regions2Pillow.push_back(current_region);
			}
                }
        }

	// get the set of fakefaces (the skin) necessary to pass as parameter to the pillowing
	std::map<gmds::FakeFace::FaceID,int> facesOccurence;
	std::map<gmds::FakeFace::FaceID,gmds::Region> faces2Region;
	for(int iRegion=0; iRegion<regions2Pillow.size(); iRegion++) {
		gmds::Region current_region = regions2Pillow[iRegion];

		std::vector<std::vector<Node> > nodesFaces = current_region.getOrderedNodesFaces();
		
		for(int iFace=0; iFace<nodesFaces.size(); iFace++) {
			std::vector<gmds::TCellID> nodesIDs;

			for(int iNode=0; iNode<nodesFaces[iFace].size(); iNode++) {
				nodesIDs.push_back(nodesFaces[iFace][iNode].getID());
			}
			FakeFace f(nodesIDs);

			if(facesOccurence.find(f.getID()) == facesOccurence.end()) {
				facesOccurence[f.getID()] = 1;
			} else {
				facesOccurence[f.getID()]++;
			}

			// it is not strictly consistent but as we plan only on creating fakefaces that appear only once,
			// it is OK.
			faces2Region[f.getID()] = current_region;
		}
	}
	
	// now we will build the faces that do not already exist
	// R2F F2R adjacencies have to be consistent

	// first get all the existing faces
	std::map<gmds::FakeFace::FaceID, gmds::Face> fakeFace2Face;
	gmds::IGMesh::face_iterator itf  = this->mesh_.faces_begin();
        for(;!itf.isDone();itf.next()) {
                gmds::Face current_face = itf.value();
		std::vector<gmds::TCellID> nodesIDs = current_face.getIDs<gmds::Node> ();

		FakeFace f(nodesIDs);
		fakeFace2Face[f.getID()] = current_face;
	}

	// create the relevant faces
	std::map<gmds::FakeFace::FaceID,int>::iterator itff = facesOccurence.begin();
	for(; itff != facesOccurence.end(); itff++) {

		gmds::FakeFace::FaceID fID = itff->first;
		if(facesOccurence[fID] == 1) {
			if(fakeFace2Face.find(fID) == fakeFace2Face.end()) {
				gmds::Face newF = this->mesh_.newFace(fID.getNodeIDs());
				faces2Pillow.push_back(newF);

				faces2Region[fID].add(newF);
				newF.add(faces2Region[fID]);
			} else {
				faces2Pillow.push_back(fakeFace2Face[fID]);
			}
		}
	}

std::cout<<"POYOP "<<regions2Pillow.size()<<std::endl;
std::cout<<"POYOP "<<faces2Pillow.size()<<std::endl;

	// pillow	
	gmds::SheetOperator sheetOp(this->mesh_);
        sheetOp.pillow(regions2Pillow,faces2Pillow,newRegions,true);

        // associate the newly created external nodes after pillowing
        std::set<gmds::Node> newCreatedNodes = sheetOp.getNewInternalNodes();
        for(std::set<gmds::Node>::iterator it = newCreatedNodes.begin(); it != newCreatedNodes.end(); it++) {
                (*nodesClassification)[(*it).getID()] = NULL;
        }

        // put the regions2Pillow and newRegions in variables in order to debug
        Variable<int>* regions2PillowVar = this->mesh_.newVariable<int>(GMDS_REGION,"regions2PillowEXT");
        Variable<int>* newRegionsPillowVar = this->mesh_.newVariable<int>(GMDS_REGION,"newRegionsPillowEXT");
        for(int iRegion=0; iRegion<regions2Pillow.size(); iRegion++) {
                (*regions2PillowVar)[regions2Pillow[iRegion].getID()] = 1;
        }

        for(int iRegion=0; iRegion<newRegions.size(); iRegion++) {
                (*newRegionsPillowVar)[newRegions[iRegion].getID()] = 1;
        }

}
/*----------------------------------------------------------------------------*/
void
MeshInsertDetail::pillowIntWithoutVerticesCurves_bis(gmds::geom::GeomVolume* AVol)
{
        std::vector<gmds::geom::GeomSurface* > surfaces;
        AVol->get(surfaces);

	Variable<geom::GeomEntity* >* nodesClassification= this->mesh_.getGeometricClassification(0);
        Variable<geom::GeomEntity* >* surfaceClassification= this->mesh_.getGeometricClassification(2);
        Variable<geom::GeomEntity* >* volumeClassification= this->mesh_.getGeometricClassification(3);

        std::vector<gmds::Region> regions2Pillow;
        std::vector<gmds::Face> faces2Pillow;
        std::vector<gmds::Region> newRegions;

	// first get the shrink set of regions for the pillowing
        gmds::IGMesh::region_iterator itr  = this->mesh_.regions_begin();
        for(;!itr.isDone();itr.next()) {
                gmds::Region current_region = itr.value();

                if((*volumeClassification)[current_region.getID()] == AVol) {

			// do not keep the regions adjacent to a vertex or a curve
			std::vector<gmds::TCellID> nodesIDs = current_region.getIDs<gmds::Node> ();
			bool toKeep = true;
			for(int iNode=0; iNode<nodesIDs.size(); iNode++) {
				if(((*nodesClassification)[nodesIDs[iNode]] != NULL) && ((*nodesClassification)[nodesIDs[iNode]])->getDim() < 2) {
					toKeep = false;
					break;
				}
			}

			// but keep it when region has two or more faces on the same surface.
			// Since they are associated to a surface, they actually exist
			if (!toKeep) {
				std::vector<gmds::TCellID> facesIDs = current_region.getIDs<gmds::Face> ();
				std::set<gmds::geom::GeomEntity*> surfacesFound;
				for(int iFace=0; iFace<facesIDs.size(); iFace++) {
					if((*surfaceClassification)[facesIDs[iFace]] != NULL) {
						if(surfacesFound.find((*surfaceClassification)[facesIDs[iFace]]) != surfacesFound.end()) {
							toKeep = true;
							break;
						} else {
							surfacesFound.insert((*surfaceClassification)[facesIDs[iFace]]);
						}
					}
				}
			}
			
			if(toKeep) {
                        	regions2Pillow.push_back(current_region);
			}
                }
        }


std::cout<<"POYOP "<<regions2Pillow.size()<<std::endl;
std::cout<<"POYOP "<<faces2Pillow.size()<<std::endl;

	// pillow
	gmds::Pillowing op(this->mesh_);
        op.pillow(faces2Pillow,regions2Pillow,newRegions,true);

        // associate the newly created external nodes after pillowing
        std::set<gmds::Node> newCreatedNodes = op.getNewInternalNodes();
        for(std::set<gmds::Node>::iterator it = newCreatedNodes.begin(); it != newCreatedNodes.end(); it++) {
                (*nodesClassification)[(*it).getID()] = NULL;
        }

        // put the regions2Pillow and newRegions in variables in order to debug
        Variable<int>* regions2PillowVar = this->mesh_.newVariable<int>(GMDS_REGION,"regions2PillowINT");
        Variable<int>* newRegionsPillowVar = this->mesh_.newVariable<int>(GMDS_REGION,"newRegionsPillowINT");
        for(int iRegion=0; iRegion<regions2Pillow.size(); iRegion++) {
                (*regions2PillowVar)[regions2Pillow[iRegion].getID()] = 1;
        }

        for(int iRegion=0; iRegion<newRegions.size(); iRegion++) {
                (*newRegionsPillowVar)[newRegions[iRegion].getID()] = 1;
        }
}
/*----------------------------------------------------------------------------*/
void
MeshInsertDetail::pillowExtWithoutVerticesCurves_bis(gmds::geom::GeomVolume* AVol)
{
        std::vector<gmds::geom::GeomSurface* > surfaces;
        AVol->get(surfaces);

	Variable<geom::GeomEntity* >* nodesClassification= this->mesh_.getGeometricClassification(0);
        Variable<geom::GeomEntity* >* surfaceClassification= this->mesh_.getGeometricClassification(2);
        Variable<geom::GeomEntity* >* volumeClassification= this->mesh_.getGeometricClassification(3);

        std::vector<gmds::Region> regions2Pillow;
        std::vector<gmds::Face> faces2Pillow;
        std::vector<gmds::Region> newRegions;

	// first get the shrink set of regions for the pillowing
        gmds::IGMesh::region_iterator itr  = this->mesh_.regions_begin();
        for(;!itr.isDone();itr.next()) {
                gmds::Region current_region = itr.value();

                if((*volumeClassification)[current_region.getID()] != AVol) {

			// do not keep the regions adjacent to a vertex or a curve
			std::vector<gmds::TCellID> nodesIDs = current_region.getIDs<gmds::Node> ();
			bool toKeep = true;
			for(int iNode=0; iNode<nodesIDs.size(); iNode++) {
				if(((*nodesClassification)[nodesIDs[iNode]] != NULL) && ((*nodesClassification)[nodesIDs[iNode]])->getDim() < 2) {
					toKeep = false;
					break;
				}
			}

			// but keep it when region has two or more faces on the same surface.
			// Since they are associated to a surface, they actually exist
			if (!toKeep) {
				std::vector<gmds::TCellID> facesIDs = current_region.getIDs<gmds::Face> ();
				std::set<gmds::geom::GeomEntity*> surfacesFound;
				for(int iFace=0; iFace<facesIDs.size(); iFace++) {
					if((*surfaceClassification)[facesIDs[iFace]] != NULL) {
						if(surfacesFound.find((*surfaceClassification)[facesIDs[iFace]]) != surfacesFound.end()) {
							toKeep = true;
							break;
						} else {
							surfacesFound.insert((*surfaceClassification)[facesIDs[iFace]]);
						}
					}
				}
			}
			
			if(toKeep) {
                        	regions2Pillow.push_back(current_region);
			}
                }
        }


std::cout<<"POYOP "<<regions2Pillow.size()<<std::endl;
std::cout<<"POYOP "<<faces2Pillow.size()<<std::endl;

	// pillow	
	gmds::Pillowing op(this->mesh_);
        op.pillow(faces2Pillow,regions2Pillow,newRegions,true);

        // associate the newly created external nodes after pillowing
        std::set<gmds::Node> newCreatedNodes = op.getNewInternalNodes();
        for(std::set<gmds::Node>::iterator it = newCreatedNodes.begin(); it != newCreatedNodes.end(); it++) {
                (*nodesClassification)[(*it).getID()] = NULL;
        }

        // put the regions2Pillow and newRegions in variables in order to debug
        Variable<int>* regions2PillowVar = this->mesh_.newVariable<int>(GMDS_REGION,"regions2PillowEXT");
        Variable<int>* newRegionsPillowVar = this->mesh_.newVariable<int>(GMDS_REGION,"newRegionsPillowEXT");
        for(int iRegion=0; iRegion<regions2Pillow.size(); iRegion++) {
                (*regions2PillowVar)[regions2Pillow[iRegion].getID()] = 1;
        }

        for(int iRegion=0; iRegion<newRegions.size(); iRegion++) {
                (*newRegionsPillowVar)[newRegions[iRegion].getID()] = 1;
        }

}
/*----------------------------------------------------------------------------*/
void
MeshInsertDetail::pillowIntSolitaryRegions(gmds::geom::GeomVolume* AVol)
{
        std::vector<gmds::geom::GeomSurface* > surfaces;
        AVol->get(surfaces);

	Variable<geom::GeomEntity* >* nodesClassification= this->mesh_.getGeometricClassification(0);
        Variable<geom::GeomEntity* >* surfaceClassification= this->mesh_.getGeometricClassification(2);
        Variable<geom::GeomEntity* >* volumeClassification= this->mesh_.getGeometricClassification(3);

        std::vector<gmds::Region> regions2Pillow;
        std::vector<gmds::Face> faces2Pillow;
        std::vector<gmds::Region> newRegions;

	// first get the shrink set of regions for the pillowing
        gmds::IGMesh::region_iterator itr  = this->mesh_.regions_begin();
        for(;!itr.isDone();itr.next()) {
                gmds::Region current_region = itr.value();

                if((*volumeClassification)[current_region.getID()] == AVol) {

			// keep it when region has two or more faces on the same surface.
			// Since they are associated to a surface, they actually exist
		  bool toKeep = false;
				std::vector<gmds::TCellID> facesIDs = current_region.getIDs<gmds::Face> ();
				std::set<gmds::geom::GeomEntity*> surfacesFound;
				for(int iFace=0; iFace<facesIDs.size(); iFace++) {
					if((*surfaceClassification)[facesIDs[iFace]] != NULL) {
						if(surfacesFound.find((*surfaceClassification)[facesIDs[iFace]]) != surfacesFound.end()) {
							toKeep = true;
							break;
						} else {
							surfacesFound.insert((*surfaceClassification)[facesIDs[iFace]]);
						}
					}
				}
			
			
			if(toKeep) {
                        	regions2Pillow.push_back(current_region);
			}
                }
        }

std::cout<<"POYOP "<<regions2Pillow.size()<<std::endl;
std::cout<<"POYOP "<<faces2Pillow.size()<<std::endl;

	// pillow
	gmds::Pillowing op(this->mesh_);
        op.pillow(faces2Pillow,regions2Pillow,newRegions,true);

        // associate the newly created external nodes after pillowing
        std::set<gmds::Node> newCreatedNodes = op.getNewInternalNodes();
        for(std::set<gmds::Node>::iterator it = newCreatedNodes.begin(); it != newCreatedNodes.end(); it++) {
                (*nodesClassification)[(*it).getID()] = NULL;
        }

        // put the regions2Pillow and newRegions in variables in order to debug
        Variable<int>* regions2PillowVar = this->mesh_.newVariable<int>(GMDS_REGION,"regions2PillowINT");
        Variable<int>* newRegionsPillowVar = this->mesh_.newVariable<int>(GMDS_REGION,"newRegionsPillowINT");
        for(int iRegion=0; iRegion<regions2Pillow.size(); iRegion++) {
                (*regions2PillowVar)[regions2Pillow[iRegion].getID()] = 1;
        }

        for(int iRegion=0; iRegion<newRegions.size(); iRegion++) {
                (*newRegionsPillowVar)[newRegions[iRegion].getID()] = 1;
        }
}
/*----------------------------------------------------------------------------*/
void
MeshInsertDetail::pillowExtSolitaryRegions(gmds::geom::GeomVolume* AVol)
{
        std::vector<gmds::geom::GeomSurface* > surfaces;
        AVol->get(surfaces);

	Variable<geom::GeomEntity* >* nodesClassification= this->mesh_.getGeometricClassification(0);
        Variable<geom::GeomEntity* >* surfaceClassification= this->mesh_.getGeometricClassification(2);
        Variable<geom::GeomEntity* >* volumeClassification= this->mesh_.getGeometricClassification(3);

        std::vector<gmds::Region> regions2Pillow;
        std::vector<gmds::Face> faces2Pillow;
        std::vector<gmds::Region> newRegions;

	// first get the shrink set of regions for the pillowing
        gmds::IGMesh::region_iterator itr  = this->mesh_.regions_begin();
        for(;!itr.isDone();itr.next()) {
                gmds::Region current_region = itr.value();

                if((*volumeClassification)[current_region.getID()] != AVol) {

			// keep it when region has two or more faces on the same surface.
			// Since they are associated to a surface, they actually exist
		  bool toKeep = false;
				std::vector<gmds::TCellID> facesIDs = current_region.getIDs<gmds::Face> ();
				std::set<gmds::geom::GeomEntity*> surfacesFound;
				for(int iFace=0; iFace<facesIDs.size(); iFace++) {
					if((*surfaceClassification)[facesIDs[iFace]] != NULL) {
						if(surfacesFound.find((*surfaceClassification)[facesIDs[iFace]]) != surfacesFound.end()) {
							toKeep = true;
							break;
						} else {
							surfacesFound.insert((*surfaceClassification)[facesIDs[iFace]]);
						}
					}
				}
			
			
			if(toKeep) {
                        	regions2Pillow.push_back(current_region);
			}
                }
        }

std::cout<<"POYOP "<<regions2Pillow.size()<<std::endl;
std::cout<<"POYOP "<<faces2Pillow.size()<<std::endl;

	// pillow
	gmds::Pillowing op(this->mesh_);
        op.pillow(faces2Pillow,regions2Pillow,newRegions,true);

        // associate the newly created external nodes after pillowing
        std::set<gmds::Node> newCreatedNodes = op.getNewInternalNodes();
        for(std::set<gmds::Node>::iterator it = newCreatedNodes.begin(); it != newCreatedNodes.end(); it++) {
                (*nodesClassification)[(*it).getID()] = NULL;
        }

        // put the regions2Pillow and newRegions in variables in order to debug
        Variable<int>* regions2PillowVar = this->mesh_.newVariable<int>(GMDS_REGION,"regions2PillowEXT");
        Variable<int>* newRegionsPillowVar = this->mesh_.newVariable<int>(GMDS_REGION,"newRegionsPillowEXT");
        for(int iRegion=0; iRegion<regions2Pillow.size(); iRegion++) {
                (*regions2PillowVar)[regions2Pillow[iRegion].getID()] = 1;
        }

        for(int iRegion=0; iRegion<newRegions.size(); iRegion++) {
                (*newRegionsPillowVar)[newRegions[iRegion].getID()] = 1;
        }
}
/*----------------------------------------------------------------------------*/
void
MeshInsertDetail::removeUnassociatedFaces()
{
	// first list all the faces to remove
	Variable<geom::GeomEntity* >* surfaceClassification= this->mesh_.getGeometricClassification(2);	

	std::set<gmds::Face> faces2Remove;
	gmds::IGMesh::face_iterator itf  = this->mesh_.faces_begin();
        for(;!itf.isDone();itf.next()) {
        	gmds::Face current_face = itf.value();
		
		if((*surfaceClassification)[current_face.getID()] == NULL) {
			faces2Remove.insert(current_face);
		}
	}


	// remove adjacencies that contain any of those faces
	if(this->mesh_.getModel().has(R2F)) {
		gmds::IGMesh::region_iterator itr  = this->mesh_.regions_begin();
	        for(;!itr.isDone();itr.next()) {
        	        gmds::Region current_region = itr.value();
			
			std::vector<gmds::Face> faces = current_region.get<gmds::Face>();
			for(int iFace=0; iFace<faces.size(); iFace++) {
				if(faces2Remove.find(faces[iFace]) != faces2Remove.end()) {
					current_region.remove<gmds::Face>(faces[iFace]);
				}
			}
		}
	}	
	
	if(this->mesh_.getModel().has(F2F)) {
                gmds::IGMesh::face_iterator itf  = this->mesh_.faces_begin();
                for(;!itf.isDone();itf.next()) {
                        gmds::Face current_face = itf.value();

                        std::vector<gmds::Face> faces = current_face.get<gmds::Face>();
                        for(int iFace=0; iFace<faces.size(); iFace++) {
                                if(faces2Remove.find(faces[iFace]) != faces2Remove.end()) {
                                        current_face.remove<gmds::Face>(faces[iFace]);
                                }
                        }
                }
        }

	if(this->mesh_.getModel().has(E2F)) {
                gmds::IGMesh::edge_iterator ite  = this->mesh_.edges_begin();
                for(;!ite.isDone();ite.next()) {
                        gmds::Edge current_edge = ite.value();

                        std::vector<gmds::Face> faces = current_edge.get<gmds::Face>();
                        for(int iFace=0; iFace<faces.size(); iFace++) {
                                if(faces2Remove.find(faces[iFace]) != faces2Remove.end()) {
                                        current_edge.remove<gmds::Face>(faces[iFace]);
                                }
                        }
                }
        }	

	if(this->mesh_.getModel().has(N2F)) {
                gmds::IGMesh::node_iterator itn  = this->mesh_.nodes_begin();
                for(;!itn.isDone();itn.next()) {
                        gmds::Node current_node = itn.value();

                        std::vector<gmds::Face> faces = current_node.get<gmds::Face>();
                        for(int iFace=0; iFace<faces.size(); iFace++) {
                                if(faces2Remove.find(faces[iFace]) != faces2Remove.end()) {
                                        current_node.remove<gmds::Face>(faces[iFace]);
                                }
                        }
                }
        }

	// finally delete the faces
	std::set<gmds::Face>::iterator it = faces2Remove.begin();
	for(; it!=faces2Remove.end(); it++) {
		this->mesh_.deleteFace(*it);
	}
}
/*----------------------------------------------------------------------------*/
void
MeshInsertDetail::displayMeshQuality()
{
	IGMeshQualityEvaluation qualEval;

	// scaled jacobian
	Variable<double>* scaledJacobianRegionVar;
	if(this->mesh_.doesVariableExist(GMDS_REGION,"scaledJacobianRegionVar")) {
		scaledJacobianRegionVar = this->mesh_.getVariable<double>(GMDS_REGION,"scaledJacobianRegionVar");
	} else {
		scaledJacobianRegionVar = this->mesh_.newVariable<double>(GMDS_REGION,"scaledJacobianRegionVar");
	}
	
	double minScaledJacobian =  HUGE_VALF;
	double maxScaledJacobian = -HUGE_VALF;
	double meanScaledJacobian = 0;

        gmds::IGMesh::region_iterator itr  = this->mesh_.regions_begin();
        for(;!itr.isDone();itr.next()) {
        	gmds::Region current_region = itr.value();

		double scaledJacobian = current_region.computeScaledJacobian();
		//double scaledJacobian = qualEval.scaledJacobian(current_region);
		if(scaledJacobian < minScaledJacobian) minScaledJacobian = scaledJacobian;
		if(scaledJacobian > maxScaledJacobian) maxScaledJacobian = scaledJacobian;
		meanScaledJacobian += scaledJacobian;

		(*scaledJacobianRegionVar)[current_region.getID()] = scaledJacobian;
	}
	meanScaledJacobian = meanScaledJacobian / this->mesh_.getNbRegions();

	std::cout<<"MeshQuality scaled jacobian min "<<minScaledJacobian<<" max "<<maxScaledJacobian<<" mean "<<meanScaledJacobian<<std::endl;


}
/*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/
