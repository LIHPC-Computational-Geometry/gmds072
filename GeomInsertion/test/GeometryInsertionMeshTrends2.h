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
 * GeometryInsertionMeshTrends2.h
 *
 *  Created on: 14 april 2015
 *      Author: legoff
 */
/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
/*----------------------------------------------------------------------------*/
#include <GMDS/MeshInsertDetailInOut.h>
#include <GMDS/CAD/FACFacetedGeomReadAndWrite.h>
#include <GMDS/IO/VTKFacetedGeomReadAndWrite.h>
#include <GMDS/GeomMeshIntersectionService.h>
#include <GMDS/FacetedMeshIntersectionService.h>
#include <GMDS/LaplacianSmoothing.h>
#include <GMDS/MeshSplitter3Refinement.h>
#include <GMDS/SmartDampenedLaplacianSmoothing.h>
#include <GMDS/SmartLaplacianSmoothing.h>
#include <GMDS/OrderedSmartLaplacianSmoothing.h>
#include <GMDS/NormalLaplacianSmoothing.h>
#include <GMDS/MeshModelAlgo.h>

#include <GMDS/Algo/SheetOperator.h>
#include <GMDS/IO/VTKWriter.h>
#include <GMDS/IO/VTKReader.h>
/*----------------------------------------------------------------------------*/
class GeometryInsertionMeshTrends2: public ::testing::Test {

  protected:
	GeometryInsertionMeshTrends2(){;}
    virtual ~GeometryInsertionMeshTrends2(){;}
};
/*----------------------------------------------------------------------------*/
TEST_F(GeometryInsertionMeshTrends2,DISABLED_cas_pyramid) {

    std::cout<<"=================================================="<<std::endl;
    std::cout<<"GeometryInsertionSheetTest::cas_pyramid"<<std::endl;

    std::string case_name("cas_pyramid");

    std::vector<gmds::geom::GeomVolume* > vols;
    
    gmds::MeshModel mod = gmds::DIM3|gmds::N|gmds::F|gmds::R|gmds::F2N|gmds::R2N|gmds::R2F|gmds::F2R|
      		gmds::E|gmds::E2N|gmds::F2E;
    gmds::IGMesh mesh(mod);

    gmds::geom::FacetedGeomManager model;
    gmds::VTKFacetedGeomReadAndWrite vtkReader;
    vtkReader.import(model,"/home/legoff/travail/GMDS/gmds_workingcopy/GeomInsertion/test/Samples/pyramid.vtk");

    model.getVolumes(vols);
    gmds::geom::FacetedVolume* vol = dynamic_cast<gmds::geom::FacetedVolume*> (vols[0]);
    vol->reorient();

    // generate initial grid
    {
    	const int nx = 20;
    	const int ny = 20;
    	const int nz = 20;

    	gmds::Node nodes[nx+1][ny+1][nz+1];

        double minXYZ[3], maxXYZ[3];
        vol->computeBoundingBox(minXYZ,maxXYZ);
/*
    	double xmin = minXYZ[0] - (maxXYZ[0] - minXYZ[0])/(nx-3);
    	double ymin = minXYZ[1] - (maxXYZ[1] - minXYZ[1])/(ny-3);
    	double zmin = minXYZ[2] - (maxXYZ[2] - minXYZ[2])/(nz-3);
    	double xmax = maxXYZ[0] + (maxXYZ[0] - minXYZ[0])/(nx-3);
    	double ymax = maxXYZ[1] + (maxXYZ[1] - minXYZ[1])/(ny-3);
    	double zmax = maxXYZ[2] + (maxXYZ[2] - minXYZ[2])/(nz-3);
*/
	double xmin = minXYZ[0] - (maxXYZ[0] - minXYZ[0])/(nx-9);
        double ymin = minXYZ[1] - (maxXYZ[1] - minXYZ[1])/(ny-9);
        double zmin = minXYZ[2] - (maxXYZ[2] - minXYZ[2])/(nz-9);
        double xmax = maxXYZ[0] + (maxXYZ[0] - minXYZ[0])/(nx-9);
        double ymax = maxXYZ[1] + (maxXYZ[1] - minXYZ[1])/(ny-9);
        double zmax = maxXYZ[2] + (maxXYZ[2] - minXYZ[2])/(nz-9);

    	double dx = (xmax - xmin) / nx;
    	double dy = (ymax - ymin) / ny;
    	double dz = (zmax - zmin) / nz;

    	for(int i=0; i<=nx; i++) {
    		for(int j=0; j<=ny; j++) {
    			for(int k=0; k<=nz; k++) {
    				nodes[i][j][k] = mesh.newNode(xmin+i*dx,ymin+j*dy,zmin+k*dz);
    			}
    		}
    	}

    	for(int i=0; i<nx; i++) {
    		for(int j=0; j<ny; j++) {
    			for(int k=0; k<nz; k++) {
    				mesh.newHex(
    						nodes[i][j][k+1],
    						nodes[i][j+1][k+1],
    						nodes[i+1][j+1][k+1],
    						nodes[i+1][j][k+1],
    						nodes[i][j][k],
    						nodes[i][j+1][k],
    						nodes[i+1][j+1][k],
    						nodes[i+1][j][k]
    				);
    			}
    		}
    	}

    }

    gmds::geom::GeomMeshIntersectionService*  service =
                                new gmds::geom::FacetedMeshIntersectionService();

    std::vector<gmds::geom::GeomSurface* > surfs;
    model.getSurfaces(surfs);
    service->initialization(surfs);
    service->buildAABBSurfacesTriangulationTree(surfs);
//    gmds::geom::FacetedVolume* vol = dynamic_cast<gmds::geom::FacetedVolume*> (vols[0]);

    gmds::MeshSplitter3Refinement splitter(mesh,model,*service);
    splitter.splitMesh();
    splitter.exportVTK(case_name+std::string("splitMesh"));

    // change model
    gmds::MeshModel mod_bis = gmds::DIM3|gmds::N|gmds::F|gmds::R|gmds::F2N|gmds::R2N|gmds::R2F|gmds::F2R|gmds::N2F|gmds::N2R;
    mesh.changeModel(mod_bis,false);

    gmds::IGMeshDoctor doc(&mesh);
    doc.buildFacesAndR2F();
//    doc.buildEdgesAndX2E();
    doc.updateUpwardConnectivity();

    mesh.initializeGeometryClassification();

    gmds::MeshInsertDetailInOut insertion(mesh,model,*service);

    insertion.exportModelVTK(case_name+std::string("_model"));

    insertion.autoInsert(vols[0]);

    insertion.exportMeshVTK(case_name,gmds::N|gmds::E|gmds::F|gmds::R);
    insertion.exportCurvesEdgesVTK(case_name+std::string("_curves"));
    insertion.project();
    insertion.exportCurvesEdgesVTK(case_name+std::string("_projectCurves"));
	insertion.displayMeshQuality();

    insertion.exportMeshVTK(case_name+std::string("_project"),gmds::N|gmds::E|gmds::F|gmds::R);

    gmds::MeshModelAlgo meshModelAlgo(mesh,model);
    int markBoundaryNodes = mesh.getNewMark<gmds::Node>();
    meshModelAlgo.markBoundaryNodes(markBoundaryNodes);    

    //gmds::OrderedSmartLaplacianSmoothing* smooth = new gmds::OrderedSmartLaplacianSmoothing(mesh,model,*service);
    gmds::SmartLaplacianSmoothing* smooth = new gmds::SmartLaplacianSmoothing(mesh,model,*service);
    //gmds::LaplacianSmoothing* smooth = new gmds::LaplacianSmoothing(mesh,model,*service);
    smooth->initNodesAdjacencies();
    for(int i=0; i<10; i++) {
        smooth->exec(markBoundaryNodes);
    }
    insertion.exportMeshVTK(case_name+std::string("_smartsmooth"),gmds::N|gmds::E|gmds::F|gmds::R);

    std::cout<<"mesh.getNbFaces() "<<mesh.getNbFaces()<<std::endl;
    gmds::MeshModel mod_AAAAAA = gmds::DIM3|gmds::N|gmds::F|gmds::R|gmds::F2N|gmds::R2N|gmds::R2F|gmds::F2R|
                gmds::E|gmds::E2N|gmds::F2E;
    mesh.changeModel(mod_AAAAAA,false);
	std::cout<<"mesh.getNbEdges() "<<mesh.getNbEdges()<<std::endl;
    //gmds::Node node2Refine = mesh.get<gmds::Node>(574);
    //gmds::Face face2Refine = mesh.get<gmds::Face>(1292);
    //splitter.addEdges2Node(node2Refine,face2Refine);
    //splitter.refineQuads2EdgesOnCurve();
    splitter.refineQuads2EdgesOnCurveImproved();
    //gmds::MeshModelAlgo meshModelAlgo(mesh,model);
    meshModelAlgo.associateNodes();
    mesh.changeModel(mod_bis,false);
    doc.buildFacesAndR2F();  
    doc.updateUpwardConnectivity();
 
    std::cout<<"mesh.getNbFaces() "<<mesh.getNbFaces()<<std::endl;
    insertion.exportMeshVTK(case_name+std::string("_projectrefine"),gmds::N|gmds::E|gmds::F|gmds::R);
    mesh.changeModel(mod_bis,false);

//    insertion.exportMeshVTK("cas_test_cas3_project.mli",gmds::N|gmds::E|gmds::F|gmds::R);
/*
	gmds::IGMesh::node_iterator itn  = mesh.nodes_begin();
	for(;!itn.isDone();itn.next()) {
		gmds::Node current_node = itn.value();
		if(current_node.getID() != 5687) mesh.mark(current_node,markBoundaryNodes);
	}
*/
    smooth->initNodesAdjacencies();
    for(int i=0; i<10; i++) {
    	smooth->exec(markBoundaryNodes);
    }
	insertion.displayMeshQuality();
    insertion.exportMeshVTK(case_name+std::string("_projectrefinesmartsmooth"),gmds::N|gmds::E|gmds::F|gmds::R);

    // pillowing
    gmds::MeshModel mod_ter = gmds::DIM3|gmds::N|gmds::F|gmds::R|gmds::F2N|gmds::R2N|gmds::R2F|gmds::F2R|gmds::N2R;
    mesh.changeModel(mod_ter,false);

    insertion.removeUnassociatedFaces();
    insertion.pillow(vols[0]);
    insertion.removeUnassociatedFaces();
    insertion.pillowExt(vols[0]);
	insertion.displayMeshQuality();
    insertion.exportMeshVTK(case_name+std::string("_pillow"),gmds::N|gmds::E|gmds::F|gmds::R);

    smooth->initNodesAdjacencies();
    for(int i=0; i<10; i++) {
        smooth->exec(markBoundaryNodes);
    }

//    gmds::Node node2Refine = mesh.get<gmds::Node>(1633);
//    gmds::Face face2Refine = mesh.get<gmds::Face>(2007);
//    splitter.addEdges2Node(node2Refine,face2Refine);

    insertion.displayMeshQuality();

    insertion.exportMeshVTK(case_name+std::string("_pillowsmooth"),gmds::N|gmds::F|gmds::R);
   
//    insertion.displayMeshQuality();
 
    EXPECT_EQ(0.0, 1.0);
}
/*----------------------------------------------------------------------------*/
TEST_F(GeometryInsertionMeshTrends2,DISABLED_cas_plate) {

    std::cout<<"=================================================="<<std::endl;
    std::cout<<"GeometryInsertionSheetTest::cas_plate"<<std::endl;

    std::string case_name("cas_plate");

    std::vector<gmds::geom::GeomVolume* > vols;
    
    gmds::MeshModel mod = gmds::DIM3|gmds::N|gmds::F|gmds::R|gmds::F2N|gmds::R2N|gmds::R2F|gmds::F2R|
      		gmds::E|gmds::E2N|gmds::F2E;
    gmds::IGMesh mesh(mod);

    gmds::geom::FacetedGeomManager model;
    gmds::VTKFacetedGeomReadAndWrite vtkReader;
    vtkReader.import(model,"/home/legoff/travail/GMDS/gmds_workingcopy/GeomInsertion/test/Samples/Slotted_Angle_Plate.vtk");

    model.getVolumes(vols);
    gmds::geom::FacetedVolume* vol = dynamic_cast<gmds::geom::FacetedVolume*> (vols[0]);
    vol->reorient();

    // generate initial grid
    {
    	const int nx = 30;
    	const int ny = 30;
    	const int nz = 30;

    	gmds::Node nodes[nx+1][ny+1][nz+1];

        double minXYZ[3], maxXYZ[3];
        vol->computeBoundingBox(minXYZ,maxXYZ);
/*
    	double xmin = minXYZ[0] - (maxXYZ[0] - minXYZ[0])/(nx-3);
    	double ymin = minXYZ[1] - (maxXYZ[1] - minXYZ[1])/(ny-3);
    	double zmin = minXYZ[2] - (maxXYZ[2] - minXYZ[2])/(nz-3);
    	double xmax = maxXYZ[0] + (maxXYZ[0] - minXYZ[0])/(nx-3);
    	double ymax = maxXYZ[1] + (maxXYZ[1] - minXYZ[1])/(ny-3);
    	double zmax = maxXYZ[2] + (maxXYZ[2] - minXYZ[2])/(nz-3);
*/
	double xmin = minXYZ[0] - (maxXYZ[0] - minXYZ[0])/(nx-9);
        double ymin = minXYZ[1] - (maxXYZ[1] - minXYZ[1])/(ny-9);
        double zmin = minXYZ[2] - (maxXYZ[2] - minXYZ[2])/(nz-9);
        double xmax = maxXYZ[0] + (maxXYZ[0] - minXYZ[0])/(nx-9);
        double ymax = maxXYZ[1] + (maxXYZ[1] - minXYZ[1])/(ny-9);
        double zmax = maxXYZ[2] + (maxXYZ[2] - minXYZ[2])/(nz-9);

    	double dx = (xmax - xmin) / nx;
    	double dy = (ymax - ymin) / ny;
    	double dz = (zmax - zmin) / nz;

    	for(int i=0; i<=nx; i++) {
    		for(int j=0; j<=ny; j++) {
    			for(int k=0; k<=nz; k++) {
    				nodes[i][j][k] = mesh.newNode(xmin+i*dx,ymin+j*dy,zmin+k*dz);
    			}
    		}
    	}

    	for(int i=0; i<nx; i++) {
    		for(int j=0; j<ny; j++) {
    			for(int k=0; k<nz; k++) {
    				mesh.newHex(
    						nodes[i][j][k+1],
    						nodes[i][j+1][k+1],
    						nodes[i+1][j+1][k+1],
    						nodes[i+1][j][k+1],
    						nodes[i][j][k],
    						nodes[i][j+1][k],
    						nodes[i+1][j+1][k],
    						nodes[i+1][j][k]
    				);
    			}
    		}
    	}

    }

    gmds::geom::GeomMeshIntersectionService*  service =
                                new gmds::geom::FacetedMeshIntersectionService();

    std::vector<gmds::geom::GeomSurface* > surfs;
    model.getSurfaces(surfs);
    service->initialization(surfs);
    service->buildAABBSurfacesTriangulationTree(surfs);
//    gmds::geom::FacetedVolume* vol = dynamic_cast<gmds::geom::FacetedVolume*> (vols[0]);

    gmds::MeshSplitter3Refinement splitter(mesh,model,*service);
    //splitter.splitMesh();
    //splitter.exportVTK(case_name+std::string("splitMesh"));

    // change model
    gmds::MeshModel mod_bis = gmds::DIM3|gmds::N|gmds::F|gmds::R|gmds::F2N|gmds::R2N|gmds::R2F|gmds::F2R|gmds::N2F|gmds::N2R;
    mesh.changeModel(mod_bis,false);

    gmds::IGMeshDoctor doc(&mesh);
    doc.buildFacesAndR2F();
//    doc.buildEdgesAndX2E();
    doc.updateUpwardConnectivity();

    mesh.initializeGeometryClassification();

    gmds::MeshInsertDetailInOut insertion(mesh,model,*service);

    insertion.exportModelVTK(case_name+std::string("_model"));

    insertion.autoInsert(vols[0]);

    insertion.exportMeshVTK(case_name,gmds::N|gmds::E|gmds::F|gmds::R);
    insertion.exportCurvesEdgesVTK(case_name+std::string("_curves"));
    insertion.project();
    insertion.exportCurvesEdgesVTK(case_name+std::string("_projectCurves"));

	insertion.displayMeshQuality();

    insertion.exportMeshVTK(case_name+std::string("_project"),gmds::N|gmds::E|gmds::F|gmds::R);

    gmds::MeshModelAlgo meshModelAlgo(mesh,model);
    int markBoundaryNodes = mesh.getNewMark<gmds::Node>();
    meshModelAlgo.markBoundaryNodes(markBoundaryNodes);    

    //gmds::OrderedSmartLaplacianSmoothing* smooth = new gmds::OrderedSmartLaplacianSmoothing(mesh,model,*service);
    gmds::SmartLaplacianSmoothing* smooth = new gmds::SmartLaplacianSmoothing(mesh,model,*service);
    //gmds::LaplacianSmoothing* smooth = new gmds::LaplacianSmoothing(mesh,model,*service);
    smooth->initNodesAdjacencies();
    for(int i=0; i<10; i++) {
        smooth->exec(markBoundaryNodes);
    }
    insertion.exportMeshVTK(case_name+std::string("_smartsmooth"),gmds::N|gmds::E|gmds::F|gmds::R);
    /*
    std::cout<<"mesh.getNbFaces() "<<mesh.getNbFaces()<<std::endl;
    gmds::MeshModel mod_AAAAAA = gmds::DIM3|gmds::N|gmds::F|gmds::R|gmds::F2N|gmds::R2N|gmds::R2F|gmds::F2R|
                gmds::E|gmds::E2N|gmds::F2E;
    mesh.changeModel(mod_AAAAAA,false);
	std::cout<<"mesh.getNbEdges() "<<mesh.getNbEdges()<<std::endl;
    //gmds::Node node2Refine = mesh.get<gmds::Node>(574);
    //gmds::Face face2Refine = mesh.get<gmds::Face>(1292);
    //splitter.addEdges2Node(node2Refine,face2Refine);
    //splitter.refineQuads2EdgesOnCurve();
    //splitter.refineQuads2EdgesOnCurveImproved();
    //gmds::MeshModelAlgo meshModelAlgo(mesh,model);
    meshModelAlgo.associateNodes();
    mesh.changeModel(mod_bis,false);
    doc.buildFacesAndR2F();  
    doc.updateUpwardConnectivity();
 
    std::cout<<"mesh.getNbFaces() "<<mesh.getNbFaces()<<std::endl;
    insertion.exportMeshVTK(case_name+std::string("_projectrefine"),gmds::N|gmds::E|gmds::F|gmds::R);
    mesh.changeModel(mod_bis,false);

//    insertion.exportMeshVTK("cas_test_cas3_project.mli",gmds::N|gmds::E|gmds::F|gmds::R);

    smooth->initNodesAdjacencies();
    for(int i=0; i<10; i++) {
    	smooth->exec(markBoundaryNodes);
    }
    insertion.exportMeshVTK(case_name+std::string("_projectrefinesmartsmooth"),gmds::N|gmds::E|gmds::F|gmds::R);
    
    // pillowing
    gmds::MeshModel mod_ter = gmds::DIM3|gmds::N|gmds::F|gmds::R|gmds::F2N|gmds::R2N|gmds::R2F|gmds::F2R|gmds::N2R;
    mesh.changeModel(mod_ter,false);

    insertion.removeUnassociatedFaces();
    insertion.pillow(vols[0]);
    //insertion.removeUnassociatedFaces();
    //insertion.pillowExt(vols[0]);
    insertion.exportMeshVTK(case_name+std::string("_pillow"),gmds::N|gmds::E|gmds::F|gmds::R);

    smooth->initNodesAdjacencies();
    for(int i=0; i<10; i++) {
        smooth->exec(markBoundaryNodes);
    }

//    gmds::Node node2Refine = mesh.get<gmds::Node>(1633);
//    gmds::Face face2Refine = mesh.get<gmds::Face>(2007);
//    splitter.addEdges2Node(node2Refine,face2Refine);

    insertion.displayMeshQuality();

    insertion.exportMeshVTK(case_name+std::string("_pillowsmooth"),gmds::N|gmds::F|gmds::R);
   
//    insertion.displayMeshQuality();
*/
    EXPECT_EQ(0.0, 1.0);
}
/*----------------------------------------------------------------------------*/
TEST_F(GeometryInsertionMeshTrends2,DISABLED_cas_disconnector) {

    std::cout<<"=================================================="<<std::endl;
    std::cout<<"GeometryInsertionSheetTest::cas_disconnector"<<std::endl;

    std::string case_name("cas_disconnector");

    std::vector<gmds::geom::GeomVolume* > vols;
    
    gmds::MeshModel mod = gmds::DIM3|gmds::N|gmds::F|gmds::R|gmds::F2N|gmds::R2N|gmds::R2F|gmds::F2R|
      		gmds::E|gmds::E2N|gmds::F2E;
    gmds::IGMesh mesh(mod);

    gmds::geom::FacetedGeomManager model;
    gmds::VTKFacetedGeomReadAndWrite vtkReader;
    vtkReader.import(model,"/home/legoff/travail/GMDS/gmds_workingcopy/GeomInsertion/test/Samples/M4-M16_Disconnector.vtk");

    model.getVolumes(vols);
    gmds::geom::FacetedVolume* vol = dynamic_cast<gmds::geom::FacetedVolume*> (vols[0]);
    vol->reorient();

    // generate initial grid
    {
    	const int nx = 50;
    	const int ny = 50;
    	const int nz = 30;

    	gmds::Node nodes[nx+1][ny+1][nz+1];

        double minXYZ[3], maxXYZ[3];
        vol->computeBoundingBox(minXYZ,maxXYZ);
/*
    	double xmin = minXYZ[0] - (maxXYZ[0] - minXYZ[0])/(nx-3);
    	double ymin = minXYZ[1] - (maxXYZ[1] - minXYZ[1])/(ny-3);
    	double zmin = minXYZ[2] - (maxXYZ[2] - minXYZ[2])/(nz-3);
    	double xmax = maxXYZ[0] + (maxXYZ[0] - minXYZ[0])/(nx-3);
    	double ymax = maxXYZ[1] + (maxXYZ[1] - minXYZ[1])/(ny-3);
    	double zmax = maxXYZ[2] + (maxXYZ[2] - minXYZ[2])/(nz-3);
*/
	double xmin = minXYZ[0] - (maxXYZ[0] - minXYZ[0])/(nx-9);
        double ymin = minXYZ[1] - (maxXYZ[1] - minXYZ[1])/(ny-9);
        double zmin = minXYZ[2] - (maxXYZ[2] - minXYZ[2])/(nz-9);
        double xmax = maxXYZ[0] + (maxXYZ[0] - minXYZ[0])/(nx-9);
        double ymax = maxXYZ[1] + (maxXYZ[1] - minXYZ[1])/(ny-9);
        double zmax = maxXYZ[2] + (maxXYZ[2] - minXYZ[2])/(nz-9);

    	double dx = (xmax - xmin) / nx;
    	double dy = (ymax - ymin) / ny;
    	double dz = (zmax - zmin) / nz;

    	for(int i=0; i<=nx; i++) {
    		for(int j=0; j<=ny; j++) {
    			for(int k=0; k<=nz; k++) {
    				nodes[i][j][k] = mesh.newNode(xmin+i*dx,ymin+j*dy,zmin+k*dz);
    			}
    		}
    	}

    	for(int i=0; i<nx; i++) {
    		for(int j=0; j<ny; j++) {
    			for(int k=0; k<nz; k++) {
    				mesh.newHex(
    						nodes[i][j][k+1],
    						nodes[i][j+1][k+1],
    						nodes[i+1][j+1][k+1],
    						nodes[i+1][j][k+1],
    						nodes[i][j][k],
    						nodes[i][j+1][k],
    						nodes[i+1][j+1][k],
    						nodes[i+1][j][k]
    				);
    			}
    		}
    	}

    }

    gmds::geom::GeomMeshIntersectionService*  service =
                                new gmds::geom::FacetedMeshIntersectionService();

    std::vector<gmds::geom::GeomSurface* > surfs;
    model.getSurfaces(surfs);
    service->initialization(surfs);
    service->buildAABBSurfacesTriangulationTree(surfs);
//    gmds::geom::FacetedVolume* vol = dynamic_cast<gmds::geom::FacetedVolume*> (vols[0]);

    gmds::MeshSplitter3Refinement splitter(mesh,model,*service);
    //splitter.splitMesh();
    //splitter.exportVTK(case_name+std::string("splitMesh"));

    // change model
    gmds::MeshModel mod_bis = gmds::DIM3|gmds::N|gmds::F|gmds::R|gmds::F2N|gmds::R2N|gmds::R2F|gmds::F2R|gmds::N2F|gmds::N2R;
    mesh.changeModel(mod_bis,false);

    gmds::IGMeshDoctor doc(&mesh);
    doc.buildFacesAndR2F();
//    doc.buildEdgesAndX2E();
    doc.updateUpwardConnectivity();

    mesh.initializeGeometryClassification();

    gmds::MeshInsertDetailInOut insertion(mesh,model,*service);

    insertion.exportModelVTK(case_name+std::string("_model"));

    insertion.autoInsert(vols[0]);

    insertion.exportMeshVTK(case_name,gmds::N|gmds::E|gmds::F|gmds::R);
    insertion.exportCurvesEdgesVTK(case_name+std::string("_curves"));
    insertion.project();
    insertion.exportCurvesEdgesVTK(case_name+std::string("_projectCurves"));

	insertion.displayMeshQuality();

    insertion.exportMeshVTK(case_name+std::string("_project"),gmds::N|gmds::E|gmds::F|gmds::R);

    gmds::MeshModelAlgo meshModelAlgo(mesh,model);
    int markBoundaryNodes = mesh.getNewMark<gmds::Node>();
    meshModelAlgo.markBoundaryNodes(markBoundaryNodes);    

    //gmds::OrderedSmartLaplacianSmoothing* smooth = new gmds::OrderedSmartLaplacianSmoothing(mesh,model,*service);
    gmds::SmartLaplacianSmoothing* smooth = new gmds::SmartLaplacianSmoothing(mesh,model,*service);
    //gmds::LaplacianSmoothing* smooth = new gmds::LaplacianSmoothing(mesh,model,*service);
    smooth->initNodesAdjacencies();
    for(int i=0; i<10; i++) {
        smooth->exec(markBoundaryNodes);
    }
    insertion.exportMeshVTK(case_name+std::string("_smartsmooth"),gmds::N|gmds::E|gmds::F|gmds::R);
    /*
    std::cout<<"mesh.getNbFaces() "<<mesh.getNbFaces()<<std::endl;
    gmds::MeshModel mod_AAAAAA = gmds::DIM3|gmds::N|gmds::F|gmds::R|gmds::F2N|gmds::R2N|gmds::R2F|gmds::F2R|
                gmds::E|gmds::E2N|gmds::F2E;
    mesh.changeModel(mod_AAAAAA,false);
	std::cout<<"mesh.getNbEdges() "<<mesh.getNbEdges()<<std::endl;
    //gmds::Node node2Refine = mesh.get<gmds::Node>(574);
    //gmds::Face face2Refine = mesh.get<gmds::Face>(1292);
    //splitter.addEdges2Node(node2Refine,face2Refine);
    //splitter.refineQuads2EdgesOnCurve();
    //splitter.refineQuads2EdgesOnCurveImproved();
    //gmds::MeshModelAlgo meshModelAlgo(mesh,model);
    meshModelAlgo.associateNodes();
    mesh.changeModel(mod_bis,false);
    doc.buildFacesAndR2F();  
    doc.updateUpwardConnectivity();
 
    std::cout<<"mesh.getNbFaces() "<<mesh.getNbFaces()<<std::endl;
    insertion.exportMeshVTK(case_name+std::string("_projectrefine"),gmds::N|gmds::E|gmds::F|gmds::R);
    mesh.changeModel(mod_bis,false);

//    insertion.exportMeshVTK("cas_test_cas3_project.mli",gmds::N|gmds::E|gmds::F|gmds::R);

    smooth->initNodesAdjacencies();
    for(int i=0; i<10; i++) {
    	smooth->exec(markBoundaryNodes);
    }
    insertion.exportMeshVTK(case_name+std::string("_projectrefinesmartsmooth"),gmds::N|gmds::E|gmds::F|gmds::R);
    
    // pillowing
    gmds::MeshModel mod_ter = gmds::DIM3|gmds::N|gmds::F|gmds::R|gmds::F2N|gmds::R2N|gmds::R2F|gmds::F2R|gmds::N2R;
    mesh.changeModel(mod_ter,false);

    insertion.removeUnassociatedFaces();
    insertion.pillow(vols[0]);
    //insertion.removeUnassociatedFaces();
    //insertion.pillowExt(vols[0]);
    insertion.exportMeshVTK(case_name+std::string("_pillow"),gmds::N|gmds::E|gmds::F|gmds::R);

    smooth->initNodesAdjacencies();
    for(int i=0; i<10; i++) {
        smooth->exec(markBoundaryNodes);
    }

//    gmds::Node node2Refine = mesh.get<gmds::Node>(1633);
//    gmds::Face face2Refine = mesh.get<gmds::Face>(2007);
//    splitter.addEdges2Node(node2Refine,face2Refine);

    insertion.displayMeshQuality();

    insertion.exportMeshVTK(case_name+std::string("_pillowsmooth"),gmds::N|gmds::F|gmds::R);
   
//    insertion.displayMeshQuality();
*/
    EXPECT_EQ(0.0, 1.0);
}
/*----------------------------------------------------------------------------*/
TEST_F(GeometryInsertionMeshTrends2,DISABLED_cas_ex1) {

    std::cout<<"=================================================="<<std::endl;
    std::cout<<"GeometryInsertionSheetTest::cas_ex1"<<std::endl;

    std::string case_name("cas_ex1");

    std::vector<gmds::geom::GeomVolume* > vols;
    
    gmds::MeshModel mod = gmds::DIM3|gmds::N|gmds::F|gmds::R|gmds::F2N|gmds::R2N|gmds::R2F|gmds::F2R|
      		gmds::E|gmds::E2N|gmds::F2E;
    gmds::IGMesh mesh(mod);

    gmds::geom::FacetedGeomManager model;
    gmds::VTKFacetedGeomReadAndWrite vtkReader;
    vtkReader.import(model,"/home/legoff/travail/GMDS/gmds_workingcopy/GeomInsertion/test/Samples/ex1.vtk");

    model.getVolumes(vols);
    gmds::geom::FacetedVolume* vol = dynamic_cast<gmds::geom::FacetedVolume*> (vols[0]);
    vol->reorient();

    // generate initial grid
    {
    	const int nx = 30;
    	const int ny = 30;
    	const int nz = 30;

    	gmds::Node nodes[nx+1][ny+1][nz+1];

        double minXYZ[3], maxXYZ[3];
        vol->computeBoundingBox(minXYZ,maxXYZ);
/*
    	double xmin = minXYZ[0] - (maxXYZ[0] - minXYZ[0])/(nx-3);
    	double ymin = minXYZ[1] - (maxXYZ[1] - minXYZ[1])/(ny-3);
    	double zmin = minXYZ[2] - (maxXYZ[2] - minXYZ[2])/(nz-3);
    	double xmax = maxXYZ[0] + (maxXYZ[0] - minXYZ[0])/(nx-3);
    	double ymax = maxXYZ[1] + (maxXYZ[1] - minXYZ[1])/(ny-3);
    	double zmax = maxXYZ[2] + (maxXYZ[2] - minXYZ[2])/(nz-3);
*/
	double xmin = minXYZ[0] - (maxXYZ[0] - minXYZ[0])/(nx-9);
        double ymin = minXYZ[1] - (maxXYZ[1] - minXYZ[1])/(ny-9);
        double zmin = minXYZ[2] - (maxXYZ[2] - minXYZ[2])/(nz-9);
        double xmax = maxXYZ[0] + (maxXYZ[0] - minXYZ[0])/(nx-9);
        double ymax = maxXYZ[1] + (maxXYZ[1] - minXYZ[1])/(ny-9);
        double zmax = maxXYZ[2] + (maxXYZ[2] - minXYZ[2])/(nz-9);

    	double dx = (xmax - xmin) / nx;
    	double dy = (ymax - ymin) / ny;
    	double dz = (zmax - zmin) / nz;

    	for(int i=0; i<=nx; i++) {
    		for(int j=0; j<=ny; j++) {
    			for(int k=0; k<=nz; k++) {
    				nodes[i][j][k] = mesh.newNode(xmin+i*dx,ymin+j*dy,zmin+k*dz);
    			}
    		}
    	}

    	for(int i=0; i<nx; i++) {
    		for(int j=0; j<ny; j++) {
    			for(int k=0; k<nz; k++) {
    				mesh.newHex(
    						nodes[i][j][k+1],
    						nodes[i][j+1][k+1],
    						nodes[i+1][j+1][k+1],
    						nodes[i+1][j][k+1],
    						nodes[i][j][k],
    						nodes[i][j+1][k],
    						nodes[i+1][j+1][k],
    						nodes[i+1][j][k]
    				);
    			}
    		}
    	}

    }

    gmds::geom::GeomMeshIntersectionService*  service =
                                new gmds::geom::FacetedMeshIntersectionService();

    std::vector<gmds::geom::GeomSurface* > surfs;
    model.getSurfaces(surfs);
    service->initialization(surfs);
    service->buildAABBSurfacesTriangulationTree(surfs);
//    gmds::geom::FacetedVolume* vol = dynamic_cast<gmds::geom::FacetedVolume*> (vols[0]);

    gmds::MeshSplitter3Refinement splitter(mesh,model,*service);
    //splitter.splitMesh();
    //splitter.exportVTK(case_name+std::string("splitMesh"));

    // change model
    gmds::MeshModel mod_bis = gmds::DIM3|gmds::N|gmds::F|gmds::R|gmds::F2N|gmds::R2N|gmds::R2F|gmds::F2R|gmds::N2F|gmds::N2R;
    mesh.changeModel(mod_bis,false);

    gmds::IGMeshDoctor doc(&mesh);
    doc.buildFacesAndR2F();
//    doc.buildEdgesAndX2E();
    doc.updateUpwardConnectivity();

    mesh.initializeGeometryClassification();

    gmds::MeshInsertDetailInOut insertion(mesh,model,*service);

    insertion.exportModelVTK(case_name+std::string("_model"));

    insertion.autoInsert(vols[0]);

    insertion.exportMeshVTK(case_name,gmds::N|gmds::E|gmds::F|gmds::R);
    insertion.exportCurvesEdgesVTK(case_name+std::string("_curves"));
    insertion.project();
    insertion.exportCurvesEdgesVTK(case_name+std::string("_projectCurves"));

	insertion.displayMeshQuality();

    insertion.exportMeshVTK(case_name+std::string("_project"),gmds::N|gmds::E|gmds::F|gmds::R);

    gmds::MeshModelAlgo meshModelAlgo(mesh,model);
    int markBoundaryNodes = mesh.getNewMark<gmds::Node>();
    meshModelAlgo.markBoundaryNodes(markBoundaryNodes);    

    //gmds::OrderedSmartLaplacianSmoothing* smooth = new gmds::OrderedSmartLaplacianSmoothing(mesh,model,*service);
    gmds::SmartLaplacianSmoothing* smooth = new gmds::SmartLaplacianSmoothing(mesh,model,*service);
    //gmds::LaplacianSmoothing* smooth = new gmds::LaplacianSmoothing(mesh,model,*service);
    smooth->initNodesAdjacencies();
    for(int i=0; i<10; i++) {
        smooth->exec(markBoundaryNodes);
    }
    insertion.exportMeshVTK(case_name+std::string("_smartsmooth"),gmds::N|gmds::E|gmds::F|gmds::R);
    /*
    std::cout<<"mesh.getNbFaces() "<<mesh.getNbFaces()<<std::endl;
    gmds::MeshModel mod_AAAAAA = gmds::DIM3|gmds::N|gmds::F|gmds::R|gmds::F2N|gmds::R2N|gmds::R2F|gmds::F2R|
                gmds::E|gmds::E2N|gmds::F2E;
    mesh.changeModel(mod_AAAAAA,false);
	std::cout<<"mesh.getNbEdges() "<<mesh.getNbEdges()<<std::endl;
    //gmds::Node node2Refine = mesh.get<gmds::Node>(574);
    //gmds::Face face2Refine = mesh.get<gmds::Face>(1292);
    //splitter.addEdges2Node(node2Refine,face2Refine);
    //splitter.refineQuads2EdgesOnCurve();
    //splitter.refineQuads2EdgesOnCurveImproved();
    //gmds::MeshModelAlgo meshModelAlgo(mesh,model);
    meshModelAlgo.associateNodes();
    mesh.changeModel(mod_bis,false);
    doc.buildFacesAndR2F();  
    doc.updateUpwardConnectivity();
 
    std::cout<<"mesh.getNbFaces() "<<mesh.getNbFaces()<<std::endl;
    insertion.exportMeshVTK(case_name+std::string("_projectrefine"),gmds::N|gmds::E|gmds::F|gmds::R);
    mesh.changeModel(mod_bis,false);

//    insertion.exportMeshVTK("cas_test_cas3_project.mli",gmds::N|gmds::E|gmds::F|gmds::R);

    smooth->initNodesAdjacencies();
    for(int i=0; i<10; i++) {
    	smooth->exec(markBoundaryNodes);
    }
    insertion.exportMeshVTK(case_name+std::string("_projectrefinesmartsmooth"),gmds::N|gmds::E|gmds::F|gmds::R);
    
    // pillowing
    gmds::MeshModel mod_ter = gmds::DIM3|gmds::N|gmds::F|gmds::R|gmds::F2N|gmds::R2N|gmds::R2F|gmds::F2R|gmds::N2R;
    mesh.changeModel(mod_ter,false);

    insertion.removeUnassociatedFaces();
    insertion.pillow(vols[0]);
    //insertion.removeUnassociatedFaces();
    //insertion.pillowExt(vols[0]);
    insertion.exportMeshVTK(case_name+std::string("_pillow"),gmds::N|gmds::E|gmds::F|gmds::R);

    smooth->initNodesAdjacencies();
    for(int i=0; i<10; i++) {
        smooth->exec(markBoundaryNodes);
    }

//    gmds::Node node2Refine = mesh.get<gmds::Node>(1633);
//    gmds::Face face2Refine = mesh.get<gmds::Face>(2007);
//    splitter.addEdges2Node(node2Refine,face2Refine);

    insertion.displayMeshQuality();

    insertion.exportMeshVTK(case_name+std::string("_pillowsmooth"),gmds::N|gmds::F|gmds::R);
   
//    insertion.displayMeshQuality();
*/
    EXPECT_EQ(0.0, 1.0);
}
/*----------------------------------------------------------------------------*/
TEST_F(GeometryInsertionMeshTrends2,DISABLED_cas_ex1thex) {

    std::cout<<"=================================================="<<std::endl;
    std::cout<<"GeometryInsertionSheetTest::cas_ex1thex"<<std::endl;

    std::string case_name("cas_ex1thex");

    std::vector<gmds::geom::GeomVolume* > vols;
    
    gmds::MeshModel mod = gmds::DIM3|gmds::N|gmds::F|gmds::R|gmds::F2N|gmds::R2N|gmds::R2F|gmds::F2R|
      		gmds::E|gmds::E2N|gmds::F2E;
    gmds::IGMesh mesh(mod);

    gmds::geom::FacetedGeomManager model;
    gmds::VTKFacetedGeomReadAndWrite vtkReader;
    vtkReader.import(model,"/home/legoff/travail/GMDS/gmds_workingcopy/GeomInsertion/test/Samples/ex1.vtk");

    model.getVolumes(vols);
    gmds::geom::FacetedVolume* vol = dynamic_cast<gmds::geom::FacetedVolume*> (vols[0]);
    vol->reorient();

    // generate initial grid
    {
      /*
    	int nx = 30;
    	int ny = 30;
    	int nz = 30;

    	gmds::Node nodes[nx+1][ny+1][nz+1];

        double minXYZ[3], maxXYZ[3];
        vol->computeBoundingBox(minXYZ,maxXYZ);

	double xmin = minXYZ[0] - (maxXYZ[0] - minXYZ[0])/(nx-9);
        double ymin = minXYZ[1] - (maxXYZ[1] - minXYZ[1])/(ny-9);
        double zmin = minXYZ[2] - (maxXYZ[2] - minXYZ[2])/(nz-9);
        double xmax = maxXYZ[0] + (maxXYZ[0] - minXYZ[0])/(nx-9);
        double ymax = maxXYZ[1] + (maxXYZ[1] - minXYZ[1])/(ny-9);
        double zmax = maxXYZ[2] + (maxXYZ[2] - minXYZ[2])/(nz-9);

    	double dx = (xmax - xmin) / nx;
    	double dy = (ymax - ymin) / ny;
    	double dz = (zmax - zmin) / nz;

    	for(int i=0; i<=nx; i++) {
    		for(int j=0; j<=ny; j++) {
    			for(int k=0; k<=nz; k++) {
    				nodes[i][j][k] = mesh.newNode(xmin+i*dx,ymin+j*dy,zmin+k*dz);
    			}
    		}
    	}

    	for(int i=0; i<nx; i++) {
    		for(int j=0; j<ny; j++) {
    			for(int k=0; k<nz; k++) {
    				mesh.newHex(
    						nodes[i][j][k+1],
    						nodes[i][j+1][k+1],
    						nodes[i+1][j+1][k+1],
    						nodes[i+1][j][k+1],
    						nodes[i][j][k],
    						nodes[i][j+1][k],
    						nodes[i+1][j+1][k],
    						nodes[i+1][j][k]
    				);
    			}
    		}
    	}
	*/

      gmds::VTKReader<gmds::IGMesh> reader(mesh);
      reader.read("/home/legoff/travail/GMDS/gmds_workingcopy/GeomInsertion/test/Samples/box_Stepthex.vtk");
    }

    gmds::geom::GeomMeshIntersectionService*  service =
                                new gmds::geom::FacetedMeshIntersectionService();

    std::vector<gmds::geom::GeomSurface* > surfs;
    model.getSurfaces(surfs);
    service->initialization(surfs);
    service->buildAABBSurfacesTriangulationTree(surfs);
//    gmds::geom::FacetedVolume* vol = dynamic_cast<gmds::geom::FacetedVolume*> (vols[0]);

    gmds::MeshSplitter3Refinement splitter(mesh,model,*service);
    //splitter.splitMesh();
    //splitter.exportVTK(case_name+std::string("splitMesh"));

    // change model
    gmds::MeshModel mod_bis = gmds::DIM3|gmds::N|gmds::F|gmds::R|gmds::F2N|gmds::R2N|gmds::R2F|gmds::F2R|gmds::N2F|gmds::N2R;
    mesh.changeModel(mod_bis,false);

    gmds::IGMeshDoctor doc(&mesh);
    doc.buildFacesAndR2F();
//    doc.buildEdgesAndX2E();
    doc.updateUpwardConnectivity();

    mesh.initializeGeometryClassification();

    gmds::MeshInsertDetailInOut insertion(mesh,model,*service);

    insertion.exportModelVTK(case_name+std::string("_model"));

    insertion.autoInsert(vols[0]);

    insertion.exportMeshVTK(case_name,gmds::N|gmds::E|gmds::F|gmds::R);
    insertion.exportCurvesEdgesVTK(case_name+std::string("_curves"));
    insertion.project();
    insertion.exportCurvesEdgesVTK(case_name+std::string("_projectCurves"));

	insertion.displayMeshQuality();

    insertion.exportMeshVTK(case_name+std::string("_project"),gmds::N|gmds::E|gmds::F|gmds::R);

    gmds::MeshModelAlgo meshModelAlgo(mesh,model);
    int markBoundaryNodes = mesh.getNewMark<gmds::Node>();
    meshModelAlgo.markBoundaryNodes(markBoundaryNodes);    

    //gmds::OrderedSmartLaplacianSmoothing* smooth = new gmds::OrderedSmartLaplacianSmoothing(mesh,model,*service);
    gmds::SmartLaplacianSmoothing* smooth = new gmds::SmartLaplacianSmoothing(mesh,model,*service);
    //gmds::LaplacianSmoothing* smooth = new gmds::LaplacianSmoothing(mesh,model,*service);
    smooth->initNodesAdjacencies();
    for(int i=0; i<10; i++) {
        smooth->exec(markBoundaryNodes);
    }
    insertion.exportMeshVTK(case_name+std::string("_smartsmooth"),gmds::N|gmds::E|gmds::F|gmds::R);
    /*
    std::cout<<"mesh.getNbFaces() "<<mesh.getNbFaces()<<std::endl;
    gmds::MeshModel mod_AAAAAA = gmds::DIM3|gmds::N|gmds::F|gmds::R|gmds::F2N|gmds::R2N|gmds::R2F|gmds::F2R|
                gmds::E|gmds::E2N|gmds::F2E;
    mesh.changeModel(mod_AAAAAA,false);
	std::cout<<"mesh.getNbEdges() "<<mesh.getNbEdges()<<std::endl;
    //gmds::Node node2Refine = mesh.get<gmds::Node>(574);
    //gmds::Face face2Refine = mesh.get<gmds::Face>(1292);
    //splitter.addEdges2Node(node2Refine,face2Refine);
    //splitter.refineQuads2EdgesOnCurve();
    //splitter.refineQuads2EdgesOnCurveImproved();
    //gmds::MeshModelAlgo meshModelAlgo(mesh,model);
    meshModelAlgo.associateNodes();
    mesh.changeModel(mod_bis,false);
    doc.buildFacesAndR2F();  
    doc.updateUpwardConnectivity();
 
    std::cout<<"mesh.getNbFaces() "<<mesh.getNbFaces()<<std::endl;
    insertion.exportMeshVTK(case_name+std::string("_projectrefine"),gmds::N|gmds::E|gmds::F|gmds::R);
    mesh.changeModel(mod_bis,false);

//    insertion.exportMeshVTK("cas_test_cas3_project.mli",gmds::N|gmds::E|gmds::F|gmds::R);

    smooth->initNodesAdjacencies();
    for(int i=0; i<10; i++) {
    	smooth->exec(markBoundaryNodes);
    }
    insertion.exportMeshVTK(case_name+std::string("_projectrefinesmartsmooth"),gmds::N|gmds::E|gmds::F|gmds::R);
    
    // pillowing
    gmds::MeshModel mod_ter = gmds::DIM3|gmds::N|gmds::F|gmds::R|gmds::F2N|gmds::R2N|gmds::R2F|gmds::F2R|gmds::N2R;
    mesh.changeModel(mod_ter,false);

    insertion.removeUnassociatedFaces();
    insertion.pillow(vols[0]);
    //insertion.removeUnassociatedFaces();
    //insertion.pillowExt(vols[0]);
    insertion.exportMeshVTK(case_name+std::string("_pillow"),gmds::N|gmds::E|gmds::F|gmds::R);

    smooth->initNodesAdjacencies();
    for(int i=0; i<10; i++) {
        smooth->exec(markBoundaryNodes);
    }

//    gmds::Node node2Refine = mesh.get<gmds::Node>(1633);
//    gmds::Face face2Refine = mesh.get<gmds::Face>(2007);
//    splitter.addEdges2Node(node2Refine,face2Refine);

    insertion.displayMeshQuality();

    insertion.exportMeshVTK(case_name+std::string("_pillowsmooth"),gmds::N|gmds::F|gmds::R);
   
//    insertion.displayMeshQuality();
*/
    EXPECT_EQ(0.0, 1.0);
}
/*----------------------------------------------------------------------------*/
TEST_F(GeometryInsertionMeshTrends2,DISABLED_cas_ex1pyr) {

    std::cout<<"=================================================="<<std::endl;
    std::cout<<"GeometryInsertionSheetTest::cas_ex1pyr"<<std::endl;

    std::string case_name("cas_ex1pyr");

    std::vector<gmds::geom::GeomVolume* > vols;
    
    gmds::MeshModel mod = gmds::DIM3|gmds::N|gmds::F|gmds::R|gmds::F2N|gmds::R2N|gmds::R2F|gmds::F2R|
      		gmds::E|gmds::E2N|gmds::F2E;
    gmds::IGMesh mesh(mod);

    gmds::geom::FacetedGeomManager model;
    gmds::VTKFacetedGeomReadAndWrite vtkReader;
    vtkReader.import(model,"/home/legoff/travail/GMDS/gmds_workingcopy/GeomInsertion/test/Samples/ex1.vtk");

    model.getVolumes(vols);
    gmds::geom::FacetedVolume* vol = dynamic_cast<gmds::geom::FacetedVolume*> (vols[0]);
    vol->reorient();

    // generate initial grid
    {
    	const int nx = 19;
    	const int ny = 19;
    	const int nz = 19;

    	gmds::Node nodes[nx+1][ny+1][nz+1];
	gmds::Node nodes_centers[nx][ny][nz];

        double minXYZ[3], maxXYZ[3];
        vol->computeBoundingBox(minXYZ,maxXYZ);
/*
    	double xmin = minXYZ[0] - (maxXYZ[0] - minXYZ[0])/(nx-3);
    	double ymin = minXYZ[1] - (maxXYZ[1] - minXYZ[1])/(ny-3);
    	double zmin = minXYZ[2] - (maxXYZ[2] - minXYZ[2])/(nz-3);
    	double xmax = maxXYZ[0] + (maxXYZ[0] - minXYZ[0])/(nx-3);
    	double ymax = maxXYZ[1] + (maxXYZ[1] - minXYZ[1])/(ny-3);
    	double zmax = maxXYZ[2] + (maxXYZ[2] - minXYZ[2])/(nz-3);
*/
	double xmin = minXYZ[0] - (maxXYZ[0] - minXYZ[0])/(nx-8);
        double ymin = minXYZ[1] - (maxXYZ[1] - minXYZ[1])/(ny-7);
        double zmin = minXYZ[2] - (maxXYZ[2] - minXYZ[2])/(nz-9);
        double xmax = maxXYZ[0] + (maxXYZ[0] - minXYZ[0])/(nx-9);
        double ymax = maxXYZ[1] + (maxXYZ[1] - minXYZ[1])/(ny-9);
        double zmax = maxXYZ[2] + (maxXYZ[2] - minXYZ[2])/(nz-9);

    	double dx = (xmax - xmin) / nx;
    	double dy = (ymax - ymin) / ny;
    	double dz = (zmax - zmin) / nz;

    	for(int i=0; i<=nx; i++) {
    		for(int j=0; j<=ny; j++) {
    			for(int k=0; k<=nz; k++) {
    				nodes[i][j][k] = mesh.newNode(xmin+i*dx,ymin+j*dy,zmin+k*dz);
    			}
    		}
    	}
	for(int i=0; i<nx; i++) {
	  for(int j=0; j<ny; j++) {
	    for(int k=0; k<nz; k++) {
	      nodes_centers[i][j][k] = mesh.newNode(xmin+i*dx+dx/2.1,ymin+j*dy+dy/2.2,zmin+k*dz+dz/2.15);
	    }
	  }
    	}

    	for(int i=0; i<nx; i++) {
    		for(int j=0; j<ny; j++) {
    			for(int k=0; k<nz; k++) {

				mesh.newPyramid(//top
					    nodes[i][j][k+1],
					    nodes[i][j+1][k+1],
					    nodes[i+1][j+1][k+1],
					    nodes[i+1][j][k+1],
					    nodes_centers[i][j][k]
				);
				mesh.newPyramid(//bottom
					    nodes[i][j][k],
					    nodes[i+1][j][k],
					    nodes[i+1][j+1][k],
					    nodes[i][j+1][k],
					    nodes_centers[i][j][k]
				);
				mesh.newPyramid(//left
					    nodes[i][j][k],
					    nodes[i][j+1][k],
					    nodes[i][j+1][k+1],
					    nodes[i][j][k+1],
					    nodes_centers[i][j][k]
				);
				mesh.newPyramid(//right
					    nodes[i+1][j+1][k],
					    nodes[i+1][j][k],
					    nodes[i+1][j][k+1],
					    nodes[i+1][j+1][k+1],
					    nodes_centers[i][j][k]
				);
				mesh.newPyramid(//front
					    nodes[i+1][j][k],
					    nodes[i][j][k],
					    nodes[i][j][k+1],
					    nodes[i+1][j][k+1],
					    nodes_centers[i][j][k]
				);
				mesh.newPyramid(//back
					    nodes[i][j+1][k],
					    nodes[i+1][j+1][k],
					    nodes[i+1][j+1][k+1],
					    nodes[i][j+1][k+1],
					    nodes_centers[i][j][k]
				);
				
    			}
    		}
    	}


    }

    gmds::geom::GeomMeshIntersectionService*  service =
                                new gmds::geom::FacetedMeshIntersectionService();

    std::vector<gmds::geom::GeomSurface* > surfs;
    model.getSurfaces(surfs);
    service->initialization(surfs);
    service->buildAABBSurfacesTriangulationTree(surfs);
//    gmds::geom::FacetedVolume* vol = dynamic_cast<gmds::geom::FacetedVolume*> (vols[0]);

    gmds::MeshSplitter3Refinement splitter(mesh,model,*service);
    //splitter.splitMesh();
    //splitter.exportVTK(case_name+std::string("splitMesh"));

    // change model
    gmds::MeshModel mod_bis = gmds::DIM3|gmds::N|gmds::F|gmds::R|gmds::F2N|gmds::R2N|gmds::R2F|gmds::F2R|gmds::N2F|gmds::N2R;
    mesh.changeModel(mod_bis,false);

    gmds::IGMeshDoctor doc(&mesh);
    doc.buildFacesAndR2F();
//    doc.buildEdgesAndX2E();
    doc.updateUpwardConnectivity();

    mesh.initializeGeometryClassification();

    gmds::MeshInsertDetailInOut insertion(mesh,model,*service);

    insertion.exportModelVTK(case_name+std::string("_model"));

    insertion.autoInsert(vols[0]);

    insertion.exportMeshVTK(case_name,gmds::N|gmds::E|gmds::F|gmds::R);
    insertion.exportCurvesEdgesVTK(case_name+std::string("_curves"));
    insertion.project();
    insertion.exportCurvesEdgesVTK(case_name+std::string("_projectCurves"));

	insertion.displayMeshQuality();

    insertion.exportMeshVTK(case_name+std::string("_project"),gmds::N|gmds::E|gmds::F|gmds::R);

    gmds::MeshModelAlgo meshModelAlgo(mesh,model);
    int markBoundaryNodes = mesh.getNewMark<gmds::Node>();
    meshModelAlgo.markBoundaryNodes(markBoundaryNodes);    

    //gmds::OrderedSmartLaplacianSmoothing* smooth = new gmds::OrderedSmartLaplacianSmoothing(mesh,model,*service);
    gmds::SmartLaplacianSmoothing* smooth = new gmds::SmartLaplacianSmoothing(mesh,model,*service);
    //gmds::LaplacianSmoothing* smooth = new gmds::LaplacianSmoothing(mesh,model,*service);
    smooth->initNodesAdjacencies();
    for(int i=0; i<10; i++) {
        smooth->exec(markBoundaryNodes);
    }
    insertion.exportMeshVTK(case_name+std::string("_smartsmooth"),gmds::N|gmds::E|gmds::F|gmds::R);
    /*
    std::cout<<"mesh.getNbFaces() "<<mesh.getNbFaces()<<std::endl;
    gmds::MeshModel mod_AAAAAA = gmds::DIM3|gmds::N|gmds::F|gmds::R|gmds::F2N|gmds::R2N|gmds::R2F|gmds::F2R|
                gmds::E|gmds::E2N|gmds::F2E;
    mesh.changeModel(mod_AAAAAA,false);
	std::cout<<"mesh.getNbEdges() "<<mesh.getNbEdges()<<std::endl;
    //gmds::Node node2Refine = mesh.get<gmds::Node>(574);
    //gmds::Face face2Refine = mesh.get<gmds::Face>(1292);
    //splitter.addEdges2Node(node2Refine,face2Refine);
    //splitter.refineQuads2EdgesOnCurve();
    //splitter.refineQuads2EdgesOnCurveImproved();
    //gmds::MeshModelAlgo meshModelAlgo(mesh,model);
    meshModelAlgo.associateNodes();
    mesh.changeModel(mod_bis,false);
    doc.buildFacesAndR2F();  
    doc.updateUpwardConnectivity();
 
    std::cout<<"mesh.getNbFaces() "<<mesh.getNbFaces()<<std::endl;
    insertion.exportMeshVTK(case_name+std::string("_projectrefine"),gmds::N|gmds::E|gmds::F|gmds::R);
    mesh.changeModel(mod_bis,false);

//    insertion.exportMeshVTK("cas_test_cas3_project.mli",gmds::N|gmds::E|gmds::F|gmds::R);

    smooth->initNodesAdjacencies();
    for(int i=0; i<10; i++) {
    	smooth->exec(markBoundaryNodes);
    }
    insertion.exportMeshVTK(case_name+std::string("_projectrefinesmartsmooth"),gmds::N|gmds::E|gmds::F|gmds::R);
    
    // pillowing
    gmds::MeshModel mod_ter = gmds::DIM3|gmds::N|gmds::F|gmds::R|gmds::F2N|gmds::R2N|gmds::R2F|gmds::F2R|gmds::N2R;
    mesh.changeModel(mod_ter,false);

    insertion.removeUnassociatedFaces();
    insertion.pillow(vols[0]);
    //insertion.removeUnassociatedFaces();
    //insertion.pillowExt(vols[0]);
    insertion.exportMeshVTK(case_name+std::string("_pillow"),gmds::N|gmds::E|gmds::F|gmds::R);

    smooth->initNodesAdjacencies();
    for(int i=0; i<10; i++) {
        smooth->exec(markBoundaryNodes);
    }

//    gmds::Node node2Refine = mesh.get<gmds::Node>(1633);
//    gmds::Face face2Refine = mesh.get<gmds::Face>(2007);
//    splitter.addEdges2Node(node2Refine,face2Refine);

    insertion.displayMeshQuality();

    insertion.exportMeshVTK(case_name+std::string("_pillowsmooth"),gmds::N|gmds::F|gmds::R);
   
//    insertion.displayMeshQuality();
*/
    EXPECT_EQ(0.0, 1.0);
}
/*----------------------------------------------------------------------------*/
TEST_F(GeometryInsertionMeshTrends2,DISABLED_cas_bunny) {

    std::cout<<"=================================================="<<std::endl;
    std::cout<<"GeometryInsertionSheetTest::cas_bunny"<<std::endl;

    std::string case_name("cas_bunny");

    std::vector<gmds::geom::GeomVolume* > vols;
    
    gmds::MeshModel mod = gmds::DIM3|gmds::N|gmds::F|gmds::R|gmds::F2N|gmds::R2N|gmds::R2F|gmds::F2R|
      		gmds::E|gmds::E2N|gmds::F2E;
    gmds::IGMesh mesh(mod);

    gmds::geom::FacetedGeomManager model;
    gmds::VTKFacetedGeomReadAndWrite vtkReader;
    vtkReader.build(model,"/home/legoff/travail/GMDS/gmds_workingcopy/GeomInsertion/test/Samples/bunny.vtk",true);

    model.getVolumes(vols);
    gmds::geom::FacetedVolume* vol = dynamic_cast<gmds::geom::FacetedVolume*> (vols[0]);
    //vol->reorient();

    // generate initial grid
    {
    	const int nx = 50;
    	const int ny = 50;
    	const int nz = 50;

    	gmds::Node nodes[nx+1][ny+1][nz+1];

        double minXYZ[3], maxXYZ[3];
        vol->computeBoundingBox(minXYZ,maxXYZ);
/*
    	double xmin = minXYZ[0] - (maxXYZ[0] - minXYZ[0])/(nx-3);
    	double ymin = minXYZ[1] - (maxXYZ[1] - minXYZ[1])/(ny-3);
    	double zmin = minXYZ[2] - (maxXYZ[2] - minXYZ[2])/(nz-3);
    	double xmax = maxXYZ[0] + (maxXYZ[0] - minXYZ[0])/(nx-3);
    	double ymax = maxXYZ[1] + (maxXYZ[1] - minXYZ[1])/(ny-3);
    	double zmax = maxXYZ[2] + (maxXYZ[2] - minXYZ[2])/(nz-3);
*/
	double xmin = minXYZ[0] - (maxXYZ[0] - minXYZ[0])/(nx-9);
        double ymin = minXYZ[1] - (maxXYZ[1] - minXYZ[1])/(ny-9);
        double zmin = minXYZ[2] - (maxXYZ[2] - minXYZ[2])/(nz-9);
        double xmax = maxXYZ[0] + (maxXYZ[0] - minXYZ[0])/(nx-9);
        double ymax = maxXYZ[1] + (maxXYZ[1] - minXYZ[1])/(ny-9);
        double zmax = maxXYZ[2] + (maxXYZ[2] - minXYZ[2])/(nz-9);

    	double dx = (xmax - xmin) / nx;
    	double dy = (ymax - ymin) / ny;
    	double dz = (zmax - zmin) / nz;

    	for(int i=0; i<=nx; i++) {
    		for(int j=0; j<=ny; j++) {
    			for(int k=0; k<=nz; k++) {
    				nodes[i][j][k] = mesh.newNode(xmin+i*dx,ymin+j*dy,zmin+k*dz);
    			}
    		}
    	}

    	for(int i=0; i<nx; i++) {
    		for(int j=0; j<ny; j++) {
    			for(int k=0; k<nz; k++) {
    				mesh.newHex(
    						nodes[i][j][k+1],
    						nodes[i][j+1][k+1],
    						nodes[i+1][j+1][k+1],
    						nodes[i+1][j][k+1],
    						nodes[i][j][k],
    						nodes[i][j+1][k],
    						nodes[i+1][j+1][k],
    						nodes[i+1][j][k]
    				);
    			}
    		}
    	}

    }

    gmds::geom::GeomMeshIntersectionService*  service =
                                new gmds::geom::FacetedMeshIntersectionService();

    std::vector<gmds::geom::GeomSurface* > surfs;
    model.getSurfaces(surfs);
    service->initialization(surfs);
    service->buildAABBSurfacesTriangulationTree(surfs);
//    gmds::geom::FacetedVolume* vol = dynamic_cast<gmds::geom::FacetedVolume*> (vols[0]);

    gmds::MeshSplitter3Refinement splitter(mesh,model,*service);
    //splitter.splitMesh();
    //splitter.exportVTK(case_name+std::string("splitMesh"));

    // change model
    gmds::MeshModel mod_bis = gmds::DIM3|gmds::N|gmds::F|gmds::R|gmds::F2N|gmds::R2N|gmds::R2F|gmds::F2R|gmds::N2F|gmds::N2R;
    mesh.changeModel(mod_bis,false);

    gmds::IGMeshDoctor doc(&mesh);
    doc.buildFacesAndR2F();
//    doc.buildEdgesAndX2E();
    doc.updateUpwardConnectivity();

    mesh.initializeGeometryClassification();

    gmds::MeshInsertDetailInOut insertion(mesh,model,*service);

    insertion.exportModelVTK(case_name+std::string("_model"));

    insertion.autoInsert(vols[0]);

    insertion.exportMeshVTK(case_name,gmds::N|gmds::E|gmds::F|gmds::R);
    insertion.exportCurvesEdgesVTK(case_name+std::string("_curves"));
    insertion.project();
    insertion.exportCurvesEdgesVTK(case_name+std::string("_projectCurves"));

	insertion.displayMeshQuality();

    insertion.exportMeshVTK(case_name+std::string("_project"),gmds::N|gmds::E|gmds::F|gmds::R);

    gmds::MeshModelAlgo meshModelAlgo(mesh,model);
    int markBoundaryNodes = mesh.getNewMark<gmds::Node>();
    meshModelAlgo.markBoundaryNodes(markBoundaryNodes);    

    //gmds::OrderedSmartLaplacianSmoothing* smooth = new gmds::OrderedSmartLaplacianSmoothing(mesh,model,*service);
    gmds::SmartLaplacianSmoothing* smooth = new gmds::SmartLaplacianSmoothing(mesh,model,*service);
    //gmds::LaplacianSmoothing* smooth = new gmds::LaplacianSmoothing(mesh,model,*service);
    smooth->initNodesAdjacencies();
    for(int i=0; i<10; i++) {
        smooth->exec(markBoundaryNodes);
    }
    insertion.exportMeshVTK(case_name+std::string("_smartsmooth"),gmds::N|gmds::E|gmds::F|gmds::R);
    /*
    std::cout<<"mesh.getNbFaces() "<<mesh.getNbFaces()<<std::endl;
    gmds::MeshModel mod_AAAAAA = gmds::DIM3|gmds::N|gmds::F|gmds::R|gmds::F2N|gmds::R2N|gmds::R2F|gmds::F2R|
                gmds::E|gmds::E2N|gmds::F2E;
    mesh.changeModel(mod_AAAAAA,false);
	std::cout<<"mesh.getNbEdges() "<<mesh.getNbEdges()<<std::endl;
    //gmds::Node node2Refine = mesh.get<gmds::Node>(574);
    //gmds::Face face2Refine = mesh.get<gmds::Face>(1292);
    //splitter.addEdges2Node(node2Refine,face2Refine);
    //splitter.refineQuads2EdgesOnCurve();
    //splitter.refineQuads2EdgesOnCurveImproved();
    //gmds::MeshModelAlgo meshModelAlgo(mesh,model);
    meshModelAlgo.associateNodes();
    mesh.changeModel(mod_bis,false);
    doc.buildFacesAndR2F();  
    doc.updateUpwardConnectivity();
 
    std::cout<<"mesh.getNbFaces() "<<mesh.getNbFaces()<<std::endl;
    insertion.exportMeshVTK(case_name+std::string("_projectrefine"),gmds::N|gmds::E|gmds::F|gmds::R);
    mesh.changeModel(mod_bis,false);

//    insertion.exportMeshVTK("cas_test_cas3_project.mli",gmds::N|gmds::E|gmds::F|gmds::R);

    smooth->initNodesAdjacencies();
    for(int i=0; i<10; i++) {
    	smooth->exec(markBoundaryNodes);
    }
    insertion.exportMeshVTK(case_name+std::string("_projectrefinesmartsmooth"),gmds::N|gmds::E|gmds::F|gmds::R);
    
    // pillowing
    gmds::MeshModel mod_ter = gmds::DIM3|gmds::N|gmds::F|gmds::R|gmds::F2N|gmds::R2N|gmds::R2F|gmds::F2R|gmds::N2R;
    mesh.changeModel(mod_ter,false);

    insertion.removeUnassociatedFaces();
    insertion.pillow(vols[0]);
    //insertion.removeUnassociatedFaces();
    //insertion.pillowExt(vols[0]);
    insertion.exportMeshVTK(case_name+std::string("_pillow"),gmds::N|gmds::E|gmds::F|gmds::R);

    smooth->initNodesAdjacencies();
    for(int i=0; i<10; i++) {
        smooth->exec(markBoundaryNodes);
    }

//    gmds::Node node2Refine = mesh.get<gmds::Node>(1633);
//    gmds::Face face2Refine = mesh.get<gmds::Face>(2007);
//    splitter.addEdges2Node(node2Refine,face2Refine);

    insertion.displayMeshQuality();

    insertion.exportMeshVTK(case_name+std::string("_pillowsmooth"),gmds::N|gmds::F|gmds::R);
   
//    insertion.displayMeshQuality();
*/
    EXPECT_EQ(0.0, 1.0);
}
/*----------------------------------------------------------------------------*/
TEST_F(GeometryInsertionMeshTrends2,DISABLED_cas_c2) {

    std::cout<<"=================================================="<<std::endl;
    std::cout<<"GeometryInsertionSheetTest::cas_c2"<<std::endl;

    std::string case_name("cas_c2");

    std::vector<gmds::geom::GeomVolume* > vols;
    
    gmds::MeshModel mod = gmds::DIM3|gmds::N|gmds::F|gmds::R|gmds::F2N|gmds::R2N|gmds::R2F|gmds::F2R|
      		gmds::E|gmds::E2N|gmds::F2E;
    gmds::IGMesh mesh(mod);

    gmds::geom::FacetedGeomManager model;
    gmds::VTKFacetedGeomReadAndWrite vtkReader;
    vtkReader.import(model,"/home/legoff/travail/GMDS/gmds_workingcopy/GeomInsertion/test/Samples/cas2.vtk");

    model.getVolumes(vols);
    gmds::geom::FacetedVolume* vol = dynamic_cast<gmds::geom::FacetedVolume*> (vols[0]);
    vol->reorient();

    // generate initial grid
    {
    	const int nx = 10;
    	const int ny = 10;
    	const int nz = 10;

    	gmds::Node nodes[nx+1][ny+1][nz+1];

        double minXYZ[3], maxXYZ[3];
        vol->computeBoundingBox(minXYZ,maxXYZ);

    	double xmin = minXYZ[0] - (maxXYZ[0] - minXYZ[0])/(nx-3);
    	double ymin = minXYZ[1] - (maxXYZ[1] - minXYZ[1])/(ny-3);
    	double zmin = minXYZ[2] - (maxXYZ[2] - minXYZ[2])/(nz-3);
    	double xmax = maxXYZ[0] + (maxXYZ[0] - minXYZ[0])/(nx-3);
    	double ymax = maxXYZ[1] + (maxXYZ[1] - minXYZ[1])/(ny-3);
    	double zmax = maxXYZ[2] + (maxXYZ[2] - minXYZ[2])/(nz-3);
	/*
	double xmin = minXYZ[0] - (maxXYZ[0] - minXYZ[0])/(nx-9);
        double ymin = minXYZ[1] - (maxXYZ[1] - minXYZ[1])/(ny-9);
        double zmin = minXYZ[2] - (maxXYZ[2] - minXYZ[2])/(nz-9);
        double xmax = maxXYZ[0] + (maxXYZ[0] - minXYZ[0])/(nx-9);
        double ymax = maxXYZ[1] + (maxXYZ[1] - minXYZ[1])/(ny-9);
        double zmax = maxXYZ[2] + (maxXYZ[2] - minXYZ[2])/(nz-9);
	*/
    	double dx = (xmax - xmin) / nx;
    	double dy = (ymax - ymin) / ny;
    	double dz = (zmax - zmin) / nz;

    	for(int i=0; i<=nx; i++) {
    		for(int j=0; j<=ny; j++) {
    			for(int k=0; k<=nz; k++) {
    				nodes[i][j][k] = mesh.newNode(xmin+i*dx,ymin+j*dy,zmin+k*dz);
    			}
    		}
    	}

    	for(int i=0; i<nx; i++) {
    		for(int j=0; j<ny; j++) {
    			for(int k=0; k<nz; k++) {
    				mesh.newHex(
    						nodes[i][j][k+1],
    						nodes[i][j+1][k+1],
    						nodes[i+1][j+1][k+1],
    						nodes[i+1][j][k+1],
    						nodes[i][j][k],
    						nodes[i][j+1][k],
    						nodes[i+1][j+1][k],
    						nodes[i+1][j][k]
    				);
    			}
    		}
    	}

    }

    gmds::geom::GeomMeshIntersectionService*  service =
                                new gmds::geom::FacetedMeshIntersectionService();

    std::vector<gmds::geom::GeomSurface* > surfs;
    model.getSurfaces(surfs);
    service->initialization(surfs);
    service->buildAABBSurfacesTriangulationTree(surfs);
//    gmds::geom::FacetedVolume* vol = dynamic_cast<gmds::geom::FacetedVolume*> (vols[0]);

    gmds::MeshSplitter3Refinement splitter(mesh,model,*service);
    splitter.splitMesh();
    splitter.exportVTK(case_name+std::string("splitMesh"));

    // change model
    gmds::MeshModel mod_bis = gmds::DIM3|gmds::N|gmds::F|gmds::R|gmds::F2N|gmds::R2N|gmds::R2F|gmds::F2R|gmds::N2F|gmds::N2R;
    mesh.changeModel(mod_bis,false);

    gmds::IGMeshDoctor doc(&mesh);
    doc.buildFacesAndR2F();
//    doc.buildEdgesAndX2E();
    doc.updateUpwardConnectivity();

    mesh.initializeGeometryClassification();

    gmds::MeshInsertDetailInOut insertion(mesh,model,*service);

    insertion.exportModelVTK(case_name+std::string("_model"));

    insertion.autoInsert(vols[0]);

    insertion.exportMeshVTK(case_name,gmds::N|gmds::E|gmds::F|gmds::R);
    insertion.exportCurvesEdgesVTK(case_name+std::string("_curves"));
    insertion.project();
    insertion.exportCurvesEdgesVTK(case_name+std::string("_projectCurves"));

	insertion.displayMeshQuality();

    insertion.exportMeshVTK(case_name+std::string("_project"),gmds::N|gmds::E|gmds::F|gmds::R);

    gmds::MeshModelAlgo meshModelAlgo(mesh,model);
    int markBoundaryNodes = mesh.getNewMark<gmds::Node>();
    meshModelAlgo.markBoundaryNodes(markBoundaryNodes);    

    //gmds::OrderedSmartLaplacianSmoothing* smooth = new gmds::OrderedSmartLaplacianSmoothing(mesh,model,*service);
    gmds::SmartLaplacianSmoothing* smooth = new gmds::SmartLaplacianSmoothing(mesh,model,*service);
    //gmds::LaplacianSmoothing* smooth = new gmds::LaplacianSmoothing(mesh,model,*service);
    smooth->initNodesAdjacencies();
    for(int i=0; i<10; i++) {
        smooth->exec(markBoundaryNodes);
    }
    insertion.exportMeshVTK(case_name+std::string("_smartsmooth"),gmds::N|gmds::E|gmds::F|gmds::R);
    /*
    std::cout<<"mesh.getNbFaces() "<<mesh.getNbFaces()<<std::endl;
    gmds::MeshModel mod_AAAAAA = gmds::DIM3|gmds::N|gmds::F|gmds::R|gmds::F2N|gmds::R2N|gmds::R2F|gmds::F2R|
                gmds::E|gmds::E2N|gmds::F2E;
    mesh.changeModel(mod_AAAAAA,false);
	std::cout<<"mesh.getNbEdges() "<<mesh.getNbEdges()<<std::endl;
    //gmds::Node node2Refine = mesh.get<gmds::Node>(574);
    //gmds::Face face2Refine = mesh.get<gmds::Face>(1292);
    //splitter.addEdges2Node(node2Refine,face2Refine);
    //splitter.refineQuads2EdgesOnCurve();
    //splitter.refineQuads2EdgesOnCurveImproved();
    //gmds::MeshModelAlgo meshModelAlgo(mesh,model);
    meshModelAlgo.associateNodes();
    mesh.changeModel(mod_bis,false);
    doc.buildFacesAndR2F();  
    doc.updateUpwardConnectivity();
 
    std::cout<<"mesh.getNbFaces() "<<mesh.getNbFaces()<<std::endl;
    insertion.exportMeshVTK(case_name+std::string("_projectrefine"),gmds::N|gmds::E|gmds::F|gmds::R);
    mesh.changeModel(mod_bis,false);

//    insertion.exportMeshVTK("cas_test_cas3_project.mli",gmds::N|gmds::E|gmds::F|gmds::R);

    smooth->initNodesAdjacencies();
    for(int i=0; i<10; i++) {
    	smooth->exec(markBoundaryNodes);
    }
    insertion.exportMeshVTK(case_name+std::string("_projectrefinesmartsmooth"),gmds::N|gmds::E|gmds::F|gmds::R);
    
    // pillowing
    gmds::MeshModel mod_ter = gmds::DIM3|gmds::N|gmds::F|gmds::R|gmds::F2N|gmds::R2N|gmds::R2F|gmds::F2R|gmds::N2R;
    mesh.changeModel(mod_ter,false);

    insertion.removeUnassociatedFaces();
    insertion.pillow(vols[0]);
    //insertion.removeUnassociatedFaces();
    //insertion.pillowExt(vols[0]);
    insertion.exportMeshVTK(case_name+std::string("_pillow"),gmds::N|gmds::E|gmds::F|gmds::R);

    smooth->initNodesAdjacencies();
    for(int i=0; i<10; i++) {
        smooth->exec(markBoundaryNodes);
    }

//    gmds::Node node2Refine = mesh.get<gmds::Node>(1633);
//    gmds::Face face2Refine = mesh.get<gmds::Face>(2007);
//    splitter.addEdges2Node(node2Refine,face2Refine);

    insertion.displayMeshQuality();

    insertion.exportMeshVTK(case_name+std::string("_pillowsmooth"),gmds::N|gmds::F|gmds::R);
   
//    insertion.displayMeshQuality();
*/
    EXPECT_EQ(0.0, 1.0);
}
/*----------------------------------------------------------------------------*/
TEST_F(GeometryInsertionMeshTrends2,DISABLED_cas_fan) {

    std::cout<<"=================================================="<<std::endl;
    std::cout<<"GeometryInsertionSheetTest::cas_fan"<<std::endl;

    std::string case_name("cas_fan");

    std::vector<gmds::geom::GeomVolume* > vols;
    
    gmds::MeshModel mod = gmds::DIM3|gmds::N|gmds::F|gmds::R|gmds::F2N|gmds::R2N|gmds::R2F|gmds::F2R|
      		gmds::E|gmds::E2N|gmds::F2E;
    gmds::IGMesh mesh(mod);

    gmds::geom::FacetedGeomManager model;
    gmds::VTKFacetedGeomReadAndWrite vtkReader;
    vtkReader.import(model,"/home/legoff/travail/GMDS/gmds_workingcopy/GeomInsertion/test/Samples/fan.vtk");

    model.getVolumes(vols);
    gmds::geom::FacetedVolume* vol = dynamic_cast<gmds::geom::FacetedVolume*> (vols[0]);
    vol->reorient();

    // generate initial grid
    {
    	const int nx = 197;
    	const int ny = 197;
    	const int nz = 53;

    	gmds::Node nodes[nx+1][ny+1][nz+1];

        double minXYZ[3], maxXYZ[3];
        vol->computeBoundingBox(minXYZ,maxXYZ);

    	double xmin = minXYZ[0] - (maxXYZ[0] - minXYZ[0])/(nx-3);
    	double ymin = minXYZ[1] - (maxXYZ[1] - minXYZ[1])/(ny-3);
    	double zmin = minXYZ[2] - (maxXYZ[2] - minXYZ[2])/(nz-3);
    	double xmax = maxXYZ[0] + (maxXYZ[0] - minXYZ[0])/(nx-3);
    	double ymax = maxXYZ[1] + (maxXYZ[1] - minXYZ[1])/(ny-3);
    	double zmax = maxXYZ[2] + (maxXYZ[2] - minXYZ[2])/(nz-3);
	/*
	double xmin = minXYZ[0] - (maxXYZ[0] - minXYZ[0])/(nx-9);
        double ymin = minXYZ[1] - (maxXYZ[1] - minXYZ[1])/(ny-9);
        double zmin = minXYZ[2] - (maxXYZ[2] - minXYZ[2])/(nz-9);
        double xmax = maxXYZ[0] + (maxXYZ[0] - minXYZ[0])/(nx-9);
        double ymax = maxXYZ[1] + (maxXYZ[1] - minXYZ[1])/(ny-9);
        double zmax = maxXYZ[2] + (maxXYZ[2] - minXYZ[2])/(nz-9);
	*/
    	double dx = (xmax - xmin) / nx;
    	double dy = (ymax - ymin) / ny;
    	double dz = (zmax - zmin) / nz;

    	for(int i=0; i<=nx; i++) {
    		for(int j=0; j<=ny; j++) {
    			for(int k=0; k<=nz; k++) {
    				nodes[i][j][k] = mesh.newNode(xmin+i*dx,ymin+j*dy,zmin+k*dz);
    			}
    		}
    	}

    	for(int i=0; i<nx; i++) {
    		for(int j=0; j<ny; j++) {
    			for(int k=0; k<nz; k++) {
    				mesh.newHex(
    						nodes[i][j][k+1],
    						nodes[i][j+1][k+1],
    						nodes[i+1][j+1][k+1],
    						nodes[i+1][j][k+1],
    						nodes[i][j][k],
    						nodes[i][j+1][k],
    						nodes[i+1][j+1][k],
    						nodes[i+1][j][k]
    				);
    			}
    		}
    	}

    }

    gmds::geom::GeomMeshIntersectionService*  service =
                                new gmds::geom::FacetedMeshIntersectionService();

    std::vector<gmds::geom::GeomSurface* > surfs;
    model.getSurfaces(surfs);
    service->initialization(surfs);
    service->buildAABBSurfacesTriangulationTree(surfs);
//    gmds::geom::FacetedVolume* vol = dynamic_cast<gmds::geom::FacetedVolume*> (vols[0]);

    gmds::MeshSplitter3Refinement splitter(mesh,model,*service);
    //splitter.splitMesh();
    //splitter.exportVTK(case_name+std::string("splitMesh"));

    // change model
    gmds::MeshModel mod_bis = gmds::DIM3|gmds::N|gmds::F|gmds::R|gmds::F2N|gmds::R2N|gmds::R2F|gmds::F2R|gmds::N2F|gmds::N2R;
    mesh.changeModel(mod_bis,false);

    gmds::IGMeshDoctor doc(&mesh);
    doc.buildFacesAndR2F();
//    doc.buildEdgesAndX2E();
    doc.updateUpwardConnectivity();

    mesh.initializeGeometryClassification();

    gmds::MeshInsertDetailInOut insertion(mesh,model,*service);

    insertion.exportModelVTK(case_name+std::string("_model"));

    insertion.autoInsert(vols[0]);

    insertion.exportMeshVTK(case_name,gmds::N|gmds::E|gmds::F|gmds::R);
    insertion.exportCurvesEdgesVTK(case_name+std::string("_curves"));
    insertion.project();
    insertion.exportCurvesEdgesVTK(case_name+std::string("_projectCurves"));
/*
	insertion.displayMeshQuality();

    insertion.exportMeshVTK(case_name+std::string("_project"),gmds::N|gmds::E|gmds::F|gmds::R);

    gmds::MeshModelAlgo meshModelAlgo(mesh,model);
    int markBoundaryNodes = mesh.getNewMark<gmds::Node>();
    meshModelAlgo.markBoundaryNodes(markBoundaryNodes);    

    //gmds::OrderedSmartLaplacianSmoothing* smooth = new gmds::OrderedSmartLaplacianSmoothing(mesh,model,*service);
    gmds::SmartLaplacianSmoothing* smooth = new gmds::SmartLaplacianSmoothing(mesh,model,*service);
    //gmds::LaplacianSmoothing* smooth = new gmds::LaplacianSmoothing(mesh,model,*service);
    smooth->initNodesAdjacencies();
    for(int i=0; i<10; i++) {
        smooth->exec(markBoundaryNodes);
    }
    insertion.exportMeshVTK(case_name+std::string("_smartsmooth"),gmds::N|gmds::E|gmds::F|gmds::R);
    
    std::cout<<"mesh.getNbFaces() "<<mesh.getNbFaces()<<std::endl;
    gmds::MeshModel mod_AAAAAA = gmds::DIM3|gmds::N|gmds::F|gmds::R|gmds::F2N|gmds::R2N|gmds::R2F|gmds::F2R|
                gmds::E|gmds::E2N|gmds::F2E;
    mesh.changeModel(mod_AAAAAA,false);
	std::cout<<"mesh.getNbEdges() "<<mesh.getNbEdges()<<std::endl;
    //gmds::Node node2Refine = mesh.get<gmds::Node>(574);
    //gmds::Face face2Refine = mesh.get<gmds::Face>(1292);
    //splitter.addEdges2Node(node2Refine,face2Refine);
    //splitter.refineQuads2EdgesOnCurve();
    //splitter.refineQuads2EdgesOnCurveImproved();
    //gmds::MeshModelAlgo meshModelAlgo(mesh,model);
    meshModelAlgo.associateNodes();
    mesh.changeModel(mod_bis,false);
    doc.buildFacesAndR2F();  
    doc.updateUpwardConnectivity();
 
    std::cout<<"mesh.getNbFaces() "<<mesh.getNbFaces()<<std::endl;
    insertion.exportMeshVTK(case_name+std::string("_projectrefine"),gmds::N|gmds::E|gmds::F|gmds::R);
    mesh.changeModel(mod_bis,false);

//    insertion.exportMeshVTK("cas_test_cas3_project.mli",gmds::N|gmds::E|gmds::F|gmds::R);

    smooth->initNodesAdjacencies();
    for(int i=0; i<10; i++) {
    	smooth->exec(markBoundaryNodes);
    }
    insertion.exportMeshVTK(case_name+std::string("_projectrefinesmartsmooth"),gmds::N|gmds::E|gmds::F|gmds::R);
    
    // pillowing
    gmds::MeshModel mod_ter = gmds::DIM3|gmds::N|gmds::F|gmds::R|gmds::F2N|gmds::R2N|gmds::R2F|gmds::F2R|gmds::N2R;
    mesh.changeModel(mod_ter,false);

    insertion.removeUnassociatedFaces();
    insertion.pillow(vols[0]);
    //insertion.removeUnassociatedFaces();
    //insertion.pillowExt(vols[0]);
    insertion.exportMeshVTK(case_name+std::string("_pillow"),gmds::N|gmds::E|gmds::F|gmds::R);

    smooth->initNodesAdjacencies();
    for(int i=0; i<10; i++) {
        smooth->exec(markBoundaryNodes);
    }

//    gmds::Node node2Refine = mesh.get<gmds::Node>(1633);
//    gmds::Face face2Refine = mesh.get<gmds::Face>(2007);
//    splitter.addEdges2Node(node2Refine,face2Refine);

    insertion.displayMeshQuality();

    insertion.exportMeshVTK(case_name+std::string("_pillowsmooth"),gmds::N|gmds::F|gmds::R);
   
//    insertion.displayMeshQuality();
*/
    EXPECT_EQ(0.0, 1.0);
}
/*----------------------------------------------------------------------------*/
TEST_F(GeometryInsertionMeshTrends2,DISABLED_cas_plop) {

    std::cout<<"=================================================="<<std::endl;
    std::cout<<"GeometryInsertionSheetTest::cas_ex1"<<std::endl;

    std::vector<gmds::geom::GeomVolume* > vols;

    gmds::MeshModel mod = gmds::DIM3|gmds::N|gmds::F|gmds::R|gmds::F2N|gmds::R2N|gmds::R2F|gmds::F2R|
                gmds::E|gmds::E2N|gmds::F2E;
    gmds::IGMesh mesh(mod);

    gmds::geom::FacetedGeomManager model;
    gmds::VTKFacetedGeomReadAndWrite vtkReader;
    vtkReader.import(model,"/home/legoff/travail/GMDS/gmds_workingcopy/GeomInsertion/test/Samples/fan.vtk");

    model.getVolumes(vols);
    gmds::geom::FacetedVolume* vol = dynamic_cast<gmds::geom::FacetedVolume*> (vols[0]);
    vol->reorient();

    // generate initial grid
    {
    	const int nx = 197;
    	const int ny = 197;
    	const int nz = 53;

    	gmds::Node nodes[nx+1][ny+1][nz+1];

        double minXYZ[3], maxXYZ[3];
        vol->computeBoundingBox(minXYZ,maxXYZ);
/*
    	double xmin = minXYZ[0] - (maxXYZ[0] - minXYZ[0])/(nx-3);
    	double ymin = minXYZ[1] - (maxXYZ[1] - minXYZ[1])/(ny-3);
    	double zmin = minXYZ[2] - (maxXYZ[2] - minXYZ[2])/(nz-3);
    	double xmax = maxXYZ[0] + (maxXYZ[0] - minXYZ[0])/(nx-3);
    	double ymax = maxXYZ[1] + (maxXYZ[1] - minXYZ[1])/(ny-3);
    	double zmax = maxXYZ[2] + (maxXYZ[2] - minXYZ[2])/(nz-3);
*/
	double xmin = minXYZ[0] - (maxXYZ[0] - minXYZ[0])/(nx-9);
        double ymin = minXYZ[1] - (maxXYZ[1] - minXYZ[1])/(ny-9);
        double zmin = minXYZ[2] - (maxXYZ[2] - minXYZ[2])/(nz-9);
        double xmax = maxXYZ[0] + (maxXYZ[0] - minXYZ[0])/(nx-9);
        double ymax = maxXYZ[1] + (maxXYZ[1] - minXYZ[1])/(ny-9);
        double zmax = maxXYZ[2] + (maxXYZ[2] - minXYZ[2])/(nz-9);

    	double dx = (xmax - xmin) / nx;
    	double dy = (ymax - ymin) / ny;
    	double dz = (zmax - zmin) / nz;

    	for(int i=0; i<=nx; i++) {
    		for(int j=0; j<=ny; j++) {
    			for(int k=0; k<=nz; k++) {
    				nodes[i][j][k] = mesh.newNode(xmin+i*dx,ymin+j*dy,zmin+k*dz);
    			}
    		}
    	}

    	for(int i=0; i<nx; i++) {
    		for(int j=0; j<ny; j++) {
    			for(int k=0; k<nz; k++) {
    				mesh.newHex(
    						nodes[i][j][k+1],
    						nodes[i][j+1][k+1],
    						nodes[i+1][j+1][k+1],
    						nodes[i+1][j][k+1],
    						nodes[i][j][k],
    						nodes[i][j+1][k],
    						nodes[i+1][j+1][k],
    						nodes[i+1][j][k]
    				);
    			}
    		}
    	}

    }

    gmds::geom::GeomMeshIntersectionService*  service =
                                new gmds::geom::FacetedMeshIntersectionService();

    std::vector<gmds::geom::GeomSurface* > surfs;
    model.getSurfaces(surfs);
    service->initialization(surfs);
    service->buildAABBSurfacesTriangulationTree(surfs);
//    gmds::geom::FacetedVolume* vol = dynamic_cast<gmds::geom::FacetedVolume*> (vols[0]);

    gmds::MeshSplitter3Refinement splitter(mesh,model,*service);
    //splitter.splitMesh();
    //splitter.exportVTK("splitMesh.mli");

    // change model
    gmds::MeshModel mod_bis = gmds::DIM3|gmds::N|gmds::F|gmds::R|gmds::F2N|gmds::R2N|gmds::R2F|gmds::F2R|gmds::N2F|gmds::N2R;
    mesh.changeModel(mod_bis,false);

    gmds::IGMeshDoctor doc(&mesh);
    doc.buildFacesAndR2F();
//    doc.buildEdgesAndX2E();
    doc.updateUpwardConnectivity();

    mesh.initializeGeometryClassification();

    gmds::MeshInsertDetailInOut insertion(mesh,model,*service);

    insertion.exportModelVTK("cas_test_bunny");


    insertion.autoInsert(vols[0]);

    insertion.exportMeshVTK("cas_test_cas3.mli",gmds::N|gmds::E|gmds::F|gmds::R);

    insertion.project();

    insertion.displayMeshQuality();

    insertion.exportMeshVTK("cas_test_cas3_project.mli",gmds::N|gmds::E|gmds::F|gmds::R);


    gmds::MeshModelAlgo meshModelAlgo(mesh,model);
    int markBoundaryNodes = mesh.getNewMark<gmds::Node>();
    meshModelAlgo.markBoundaryNodes(markBoundaryNodes);    

    //gmds::OrderedSmartLaplacianSmoothing* smooth = new gmds::OrderedSmartLaplacianSmoothing(mesh,model,*service);
    gmds::SmartLaplacianSmoothing* smooth = new gmds::SmartLaplacianSmoothing(mesh,model,*service);
    //gmds::LaplacianSmoothing* smooth = new gmds::LaplacianSmoothing(mesh,model,*service);
    smooth->initNodesAdjacencies();
    for(int i=0; i<10; i++) {
        smooth->exec(markBoundaryNodes);
    }
    insertion.exportMeshVTK("cas_test_cas3_smartsmooth.mli",gmds::N|gmds::E|gmds::F|gmds::R);

    gmds::MeshModel mod_ter = gmds::DIM3|gmds::N|gmds::F|gmds::R|gmds::F2N|gmds::R2N|gmds::R2F|gmds::F2R|gmds::N2R;
    mesh.changeModel(mod_ter,false);

    insertion.removeUnassociatedFaces();
    insertion.pillow(vols[0]);
    //insertion.removeUnassociatedFaces();
    //insertion.pillowExt(vols[0]);
    insertion.exportMeshVTK("cas_test_cas3_pillow.mli",gmds::N|gmds::E|gmds::F|gmds::R);

    smooth->initNodesAdjacencies();
    for(int i=0; i<10; i++) {
        smooth->exec(markBoundaryNodes);
    }

//    gmds::Node node2Refine = mesh.get<gmds::Node>(1633);
//    gmds::Face face2Refine = mesh.get<gmds::Face>(2007);
//    splitter.addEdges2Node(node2Refine,face2Refine);

    insertion.displayMeshQuality();

    insertion.exportMeshVTK("cas_test_cas3_pillowsmooth.mli",gmds::N|gmds::F|gmds::R);

    EXPECT_EQ(0.0, 1.0);

}
/*----------------------------------------------------------------------------*/
TEST_F(GeometryInsertionMeshTrends2,cas_pyramid_EX) {

    std::cout<<"=================================================="<<std::endl;
    std::cout<<"GeometryInsertionSheetTest::cas_pyramid_EX"<<std::endl;

    std::string case_name("cas_pyramid");

    std::vector<gmds::geom::GeomVolume* > vols;
    
    gmds::MeshModel mod = gmds::DIM3|gmds::N|gmds::F|gmds::R|gmds::F2N|gmds::R2N|gmds::R2F|gmds::F2R|
      		gmds::E|gmds::E2N|gmds::F2E;
    gmds::IGMesh mesh(mod);

    gmds::geom::FacetedGeomManager model;
    gmds::VTKFacetedGeomReadAndWrite vtkReader;
    vtkReader.import(model,"/home/legoff/travail/GMDS/gmds_workingcopy/GeomInsertion/test/Samples/pyramid.vtk");

    model.getVolumes(vols);
    gmds::geom::FacetedVolume* vol = dynamic_cast<gmds::geom::FacetedVolume*> (vols[0]);
    vol->reorient();

    // generate initial grid
    {
    	const int nx = 20;
    	const int ny = 20;
    	const int nz = 20;

    	gmds::Node nodes[nx+1][ny+1][nz+1];

        double minXYZ[3], maxXYZ[3];
        vol->computeBoundingBox(minXYZ,maxXYZ);
/*
    	double xmin = minXYZ[0] - (maxXYZ[0] - minXYZ[0])/(nx-3);
    	double ymin = minXYZ[1] - (maxXYZ[1] - minXYZ[1])/(ny-3);
    	double zmin = minXYZ[2] - (maxXYZ[2] - minXYZ[2])/(nz-3);
    	double xmax = maxXYZ[0] + (maxXYZ[0] - minXYZ[0])/(nx-3);
    	double ymax = maxXYZ[1] + (maxXYZ[1] - minXYZ[1])/(ny-3);
    	double zmax = maxXYZ[2] + (maxXYZ[2] - minXYZ[2])/(nz-3);
*/
	double xmin = minXYZ[0] - (maxXYZ[0] - minXYZ[0])/(nx-9);
        double ymin = minXYZ[1] - (maxXYZ[1] - minXYZ[1])/(ny-9);
        double zmin = minXYZ[2] - (maxXYZ[2] - minXYZ[2])/(nz-9);
        double xmax = maxXYZ[0] + (maxXYZ[0] - minXYZ[0])/(nx-9);
        double ymax = maxXYZ[1] + (maxXYZ[1] - minXYZ[1])/(ny-9);
        double zmax = maxXYZ[2] + (maxXYZ[2] - minXYZ[2])/(nz-9);

    	double dx = (xmax - xmin) / nx;
    	double dy = (ymax - ymin) / ny;
    	double dz = (zmax - zmin) / nz;

    	for(int i=0; i<=nx; i++) {
    		for(int j=0; j<=ny; j++) {
    			for(int k=0; k<=nz; k++) {
    				nodes[i][j][k] = mesh.newNode(xmin+i*dx,ymin+j*dy,zmin+k*dz);
    			}
    		}
    	}

    	for(int i=0; i<nx; i++) {
    		for(int j=0; j<ny; j++) {
    			for(int k=0; k<nz; k++) {
    				mesh.newHex(
    						nodes[i][j][k+1],
    						nodes[i][j+1][k+1],
    						nodes[i+1][j+1][k+1],
    						nodes[i+1][j][k+1],
    						nodes[i][j][k],
    						nodes[i][j+1][k],
    						nodes[i+1][j+1][k],
    						nodes[i+1][j][k]
    				);
    			}
    		}
    	}

    }

    gmds::geom::GeomMeshIntersectionService*  service =
                                new gmds::geom::FacetedMeshIntersectionService();

    std::vector<gmds::geom::GeomSurface* > surfs;
    model.getSurfaces(surfs);
    service->initialization(surfs);
    service->buildAABBSurfacesTriangulationTree(surfs);
//    gmds::geom::FacetedVolume* vol = dynamic_cast<gmds::geom::FacetedVolume*> (vols[0]);

    gmds::MeshSplitter3Refinement splitter(mesh,model,*service);
    splitter.splitMesh();
    splitter.exportVTK(case_name+std::string("splitMesh"));

    // change model
    gmds::MeshModel mod_bis = gmds::DIM3|gmds::N|gmds::F|gmds::R|gmds::F2N|gmds::R2N|gmds::R2F|gmds::F2R|gmds::N2F|gmds::N2R;
    mesh.changeModel(mod_bis,false);

    gmds::IGMeshDoctor doc(&mesh);
    doc.buildFacesAndR2F();
//    doc.buildEdgesAndX2E();
    doc.updateUpwardConnectivity();

    mesh.initializeGeometryClassification();

    gmds::MeshInsertDetailInOut insertion(mesh,model,*service);

    insertion.exportModelVTK(case_name+std::string("_model"));

    insertion.autoInsert(vols[0]);

    insertion.exportMeshVTK(case_name,gmds::N|gmds::E|gmds::F|gmds::R);
    insertion.exportCurvesEdgesVTK(case_name+std::string("_curves"));
    insertion.project();
    insertion.exportCurvesEdgesVTK(case_name+std::string("_projectCurves"));
	insertion.displayMeshQuality();

    insertion.exportMeshVTK(case_name+std::string("_project"),gmds::N|gmds::E|gmds::F|gmds::R);

    gmds::MeshModelAlgo meshModelAlgo(mesh,model);
    int markBoundaryNodes = mesh.getNewMark<gmds::Node>();
    meshModelAlgo.markBoundaryNodes(markBoundaryNodes);    

    //gmds::OrderedSmartLaplacianSmoothing* smooth = new gmds::OrderedSmartLaplacianSmoothing(mesh,model,*service);
    gmds::SmartLaplacianSmoothing* smooth = new gmds::SmartLaplacianSmoothing(mesh,model,*service);
    //gmds::SmartDampenedLaplacianSmoothing* smooth = new gmds::SmartDampenedLaplacianSmoothing(mesh,model,*service);
    //gmds::NormalLaplacianSmoothing* smooth = new gmds::NormalLaplacianSmoothing(mesh,model,*service);
    //gmds::LaplacianSmoothing* smooth = new gmds::LaplacianSmoothing(mesh,model,*service);
    smooth->initNodesAdjacencies();
    for(int i=0; i<10; i++) {
        smooth->exec(markBoundaryNodes);
    }
    insertion.exportMeshVTK(case_name+std::string("_smartsmooth"),gmds::N|gmds::E|gmds::F|gmds::R);

    std::cout<<"mesh.getNbFaces() "<<mesh.getNbFaces()<<std::endl;
    gmds::MeshModel mod_AAAAAA = gmds::DIM3|gmds::N|gmds::F|gmds::R|gmds::F2N|gmds::R2N|gmds::R2F|gmds::F2R|
                gmds::E|gmds::E2N|gmds::F2E;
    mesh.changeModel(mod_AAAAAA,false);
	std::cout<<"mesh.getNbEdges() "<<mesh.getNbEdges()<<std::endl;
    //gmds::Node node2Refine = mesh.get<gmds::Node>(574);
    //gmds::Face face2Refine = mesh.get<gmds::Face>(1292);
    //splitter.addEdges2Node(node2Refine,face2Refine);
    //splitter.refineQuads2EdgesOnCurve();
    splitter.refineQuads2EdgesOnCurveImproved();
    //gmds::MeshModelAlgo meshModelAlgo(mesh,model);
    meshModelAlgo.associateNodes();
    mesh.changeModel(mod_bis,false);
    doc.buildFacesAndR2F();  
    doc.updateUpwardConnectivity();
 
    std::cout<<"mesh.getNbFaces() "<<mesh.getNbFaces()<<std::endl;
    insertion.exportMeshVTK(case_name+std::string("_projectrefine"),gmds::N|gmds::E|gmds::F|gmds::R);
    mesh.changeModel(mod_bis,false);

//    insertion.exportMeshVTK("cas_test_cas3_project.mli",gmds::N|gmds::E|gmds::F|gmds::R);
/*
	gmds::IGMesh::node_iterator itn  = mesh.nodes_begin();
	for(;!itn.isDone();itn.next()) {
		gmds::Node current_node = itn.value();
		if(current_node.getID() != 5687) mesh.mark(current_node,markBoundaryNodes);
	}
*/
    smooth->initNodesAdjacencies();
    for(int i=0; i<10; i++) {
       	smooth->exec(markBoundaryNodes);
    }
	insertion.displayMeshQuality();
    insertion.exportMeshVTK(case_name+std::string("_projectrefinesmartsmooth"),gmds::N|gmds::E|gmds::F|gmds::R);

    // pillowing
    gmds::MeshModel mod_ter = gmds::DIM3|gmds::N|gmds::F|gmds::R|gmds::F2N|gmds::R2N|gmds::R2F|gmds::F2R|gmds::N2R;
    mesh.changeModel(mod_ter,false);

    insertion.removeUnassociatedFaces();
    insertion.pillow(vols[0]);
	std::cout<<"pillowINT"<<std::endl;
    insertion.removeUnassociatedFaces();
    insertion.pillowExt(vols[0]);
	std::cout<<"pillowEXT"<<std::endl;
	insertion.displayMeshQuality();
    insertion.exportMeshVTK(case_name+std::string("_pillow"),gmds::N|gmds::E|gmds::F|gmds::R);
	
//	insertion.project();
	insertion.displayMeshQuality();
	insertion.exportMeshVTK(case_name+std::string("_poyop"),gmds::N|gmds::E|gmds::F|gmds::R);
    smooth->initNodesAdjacencies();
    for(int i=0; i<20; i++) {
	insertion.displayMeshQuality();
        smooth->exec(markBoundaryNodes);
    }

//	insertion.project();

//    gmds::Node node2Refine = mesh.get<gmds::Node>(1633);
//    gmds::Face face2Refine = mesh.get<gmds::Face>(2007);
//    splitter.addEdges2Node(node2Refine,face2Refine);

    insertion.displayMeshQuality();

    insertion.exportMeshVTK(case_name+std::string("_pillowsmooth"),gmds::N|gmds::F|gmds::R);

	//mesh.changeModel(mod_AAAAAA,false); 
	//meshModelAlgo.associateNodes();
	gmds::NormalLaplacianSmoothing* smooth2 = new gmds::NormalLaplacianSmoothing(mesh,model,*service);
	smooth2->initNodesAdjacencies();
	for(int i=0; i<20; i++) {
		insertion.displayMeshQuality();
		smooth2->exec(markBoundaryNodes,vol);
	}
 
    insertion.displayMeshQuality();
    insertion.exportMeshVTK(case_name+std::string("_pillowsmoothNormal"),gmds::N|gmds::F|gmds::R);
 
    EXPECT_EQ(0.0, 1.0);
}
/*----------------------------------------------------------------------------*/
